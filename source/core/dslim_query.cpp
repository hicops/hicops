/*
 * Copyright (C) 2019  Muhammad Haseeb, and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#include "dslim_fileout.h"
#include "msquery.h"
#include "dslim.h"
#include "lwqueue.h"
#include "lwbuff.h"
#include "scheduler.h"
#include "hicops_instr.hpp"

using namespace std;

extern gParams   params;
extern BYICount  *Score;
extern vector<string_t> queryfiles;

/* Global variables */
float_t *hyperscores         = NULL;
uchar_t *sCArr               = NULL;
BOOL   ExitSignal          = false;

#ifdef USE_MPI
DSLIM_Comm *CommHandle    = NULL;
hCell      *CandidatePSMS  = NULL;
#endif /* USE_MPI */

Scheduler  *SchedHandle    = NULL;
expeRT     *ePtrs          = NULL;

ebuffer    *iBuff          = NULL;
int_t         ciBuff         = -1;
lock_t       writer;

/* Lock for query file vector and thread manager */
lock_t qfilelock;
lwqueue<MSQuery *> *qfPtrs = NULL;
MSQuery **ptrs             = NULL;
int_t spectrumID             = 0;
int_t nBatches               = 0;
int_t dssize                 = 0;

/* Expt spectra data buffer */
lwbuff<Queries> *qPtrs     = NULL;
Queries *workPtr           = NULL;

#ifdef USE_MPI
VOID *DSLIM_FOut_Thread_Entry(VOID *argv);
#endif

/* A queue containing I/O thread state when preempted */
lwqueue<MSQuery *> *ioQ = NULL;
lock_t ioQlock;
std::atomic<bool> scheduler_init(false);

// DEBUG ONLY - Remove me
std::ostream& operator<<(std::ostream &out, const _heapEntry &c)
{
    out << "shd_pk: " << (int)c.sharedions<< ", idxoff: " << (int)c.idxoffset;
    out << ", total: " << c.totalions << ", psid: " << c.psid;
    out << ", pmass: " << c.pmass << ", hyp: " << c.hyperscore;
    return out;
}

//
// ----------------------------------------------------------------------------------
//

VOID *DSLIM_IO_Threads_Entry(VOID *argv);

/* Static function */
static BOOL DSLIM_BinarySearch(Index *, float_t, int_t&, int_t&);
static int_t  DSLIM_BinFindMin(pepEntry *entries, float_t pmass1, int_t min, int_t max);
static int_t  DSLIM_BinFindMax(pepEntry *entries, float_t pmass2, int_t min, int_t max);
static inline status_t DSLIM_WaitFor_IO(int_t &);
static inline status_t DSLIM_Deinit_IO();

/* FUNCTION: DSLIM_WaitFor_IO
 *
 * DESCRIPTION:
 *
 * INPUT:
 * @none
 *
 * OUTPUT:
 * @status: status of execution
 *
 */
static inline status_t DSLIM_WaitFor_IO(int_t &batchsize)
{
    status_t status;

    batchsize = 0;

    /* Wait for a I/O request */
    status = qPtrs->lockr_();

    /* Check for either a buffer or stopSignal */
    while (qPtrs->isEmptyReadyQ())
    {
        /* Check if endSignal has been received */
        if (SchedHandle->checkSignal())
        {
            status = qPtrs->unlockr_();
            status = ENDSIGNAL;
            break;
        }

        /* If both conditions fail,
         * the I/O threads still working */
        status = qPtrs->unlockr_();

        sleep(0.1);

        status = qPtrs->lockr_();
    }

    if (status == SLM_SUCCESS)
    {
        /* Get the I/O ptr from the wait queue */
        workPtr = qPtrs->getWorkPtr();

        if (workPtr == NULL)
        {
            status = ERR_INVLD_PTR;
        }

        status = qPtrs->unlockr_();

        batchsize = workPtr->numSpecs;
    }

    return status;
}

/* FUNCTION: DSLIM_SearchManager
 *
 * DESCRIPTION: Manages and performs the Peptide Search
 *
 * INPUT:
 * @slm_index : Pointer to the SLM_Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
static status_t DSLIM_InitializeMS2Data()
{
    status_t status = SLM_SUCCESS;

    int_t nfiles = queryfiles.size();
    ptrs = new MSQuery*[nfiles];

    /* Initialize the queue with already created nfiles */
    qfPtrs = new lwqueue<MSQuery*>(nfiles, false);

#ifdef USE_OMP
#pragma omp parallel for schedule (dynamic, 1)
#endif/* _OPENMP */
    for (auto fid = 0; fid < nfiles; fid++)
    {
        ptrs[fid] = new MSQuery;
        ptrs[fid]->InitQueryFile(&queryfiles[fid], fid);
    }

    /* Push zeroth as is */
    qfPtrs->push(ptrs[0]);
    dssize += ptrs[0]->QAcount;

    /* Update batch numbers */
    for (auto fid = 1; fid < nfiles; fid++)
    {
        ptrs[fid]->curr_chunk = ptrs[fid - 1]->curr_chunk + ptrs[fid - 1]->nqchunks;
        qfPtrs->push(ptrs[fid]);
        dssize += ptrs[fid]->QAcount;
    }

    /* Compute the total number of batches in the dataset */
    nBatches = ptrs[nfiles-1]->curr_chunk + ptrs[nfiles-1]->nqchunks;

    return status;
}


/* FUNCTION: DSLIM_SearchManager
 *
 * DESCRIPTION: Manages and performs the Peptide Search
 *
 * INPUT:
 * @slm_index : Pointer to the SLM_Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_SearchManager(Index *index)
{
    status_t status = SLM_SUCCESS;
    int_t batchsize = 0;

    double qtime = 0;

    int_t maxlen = params.max_len;
    int_t minlen = params.min_len;

#ifdef USE_MPI
    thread_t *wthread = new thread_t;
#endif /* USE_MPI */

    /* The mutex for queryfile vector */
    if (status == SLM_SUCCESS)
    {
        status = sem_init(&qfilelock, 0, 1);

        status = DSLIM_InitializeMS2Data();
    }

    /* Initialize the lw double buffer queues with
     * capacity, min and max thresholds */
    if (status == SLM_SUCCESS)
    {
        qPtrs = new lwbuff<Queries>(20, 5, 15); // cap, th1, th2
    }

    /* Initialize the ePtrs */
    if (status == SLM_SUCCESS)
    {
        ePtrs = new expeRT[params.threads];
    }

    /* Create queries buffers and push them to the lwbuff */
    if (status == SLM_SUCCESS)
    {
        /* Create new Queries */
        for (int_t wq = 0; wq < qPtrs->len(); wq++)
        {
            Queries *nPtr = new Queries;

            /* Initialize the query buffer */
            nPtr->init();

            /* Add them to the buffer */
            qPtrs->Add(nPtr);
        }
    }

    if (status == SLM_SUCCESS)
    {
        /* Let's do a queue of 10 MSQuery elements -
         * should be more than enough */

        ioQ = new lwqueue<MSQuery*>(10);

        /* Check for correct allocation */
        if (ioQ == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }

        /* Initialize the ioQlock */
        status = sem_init(&ioQlock, 0, 1);
    }

    if (status == SLM_SUCCESS && params.nodes == 1)
    {
        status = DFile_InitFiles();
    }
#ifdef USE_MPI
    else if (status == SLM_SUCCESS && params.nodes > 1)
    {
        iBuff = new ebuffer[NIBUFFS];

        status = sem_init(&writer, 0, 0);

        if (wthread != NULL && status == SLM_SUCCESS)
        {
            /* Pass the reference to thread block as argument */
            status = pthread_create(wthread, NULL, &DSLIM_FOut_Thread_Entry, (VOID *)NULL);
        }
    }
#endif /* USE_MPI */

    /* Initialize the Comm module */
#ifdef USE_MPI

    /* Only required if nodes > 1 */
    if (params.nodes > 1)
    {
        /* Allocate a new DSLIM Comm handle */
        if (status == SLM_SUCCESS)
        {
            CommHandle = new DSLIM_Comm(nBatches);

            if (CommHandle == NULL)
            {
                status = ERR_BAD_MEM_ALLOC;
            }
        }

        if (status == SLM_SUCCESS)
        {
            CandidatePSMS = new hCell[dssize];

            if (CandidatePSMS == NULL)
            {
                status = ERR_BAD_MEM_ALLOC;
            }
        }
    }

#endif /* USE_MPI */

    /* Create a new Scheduler handle */
    if (status == SLM_SUCCESS)
    {
        SchedHandle = new Scheduler;

        /* Check for correct allocation */
        if (SchedHandle == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
        else
        {
            scheduler_init = true;
        }
    }

    //
    // ----------------------------------------------------------------------------------
    //
    /* The main query loop starts here */
    while (status == SLM_SUCCESS)
    {
#if (USE_TIMEMORY)
        static wall_tuple_t sched_penalty("DAG_penalty", false);
        sched_penalty.start();
#endif
        /* Start computing penalty */
        MARK_START(penal);

        status = DSLIM_WaitFor_IO(batchsize);
        
#if (USE_TIMEMORY)
        sched_penalty.stop();
#endif
        /* Compute the penalty */
        MARK_END(penal);

        /* Check if endsignal */
        if (status == ENDSIGNAL)
        {
            break;
        }

        auto penalty = ELAPSED_SECONDS(penal);

#ifndef DIAGNOSE
        if (params.myid == 0)
        {
            std::cout << "PENALTY: " << penalty << 's'<< std::endl;
        }
#endif /* DIAGNOSE */

        /* Check the status of buffer queues */
        qPtrs->lockr_();
        int_t dec = qPtrs->readyQStatus();
        qPtrs->unlockr_();

        /* Run the Scheduler to manage thread between compute and I/O */
        SchedHandle->runManager(penalty, dec);

#ifndef DIAGNOSE
        if (params.myid == 0)
        {
            std::cout << "Querying: \n" << endl;
        }
#endif /* DIAGNOSE */

        MARK_START(sch_time);

        if (status == SLM_SUCCESS)
        {
            /* Query the chunk */
            status = DSLIM_QuerySpectrum(workPtr, index, (maxlen - minlen + 1));
        }

#ifdef USE_MPI
        /* Transfer my partial results to others */
        if (status == SLM_SUCCESS && params.nodes > 1)
        {
            status = sem_post(&writer);
        }
#endif /* USE_MPI */

        status = qPtrs->lockw_();

        /* Request next I/O chunk */
        qPtrs->Replenish(workPtr);

        status = qPtrs->unlockw_();

        MARK_END(sch_time);

        /* Compute Duration */
        qtime +=  ELAPSED_SECONDS(sch_time);

#ifndef DIAGNOSE
        if (params.myid == 0)
        {
            std::cout << "\nQuery Time: " << qtime << "s" << endl;
            std::cout << "DONE: Querying\tstatus: " << status << endl << endl;
        }
#endif /* DIAGNOSE */
    }

    /* Deinitialize the IO module */
    status = DSLIM_Deinit_IO();

    /* Delete the scheduler object */
    if (SchedHandle != NULL)
    {
        /* Deallocate the scheduler module */
        delete SchedHandle;

        SchedHandle = NULL;
    }

#ifdef USE_MPI
    /* Deinitialize the Communication module */
    if (params.nodes > 1)
    {
        ciBuff ++;
        ebuffer *liBuff = iBuff + (ciBuff % NIBUFFS);

        /* Wait for FOut thread to take out iBuff */
        for (; liBuff->isDone == false; usleep(1000));

        status = sem_post(&writer);

        VOID *ptr = NULL;

#if defined (USE_TIMEMORY)
        static wall_tuple_t comm_penalty("DAG_penalty");
        comm_penalty.start();
#endif // USE_TIMEMORY

        pthread_join(*wthread, &ptr);

#if defined (USE_TIMEMORY)
        comm_penalty.stop();
#endif // USE_TIMEMORY

        delete wthread;

        sem_destroy(&writer);

        if (iBuff != NULL)
        {
            delete[] iBuff;
            iBuff = NULL;
        }

#ifdef DIAGNOSE
        std::cout << "ExitSignal: " << params.myid << endl;
#endif /* DIAGNOSE */

        /* Wait for everyone to synchronize */
        status = MPI_Barrier(MPI_COMM_WORLD);

        /* Carry forward the data to the distributed scoring module */
        status = DSLIM_CarryForward(index, CommHandle, ePtrs, CandidatePSMS, spectrumID);

        /* Delete the instance of CommHandle */
        if (CommHandle != NULL)
        {
            delete CommHandle;
            CommHandle = NULL;
        }
    }
#endif /* USE_MPI */

    if (status == SLM_SUCCESS && params.nodes == 1)
    {
        status = DFile_DeinitFiles();

        delete[] ePtrs;
        ePtrs = NULL;
    }

    /* Delete ptrs */
    if (ptrs)
    {
        delete[] ptrs;
    }

    /* Return the status of execution */
    return status;
}

/* FUNCTION: DSLIM_QuerySpectrum
 *
 * DESCRIPTION: Query the DSLIM for all query peaks
 *              and count the number of hits per chunk
 *
 * INPUT:
 * @QA     : Query Spectra Array
 * @len    : Number of spectra in the array
 * @Matches: Array to fill in the number of hits per chunk
 * @threads: Number of parallel threads to launch
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_QuerySpectrum(Queries *ss, Index *index, uint_t idxchunk)
{
    status_t status = SLM_SUCCESS;
    uint_t maxz = params.maxz;
    uint_t dF = params.dF;
    int_t threads = (int_t) params.threads - (int_t) SchedHandle->getNumActivThds();
    uint_t scale = params.scale;
    double_t maxmass = params.max_mass;
    ebuffer *liBuff = NULL;
    partRes *txArray = NULL;

#if defined (USE_TIMEMORY)
    static search_tuple_t search_inst("theSrch");
    search_inst.start();

    // PAPI is only available on Windows
#   if defined (_UNIX)
    static hw_counters_t search_cntr ("theSrch");
    search_cntr.start();
#   endif // _UNIX
#endif // USE_TIMEMORY

    if (params.nodes > 1)
    {
#if defined (USE_TIMEMORY)
        static wall_tuple_t comm_penalty("DAG_penalty");
        comm_penalty.start();
#endif // USE_TIMEMORY

        ciBuff ++;
        liBuff = iBuff + (ciBuff % NIBUFFS);

        /* Wait for FOut thread to take out iBuff */
        for (; liBuff->isDone == false; usleep(10000));

        txArray = liBuff->packs;
        liBuff->isDone = false;
        liBuff->batchNum = ss->batchNum;

#if defined (USE_TIMEMORY)
        comm_penalty.stop();
#endif // USE_TIMEMORY
    }
    else
    {
        LBE_UNUSED_PARAM(liBuff);
    }

#ifndef _OPENMP
    LBE_UNUSED_PARAM(threads);
#endif /* USE_OMP */

    /* Sanity checks */
    if (Score == NULL || (txArray == NULL && params.nodes > 1))
    {
        status = ERR_INVLD_MEMORY;
    }

    if (status == SLM_SUCCESS)
    {
        /* Should at least be 1 and min 75% */
        int_t minthreads = MAX(1, (params.threads * 3)/4);

        threads = MAX(threads, minthreads);

#ifndef DIAGNOSE
        /* Print how many threads are we using here */
        if (params.myid == 0)
        {
            /* Print the number of query threads */
            std::cout << "\n#QThds: " << threads << endl;
        }
#endif /* DIAGNOSE */

        /* Process all the queries in the chunk.
         * Setting chunk size to 4 to avoid false sharing
         */
#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 4)
#endif /* USE_OMP */
        for (int_t queries = 0; queries < ss->numSpecs; queries++)
        {
            /* Pointer to each query spectrum */
            uint_t *QAPtr = ss->moz + ss->idx[queries];
            float_t pmass = ss->precurse[queries];
            uint_t *iPtr = ss->intensity + ss->idx[queries];
            uint_t qspeclen = ss->idx[queries + 1] - ss->idx[queries];
            uint_t thno = omp_get_thread_num();

            BYC *bycPtr = Score[thno].byc;
            Results *resPtr = &Score[thno].res;
            expeRT  *expPtr = ePtrs + thno;
            ebuffer *inBuff = inBuff + thno;

#ifndef DIAGNOSE
            if (thno == 0 && params.myid == 0)
            {
                std::cout << "\rDONE: " << (queries * 100) /ss->numSpecs << "%";
            }
#endif /* DIAGNOSE */

            for (uint_t ixx = 0; ixx < idxchunk; ixx++)
            {
                uint_t speclen = (index[ixx].pepIndex.peplen - 1) * maxz * iSERIES;
                uint_t halfspeclen = speclen / 2;

                for (uint_t chno = 0; chno < index[ixx].nChunks; chno++)
                {
                    /* Query each chunk in parallel */
                    uint_t *bAPtr = index[ixx].ionIndex[chno].bA;
                    uint_t *iAPtr = index[ixx].ionIndex[chno].iA;

                    int_t minlimit = 0;
                    int_t maxlimit = 0;

                    BOOL val = DSLIM_BinarySearch(index + ixx, ss->precurse[queries], minlimit, maxlimit);

                    /* Spectrum violates limits */
                    if (val == false || (maxlimit < minlimit))
                    {
                        continue;
                    }

                    /* Query all fragments in each spectrum */
                    for (uint_t k = 0; k < qspeclen; k++)
                    {
                        /* Do this to save mem boundedness */
                        auto qion = QAPtr[k];
                        uint_t intn = iPtr[k];

                        /* Check for any zeros
                         * Zero = Trivial query */
                        if (qion > dF && qion < ((maxmass * scale) - 1 - dF))
                        {
                            for (auto bin = qion - dF; bin < qion + 1 + dF; bin++)
                            {
                                /* Locate iAPtr start and end */
                                uint_t start = bAPtr[bin];
                                uint_t end = bAPtr[bin + 1];

                                /* If no ions in the bin */
                                if (end - start < 1)
                                {
                                    continue;
                                }

                                auto ptr = std::lower_bound(iAPtr + start, iAPtr + end, minlimit * speclen);
                                int_t stt = start + std::distance(iAPtr + start, ptr);

                                ptr = std::upper_bound(iAPtr + stt, iAPtr + end, (((maxlimit + 1) * speclen) - 1));
                                int_t ends = stt + std::distance(iAPtr + stt, ptr) - 1;

                                /* Loop through located iAions */
                                for (auto ion = stt; ion <= ends; ion++)
                                {
                                    uint_t raw = iAPtr[ion];

                                    /* Calculate parent peptide ID */
                                    int_t ppid = (raw / speclen);

                                    /* Calculate the residue */
                                    int_t residue = (raw % speclen);

                                    /* Either 0 or 1 */
                                    int_t isY = residue / halfspeclen;
                                    int_t isB = 1 - isY;

                                    while (isY < 0 || isY > 1);

                                    /* Get the map element */
                                    BYC *elmnt = bycPtr + ppid;

                                    /* Update */
                                    elmnt->bc += isB;
                                    elmnt->ibc += intn * isB;

                                    elmnt->yc += isY;
                                    elmnt->iyc += intn * isY;
                                }
                            }
                        }
                    }

                    /* Compute the chunksize to look further into */
                    int_t csize = maxlimit - minlimit + 1;

                    /* Look for candidate PSMs */
                    for (int_t it = minlimit; it <= maxlimit; it++)
                    {
                        ushort_t bcc = bycPtr[it].bc;
                        ushort_t ycc = bycPtr[it].yc;
                        ushort_t shpk = bcc + ycc;

                        /* Filter by the min shared peaks */
                        if (shpk >= params.min_shp)
                        {                            
                            /* Create a heap cell */
                            hCell cell;

                            ull_t pp = UTILS_Factorial(bcc) *
                                    UTILS_Factorial(ycc);

                            /* Fill in the information */
                            cell.hyperscore = 0.001 + pp * bycPtr[it].ibc * bycPtr[it].iyc;

                            cell.hyperscore = log10(cell.hyperscore) - 6;

                            /* hyperscore < 0 means either b- or y- ions were not matched */
                            if (cell.hyperscore > 0)
                            {
                                cell.idxoffset = ixx;
                                cell.psid = it;
                                cell.sharedions = shpk;
                                cell.totalions = speclen;
                                cell.pmass = pmass;

                                /* Insert the cell in the heap dst */
                                resPtr->topK.insert(cell);

                                /* Increase the N */
                                resPtr->cpsms += 1;

                                /* Update the histogram */
                                resPtr->survival[(int_t) (cell.hyperscore * 10 + 0.5)] += 1;
                            }
                        }
                    }

                    /* Clear the scorecard */
                    std::memset(bycPtr + minlimit, 0x0, sizeof(BYC) * csize);
                }
            }

#ifdef USE_MPI
            /* Distributed memory mode - Model partial Gumbel
             * and transmit parameters to rx machine */
            if (params.nodes > 1)
            {
                /* Set the params.min_cpsm in dist mem mode to 1 */
                if (resPtr->cpsms >= 1)
                {
                    /* Extract the top PSM */
                    hCell psm = resPtr->topK.getMax();

                    /* Put it in the list */
                    CandidatePSMS[spectrumID + queries] = psm;

                    resPtr->maxhypscore = (psm.hyperscore * 10 + 0.5);

                    status = expPtr->StoreIResults(resPtr, queries, liBuff);

                    /* Fill in the Tx array cells */
                    txArray[queries].min  = resPtr->minhypscore;
                    txArray[queries].max2 = resPtr->nexthypscore;
                    txArray[queries].max  = resPtr->maxhypscore;
                    txArray[queries].N    = resPtr->cpsms;
                    txArray[queries].qID  = spectrumID + queries;
                }
                else
                {
                    /* No need to memset as there are apt checks in dslim_score.cpp
                    memset(liBuff->ibuff + (queries * psize * sizeof(ushort_t)), 0x0, psize * sizeof(ushort_t));*/

                    /* Extract the top result
                     * and put it in the list */
                    CandidatePSMS[spectrumID + queries] = 0;

                    /* Get the handle to the txArr
                     * Fill it up and move on */
                    txArray[queries] = 0;
                    txArray[queries].qID  = spectrumID + queries;
                }
            }

            /* Shared memory mode - Do complete
             * modeling and print results */
            else
#endif /* USE_MPI */
            {
                /* Check for minimum number of PSMs */
                if (resPtr->cpsms >= params.min_cpsm)
                {
                    /* Extract the top PSM */
                    hCell psm = resPtr->topK.getMax();

                    resPtr->maxhypscore = (psm.hyperscore * 10 + 0.5);

                    /* Compute expect score if there
                     * are any candidate PSMs */
#ifdef TAILFIT
                    status = expPtr->ModelTailFit(resPtr);

                    /* Linear Regression Parameters */
                    double_t w = resPtr->mu;
                    double_t b = resPtr->beta;

                    w /= 1e6;
                    b /= 1e6;

                    /* Estimate the log (s(x)); x = log(hyperscore) */
                    double_t lgs_x = (w * resPtr->maxhypscore) + b;

                    /* Compute the s(x) */
                    double_t e_x = pow(10, lgs_x);

                    /* e(x) = n * s(x) */
                    e_x *= resPtr->cpsms;

#else
                    status = expPtr->ModelSurvivalFunction(resPtr);

                    /* Extract e(x) = n * s(x) = mu * 1e6 */
                    double_t e_x = resPtr->mu;

                    e_x /= 1e6;

#endif /* TAILFIT */

                    /* Do not print any scores just yet */
                    if (e_x < params.expect_max)
                    {
                        /* Printing the scores in OpenMP mode */
                        status = DFile_PrintScore(index, spectrumID + queries, pmass, &psm, e_x, resPtr->cpsms);
                    }
                }
            }

            /* Reset the results */
            resPtr->reset();
        }

        if (params.nodes > 1)
        {
            liBuff->currptr = ss->numSpecs * psize * sizeof(ushort_t);
        }

        /* Update the number of queried spectra */
        spectrumID += ss->numSpecs;
    }

#if defined (USE_TIMEMORY)
    search_inst.stop();

    // PAPI is only available on Windows
#   if defined (_UNIX)
    search_cntr.stop();
#   endif // _UNIX
#endif // USE_TIMEMORY

#ifndef DIAGNOSE
    if (params.myid == 0)
    {
        std::cout << "\nQueried Spectra:\t\t" << workPtr->numSpecs << endl;
    }
#endif /* DIAGNOSE */

    return status;
}

/*
 * FUNCTION: DSLIM_BinarySearch
 *
 * DESCRIPTION: The Binary Search Algorithm
 *
 * INPUT:
 *
 * OUTPUT
 * none
 */
static BOOL DSLIM_BinarySearch(Index *index, float_t precmass, int_t &minlimit, int_t &maxlimit)
{
    /* Get the float_t precursor mass */
    float_t pmass1 = precmass - params.dM;
    float_t pmass2 = precmass + params.dM;
    pepEntry *entries = index->pepEntries;

    BOOL rv = false;

    uint_t min = 0;
    uint_t max = index->lcltotCnt - 1;

    if (params.dM < 0.0)
    {
        minlimit = min;
        maxlimit = max;

        return rv;
    }

    /* Check for base case */
    if (pmass1 < entries[min].Mass)
    {
        minlimit = min;
    }
    else if (pmass1 > entries[max].Mass)
    {
        minlimit = max;
        maxlimit = max;
        return rv;
    }
    else
    {
        /* Find the minlimit here */
        minlimit = DSLIM_BinFindMin(entries, pmass1, min, max);
    }

    min = 0;
    max = index->lcltotCnt - 1;


    /* Check for base case */
    if (pmass2 > entries[max].Mass)
    {
        maxlimit = max;
    }
    else if (pmass2 < entries[min].Mass)
    {
        minlimit = min;
        maxlimit = min;
        return rv;
    }
    else
    {
        /* Find the maxlimit here */
        maxlimit = DSLIM_BinFindMax(entries, pmass2, min, max);
    }

    if (entries[maxlimit].Mass <= pmass2 && entries[minlimit].Mass >= pmass1)
    {
        rv = true;
    }

    return rv;
}


static int_t DSLIM_BinFindMin(pepEntry *entries, float_t pmass1, int_t min, int_t max)
{
    int_t half = (min + max)/2;

    if (max - min < 20)
    {
        int_t current = min;

        while (entries[current].Mass < pmass1)
        {
            current++;
        }

        return current;
    }

    if (pmass1 > entries[half].Mass)
    {
        min = half;
        return DSLIM_BinFindMin(entries, pmass1, min, max);
    }
    else if (pmass1 < entries[half].Mass)
    {
        max = half;
        return DSLIM_BinFindMin(entries, pmass1, min, max);
    }

    if (pmass1 == entries[half].Mass)
    {
        while (pmass1 == entries[half].Mass)
        {
            half--;
        }

        half++;
    }

    return half;


}

static int_t DSLIM_BinFindMax(pepEntry *entries, float_t pmass2, int_t min, int_t max)
{
    int_t half = (min + max)/2;

    if (max - min < 20)
    {
        int_t current = max;

        while (entries[current].Mass > pmass2)
        {
            current--;
        }

        return current;
    }

    if (pmass2 > entries[half].Mass)
    {
        min = half;
        return DSLIM_BinFindMax(entries, pmass2, min, max);
    }
    else if (pmass2 < entries[half].Mass)
    {
        max = half;
        return DSLIM_BinFindMax(entries, pmass2, min, max);
    }

    if (pmass2 == entries[half].Mass)
    {
        half++;

        while (pmass2 == entries[half].Mass)
        {
            half++;
        }

        half--;
    }

    return half;

}
/*
 * FUNCTION: DSLIM_IO_Threads_Entry
 *
 * DESCRIPTION: Entry function for all
 *              I/O threads
 *
 * INPUT:
 * @argv: Pointer to void arguments
 *
 * OUTPUT:
 * @NULL: Nothing
 */
VOID *DSLIM_IO_Threads_Entry(VOID *argv)
{
    status_t status = SLM_SUCCESS;
    Queries *ioPtr = NULL;

    BOOL eSignal = false;

    // TODO: verify thread local performance
#if defined (USE_TIMEMORY)
    thread_local prep_tuple_t prep_inst("preprocess");
    prep_inst.start();
#endif // USE_TIMEMORY

    /* local object is fine since it will be copied
     * to the queue object at the time of preemption */
    MSQuery *Query = NULL;

    int_t rem_spec = 0;

    while (!scheduler_init) { usleep(1); }

    /* Initialize and process Query Spectra */
    for (;status == SLM_SUCCESS;)
    {
        /* Check if the Query object is not initialized */
        if (Query == NULL || Query->isDeInit())
        {
            /* Try getting the Query object from queue if present */
            status = sem_wait(&ioQlock);

            if (!ioQ->isEmpty())
            {
                Query = ioQ->front();
                status = ioQ->pop();
            }

            status = sem_post(&ioQlock);
        }

        /* If the queue is empty */
        if (Query == NULL || Query->isDeInit())
        {
            /* Otherwise, initialize the object from a file */

            /* lock the query file */
            sem_wait(&qfilelock);

            /* Check if anymore queryfiles */
            if (!qfPtrs->isEmpty())
            {
                Query = qfPtrs->front();
                qfPtrs->pop();
                rem_spec = Query->getQAcount(); // Init to 1 for first loop to run
            }
            else
            {
                /* Raise the exit signal */
                eSignal = true;
            }

            /* Unlock the query file */
            sem_post(&qfilelock);
        }

        /* If no more files, then break the inf loop */
        if (eSignal == true)
        {
            break;
        }

        /*********************************************
         * At this point, we have the data ready     *
         *********************************************/

        /* All set */
        if (status == SLM_SUCCESS)
        {
            MARK_START(prep);

            /* Wait for a I/O request */
            status = qPtrs->lockw_();

            /* Empty wait queue or Scheduler preemption signal raised  */
            if (SchedHandle->checkPreempt() || qPtrs->isEmptyWaitQ())
            {
                status = qPtrs->unlockw_();

                if (status == SLM_SUCCESS)
                {
                    status = sem_wait(&ioQlock);
                }

                if (status == SLM_SUCCESS)
                {
                    status = ioQ->push(Query);
                }

                if (status == SLM_SUCCESS)
                {
                    status = sem_post(&ioQlock);
                }

                /* Break from loop */
                break;
            }

            /* Otherwise, get the I/O ptr from the wait queue */
            ioPtr = qPtrs->getIOPtr();

            status = qPtrs->unlockw_();

            /* Reset the ioPtr */
            ioPtr->reset();

            /* Extract a chunk and return the chunksize */
            status = Query->ExtractQueryChunk(QCHUNK, ioPtr, rem_spec);

            ioPtr->batchNum = Query->curr_chunk;
            Query->curr_chunk++;

            /* Lock the ready queue */
            qPtrs->lockr_();

#ifdef USE_MPI
            if (params.nodes > 1)
            {
                /* Add an entry of the added buffer to the CommHandle */
                status = CommHandle->AddBatch(ioPtr->batchNum, ioPtr->numSpecs, Query->getQfileIndex());
            }
#endif /* USE_MPI */

            /*************************************
             * Add available data to ready queue *
             *************************************/
            qPtrs->IODone(ioPtr);

            /* Unlock the ready queue */
            qPtrs->unlockr_();

            MARK_END(prep);
            auto es = ELAPSED_SECONDS(prep);

#ifndef DIAGNOSE
            if (params.myid == 0)
            {
                std::cout << "\nExtracted Spectra :\t\t" << ioPtr->numSpecs << endl;
                std::cout << "Elapsed Time: " << es << "s" << endl << endl;
            }
#endif /* DIAGNOSE */

            /* If no more remaining spectra, then deinit */
            if (rem_spec < 1)
            {
                status = Query->DeinitQueryFile();

                if (Query != NULL)
                {
                    delete Query;
                    Query = NULL;
                }
            }
        }
        /*else
        {
            std::cout << "ALERT: IOTHD @" << params.myid << endl;
        }*/
    }

    /* Check if we ran out of files */
    if (eSignal == true)
    {
        if (Query != NULL)
        {
            delete Query;
            Query = NULL;
        }

        /* Free the main IO thread */
        SchedHandle->ioComplete();
    }

#if defined (USE_TIMEMORY)
    prep_inst.stop();
#endif // USE_TIMEMORY

    /* Request pre-emption */
    SchedHandle->takeControl(argv);

    return NULL;
}

#ifdef USE_MPI
VOID *DSLIM_FOut_Thread_Entry(VOID *argv)
{
    //status_t status = SLM_SUCCESS;
    int_t batchSize = 0;
    int_t clbuff = -1;

#if defined (USE_TIMEMORY)
    thread_local comm_tuple_t comm_inst("result_comm");
#endif // USE_TIMEMORY
    for (;;)
    {
        sem_wait(&writer);

        clbuff += 1;
        ebuffer *lbuff = iBuff + (clbuff % NIBUFFS);

        if (lbuff->isDone == true)
        {
            //cout << "FOut_Thread_Exiting @: " << params.myid <<endl;
            break;
        }

        ofstream *fh = new ofstream;
        string_t fn = params.datapath + "/" +
                    std::to_string(lbuff->batchNum) +
                    "_" + std::to_string(params.myid) + ".dat";

        batchSize = lbuff->currptr / (psize * sizeof(ushort_t));

        fh->open(fn, ios::out | ios::binary);

        fh->write((char_t *)lbuff->packs, batchSize * sizeof(partRes));
        fh->write(lbuff->ibuff, lbuff->currptr * sizeof (char_t));

        fh->close();

        lbuff->isDone = true;
    }

#if defined (USE_TIMEMORY)
    comm_inst.stop();
#endif

    return argv;
}
#endif /* USE_MPI */

static inline status_t DSLIM_Deinit_IO()
{
    status_t status = SLM_SUCCESS;

    Queries *ptr = NULL;

    while (!qPtrs->isEmptyReadyQ())
    {
        ptr = qPtrs->getWorkPtr();

        if (ptr != NULL)
        {
            delete ptr;
            ptr = NULL;
        }
    }

    while (!qPtrs->isEmptyWaitQ())
    {
        ptr = qPtrs->getIOPtr();

        if (ptr != NULL)
        {
            delete ptr;
            ptr = NULL;
        }
    }

    /* Delete the qPtrs buffer handle */
    delete qPtrs;

    qPtrs = NULL;

    /* Deallocate the I/O queues */
    delete ioQ;
    ioQ = NULL;

    /* Destroy the ioQ lock semaphore */
    status = sem_destroy(&ioQlock);

    return status;
}
