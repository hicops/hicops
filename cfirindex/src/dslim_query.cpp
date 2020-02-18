/*
 * This file is part of PCDSFrame software
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
#include "dslim_comm.h"
#include "lwqueue.h"
#include "lwbuff.h"
#include "scheduler.h"

using namespace std;

extern gParams   params;
extern BYICount  *Score;
extern vector<STRING> queryfiles;

/* Global variables */
FLOAT *hyperscores = NULL;
UCHAR *sCArr = NULL;
BOOL ExitSignal = false;

DSLIM_Comm *CommHandle = NULL;
Scheduler *SchedHandle = NULL;
lwqueue<commRqst> *requestQ = NULL;

/* Lock for query file vector and thread manager */
LOCK qfilelock;

INT qfid = 0;

/****************************************************************/

/* Expt spectra data buffer */
lwbuff<Queries> *qPtrs = NULL;
Queries *workPtr = NULL;

#ifdef DISTMEM
VOID *DSLIM_Comm_Thread_Entry(VOID *argv);
#endif

/* A queue containing I/O thread state when preempted */
lwqueue<MSQuery *> *ioQ = NULL;
LOCK ioQlock;
/****************************************************************/


#ifdef BENCHMARK
static DOUBLE duration = 0;
extern DOUBLE compute;
extern DOUBLE fileio;
extern DOUBLE memory;
#endif /* BENCHMARK */

#ifdef DISTMEM
VOID *DSLIM_Comm_Thread_Entry(VOID *argv);
#endif /* DISTMEM */

VOID *DSLIM_IO_Threads_Entry(VOID *argv);

static VOID DSLIM_BinarySearch(Index *, FLOAT, INT&, INT&);
static INT  DSLIM_BinFindMin(pepEntry *entries, FLOAT pmass1, INT min, INT max);
static INT  DSLIM_BinFindMax(pepEntry *entries, FLOAT pmass2, INT min, INT max);
static inline STATUS DSLIM_WaitFor_IO(INT &);
static inline STATUS DSLIM_Deinit_IO();

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
static inline STATUS DSLIM_WaitFor_IO(INT &batchsize)
{
    STATUS status;

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

STATUS DSLIM_Process_RxData()
{
    //partRes *rxData = CommHandle->getRxArray()

    return SLM_SUCCESS;
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
STATUS DSLIM_SearchManager(Index *index)
{
    STATUS status = SLM_SUCCESS;

    INT batchnum = 0;
    INT batchsize = 0;
    INT buffernum = 0;
    partRes *buffers = NULL;

    auto start = chrono::system_clock::now();
    auto end   = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    chrono::duration<double> qtime = end - start;
    INT maxlen = params.max_len;
    INT minlen = params.min_len;

    /* The mutex for queryfile vector */
    if (status == SLM_SUCCESS)
    {
        status = sem_init(&qfilelock, 0, 1);
    }

    /* Initialize the lw double buffer queues with
     * capacity, min and max thresholds */
    if (status == SLM_SUCCESS)
    {
        qPtrs = new lwbuff<Queries>(20, 5, 15); // cap, th1, th2
    }

    /* Create queries buffers and push them to the lwbuff */
    if (status == SLM_SUCCESS)
    {
        /* Create new Queries */
        for (INT wq = 0; wq < qPtrs->len(); wq++)
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

    /* Initialize the Comm module */
#ifdef DISTMEM

    /* Only required if nodes > 1 */
    if (params.nodes > 1)
    {
        /* Construct the Request Queue */
        requestQ = new lwqueue<commRqst>(3);

        /* Allocate a new DSLIM Comm handle */
        if (status == SLM_SUCCESS)
        {
            CommHandle = new DSLIM_Comm;

            if (CommHandle == NULL)
            {
                status = ERR_BAD_MEM_ALLOC;
            }
        }
    }
#endif /* DISTMEM */

    /* Create a new Scheduler handle */
    if (status == SLM_SUCCESS)
    {
        SchedHandle = new Scheduler;

        /* Check for correct allocation */
        if (SchedHandle == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
    }

    /* Initialize the file handles */
    if (status == SLM_SUCCESS)
    {
        status = DFile_InitFiles();
    }

    /**************************************************************************/
    /* The main query loop starts here */
    while (status == SLM_SUCCESS)
    {
        /* Start computing penalty */
        auto spen = chrono::system_clock::now();

        status = DSLIM_WaitFor_IO(batchsize);

        /* Check if endsignal */
        if (status == ENDSIGNAL)
        {
            break;
        }

        /* Compute the penalty */
        chrono::duration<double> penalty = chrono::system_clock::now() - spen;

        cout << "PENALTY: " << penalty.count() << endl;

        /* Check the status of buffer queues */
        qPtrs->lockr_();
        INT dec = qPtrs->readyQStatus();
        qPtrs->unlockr_();

        /* Run the Scheduler to manage thread between compute and I/O */
        SchedHandle->runManager(penalty.count(), dec);

#ifdef DISTMEM
        if (params.nodes > 1)
        {
            /* Get the Tx Buffer and start Rx */
            buffers = CommHandle->getTxBuffer(batchnum, batchsize, buffernum);

        }
#endif /* DISTMEM */

        if (params.myid == 0)
        {
            cout << "Querying: \n" << endl;
        }

        start = chrono::system_clock::now();

        if (status == SLM_SUCCESS)
        {
            /* Query the chunk */
            status = DSLIM_QuerySpectrum(workPtr, index, (maxlen - minlen + 1), buffers);
        }

#ifdef DISTMEM
        /* Transfer my partial results to others */
        if (status == SLM_SUCCESS && params.nodes > 1)
        {
            commRqst txR;
            txR.bsize = batchsize;
            txR.btag = batchnum;
            txR.buff = buffernum;

            requestQ->push(txR);

            /* Signal the thread about the Tx/Rx array */
            status = CommHandle->SignalWakeup();

            batchnum++;
        }
#endif /* DISTMEM */

        status = qPtrs->lockw_();

        /* Request next I/O chunk */
        qPtrs->Replenish(workPtr);

        status = qPtrs->unlockw_();

        end = chrono::system_clock::now();

        /* Compute Duration */
        qtime += end - start;

        /* FIXME: Uncomment this thing */
        //if (params.myid == 0)
        {
            /* Compute Duration */
            cout << "Query Time: " << qtime.count() << "s" << endl;
            cout << "Queried with status:\t\t" << status << endl << endl;
        }

        end = chrono::system_clock::now();
    }

    /* Deinitialize the IO module */
    status = DSLIM_Deinit_IO();

#ifdef DISTMEM
    /* Deinitialize the Communication module */
    if (params.nodes > 1)
    {
        /* TODO: Wait for CommHandle to complete its work */

        /* Post the exit signal to Comm Handle */
        CommHandle->SignalExit();

        /* Delete the instance of CommHandle */
        delete CommHandle;

        CommHandle = NULL;

        delete requestQ;
        requestQ = NULL;
    }
#endif /* DISTMEM */

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
STATUS DSLIM_QuerySpectrum(Queries *ss, Index *index, UINT idxchunk, partRes *txArray)
{
    STATUS status = SLM_SUCCESS;
    UINT maxz = params.maxz;
    UINT dF = params.dF;
    UINT threads = params.threads - SchedHandle->getNumActivThds();
    UINT scale = params.scale;
    DOUBLE maxmass = params.max_mass;

#ifdef BENCHMARK
    DOUBLE tcons[threads];
    std::memset(tcons, 0x0, threads * sizeof(DOUBLE));
#endif /* BENCHMARK */

#ifndef _OPENMP
    LBE_UNUSED_PARAM(threads);
#endif /* _OPENMP */

#ifdef BENCHMARK
    duration = omp_get_wtime();
#endif /* BENCHMARK */

    if (Score == NULL || txArray == NULL)
    {
        status = ERR_INVLD_MEMORY;
    }

    if (status == SLM_SUCCESS)
    {
        /* Print the number of query threads */
        cout << "\n#QThreads: " << threads << " @node: " << params.myid << endl;

        /* Process all the queries in the chunk */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif /* _OPENMP */
        for (INT queries = 0; queries < ss->numSpecs; queries++)
        {
#ifdef BENCHMARK
            DOUBLE stime = omp_get_wtime();
#endif
            /* Pointer to each query spectrum */
            UINT *QAPtr = ss->moz + ss->idx[queries];
            FLOAT pmass = ss->precurse[queries];
            UINT *iPtr = ss->intensity + ss->idx[queries];
            UINT qspeclen = ss->idx[queries + 1] - ss->idx[queries];
            UINT thno = omp_get_thread_num();

            BYC *bycPtr = Score[thno].byc;
            iBYC *ibycPtr = Score[thno].ibyc;
            Results *resPtr = &Score[thno].res;

            if (thno == 0 && params.myid == 0)
            {
                std::cout << "\rDONE: " << (queries * 100) /ss->numSpecs << "%";
            }

            for (UINT ixx = 0; ixx < idxchunk; ixx++)
            {
                UINT speclen = (index[ixx].pepIndex.peplen - 1) * maxz * iSERIES;

                for (UINT chno = 0; chno < index[ixx].nChunks; chno++)
                {
                    /* Query each chunk in parallel */
                    UINT *bAPtr = index[ixx].ionIndex[chno].bA;
                    UINT *iAPtr = index[ixx].ionIndex[chno].iA;

                    INT minlimit = 0;
                    INT maxlimit = 0;

                    DSLIM_BinarySearch(index + ixx, ss->precurse[queries], minlimit, maxlimit);

                    /* Spectrum violates limits */
                    if ((maxlimit - minlimit) < 1)
                    {
                        continue;
                    }

                    /* Query all fragments in each spectrum */
                    for (UINT k = 0; k < qspeclen; k++)
                    {
                        /* Check for any zeros
                         * Zero = Trivial query */
                        if (QAPtr[k] > dF && QAPtr[k] < ((maxmass * scale) - 1 - dF))
                        {
                            /* Locate iAPtr start and end */
                            UINT start = bAPtr[QAPtr[k] - dF];
                            UINT end = bAPtr[QAPtr[k] + 1 + dF];

                            /* Loop through located iAions */
                            for (UINT ion = start; ion < end; ion++)
                            {
                                UINT raw = iAPtr[ion];

                                /* Calculate parent peptide ID */
                                INT ppid = (raw / speclen);

                                if (ppid >= minlimit && ppid <= maxlimit)
                                {
                                    /* Update corresponding scorecard entries */
                                    if ((raw % speclen) < speclen / 2)
                                    {
                                        bycPtr[ppid].bc += 1;
                                        ibycPtr[ppid].ibc += iPtr[k];
                                    }
                                    else
                                    {
                                        bycPtr[ppid].yc += 1;
                                        ibycPtr[ppid].iyc += iPtr[k];
                                    }
                                }

                            }
                        }
                    }

                    INT csize = maxlimit - minlimit;

                    for (INT it = minlimit; it < maxlimit; it++)
                    {
                        if (bycPtr[it].bc + bycPtr[it].yc >= params.min_shp)
                        {
                            /* Create a heap cell */
                            hCell cell;

                            /* Fill in the information */
                            cell.hyperscore = log10(0.001 +
                                                    DFile_Factorial(ULONGLONG(bycPtr[it].bc)) *
                                                    DFile_Factorial(ULONGLONG(bycPtr[it].yc)) *
                                                    ibycPtr[it].ibc *
                                                    ibycPtr[it].iyc);

                            /* hyperscore < 0 means either b- or y- ions were not matched */
                            if (cell.hyperscore > 0)
                            {
                                cell.idxoffset = ixx;
                                cell.psid = it;
                                cell.sharedions = bycPtr[it].bc + bycPtr[it].yc;
                                cell.totalions = speclen;

                                /* Insert the cell in the heap dst */
                                resPtr->topK.insert(cell);
                                /* Increase the N */
                                resPtr->cpsms += 1;
                                /* Update the histogram */
                                resPtr->survival[(INT) (cell.hyperscore * 10 + 0.5)] += 1;
                            }
                        }
                    }

                    /* Clear the scorecard */
                    std::memset(bycPtr + minlimit, 0x0, sizeof(BYC) * csize);
                    std::memset(ibycPtr + minlimit, 0x0, sizeof(iBYC) * csize);
                }
            }

            /* Compute expect score if there are any candidate PSMs */
            if (resPtr->cpsms > params.min_cpsm)
            {
                status = DSLIM_ModelSurvivalFunction(resPtr);

                hCell psm = resPtr->topK.getMax();

                /* Fill in the Tx array cells */
//                txArray[queries].b = resPtr->bias;
//                txArray[queries].m = resPtr->weight;
//                txArray[queries].min = resPtr->minhypscore;
//                txArray[queries].max = resPtr->nexthypscore;

                /* Estimate the log (s(x)); x = log(hyperscore) */
                DOUBLE lgs_x = resPtr->weight * (psm.hyperscore * 10 + 0.5) + resPtr->bias;

                /* Compute the s(x) */
                DOUBLE s_x = pow(10, lgs_x);

                /* e(x) = n * s(x) */
                DOUBLE e_x = resPtr->cpsms * s_x;

                if (e_x < params.expect_max)
                {
#ifndef ANALYSIS

                    /* Printing the scores in OpenMP mode */
                    status = DFile_PrintScore(index, queries, pmass, &psm, e_x, resPtr->cpsms);
#else
                    status = DFile_PrintPartials(queries, resPtr);
#endif /* ANALYSIS */
                }
            }
            else
            {
                /* Get the handle to the txArr - Fill it up and move on */
            }

            /* Reset the results */
            resPtr->reset();

#ifdef BENCHMARK
            tcons[thno] += omp_get_wtime() - stime;
#endif
        }
    }

    std::cout << '\n';

#ifdef BENCHMARK
    compute += omp_get_wtime() - duration;

    for (unsigned int thd = 0; thd < params.threads; thd++)
    {
        std::cout << "Thread #: " << thd << "\t" << tcons[thd] << std::endl;
    }
#endif

    cout << "Queried Spectra:\t\t" << workPtr->numSpecs << endl;

    return status;
}

/*
 * FUNCTION: DSLIM_DeallocateSC
 *
 * DESCRIPTION:
 *
 * INPUT:
 * none
 *
 * OUTPUT:
 * @status: status of execution
 */
STATUS DSLIM_DeallocateSC()
{
    /* Free the Scorecard memory */
    if (Score != NULL)
    {
        for (UINT thd = 0; thd < params.threads; thd++)
        {
            delete[] Score[thd].byc;
            delete[] Score[thd].ibyc;
            delete[] Score[thd].res.survival;
            delete[] Score[thd].res.xaxis;

            Score[thd].byc = NULL;
            Score[thd].ibyc = NULL;
            Score[thd].res.survival = NULL;
            Score[thd].res.xaxis = NULL;
        }

        delete[] Score;
        Score = NULL;

    }

    return SLM_SUCCESS;
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
static VOID DSLIM_BinarySearch(Index *index, FLOAT precmass, INT &minlimit, INT &maxlimit)
{
    /* Get the FLOAT precursor mass */
    FLOAT pmass1 = precmass - params.dM;
    FLOAT pmass2 = precmass + params.dM;
    pepEntry *entries = index->pepEntries;

    UINT min = 0;
    UINT max = index->lcltotCnt - 1;

    if (params.dM < 0.0)
    {
        minlimit = min;
        maxlimit = max;

        return;
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
        return;
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
        return;
    }
    else
    {
        /* Find the maxlimit here */
        maxlimit = DSLIM_BinFindMax(entries, pmass2, min, max);
    }

}


static INT DSLIM_BinFindMin(pepEntry *entries, FLOAT pmass1, INT min, INT max)
{
    INT half = (min + max)/2;

    if (max - min < 500)
    {
        INT current = min;

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


static INT DSLIM_BinFindMax(pepEntry *entries, FLOAT pmass2, INT min, INT max)
{
    INT half = (min + max)/2;

    if (max - min < 500)
    {
        INT current = max;

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

STATUS DSLIM_ModelSurvivalFunction(Results *resPtr)
{
    STATUS status = SLM_SUCCESS;

    /* Total size and tailends */
    UINT N = resPtr->cpsms;

    /* Choosing the kneePt and endPt to be
     * at 70.7% and 99% respectively */
    UINT kneePt = N - (UINT)((float)N * 0.707);

#ifndef ANALYSIS
    UINT endPt = (UINT)((float)N * 0.995);
#endif /* ANALYSIS */

    /* Copy the slope and bias into local variables */
    DOUBLE slope = resPtr->weight;
    DOUBLE bias = resPtr->bias;

    /* Histogram pointer */
    DOUBLE *histogram = resPtr->survival;
    const INT histosize = 2 + (MAX_HYPERSCORE * 10);

    /* The extracted tail */
    DOUBLE *tail = NULL;
    DOUBLE *axis = NULL;
    INT tailsize = 0;

    /* Loop indexing variable */
    INT ii = 0;

    /* Construct the model */

    /* Initialize the max hyperscore */
    for (ii = histosize - 1; ii > 0; ii--)
    {
        if (histogram[ii] > 0)
        {
            resPtr->maxhypscore = ii;
            ii--;
            break;
        }
    }

    /* Initialize nexthypscore */
    for (; ii > 0; ii--)
    {
        if (histogram[ii] > 0)
        {
            resPtr->nexthypscore = ii;
            break;
        }
    }

    /* Initialize minhypscore */
    for (ii = 0; ii < resPtr->nexthypscore; ii++)
    {
        if (histogram[ii] > 0)
        {
            resPtr->minhypscore = ii;
            break;
        }
    }

    /* Set the minhypscore at kneepoint: ~ 70.7% */
    UINT cumulative = 0;

    for (ii = resPtr->minhypscore; ii < resPtr->nexthypscore; ii++)
    {
        cumulative += histogram[ii];

        /* Mark the lower end of the tail */
        if (cumulative >= kneePt)
        {
            if (ii > resPtr->minhypscore)
            {
                resPtr->minhypscore = ii;
            }

            break;
        }
    }

#ifndef ANALYSIS
    /* Set the nexthypscore at endPt: ~99.5% */
    cumulative = N;

    for (ii = resPtr->maxhypscore; ii >= resPtr->minhypscore; ii--)
    {
        cumulative -= histogram[ii];

        /* Mark the upper end of tail */
        if (cumulative <= endPt)
        {
            if (ii < resPtr->nexthypscore)
            {
                resPtr->nexthypscore = ii;
            }

            break;
        }
    }
#endif /* ANALYSIS */

    /* If both ends at the same point,
     * set the upper end to maxhypscore
     */
    if (resPtr->nexthypscore <= resPtr->minhypscore)
    {
        resPtr->nexthypscore = resPtr->maxhypscore - 1;
    }

    /* Construct s(x) = 1 - CDF in [minhyp, nexthyp] */
    UINT count = histogram[resPtr->nexthypscore];

    for (ii = resPtr->nexthypscore - 1; ii >= resPtr->minhypscore; ii--)
    {
        UINT tmpcount = histogram[ii];

        cumulative -= tmpcount;
        histogram[ii] = (count + histogram[ii + 1]);
        count = tmpcount;
    }

    /* Construct log_10(s(x)) */
    for (ii = resPtr->minhypscore; ii <= resPtr->nexthypscore; ii++)
    {
        histogram[ii] = log10(histogram[ii] / N);
    }

    /* Set the tailPtr, score axis and tail size */
    tail     = histogram + resPtr->minhypscore;
    axis     = resPtr->xaxis + resPtr->minhypscore;
    tailsize = resPtr->nexthypscore - resPtr->minhypscore + 1;

    /*
     * Perform linear regression (least sq. error)
     * on tail curve and find slope (m) and bias (b)
     */
    (VOID) UTILS_LinearRegression(tailsize, axis, tail, slope, bias);

    /* Assign back from local variables */
    resPtr->weight = slope;
    resPtr->bias   =  bias;

    /* Return the status */
    return status;
}

#ifdef DISTMEM

/*
 * FUNCTION: Comm_Thread_Entry
 *
 * DESCRIPTION: Entry function for the MPI
 *              communication thread
 *
 * INPUT:
 * @argv: Pointer to void arguments
 *
 * OUTPUT:
 * @NULL: Nothing
 */
/*
 * FUNCTION: Comm_Thread_Entry
 *
 * DESCRIPTION: Entry function for the MPI
 *              communication thread
 *
 * INPUT:
 * @argv: Pointer to void arguments
 *
 * OUTPUT:
 * @NULL: Nothing
 */
VOID *DSLIM_Comm_Thread_Entry(VOID *argv)
{
    INT btag = 0;
    INT bsize = 0;
    INT buff = 0;

    STATUS status = SLM_SUCCESS;

    /* Avoid race conditions by waiting for
     * CommHandle pointer to initialize */
    while ((VOID *)argv != (VOID *)CommHandle);

    /* The forever loop */
    for (;status == SLM_SUCCESS;)
    {
        /* Wait for an wakeup event */
        status = CommHandle->Wait4Event();

        /* Check if the TxData has the exit signal */
        if (CommHandle->checkExitSignal() == true)
        {
            break;
        }

        /* Check if Wakeup Signal came from an IO Thread */
        if (CommHandle->checkWakeup())
        {
            /* Check if mismatch exists and if Rx is not ready */
            if (CommHandle->getRxReadyPermission())
            {
                status = CommHandle->Rx();
            }
            else
            {
                cout << "COMM THD: Something went wrong on node: "
                     << params.myid << endl;
            }
        }
        /* Wakeup came from Scheduler */
        else
        {
            /* Copy a txRqst from the front of queue */
            commRqst  txRqst = requestQ->front();

            /* Pop from the queue */
            requestQ->pop();

            btag = txRqst.btag;
            bsize = txRqst.bsize;
            buff = txRqst.buff;

            /* Check if its a Tx */
            if (btag % (params.nodes) != params.myid)
            {
                if (buff != -1)
                {
                    /* Tx the results need the batchsize */
                    status = CommHandle->Tx(btag, bsize, buff);
                }
                else
                {
                    cout << "\nDSLIM_Comm Thread: Buffer number -1 detected at process: " << params.myid << endl;
                    break;
                }

            }

            /* Check for Rx the results in every wakeup from scheduler */
            status = CommHandle->CheckRx();

            /* Check if next Rx is available and previous ended */
            if (CommHandle->getRxReadyPermission())
            {
                /* Initialize the next Rx */
                status = CommHandle->Rx();
            }
        }

        /* Check if everything is successful */
        if (status != SLM_SUCCESS)
        {
            /* Should never reach here */
            cout << "Status from Comm. Thread: " << status << " on node: "
                 << params.myid << endl;

            cout << "Aborting..." << endl;

            break;
        }

    }

    /* Return NULL */
    return NULL;
}
#endif /* DISTMEM */

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
    STATUS status = SLM_SUCCESS;
    auto start = chrono::system_clock::now();
    auto end   = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    chrono::duration<double> qtime = end - start;
    Queries *ioPtr = NULL;

    BOOL eSignal = false;

    /* local object is fine since it will be copied
     * to the queue object at the time of preemption */
    MSQuery *Query = NULL;

    INT rem_spec = 0;

    /* Initialize and process Query Spectra */
    for (;;)
    {
#ifdef BENCHMARK
        duration = omp_get_wtime();
#endif /* BENCHMARK */

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
            start = chrono::system_clock::now();

            /* lock the query file */
            sem_wait(&qfilelock);

            UINT qfid_lcl = qfid;

            /* Check if anymore queryfiles */
            if (qfid_lcl < queryfiles.size())
            {
                qfid++;
            }
            else
            {
                /* unlock the query file */
                sem_post(&qfilelock);

                /* Raise the exit signal */
                eSignal = true;
            }

            if (eSignal == false)
            {
                /* Unlock the query file */
                sem_post(&qfilelock);

                /* We need a new Query object here */
                Query = new MSQuery;

                /* Initialize Query MS/MS file */
                status = Query->InitQueryFile(&queryfiles[qfid_lcl]);

                if (params.myid == 0)
                {
                    cout << "Query File: " << queryfiles[qfid_lcl] << endl;
                    cout << "Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;
                }

                rem_spec = Query->getQAcount(); // Init to 1 for first loop to run
                end = chrono::system_clock::now();

                /* Compute Duration */
                elapsed_seconds = end - start;

#ifdef BENCHMARK
                fileio += omp_get_wtime() - duration;
#endif /* BENCHMARK */
            }
        }

        /* If no more files, then break the inf loop */
        if (eSignal == true)
        {
            break;
        }

        /*********************************************
         * At this point, we have the data ready     *
         *********************************************/

        /* All set - Run the DSLIM Query Algorithm */
        if (status == SLM_SUCCESS)
        {
            start = chrono::system_clock::now();
#ifdef BENCHMARK
            duration = omp_get_wtime();
#endif /* BENCHMARK */

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

            if (params.nodes > 1)
            {
                /* Add an entry of the added buffer to the CommHandle */
                status = CommHandle->AddBufferEntry(ioPtr->numSpecs);
            }

            /* Lock the ready queue */
            qPtrs->lockr_();

            /*************************************
             * Add available data to ready queue *
             *************************************/
            qPtrs->IODone(ioPtr);

            /* Unlock the ready queue */
            qPtrs->unlockr_();

#ifdef BENCHMARK
            fileio += omp_get_wtime() - duration;
#endif /* BENCHMARK */
            end = chrono::system_clock::now();

            /* Compute Duration */
            elapsed_seconds = end - start;

            if (params.myid == 0)
            {
                cout << "Extracted Spectra :\t\t" << ioPtr->numSpecs << endl;
                cout << "Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;
            }

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

    /* Request pre-emption */
    SchedHandle->takeControl(argv);

    return NULL;
}

static inline STATUS DSLIM_Deinit_IO()
{
    STATUS status = SLM_SUCCESS;

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

    /* Deallocate the scheduler module */
    delete SchedHandle;

    SchedHandle = NULL;

    return status;
}
