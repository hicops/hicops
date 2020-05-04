/*
 * Copyright (C) 2020  Muhammad Haseeb, and Fahad Saeed
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

#include "dslim_score.h"
#include "dslim_fileout.h"

#ifdef DISTMEM

using namespace std;

MPI_Datatype resultF;
STRING fn;

extern gParams params;
extern VOID *DSLIM_Score_Thread_Entry(VOID *);

DSLIM_Score::DSLIM_Score()
{
    INT nodes = params.nodes;

    threads = params.threads;

    /* Dataset size */
    nSpectra = 0;
    nBatches = 0;
    myRXsize = 0;

    /* These pointers will be borrowed */
    sizeArray = NULL;
    fileArray = NULL;
    indxArray = NULL;
    ePtr      = NULL;
    heapArray = NULL;
    index = NULL;
    resPtr = NULL;

    /* Data size that I expect to
     * receive from other processes */
    rxSizes = new INT[nodes];
    txSizes = new INT[nodes];

    /* Set all to zero */
    std::memset(rxSizes, 0x0, (sizeof(INT) * (nodes)));
    std::memset(txSizes, 0x0, (sizeof(INT) * (nodes)));

    /* key-values */
    keys     = NULL;
    TxValues = NULL;
    RxValues = NULL;

    InitDataTypes();

    /* Communication threads */
    comm_thd = new THREAD;

    /* Create the comm. thread */
    pthread_create(comm_thd, NULL, &DSLIM_Score_Thread_Entry, this);

    return;
}

DSLIM_Score::DSLIM_Score(BData *bd)
{
    INT nodes = params.nodes;

    threads = params.threads;

    /* Dataset size */
    nSpectra = bd->cPSMsize;
    nBatches = bd->nBatches;

    /* These pointers will be borrowed */
    sizeArray = bd->sizeArray;
    fileArray = bd->fileArray;
    indxArray = bd->indxArray;

    ePtr = bd->ePtr;
    heapArray = bd->heapArray;
    index = bd->index;
    resPtr = bd->resPtr;

    /* Data size that I expect to
     * receive from other processes */
    rxSizes = new INT[nodes];
    txSizes = new INT[nodes];

    if (rxSizes != NULL && txSizes != NULL)
    {
        /* Set all to zero */
        std::memset(rxSizes, 0x0, (sizeof(INT) * nodes));
        std::memset(txSizes, 0x0, (sizeof(INT) * nodes));
    }
    else
    {
        throw "FATAL: rxSizes or txSizes == NULL\n";
    }

    myRXsize = 0;

    /* Compute myRXsize */
    for (INT kk = 0; kk < nBatches; kk++)
    {
        myRXsize += sizeArray[kk];
    }

    /* Allocate for Rx */
    RxValues = new fResult[nSpectra];

    /* If no RX'ed data */
    if (myRXsize > 0)
    {
        TxValues = new fResult[myRXsize];
        keys = new INT[myRXsize];
    }
    else
    {
        TxValues = NULL;
        keys = NULL;
    }

    InitDataTypes();

    /* Communication threads */
    comm_thd = new THREAD;

    /* Create the comm thread */
    pthread_create(comm_thd, NULL, &DSLIM_Score_Thread_Entry, this);

    return;
}

DSLIM_Score::~DSLIM_Score()
{
    /* Wait for score thread to complete */
    Wait4RX();

    if (txSizes != NULL)
    {

        delete[] txSizes;
        txSizes = NULL;
    }

    if (rxSizes != NULL)
    {
        delete[] rxSizes;
        rxSizes = NULL;
    }

    if (sizeArray != NULL)
    {
        delete[] sizeArray;
        sizeArray = NULL;
    }

    if (fileArray != NULL)
    {
        delete[] fileArray;
        fileArray = NULL;
    }

    if (ePtr != NULL)
    {
        delete[] ePtr;
        ePtr = NULL;
    }

    if (heapArray != NULL)
    {
        delete[] heapArray;
        heapArray = NULL;
    }

    if (resPtr != NULL)
    {
        delete[] resPtr;
        resPtr = NULL;
    }

    if (keys != NULL && myRXsize != 0)
    {
        delete[] keys;
        keys = NULL;
    }

    if (TxValues != NULL && myRXsize != 0)
    {
        delete[] TxValues;
        TxValues = NULL;
    }

    if (RxValues != NULL)
    {
        delete[] RxValues;
        RxValues = NULL;
    }

    FreeDataTypes();

    nSpectra = 0;
    nBatches = 0;
    myRXsize = 0;

    return;
}

STATUS DSLIM_Score::CombineResults()
{
    STATUS status = SLM_SUCCESS;

    /* Each node sent its sample */
    const INT nSamples = params.nodes;
    auto startPos = 0;
    auto startSpec= 0;
    ifstream *fhs = NULL;
    ebuffer  *iBuffs = NULL;

    if (this->myRXsize > 0)
    {
        fhs = new ifstream[nSamples];
        iBuffs = new ebuffer[nSamples];
    }

    auto vbatch = params.myid;

    for (auto batchNum = 0; batchNum < this->nBatches; batchNum++, vbatch += nSamples)
    {
        auto bSize = sizeArray[batchNum];

        for (INT saa = 0; saa < nSamples; saa++)
        {
            fn = params.datapath + "/" + std::to_string(vbatch) + "_"
                    + std::to_string(saa) + ".dat";

            fhs[saa].open(fn, ios::in | ios::binary);

            if (fhs[saa].is_open())
            {
                fhs[saa].read(iBuffs[saa].ibuff, bSize * 128 * sizeof(USHORT));

                if (fhs[saa].fail())
                {
                    status = ERR_FILE_NOT_FOUND;
                    cout << "FATAL: File Read Failed" << endl;
                    exit(status);
                }
            }

            fn.clear();
        }

#ifdef _OPENMP
#pragma omp parallel for schedule (dynamic, 4) num_threads(params.threads)
#endif /* _OPENMP */
        for (INT spec = 0; spec < bSize; spec++)
        {
            INT thno = omp_get_thread_num();

            /* Results pointer to use */
            expeRT *expPtr = this->ePtr + thno;

            INT cpsms = 0;

            /* Compute the spectrum index in rxArray */
            INT specIDX = startPos + spec;

            /* Compute the stride in rxArray */
            INT stride = bSize;

            /* Record locators */
            INT key = params.nodes;
            INT maxhypscore = -1;

            /* For all samples, update the histogram */
            for (INT sno = 0; sno < nSamples; sno++)
            {
                /* Pointer to Result sample */
                partRes *sResult = resPtr + (specIDX + (stride * sno));

                /* Sanity check */
                if (specIDX + (stride * sno) >= (this->myRXsize * (nSamples)))
                {
                    cout << "FATAL: Segfault caught @: " << params.myid << endl;
                    cout << "spec: " << spec << endl;
                    cout << "currCount: " << startPos << endl;
                    cout << "myRXsize: " << this->myRXsize << endl;
                    cout << "specIDX: " << specIDX << endl;
                    cout << "batchNum: " << batchNum << endl;
                    cout << "sizeArray: " << sizeArray[batchNum] << endl;

                    exit(-11);
                }

                if (*sResult == 0)
                {
                    continue;
                }

                /* Update the number of samples */
                cpsms += sResult->N;

                /* Only take out data if present */
                if (sResult->N >= 1)
                {
                    /* Reconstruct the partial histogram */
                    expPtr->Reconstruct(&iBuffs[sResult->sno], spec, sResult);

                    /* Record the maxhypscore and its key */
                    if (sResult->max > 0 && sResult->max > maxhypscore)
                    {
                        maxhypscore = sResult->max;
                        key = sResult->sno;
                    }
                }
            }

            /* Combine the fResult */
            fResult *psm = &TxValues[startSpec + spec];

            /* Need further processing only if enough results */
            if (key < (INT) params.nodes && cpsms >= (INT) params.min_cpsm)
            {
                DOUBLE e_x = params.expect_max;

#ifdef TAILFIT
                /* Model the survival function */
                expPtr->ModelTailFit(e_x, maxhypscore);
#else
                expPtr->ModelSurvivalFunction(e_x, maxhypscore);
#endif /* TAILFIT */

                /* If the scores are good enough */
                if (e_x < params.expect_max)
                {
                    psm->eValue = e_x * 1e6;
                    psm->specID = indxArray[batchNum] + spec;
                    psm->npsms = cpsms;

                    /* Must be an atomic update */
#ifdef _OPENMP
#pragma omp atomic update
                    txSizes[key] += 1;
#else
                    txSizes[key] += 1;
#endif /* _OPENMP */

                    /* Update the key */
                    keys[startSpec + spec] = key;
                }
                else
                {
                    psm->eValue = params.expect_max * 1e6;
                    psm->specID = -1;
                    psm->npsms = 0;
                    keys[startSpec + spec] = params.nodes;
                }
            }
            else
            {
                expPtr->ResetPartialVectors();
                psm->eValue = params.expect_max * 1e6;
                psm->specID = -1;
                psm->npsms = 0;
                keys[startSpec + spec] = params.nodes;
            }
        }

        /* Update the counters */
        startSpec += sizeArray[batchNum];
        startPos += sizeArray[batchNum] * (nSamples);

        /* Remove the files when no longer needed */
        for (INT saa = 0; saa < nSamples; saa++)
        {
            fn = params.workspace + "/" + std::to_string(vbatch) + "_" + std::to_string(saa)
                    + ".dat";

            if (fhs[saa].is_open())
            {
                fhs[saa].close();
                std::remove((const CHAR *) fn.c_str());
            }

            fn.clear();
        }
    }

    if (this->myRXsize > 0)
    {
        delete[] fhs;
        fhs = NULL;

        delete[] iBuffs;
        iBuffs = NULL;
    }

    /* Check if we have RX data */
    if (myRXsize > 0)
    {
        /* Sort the TxValues by keys (mchID) */
        KeyVal_Parallel<INT, fResult>(keys, TxValues, myRXsize, params.threads);
    }
    else
    {
        /* Set all sizes to zero */
        for (UINT ky = 0; ky < params.nodes - 1; ky++)
        {
            txSizes[ky] = 0;
        }
    }


    /* return the status of execution */
    return status;
}

STATUS DSLIM_Score::ScatterScores()
{
    STATUS status = SLM_SUCCESS;

    /* Get the number of nodes - 1 */
    INT nodes   = params.nodes;
    INT cumulate = 0;

    /* MPI pointers */
    MPI_Request *txRqsts  = new MPI_Request [nodes];
    INT         *txStats  = new INT[nodes];

    /* Fill the txStats with zeros - not available */
    for (INT kk = 0; kk < nodes; kk++)
    {
        txStats[kk] = 1;
    }

    /* Check if everything is in order */
    if (txRqsts != NULL && txStats != NULL)
    {
        status = TXSizes(txRqsts, txStats);
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }

    if (status == SLM_SUCCESS)
    {
        /* Wait 500ms between each loop */
         for (INT cumulate = 0; cumulate < nodes; usleep(500000))
         {
             cumulate = 0;

             for (INT ll = 0; ll < nodes; ll ++)
             {
                 /* Only check if not already received */
                 if (txStats[ll] == 0 && ll != (INT)params.myid)
                 {
 #ifdef DIAGNOSE2
                     cout << " MPI_Test@ " << params.myid << "Stat: "
                          << txStats2[ll] << " ll: " << ll << endl;
 #endif /* DIAGNOSE2 */
                     status = MPI_Test(txRqsts + ll, &txStats[ll], MPI_STATUS_IGNORE);

                     /* Check if the results
                      * have been sent
                      */
                     if (txStats[ll])
                     {
                         cumulate++;
                     }
                 }
                 else
                 {
                     cumulate++;
                 }
             }
         }

    }

    /* Check if even I need to send anything? */
    if (status == SLM_SUCCESS && myRXsize != 0)
    {
        /* Fill the txStats with zeros - not available */
        for (INT kk = 0; kk < nodes; kk++)
        {
            txStats[kk] = 1;
        }

        status = TXResults(txRqsts, txStats);

        if (status == SLM_SUCCESS)
        {
            /* Wait 500ms between each loop */
            for (; cumulate < nodes; usleep(500000))
            {
                /* Reset cumulate to 0 */
                cumulate = 0;

                for (INT ll = 0; ll < nodes; ll++)
                {
                    /* Only check if not already received */
                    if (txStats[ll] == 0 && ll != (INT)params.myid)
                    {
#ifdef DIAGNOSE2
                        cout << " MPI_Test@ " << params.myid << "Stat: "
                             << txStats2[ll] << " ll: " << ll << endl;
#endif /* DIAGNOSE2 */
                        status = MPI_Test(txRqsts + ll, &txStats[ll], MPI_STATUS_IGNORE);

                        /* Check if the results
                         * have been sent
                         */
                        if (txStats[ll])
                        {
                            cumulate++;
                        }
                    }
                    else
                    {
                        cumulate++;
                    }
                }
            }
        }
    }

    /* Deallocate the rxRqsts */
    if (txRqsts != NULL)
    {
        delete [] txRqsts;
        txRqsts = NULL;
    }

    /* Deallocate the rxStats */
    if (txStats != NULL)
    {
        delete [] txStats;
        txStats = NULL;
    }

    return status;
}

STATUS DSLIM_Score::TXSizes(MPI_Request *txRqsts, INT *txStats)
{
    STATUS status = SLM_SUCCESS;

    for (UINT kk = 0; kk < params.nodes; kk++, usleep(5000))
    {
        /* If myself then no RX */
        if (kk == params.myid)
        {
            txStats[kk] = 1;
        }
        else
        {

#ifdef DIAGNOSE2
            cout << "TXSIZE: " << params.myid << " -> "
                 << kk << " Size: "<< txSizes[ll] << endl;
#endif /* DIAGNOSE */

            if (txSizes == NULL || txRqsts == NULL)
            {
                cout << "FATAL: txSizes = NULL @: " << params.myid << endl;
                exit(-1);
            }

            /* Send an integer to all other machines */
            status = MPI_Isend(txSizes + kk, 1, MPI_INT, kk, 0x0, MPI_COMM_WORLD, txRqsts + kk);

            txStats[kk] = 0;
        }
    }

    return status;
}

STATUS DSLIM_Score::RXSizes(MPI_Request *rxRqsts, INT *rxStats)
{
    STATUS status = SLM_SUCCESS;

    for (UINT kk = 0; kk < params.nodes; kk++)
    {
        /* If myself then no RX */
        if (kk == params.myid)
        {
            rxStats[kk] = 1;
        }
        else
        {

#ifdef DIAGNOSE2
            cout << "RXSIZE: " << params.myid << " -> "
                 << kk << " Size: "<< txSizes[ll] << endl;
#endif /* DIAGNOSE */

            if (rxSizes == NULL || rxRqsts == NULL)
            {
                cout << "FATAL: rxSizes = NULL @: " << params.myid << endl;
                exit(-1);
            }

            /* Send an integer to all other machines */
            status = MPI_Irecv(rxSizes + kk, 1, MPI_INT, kk, 0x0, MPI_COMM_WORLD, rxRqsts + kk);

            rxStats[kk] = 0;
        }
    }

    return status;
}

STATUS DSLIM_Score::TXResults(MPI_Request *txRqsts, INT* txStats)
{
    STATUS status = SLM_SUCCESS;

    INT offset = 0;

    for (UINT kk = 0; kk < params.nodes; kk++, usleep(5000))
    {
        /* If myself then no RX */
        if (kk == params.myid || txSizes[kk] == 0)
        {
            txStats[kk] = 1;
        }
        else
        {
#ifdef DIAGNOSE2
        cout << "RXSIZE: " << params.myid << " -> " << kk << " Size: "<< txSizes[ll] << endl;
#endif /* DIAGNOSE */

            if (TxValues == NULL || txRqsts == NULL || offset + txSizes[kk] > myRXsize)
            {
                cout << "FATAL: TxValues Failed @: " << params.myid << endl;
                exit(-1);
            }

            /* Send results to all other machines */
            status = MPI_Isend(TxValues + offset, txSizes[kk], resultF, kk, 0x1, MPI_COMM_WORLD, txRqsts + kk);

            txStats[kk] = 0;
        }

        offset += txSizes[kk];
    }

    return status;
}

STATUS DSLIM_Score::RXResults(MPI_Request *rxRqsts, INT *rxStats)
{
    STATUS status = SLM_SUCCESS;

    INT offset = 0;

    for (UINT kk = 0; kk < params.nodes; kk++)
    {
        /* If myself then no RX */
        if (kk == params.myid || rxSizes[kk] == 0)
        {
            rxStats[kk] = 1;
        }
        else
        {

#ifdef DIAGNOSE2
            cout << "RXSIZE: " << params.myid << " -> " << kk << " Size: "<< txSizes[ll] << endl;
#endif /* DIAGNOSE */

            if (RxValues == NULL || rxRqsts == NULL || offset + rxSizes[kk] > nSpectra)
            {
                cout << "FATAL: RxValues failed @: " << params.myid << endl;
                exit(-1);
            }

            /* Send an integer to all other machines */
            status = MPI_Irecv(RxValues + offset, rxSizes[kk], resultF, kk, 0x1, MPI_COMM_WORLD, rxRqsts + kk);

            offset += rxSizes[kk];

            rxStats[kk] = 0;
        }
    }

    return status;
}

STATUS DSLIM_Score::Wait4RX()
{
    /* Wait for score thread to complete */
    VOID *ptr;

    STATUS status = SLM_SUCCESS;

    if (comm_thd != NULL)
    {
        status = pthread_join(*comm_thd, &ptr);

        delete comm_thd;
        comm_thd = NULL;
    }

    return status;
}

STATUS DSLIM_Score::DisplayResults()
{
    STATUS status = SLM_SUCCESS;

    /* Initialize the file handles */
    status = DFile_InitFiles();

    /* Display TxArray Data */
    auto offset = 0;
    auto mysize = txSizes[params.myid];

    if (mysize > 0)
    {
        for (auto beg = 0; beg < (INT)params.myid; beg++)
        {
            offset += txSizes[beg];
        }
    }

    fResult *myPtr = TxValues + offset;

    /* Display the data */
#ifdef _OPENMP
#pragma omp parallel for num_threads(params.threads)
#endif /* _OPENMP */
    for (auto ik = myPtr; ik < myPtr + mysize; ik++)
    {
        hCell *psm = heapArray + ik->specID;

        DFile_PrintScore(this->index, ik->specID, psm->pmass, psm, ((DOUBLE)(ik->eValue))/1e6, ik->npsms);
    }

    /* Now display the RX data */
    myPtr = RxValues;
    mysize = 0;

    for (auto pt = rxSizes; pt < rxSizes + params.nodes; pt++)
    {
        mysize += *pt;
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(params.threads)
#endif /* _OPENMP */
    for (auto ik = myPtr; ik < myPtr + mysize; ik++)
    {
        hCell *psm = heapArray + ik->specID;

        DFile_PrintScore(this->index, ik->specID, psm->pmass, psm, ((DOUBLE)(ik->eValue))/1e6, ik->npsms);
    }

    /* Close the files and deallocate objects */
    if (status == SLM_SUCCESS)
    {
        status = DFile_DeinitFiles();
    }

    return status;
}

STATUS DSLIM_Score::InitDataTypes()
{
    MPI_Type_contiguous((INT)(sizeof(fResult) / sizeof(INT)),
                        MPI_INT,
                        &resultF);

    MPI_Type_commit(&resultF);

    return SLM_SUCCESS;
}

STATUS DSLIM_Score::FreeDataTypes()
{
    MPI_Type_free(&resultF);

    return SLM_SUCCESS;
}

#endif /* DISTMEM */
