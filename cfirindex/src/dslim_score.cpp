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

using namespace std;

MPI_Datatype resultF;

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
    scPtr = NULL;
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
    scPtr = bd->scPtr;
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

#ifdef DIAGNOSE2
    cout << "myRXsize = " << myRXsize << " @: " << params.myid << endl;
    cout << "nSpectra = " << nSpectra << " @: " << params.myid << endl;
    cout << "nSpectra + myRXsize = " << myRXsize << " @: " << params.myid << endl;

    fResult *ptr = TxValues + myRXsize;

    cout << (void *) TxValues << " TxV@: " << params.myid << endl;
    cout << (void *) RxValues << " RxV@: " << params.myid << endl;
    cout << (void *) ptr << " TxEnd@: " << params.myid << endl;
#endif /* DIAGNOSE2 */

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
#ifdef DIAGNOSE2
        cout << (void *) txSizes << " txS@: " << params.myid << endl;
#endif /* DIAGNOSE2 */

        delete[] txSizes;
        txSizes = NULL;
    }

    if (rxSizes != NULL)
    {
#ifdef DIAGNOSE2
        cout << (void *) rxSizes << " rxS@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] rxSizes;
        rxSizes = NULL;
    }

    if (sizeArray != NULL)
    {
#ifdef DIAGNOSE2
        cout << (void *) sizeArray << " sA@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] sizeArray;
        sizeArray = NULL;
    }

    if (fileArray != NULL)
    {
#ifdef DIAGNOSE2
        cout << (void *) fileArray << " fA@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] fileArray;
        fileArray = NULL;
    }

    if (scPtr != NULL)
    {
#ifdef DIAGNOSE2
        cout << (void *) scPtr << "scA@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] scPtr;
        scPtr = NULL;
    }

    if (heapArray != NULL)
    {
#ifdef DIAGNOSE2
        cout << (void *) heapArray << " hA@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] heapArray;
        heapArray = NULL;
    }

    if (resPtr != NULL)
    {
#ifdef DIAGNOSE2
        cout << (void *) resPtr << " resPtr@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] resPtr;
        resPtr = NULL;
    }

    if (keys != NULL && myRXsize != 0)
    {
#ifdef DIAGNOSE2
        cout << (void *) keys << " key@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] keys;
        keys = NULL;
    }

    if (TxValues != NULL && myRXsize != 0)
    {
#ifdef DIAGNOSE2
        cout << (void *) TxValues << " key@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] TxValues;
        TxValues = NULL;
    }

    if (RxValues != NULL)
    {
#ifdef DIAGNOSE2
        cout << (void *) RxValues << " key@: " << params.myid << endl;
#endif /* DIAGNOSE2 */
        delete[] RxValues;
        RxValues = NULL;
    }

    FreeDataTypes();

    nSpectra = 0;
    nBatches = 0;
    myRXsize = 0;

    return;
}

STATUS DSLIM_Score::ComputeDistScores()
{
    STATUS status = SLM_SUCCESS;

    /* TODO: Implement here */
    /* Take care of this during merging stage */

    /* FIXME Only stub implementation */
    for (UINT jj = 0; jj < params.nodes - 1; jj++)
    {
        txSizes[jj] = (this->myRXsize)/(params.nodes);
    }

#if 0
    /* Estimate the log (s(x)); x = log(hyperscore) */
    DOUBLE lgs_x = resPtr->weight * (psm.hyperscore * 10 + 0.5) + resPtr->bias;

    /* Compute the s(x) */
    DOUBLE s_x = pow(10, lgs_x);

    /* e(x) = n * s(x) */
    DOUBLE e_x = resPtr->cpsms * s_x;

    /* Do not print any scores just yet */
    if (e_x < params.expect_max)
    {
#ifndef ANALYSIS
        /* Printing the scores in OpenMP mode */
        //status = DFile_PrintScore(index, queries, pmass, &psm, e_x, resPtr->cpsms);
#else
        status = DFile_PrintPartials(queries, resPtr);
#endif /* ANALYSIS */
    }
#endif /* 0 */

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
                     cout << " MPI_Test@ " << params.myid << "Stat: " << txStats2[ll] << " ll: " << ll << endl;
 #endif /* DIAGNOSE2 */
 
                     //cout << (void *)((MPI_Request*)(txRqsts + ll)) << " txRqsts@:" << params.myid << endl;
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
                        cout << " MPI_Test@ " << params.myid << "Stat: " << txStats2[ll] << " ll: " << ll << endl;
#endif /* DIAGNOSE2 */
                        //cout << (void *)((MPI_Request*)(txRqsts + ll)) << " txRqsts@:" << params.myid << endl;
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
            cout << "TXSIZE: " << params.myid << " -> " << kk << " Size: "<< txSizes[ll] << endl;
#endif /* DIAGNOSE */

            if (txSizes == NULL || txRqsts == NULL)
            {
                cout << "FATAL: txSizes = NULL @: " << params.myid << endl;
                exit(-1);
            }

            //cout << (void *)((INT*)(txSizes + kk)) << " txSizes@:" << params.myid << endl;

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
            cout << "RXSIZE: " << params.myid << " -> " << kk << " Size: "<< txSizes[ll] << endl;
#endif /* DIAGNOSE */

            if (rxSizes == NULL || rxRqsts == NULL)
            {
                cout << "FATAL: rxSizes = NULL @: " << params.myid << endl;
                exit(-1);
            }

            //cout << (void *)((INT*)(rxSizes + kk)) << " rxSizes@:" << params.myid << endl;

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

            //cout << (void *)((fResult*)(TxValues + offset)) << " TxValues@:" << params.myid << endl;

            /* Send an integer to all other machines */
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

            //cout << (void *)((fResult*)(RxValues + offset)) << " RxValues@:" << params.myid << endl;

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
    //status = DFile_InitFiles();

    /* TODO: Implement stuff here */

    /* Close the files and deallocate objects */
    if (status == SLM_SUCCESS)
    {
        //status = DFile_DeinitFiles();
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
