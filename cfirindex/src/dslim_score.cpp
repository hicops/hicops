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

extern gParams params;
extern VOID *DSLIM_Score_Thread_Entry(VOID *);

DSLIM_Score::DSLIM_Score()
{
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
    rxSizes = new INT[threads-1];
    txSizes = new INT[threads-1];

    /* Set all to zero */
    std::memset(rxSizes, 0x0, (sizeof(rxSizes[0]) * (threads-1)));
    std::memset(rxSizes, 0x0, (sizeof(txSizes[0]) * (threads-1)));

    /* key-values */
    keys     = NULL;
    TxValues = NULL;
    RxValues = NULL;

    /* Communication threads */
    comm_thd = new THREAD;

    /* Create the comm. thread */
    pthread_create(comm_thd, NULL, &DSLIM_Score_Thread_Entry, this);

    return;
}

DSLIM_Score::DSLIM_Score(BData *bd)
{
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
    rxSizes = new INT[threads-1];
    txSizes = new INT[threads-1];

    /* Set all to zero */
    std::memset(rxSizes, 0x0, (sizeof(rxSizes[0]) * (threads-1)));
    std::memset(rxSizes, 0x0, (sizeof(txSizes[0]) * (threads-1)));

    myRXsize = 0;

    for (INT kk = 0; kk < nBatches; kk++)
    {
        myRXsize += sizeArray[kk];
    }

    /* FIXME: What if nothing RX'ed? */
    if (myRXsize >= 0)
    {
        keys = new INT[myRXsize];
        TxValues = new fResult[myRXsize + nSpectra];
        RxValues = TxValues + myRXsize;
    }
    else
    {
        cout << "ABORT: myRXsize = 0 @" << params.myid << endl;
    }

    /* Communication threads */
    comm_thd = new THREAD;

    /* Create the comm thread */
    pthread_create(comm_thd, NULL, &DSLIM_Score_Thread_Entry, this);

    return;
}

DSLIM_Score::~DSLIM_Score()
{
    /* Wait for score thread to complete */
    VOID *ptr;

    if (comm_thd != NULL)
    {
        pthread_join(*comm_thd, &ptr);
    }

    nSpectra = 0;
    nBatches = 0;
    myRXsize = 0;

    if (sizeArray != NULL)
    {
        delete[] sizeArray;
        sizeArray = NULL;
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

    if (scPtr != NULL)
    {
        delete[] scPtr;
        scPtr = NULL;
    }

    if (heapArray != NULL)
    {
        delete[] heapArray;
        heapArray = NULL;
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

    if (keys != NULL)
    {
        delete[] keys;
        keys = NULL;
    }

    if (TxValues != NULL)
    {
        delete[] TxValues;
        TxValues = NULL;
    }

    if (comm_thd != NULL)
    {
        delete comm_thd;
        comm_thd = NULL;
    }

    return;
}

VOID DSLIM_Score::Initialize(BData *bd)
{
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

    myRXsize = 0;

    for (INT kk = 0; kk < nBatches; kk++)
    {
        myRXsize += sizeArray[kk];
    }

    /* FIXME: What if nothing RX'ed? */
    if (myRXsize >= 0)
    {
        /* key-values */
        keys = new INT[myRXsize];
        TxValues = new fResult[myRXsize + nSpectra];
        RxValues = TxValues + myRXsize;
    }
    else
    {
        cout << "ABORT: myRXsize < 0 @" << params.myid << endl;
    }

    return;
}

STATUS DSLIM_Score::ComputeDistScores()
{
    STATUS status = SLM_SUCCESS;

    /* TODO: Implement here */
    /* Take care of this during merging stage */

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
    INT nodes_1   = threads - 1;
    INT cumulate = 0;

    /* MPI pointers */
    MPI_Request *txRqsts  = NULL;
    MPI_Request *txRqsts2 = NULL;
    INT         *txStats  = NULL;
    INT         *txStats2 = NULL;

    /* 2X per all other machines */
    txRqsts = new MPI_Request [2 * (nodes_1)];
    txStats = new INT[2 * (nodes_1)];

    /* Set the pointers to the appropriate offset */
    txRqsts2 = txRqsts + (nodes_1);
    txStats2 = txStats + (nodes_1);

    /* Fill the txStats with zeros - not available */
    std::memset(txStats, 0x0, (sizeof(txStats[0]) * (2 * (nodes_1))));

    /* Check if everything is in order */
    if (txRqsts != NULL && txStats != NULL)
    {
        status = TXSizes(txRqsts);
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }

    if (status == SLM_SUCCESS)
    {
        /* Wait 500ms between each loop */
        for (cumulate = 0; cumulate < (nodes_1); usleep(500000))
        {
            for (INT ll = 0; ll < nodes_1; ll ++)
            {
                /* Only check if not already received */
                if (!txStats[ll])
                {
                    status = MPI_Test(txRqsts + ll, &txStats[ll], MPI_STATUS_IGNORE);

                    /* If the TXsize has been sent,
                     * start the TX request.
                     */
                    if (txStats[ll])
                    {
                        status = TXResults(txRqsts2 + ll, ll);
                        cumulate++;
                    }
                }
            }
        }
    }

    if (status == SLM_SUCCESS)
    {
        /* Wait 500ms between each loop */
        for (cumulate = 0; cumulate < (nodes_1); usleep(500000))
        {
            for (INT ll = 0; ll < nodes_1; ll ++)
            {
                /* Only check if not already received */
                if (!txStats2[ll])
                {
                    status = MPI_Test(txRqsts2 + ll, &txStats2[ll], MPI_STATUS_IGNORE);

                    /* Check if the results
                     * have been sent
                     */
                    if (txStats2[ll])
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

STATUS DSLIM_Score::TXSizes(MPI_Request *txRqsts)
{
    STATUS status = SLM_SUCCESS;

    UINT ll = 0;

    for (UINT kk = 0; kk < params.nodes; kk++)
    {
        /* If myself then no RX */
        if (ll == params.myid)
        {
            continue;
        }

        /* Receive an integer from all other machines */
        status = MPI_Isend(txSizes, 1, MPI_INT, kk, 0x5125, MPI_COMM_WORLD, txRqsts + ll);

        /* Increment the ll variable */
        ll++;
    }

    return status;
}

STATUS DSLIM_Score::RXSizes(MPI_Request *rxRqsts)
{
    STATUS status = SLM_SUCCESS;

    UINT ll = 0;

    for (UINT kk = 0; kk < params.nodes; kk++)
    {
        /* If myself then no RX */
        if (ll == params.myid)
        {
            continue;
        }

        /* Receive an integer from all other machines */
        status = MPI_Irecv(rxSizes, 1, MPI_INT, kk, 0x5125, MPI_COMM_WORLD, rxRqsts + ll);

        /* Increment the ll */
        ll++;
    }

    return status;
}

STATUS DSLIM_Score::TXResults(MPI_Request *txRqst, INT rawID)
{
    INT mchID = rawID;

    /* Compute the mchID as:
     * mchID = rawID +1 if >= myid
     */
    if (rawID >= (INT)params.myid)
    {
        mchID = rawID + 1;
    }

    /* Transmit data to machine:
     * @mch: mchID, @size: rxSizes[rawID] * INT in fResult
     */
    return MPI_Isend((INT *)TxValues, txSizes[rawID] * (sizeof(TxValues[0])/sizeof(INT)),
                     MPI_INT, mchID, 0xDA2A, MPI_COMM_WORLD, txRqst);
}

STATUS DSLIM_Score::RXResults(MPI_Request *rxRqst, INT rawID)
{
   UINT mchID = rawID;

    /* Compute the mchID as:
     * mchID = rawID +1 if >= myid
     */
    if (rawID >= (INT)params.myid)
    {
        mchID = rawID + 1;
    }

    /* Receive data from machine:
     * @mch: mchID, @size: rxSizes[rawID] * INT in fResult
     */
    return MPI_Irecv((INT *)RxValues, rxSizes[rawID] * (sizeof(TxValues[0])/sizeof(INT)),
                     MPI_INT, mchID, 0xDA2A, MPI_COMM_WORLD, rxRqst);
}

STATUS DSLIM_Score::DisplayResults()
{
    STATUS status = SLM_SUCCESS;

    /* Initialize the file handles */
    status = DFile_InitFiles();

    /* TODO: Implement stuff here */

    /* Close the files and deallocate objects */
    if (status == SLM_SUCCESS)
    {
        status = DFile_DeinitFiles();
    }

    return status;
}

VOID DSLIM_Score::Deinitialize()
{
    /* Wait for score thread to complete */
    VOID *ptr;

    if (comm_thd != NULL)
    {
        pthread_join(*comm_thd, &ptr);
    }

    nSpectra = 0;
    nBatches = 0;
    myRXsize = 0;

    if (sizeArray != NULL)
    {
        delete[] sizeArray;
        sizeArray = NULL;
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

    if (scPtr != NULL)
    {
        delete[] scPtr;
        scPtr = NULL;
    }

    if (heapArray != NULL)
    {
        delete[] heapArray;
        heapArray = NULL;
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

    if (keys != NULL)
    {
        delete[] keys;
        keys = NULL;
    }

    if (TxValues != NULL)
    {
        delete[] TxValues;
        TxValues = NULL;
    }

    if (comm_thd != NULL)
    {
        delete comm_thd;
        comm_thd = NULL;
    }

    return;
}
