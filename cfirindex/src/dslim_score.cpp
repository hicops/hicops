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
    INT nodes_1 = params.nodes - 1;

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
    rxSizes = new INT[nodes_1];
    txSizes = new INT[nodes_1];

    /* Set all to zero */
    std::memset(rxSizes, 0x0, (sizeof(rxSizes[0]) * (nodes_1)));
    std::memset(rxSizes, 0x0, (sizeof(txSizes[0]) * (nodes_1)));

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
    INT nodes_1 = params.nodes - 1;

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
    rxSizes = new INT[nodes_1];
    txSizes = new INT[nodes_1];

    /* Set all to zero */
    std::memset(rxSizes, 0x0, (sizeof(rxSizes[0]) * (nodes_1)));
    std::memset(rxSizes, 0x0, (sizeof(txSizes[0]) * (nodes_1)));

    myRXsize = 0;

    for (INT kk = 0; kk < nBatches; kk++)
    {
        myRXsize += sizeArray[kk];
    }

    /* Allocate for Tx and Rx */
    TxValues = new fResult[myRXsize + nSpectra];
    RxValues = TxValues + myRXsize;

    /* If no RX'ed data */
    if (myRXsize > 0)
    {
        keys = new INT[myRXsize];
    }
    else
    {
        keys = NULL;
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
    Wait4RX();

    nSpectra = 0;
    nBatches = 0;
    myRXsize = 0;

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
        RxValues = NULL;
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

    /* Allocate for Tx and Rx */
    TxValues = new fResult[myRXsize + nSpectra];
    RxValues = TxValues + myRXsize;

    /* If no RX'ed data */
    if (myRXsize > 0)
    {
        keys = new INT[myRXsize];
    }
    else
    {
        keys = NULL;
    }

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
        txSizes[jj] = (this->myRXsize)/(params.nodes - 1);
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
    INT nodes_1   = params.nodes - 1;
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
    for (INT kk = 0; kk < 2 * nodes_1; kk++)
    {
        txStats[kk] = 0;
    }

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
        for (; cumulate < (nodes_1); usleep(500000))
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
                        /* Check if need to send */
                        if (TXResults(txRqsts2 + ll, ll) == ERR_INVLD_SIZE)
                        {
                            txStats2[ll] = 1;
                        }

                        cumulate++;
                    }
                }
            }
        }
    }

    if (status == SLM_SUCCESS)
    {
        /* Wait 500ms between each loop */
        for (; cumulate < (nodes_1); usleep(500000))
        {
            /* Reset cumulate to 0 */
            cumulate = 0;

            for (INT ll = 0; ll < nodes_1; ll ++)
            {
                /* Only check if not already received */
                if (!txStats2[ll])
                {
#ifdef DIAGNOSE2
                    cout << " MPI_Test@ " << params.myid << "Stat: " << txStats2[ll] << " ll: " << ll << endl;
#endif /* DIAGNOSE2 */

                    status = MPI_Test(txRqsts2 + ll, &txStats2[ll], MPI_STATUS_IGNORE);

                    /* Check if the results
                     * have been sent
                     */
                    if (txStats2[ll])
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
        if (kk == params.myid)
        {
            continue;
        }

#ifdef DIAGNOSE2
        cout << "TXSIZE: " << params.myid << " -> " << kk << " Size: "<< txSizes[ll] << endl;
#endif /* DIAGNOSE */

        /* Send an integer to all other machines */
        status = MPI_Isend((INT *)txSizes + ll, 1, MPI_INT, kk, 0x5125, MPI_COMM_WORLD, txRqsts + ll);

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
        if (kk == params.myid)
        {
            continue;
        }

#ifdef DIAGNOSE2
        cout << "RXSIZE: " << params.myid << " <- " << kk << endl;
#endif /* DIAGNOSE */

        /* Receive an integer from all other machines */
        status = MPI_Irecv((INT *)rxSizes + ll, 1, MPI_INT, kk, 0x5125, MPI_COMM_WORLD, rxRqsts + ll);

        /* Increment the ll */
        ll++;
    }

    return status;
}

STATUS DSLIM_Score::TXResults(MPI_Request *txRqst, INT rawID)
{
    STATUS status = SLM_SUCCESS;
    INT mchID = rawID;
    INT offset = 0;

    if (txSizes[rawID] < 1)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        /* Compute the mchID as:
         * mchID = rawID +1 if >= myid
         */
        if (rawID >= (INT) params.myid)
        {
            mchID = rawID + 1;
        }

        /* Compute the buffer offset for RxValues */
        for (INT kk = 0; kk < rawID; kk++)
        {
            offset += txSizes[kk];
        }

#ifdef DIAGNOSE2
        cout << "TXValues: " << params.myid << " -> " << mchID << " " << offset << endl;
#endif /* DIAGNOSE */

        /* Transmit data to machine:
         * @mch: mchID, @size: rxSizes[rawID] * INT in fResult
         */
        status = MPI_Isend((INT *)((fResult *)(TxValues + offset)),
                           txSizes[rawID] * (sizeof(fResult) / sizeof(INT)),
                           MPI_INT,
                           mchID,
                           0xDA2A,
                           MPI_COMM_WORLD,
                           txRqst);
    }

    return status;
}

STATUS DSLIM_Score::RXResults(MPI_Request *rxRqst, INT rawID)
{
    STATUS status = SLM_SUCCESS;

    UINT mchID = rawID;
    INT offset = 0;

    /* Check the size */
    if (rxSizes[rawID] < 1)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        /* Compute the mchID as:
         * mchID = rawID +1 if >= myid
         */
        if (rawID >= (INT) params.myid)
        {
            mchID = rawID + 1;
        }

        /* Compute the buffer offset for RxValues */
        for (INT kk = 0; kk < rawID; kk++)
        {
            offset += rxSizes[kk];
        }

#ifdef DIAGNOSE2
        cout << "RXValues: " << params.myid << " <- " << mchID << " " << offset << endl;
#endif /* DIAGNOSE */

        /* Receive data from machine:
         * @mch: mchID, @size: rxSizes[rawID] * INT in fResult
         */
        status = MPI_Irecv((INT *)((fResult *)(RxValues + offset)),
                           rxSizes[rawID] * (sizeof(fResult) / sizeof(INT)),
                           MPI_INT,
                           mchID,
                           0xDA2A,
                           MPI_COMM_WORLD,
                           rxRqst);
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
    Wait4RX();

    nSpectra = 0;
    nBatches = 0;
    myRXsize = 0;

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
        RxValues = NULL;
    }

    return;
}
