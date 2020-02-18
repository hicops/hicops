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
 * GNU General Public License for more detailSpectrum.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */
#include "dslim_comm.h"

using namespace std;

#ifdef DISTMEM

/* The custom MPI_Datatype for Tx/Rx */
MPI_Datatype resPart;

extern VOID *DSLIM_Comm_Thread_Entry(VOID *argv);

/* Global params */
extern gParams params;

VOID DSLIM_Comm::InitRx()
{
    /* We may need N-1 RxRequests at any time */
    RxRqsts = new MPI_Request[params.nodes - 1];
    RxStat = new INT[params.nodes - 1];

    for (UINT jj = 0; jj < params.nodes - 1; jj++)
    {
        /* True means free */
        RxStat[jj] = 1;
    }

    /* Initialize the RX lock */
    sem_init(&rxLock, 0, 1);

    /* Create a new rxArr with rxbufsize many partRes */
    rxArr = new partRes[rxbuffsize];

    /* Current number of entries in Rx is zero (offset) */
    currRxOffset = 0;

    /* Will be updated at each AddBuffer */
    RxTag = params.myid - params.nodes;

    /* Init the rxQueue */
    rxQueue = new lwqueue<UINT>(20);
}

VOID DSLIM_Comm::InitTx()
{
    TxRqsts = NULL;
    TxStat = NULL;

    TxRqsts = new MPI_Request[TXARRAYS];
    TxStat = new INT[TXARRAYS];

    if (TxRqsts != NULL && TxStat != NULL)
    {
        /* Initialize the two txArrays and
         * the two corresponding locks */
        for (INT jj = 0; jj < TXARRAYS; jj++)
        {
            txArr[jj] = new partRes[QCHUNK];

            /* True means available */
            TxStat[jj] = 1;
        }
    }

}

/* Default constructor */
DSLIM_Comm::DSLIM_Comm()
{
    /* Number of batches processed =  0 */
    nBatches = 0;

    /* Set the exitSignal to false */
    exitSignal = false;

    /* Set the signals */
    Wake4mIO = false;

    /* No mismatch in beginning */
    mismatch = 0;

    /* However may can fit in the RXBUFFERSIZE memory */
    rxbuffsize = RXBUFFERSIZE / (QCHUNK * sizeof(partRes));

    rxbuffsize *= QCHUNK;

    /* Initialize the wakeup signal to zero */
    sem_init(&wakeup, 0, 0);

    /* Initialize the control lock */
    sem_init(&control, 0, 1);

    /* Initialize the custom MPI_DataType */
    InitComm_DataTypes();

    /* Initialize Tx */
    InitTx();

    /* Initialize Rx */
    InitRx();

    /* Create the entry thread here */
    pthread_create(&commThd, NULL, &DSLIM_Comm_Thread_Entry, this);
}

/* Destructor */
DSLIM_Comm::~DSLIM_Comm()
{
    VOID *ptr = NULL;

    /* Wait for communication thread to complete */
    pthread_join(commThd, &ptr);

    rxbuffsize = 0;

    nBatches = 0;

    mismatch = 0;

    exitSignal = false;

    Wake4mIO = false;

    /* Destroy semaphores */
    sem_destroy(&wakeup);

    sem_destroy(&control);

    /* Destroy the Rx queue */
    if (rxQueue != NULL)
    {
        delete rxQueue;

        rxQueue  = NULL;
    }

    /* Deallocate the Rx array */
    if (rxArr != NULL)
    {
        delete[] rxArr;
        rxArr = NULL;
    }

    /* Deallocate the txArrays */
    for (INT jj = 0; jj < TXARRAYS; jj++)
    {
        if (txArr[jj] != NULL)
        {
            delete [] txArr[jj];
            txArr[jj] = NULL;
        }
    }

    /* Destroy Tx and Rx locks */
    Destroy_Locks();

    /* Deallocate Tx requests */
    if (TxRqsts != NULL)
    {
        delete[] TxRqsts;
        TxRqsts = NULL;
    }

    /* Deallocate Rx requests */
    if (RxRqsts != NULL)
    {
        delete[] RxRqsts;
        RxRqsts = NULL;
    }

    /* Decommit the MPI DataTypes */
    FreeComm_DataTypes();

}

STATUS DSLIM_Comm::InitComm_DataTypes()
{
    MPI_Type_contiguous((INT)(sizeof(partRes) / sizeof(FLOAT)),
                        MPI_FLOAT,
                        &resPart);

    MPI_Type_commit(&resPart);

    return SLM_SUCCESS;
}

STATUS DSLIM_Comm::FreeComm_DataTypes()
{
    MPI_Type_free(&resPart);

    return SLM_SUCCESS;
}

STATUS DSLIM_Comm::Tx(INT batchtag, INT batchsize, INT buff)
{
    STATUS status = SLM_SUCCESS;

    /* Print the diagnostic */
    cout << "\nTX: " << batchtag << " @node: " << (batchtag % params.nodes) << endl;

    status = MPI_Isend(txArr[buff],
                       batchsize,
                       resPart,
                       (batchtag % params.nodes),
                       batchtag,
                       MPI_COMM_WORLD,
                       TxRqsts + buff);

    return status;
}

STATUS DSLIM_Comm::Rx(INT batchtag, INT batchsize)
{
    STATUS status = SLM_SUCCESS;

    /* Calculate chunk size and last chunk size */
    INT mch = 0;
    INT rqst = 0;

    /* Print the diagnostics */
    cout << "RX: " << batchtag << " @node: " << params.myid << endl;

    /* Lock the rxArray */
    sem_wait(&rxLock);

    /* Check if myid is the last node */
    if (params.myid == params.nodes - 1)
    {
        /* Receive from all from 0 to params.nodes - 1 */
        for (mch = 0; mch < (INT)params.nodes - 1; mch++)
        {
            MPI_Request *request = RxRqsts + mch;

            /* The Rx size = lcsize */
            status = MPI_Irecv(rxArr + currRxOffset,
                               batchsize,
                               resPart,
                               mch,
                               batchtag,
                               MPI_COMM_WORLD,
                               request);

            /* Mark the Rx status as used */
            RxStat[mch] = 0;

            currRxOffset += batchsize;

            if (status != MPI_SUCCESS)
            {
                cout << endl << "Rx Failed at Node: " << params.myid << " for batchnum: " << batchtag << endl;
                break;
            }
        }
    }
    else
    {
        /* Receive from all but myid */
        for (mch = 0; mch < (INT) params.myid; mch++)
        {
            MPI_Request *request = RxRqsts + rqst;

            status = MPI_Irecv(rxArr + currRxOffset,
                               batchsize,
                               resPart,
                               mch,
                               batchtag,
                               MPI_COMM_WORLD,
                               request);

            /* Mark the Rx status as used */
            RxStat[mch] = 0;

            /* Update the pointer */
            currRxOffset += batchsize;

            if (status != MPI_SUCCESS)
            {
                cout << endl << "Rx Failed at Node: " << params.myid << " for batchnum: " << batchtag << endl;
                break;
            }

            rqst++;
        }

        for (mch = (INT) params.myid + 1; mch < (INT)params.nodes; mch++)
        {
            MPI_Request *request = RxRqsts + rqst;

            status = MPI_Irecv(rxArr + currRxOffset,
                               batchsize,
                               resPart,
                               mch,
                               batchtag,
                               MPI_COMM_WORLD,
                               request);

            /* Mark the Rx status as used */
            RxStat[mch] = 0;

            /* Update the pointer */
            currRxOffset += batchsize;

            if (status != MPI_SUCCESS)
            {
                cout << endl << "Rx Failed at Node: " << params.myid << " for batchnum: " << batchtag << endl;
                break;
            }

            rqst++;
        }
    }

    /* Unlock the rxArray */
    sem_post(&rxLock);

    return status;
}

/* An event is detected when the wakeup semaphore is signaled */
STATUS DSLIM_Comm::Wait4Event()
{
    return sem_wait(&wakeup);
}

STATUS DSLIM_Comm::Wait4Rx()
{
    STATUS status = MPI_SUCCESS;

    for (UINT loop = 0; loop < params.nodes - 1; loop++)
    {
        status = MPI_Wait(RxRqsts + loop, NULL);
    }

    sem_wait(&control);

    /* Reduce mismatch */
    mismatch--;

    isRxready = false;

    sem_post(&control);

    return status;
}

STATUS DSLIM_Comm::CheckRx()
{
    STATUS status = MPI_SUCCESS;

    UINT cumulate = 0;

    /* Check the Rx status */
    for (UINT loop = 0; loop < params.nodes - 1; loop++)
    {
        /* Check if used - false */
        if (!RxStat[loop])
        {
            status = MPI_Test(RxRqsts + loop, &RxStat[loop], MPI_STATUS_IGNORE);
        }

        /* Check if the request is complete */
        if (RxStat[loop])
        {
            cumulate++;
        }
    }

    /* Check if all requests completed */
    if (cumulate == params.nodes - 1)
    {
        sem_wait(&control);

        /* Reduce mismatch if > 0 */
        if (mismatch > 0)
        {
            /* Print the diagnostic */
            cout << "RX Complete: " << RxTag << " @node: " << params.myid << endl;

            mismatch--;
        }

        /* Set isRxReady to false */
        isRxready = false;

        sem_post(&control);
    }

    return status;
}

STATUS DSLIM_Comm::SignalExit()
{
    sem_wait(&control);

    /* Set the exitSignal to true */
    exitSignal = true;

    sem_post(&control);

    /*  Release Semaphore */
    SignalWakeup();

    return SLM_SUCCESS;
}

STATUS DSLIM_Comm::SignalWakeup()
{
    return sem_post(&wakeup);
}

VOID DSLIM_Comm::Destroy_Locks()
{
    (VOID) sem_destroy(&rxLock);

}

BOOL DSLIM_Comm::checkExitSignal()
{
    sem_wait(&control);

    BOOL signal = this->exitSignal;

    sem_post(&control);

    return signal;
}

BOOL DSLIM_Comm::checkWakeup()
{
    sem_wait(&control);

    BOOL state = Wake4mIO;

    if (state == true)
    {
        Wake4mIO = false;
    }

    sem_post(&control);

    return state;
}

STATUS DSLIM_Comm::Rx()
{
    STATUS status = SLM_SUCCESS;

    /* Get the batchtag and the batchsize */
    INT btag = RxTag + params.nodes;
    INT bsize = rxQueue->front();

    status = rxQueue->pop();

    /* Initialize the MPI_Recv */
    if (status == SLM_SUCCESS)
    {
        status = Rx(btag, bsize);
    }

    /* Set the flag to true if success */
    if (status == SLM_SUCCESS)
    {
        status = sem_wait(&control);

        isRxready = true;

        status = sem_post(&control);
    }

    return status;
}

BOOL DSLIM_Comm::getRxReadyPermission()
{
    sem_wait(&control);

    BOOL state = (mismatch > 0 && isRxready == false);

    sem_post(&control);

    return state;
}

BOOL DSLIM_Comm::checkMismatch()
{
    sem_wait(&control);

    BOOL state = (mismatch > 0);

    sem_post(&control);

    return state;
}

STATUS DSLIM_Comm::AddBufferEntry(INT bsize)
{
    STATUS status;

    /* Check the batch number */
    if (nBatches % params.nodes == params.myid)
    {
        /* Push the batchsize to the rxQueue */
        rxQueue->push(bsize);

        sem_wait (&control);

        if (mismatch == 0 && isRxready == false)
        {
            Wake4mIO = true;

            mismatch++;

            status = SignalWakeup();
        }
        else
        {
            mismatch++;
        }

        sem_post(&control);

        /* Increment the batch number */
        nBatches++;
    }

    return status;
}

partRes *DSLIM_Comm::getTxBuffer(INT batchtag, INT batchsize, INT &buffer)
{
    partRes *ptr = NULL;
    buffer = -1;

    /* Check if Tx or Rx */
    if (batchtag % params.nodes != params.myid)
    {
        INT buff = 0;

        /* Wait for an empty buffer */
        while (true)
        {
            /* Loop through all available buffers */
            for (buff = 0; buff < TXARRAYS; buff++)
            {
                /* Only check if it's not available by status */
                if (!TxStat[buff])
                {
                    MPI_Test(TxRqsts + buff, &TxStat[buff], MPI_STATUS_IGNORE);
                }

                /* Found an empty buffer */
                if (TxStat[buff])
                {
                    break;
                }

            }

            /* Check again if really empty */
            if (buff < TXARRAYS && TxStat[buff])
            {
                /* Set the return pointer to txArr[buff] */
                ptr = txArr[buff];

                /* Set the TxStat to false - Used henceforth */
                TxStat[buff] = false;

                buffer = buff;

                break;
            }

            /* No free buffer found - Wait 500ms */
            sleep(0.5);
        }
    }
    else
    {
        sem_wait(&rxLock);

        /* Set the offset to this */
        ptr = rxArr + currRxOffset;

        /* Update the currRxOffset */
        currRxOffset += batchsize;

        sem_post(&rxLock);
    }

    return ptr;
}

#endif /* DISTMEM */
