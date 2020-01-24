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

    /* Initialize the RX lock */
    sem_init(&rxLock, 0, 1);

    /* Create a new rxArr with rxbufsize many partRes */
    rxArr = new partRes[rxbuffsize];

    /* Current number of entries in Rx is zero (offset) */
    currRxOffset = 0;
}

VOID DSLIM_Comm::InitTx()
{
    TxRqsts = NULL;

    /* Initialize the counter */
    sem_init(&txLock, 0, TXARRAYS);

    TxRqsts = new MPI_Request[TXARRAYS];

    if (TxRqsts != NULL)
    {
        /* Initialize the two txArrays and
         * the two corresponding locks */
        for (INT jj = 0; jj < TXARRAYS; jj++)
        {
            txArr[jj] = new partRes[QCHUNK];
        }
    }

    /* Current Tx Array in use is 0th */
    currTxPtr = 0;

    lastTxPtr = 0;
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
    pthread_create(&commThd, NULL, &DSLIM_Comm_Thread_Entry, NULL);
}

/* Destructor */
DSLIM_Comm::~DSLIM_Comm()
{
    VOID *ptr = NULL;

    rxbuffsize = 0;

    nBatches = 0;

    mismatch = 0;

    exitSignal = false;

    Wake4mIO = false;

    /* Destroy semaphores */
    sem_destroy(&wakeup);

    sem_destroy(&control);

    /* Wait for communication thread to complete */
    pthread_join(commThd, &ptr);

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

STATUS DSLIM_Comm::Tx(INT batchtag, INT batchsize)
{
    STATUS status = SLM_SUCCESS;

    INT buff = (lastTxPtr + 1) % TXARRAYS;

    lastTxPtr = buff;

    status = MPI_Isend(txArr[buff],
                       batchsize,
                       resPart,
                       (batchtag % params.nodes),
                       batchtag,
                       MPI_COMM_WORLD,
                       TxRqsts + buff);

    /* Return the head to workQ */
    if (status == SLM_SUCCESS)
    {
        status = MPI_Wait(TxRqsts + buff, NULL);
    }

    /* Free a resource */
    sem_post(&txLock);

    return status;
}

STATUS DSLIM_Comm::Rx(INT batchtag, INT batchsize)
{
    STATUS status = SLM_SUCCESS;

    /* Calculate chunk size and last chunk size */
    INT csize = (batchsize / params.nodes);
    INT lcsize = batchsize - (csize * (params.nodes -1));
    INT mch = 0;
    INT rqst = 0;

    /* Do buffer management here */

    /* Check if myid is the last node */
    if (params.myid == params.nodes - 1)
    {
        /* Receive from all from 0 to myid - 1 */
        for (mch = 0; mch < (INT)params.myid; mch++)
        {
            MPI_Request *request = RxRqsts + mch;

            /* The Rx size = lcsize */
            status = MPI_Irecv(rxArr + (mch * rxbuffsize), lcsize, resPart, mch, batchtag, MPI_COMM_WORLD, request);

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
        for (mch = 0; mch < (INT)params.myid - 1; mch++)
        {
            MPI_Request *request = RxRqsts + rqst;
            status = MPI_Irecv(rxArr + (mch * rxbuffsize), csize, resPart, mch, batchtag, MPI_COMM_WORLD, request);

            if (status != MPI_SUCCESS)
            {
                cout << endl << "Rx Failed at Node: " << params.myid << " for batchnum: " << batchtag << endl;
                break;
            }

            rqst++;
        }

        for (mch = (INT)params.myid + 1; mch < (INT)params.nodes; mch++)
        {
            MPI_Request *request = RxRqsts + rqst;
            status = MPI_Irecv(rxArr + (mch * rxbuffsize), csize, resPart, mch, batchtag, MPI_COMM_WORLD, request);

            if (status != MPI_SUCCESS)
            {
                cout << endl << "Rx Failed at Node: " << params.myid << " for batchnum: " << batchtag << endl;
                break;
            }

            rqst++;
        }
    }

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

    sem_post(&control);

    return status;
}

STATUS DSLIM_Comm::Test4Rx()
{
    STATUS status = MPI_SUCCESS;

/*    for (UINT loop = 0; loop < params.nodes - 1; loop++)
   {
        status = MPI_Test(RxRqsts + loop, MPI_F_STATUS_IGNORE, NULL);
    }
*/
    sem_wait(&control);

    /* Reduce mismatch */
    mismatch--;

    sem_post(&control);

    return status;
}

STATUS DSLIM_Comm::SignalExit()
{
    sem_wait(&control);

    /* Set the exitSignal to true */
    exitSignal = true;

    /*  Release Semaphore */
    SignalWakeup();

    sem_post(&control);

    return SLM_SUCCESS;
}

STATUS DSLIM_Comm::SignalWakeup()
{
    return sem_post(&wakeup);
}

VOID DSLIM_Comm::Destroy_Locks()
{
    (VOID) sem_destroy(&rxLock);

    (VOID) sem_destroy(&txLock);

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

STATUS DSLIM_Comm::RxReady()
{
    /* TODO: Ready the Rx here: MPI_IRecv */

    sem_wait(&control);

    isRxready = true;

    sem_post(&control);

    return SLM_SUCCESS;
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

STATUS DSLIM_Comm::AddBufferEntry()
{
    STATUS status;

    sem_wait (&control);

    /* Check the batch number */
    if (nBatches % params.nodes == params.myid)
    {
        mismatch++;

        if (isRxready == false)
        {
            Wake4mIO = true;

            status = SignalWakeup();
        }
    }

    /* Increment the batch number */
    nBatches++;

    sem_post(&control);

    return status;
}

partRes *DSLIM_Comm::getTxBuffer(INT batchtag, INT batchsize)
{
    partRes *ptr = NULL;

    /* Check if Tx or Rx */
    if (batchtag % params.nodes != params.myid)
    {
        /* Acquire a Tx resource */
        sem_wait(&txLock);

        /* Update the currTxPtr */
        this->currTxPtr = (currTxPtr+1) % TXARRAYS;

        /* Set the return pointer to txArr + currTxPtr */
        ptr = txArr[currTxPtr];
    }
    else
    {
        /* Set the offset to this */
        ptr = rxArr + currRxOffset;

        /* TODO: Update the currRxOffset */
        currRxOffset += batchsize;
    }

    return ptr;
}

#endif /* DISTMEM */
