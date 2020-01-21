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

/* Entry function for DSLIM Comm */
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
    TxRqst = NULL;

    TxRqst = new MPI_Request;

    if (TxRqst != NULL)
    {
        /* Initialize the two txArrays and
         * the two corresponding locks */
        for (INT jj = 0; jj < TXARRAYS; jj++)
        {
            txArr[jj] = new partRes[QCHUNK];
            sem_init(txLock + jj, 0, 1);
        }
    }

    /* Current Tx Array in use is 0th */
    currTxPtr = 0;
}

/* Default constructor */
DSLIM_Comm::DSLIM_Comm()
{
    /* Number of batches processed =  0 */
    nBatches = 0;

    /* Set the exitSignal to false */
    exitSignal = false;

    /* However may can fit in the RXBUFFERSIZE memory */
    rxbuffsize = RXBUFFERSIZE / (QCHUNK * sizeof(partRes));

    /* Initialize the custom MPI_DataType */
    InitComm_DataTypes();

    /* Initialize the Tx and Rx */
    InitTx();

    InitRx();

    /* Set up the entry thread here */
    pthread_create(&commThd, NULL, &DSLIM_Comm_Thread_Entry, NULL);
}

/* Destructor */
DSLIM_Comm::~DSLIM_Comm()
{
    VOID *ptr = NULL;

    /* Wait for communication thread to complete */
    pthread_join(commThd, &ptr);

    if (rxArr != NULL)
    {
        delete[] rxArr;
        rxArr = NULL;
    }

    rxbuffsize = 0;

    Destroy_Locks();

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

partRes * DSLIM_Comm::getCurr_WorkArr()
{
    return NULL;
}

partRes * DSLIM_Comm::getCurr_CommArr()
{
    return NULL;
}

STATUS DSLIM_Comm::Tx(INT batchtag)
{
    STATUS status = SLM_SUCCESS;
    /* Calculate chunk size and last chunk size */
//    INT csize = (specs / params.nodes);
//    INT lcsize = specs - (csize * (params.nodes -1));
    //FIXME partRes *head = NULL;

    if (1)//commQ.size() < 1)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        //fixme head = commQ.front();
        INT mch = 0;

        for (UINT loop = 0; loop < params.nodes && status == SLM_SUCCESS; loop++)
        {
            if (loop == params.myid)
            {
                continue;
            }

//            INT txsize = csize;

            if (loop == params.nodes - 1)
            {
//                txsize = lcsize;
            }

            //FIXME status = MPI_Isend(head + (loop * csize), txsize, resPart, loop, 0, MPI_COMM_WORLD, TxRqsts + mch);
            mch++;
        }

        /* Wait for all requests to complete */
        for (UINT loop = 0; loop < params.nodes - 1; loop++)
        {
            //FIXME status = MPI_Wait(TxRqsts + loop, NULL);
        }
    }

    /* Return the head to workQ */
    if (status == SLM_SUCCESS)
    {
        //status = this->Sendto_workQ();
    }

    return status;
}

STATUS DSLIM_Comm::Rx(UINT specs, INT batchtag)
{
    STATUS status = SLM_SUCCESS;

    /* Calculate chunk size and last chunk size */
    INT csize = (specs / params.nodes);
    INT lcsize = specs - (csize * (params.nodes -1));
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

STATUS DSLIM_Comm::Waitfor_TxData()
{
    //return sem_wait(&comm_lock);
}

STATUS DSLIM_Comm::WaitFor_RxData()
{
    STATUS status = MPI_SUCCESS;

    for (UINT loop = 0; loop < params.nodes - 1; loop++)
    {
        status = MPI_Wait(RxRqsts + loop, NULL);
    }

    return status;
}

STATUS DSLIM_Comm::SignalExit(BOOL &signal)
{

}

STATUS DSLIM_Comm::SignalTx()
{
    STATUS status = ERR_INVLD_PARAM;

    //status = this->Sendto_commQ();

    if (status == SLM_SUCCESS)
    {
//        status = sem_post(&comm_lock);
    }

    return status;
}

VOID DSLIM_Comm::Destroy_Locks()
{
    (VOID) sem_destroy(&rxLock);

    for (INT jj = 0; jj < TXARRAYS; jj++)
        (VOID) sem_destroy(&txLock[jj]);
}

#endif /* DISTMEM */
