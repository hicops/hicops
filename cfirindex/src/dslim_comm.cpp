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


#ifdef DISTMEM

/* We would need a custom MPI
 * Datatype for Tx/Rx */
MPI_Datatype resPart;

/* Work and comm queues */
std::queue <partRes *> workQ;
std::queue <partRes *> commQ;

using namespace std;

extern gParams params;

VOID DSLIM_Comm::InitRx()
{
    RxRqsts = new MPI_Request[params.nodes - 1];

    rxArr = new partRes[(rxbuffsize * params.nodes)];
}

VOID DSLIM_Comm::InitTx()
{
    nworkQ = 0;
    TxRqsts = NULL;

    TxRqsts = new MPI_Request[params.nodes - 1];

    if (TxRqsts != NULL)
    {
        Addto_workQ();
    }
}

/* Default constructor */
DSLIM_Comm::DSLIM_Comm()
{
    INT tmp = (QCHUNK / params.nodes);

    /* Initialize the max buffer sizes */
    if (params.myid < params.nodes - 1)
    {
        rxbuffsize = tmp;
    }
    else
    {
        rxbuffsize = QCHUNK - (tmp * (params.nodes - 1));
    }

    Init_Locks();
    InitComm_DataTypes();
    InitTx();
    InitRx();
}

/* Destructor */
DSLIM_Comm::~DSLIM_Comm()
{
    if (rxArr != NULL)
    {
        delete[] rxArr;
        rxArr = NULL;
    }

    rxbuffsize = 0;

    sem_wait(&work_lock);

    for (UINT ii = 0; ii < workQ.size(); ii++)
    {
        partRes *txArray = workQ.front();

        workQ.pop();

        delete[] txArray;
    }

    sem_post(&work_lock);

    Destroy_Locks();

    FreeComm_DataTypes();

}

STATUS DSLIM_Comm::InitComm_DataTypes()
{
    MPI_Type_contiguous((INT)(sizeof(partRes) / sizeof (FLOAT)),
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
    partRes *ret = NULL;

    while (true)
    {
        sem_wait(&work_lock);

        INT size = workQ.size();

        sem_post(&work_lock);

        if (size > 0)
        {
            ret = workQ.front();
            break;
        }
    }

    return ret;
}

partRes * DSLIM_Comm::getCurr_CommArr()
{
    partRes *ret = NULL;

    if (commQ.size() > 0)
    {
        ret = commQ.front();
    }

    return ret;
}

STATUS DSLIM_Comm::Addto_workQ()
{
    STATUS status = SLM_SUCCESS;

    partRes *res = new partRes[QCHUNK];

    if (res != NULL)
    {
        nworkQ ++;

        sem_wait(&work_lock);

        workQ.push(res);

        sem_post(&work_lock);
    }
    else
    {
        status = ERR_BAD_MEM_ALLOC;
    }

    return status;
}

STATUS DSLIM_Comm::Sendto_workQ()
{
    STATUS status = SLM_SUCCESS;

    /* Pop from the comm queue */
    partRes *head = commQ.front();

    commQ.pop();

    /* Push to work queue */
    sem_wait(&work_lock);

    workQ.push(head);

    sem_post(&work_lock);

    return status;
}

STATUS DSLIM_Comm::Sendto_commQ()
{
    STATUS status = SLM_SUCCESS;

    if (workQ.size() < 1)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        status = sem_wait(&work_lock);

        partRes *head = workQ.front();

        workQ.pop();

        status = sem_post(&work_lock);

        /* Add to comm queue */
        commQ.push(head);
    }

    return status;
}

STATUS DSLIM_Comm::Tx(INT batchnum)
{
    STATUS status = SLM_SUCCESS;
    /* Calculate chunk size and last chunk size */
//    INT csize = (specs / params.nodes);
//    INT lcsize = specs - (csize * (params.nodes -1));
    //fixme partRes *head = NULL;

    if (commQ.size() < 1)
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
        status = this->Sendto_workQ();
    }

    return status;
}

STATUS DSLIM_Comm::Rx(UINT specs, INT batchnum)
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
            status = MPI_Irecv(rxArr + (mch * rxbuffsize), lcsize, resPart, mch, batchnum, MPI_COMM_WORLD, request);

            if (status != MPI_SUCCESS)
            {
                cout << endl << "Rx Failed at Node: " << params.myid << " for batchnum: " << batchnum << endl;
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
            status = MPI_Irecv(rxArr + (mch * rxbuffsize), csize, resPart, mch, batchnum, MPI_COMM_WORLD, request);

            if (status != MPI_SUCCESS)
            {
                cout << endl << "Rx Failed at Node: " << params.myid << " for batchnum: " << batchnum << endl;
                break;
            }

            rqst++;
        }

        for (mch = (INT)params.myid + 1; mch < (INT)params.nodes; mch++)
        {
            MPI_Request *request = RxRqsts + rqst;
            status = MPI_Irecv(rxArr + (mch * rxbuffsize), csize, resPart, mch, batchnum, MPI_COMM_WORLD, request);

            if (status != MPI_SUCCESS)
            {
                cout << endl << "Rx Failed at Node: " << params.myid << " for batchnum: " << batchnum << endl;
                break;
            }

            rqst++;
        }
    }

    return status;
}

STATUS DSLIM_Comm::Waitfor_TxData()
{
    return sem_wait(&comm_lock);
}

STATUS DSLIM_Comm::WaitFor_RxData()
{
    STATUS status;

    for (UINT loop = 0; loop < params.nodes - 1; loop++)
    {
        status = MPI_Wait(RxRqsts + loop, NULL);
    }

    return status;
}

STATUS DSLIM_Comm::SignalExit(BOOL &signal)
{
    signal = true;

    return sem_post(&comm_lock);
}

STATUS DSLIM_Comm::SignalTx()
{
    STATUS status = ERR_INVLD_PARAM;

    status = this->Sendto_commQ();

    if (status == SLM_SUCCESS)
    {
        status = sem_post(&comm_lock);
    }

    return status;
}

VOID DSLIM_Comm::Init_Locks()
{
    (VOID) sem_init(&comm_lock, 0, 0);
    (VOID) sem_init(&work_lock, 0, 1);
}

VOID DSLIM_Comm::Destroy_Locks()
{
    (VOID) sem_destroy(&comm_lock);
    (VOID) sem_destroy(&work_lock);
}

#endif /* DISTMEM */
