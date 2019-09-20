/*
 * Copyright (C) 2019 Muhammad Haseeb, Fahad Saeed
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

#ifndef INCLUDE_DSLIM_COMM_H_
#define INCLUDE_DSLIM_COMM_H_

#include "config.h"

#ifdef DISTMEM
/* Include headers */
#include <mpi.h>
#endif /* DISTMEM */

#include <semaphore.h>
#include <unistd.h>
#include "common.h"
#include "slm_dsts.h"
#include "slmerr.h"
#include "utils.h"
#include "msquery.h"
#include "config.h"
#include "lwqueue.h"

/* Class for DSLIM MPI Communication */
class DSLIM_Comm
{
private:

    /* MPI Communication thread */
    THREAD commThd;

    partRes *rxArr;

#ifdef DISTMEM
    /* Handle for Tx requests */
    MPI_Request *TxRqsts;

    /* Handle for Rx requests */
    MPI_Request *RxRqsts;
#endif /* DISTMEM */

    /* Number of allocated work queues */
    UINT    nworkQ;

    /* Semaphores for both queues */
    LOCK work_lock;
    LOCK comm_lock;

    INT rxbuffsize;

    VOID Init_Locks();

    VOID InitTx();

    VOID InitRx();

    VOID Destroy_Locks();

    STATUS Addto_workQ();

    STATUS FlushRxBuffer();

    STATUS DestroyRxBuffers();

    STATUS InitComm_DataTypes();

    STATUS FreeComm_DataTypes();

public:

    /* Default constructor */
    DSLIM_Comm();

    /* Destructor */
    virtual ~DSLIM_Comm();

    partRes *getCurr_WorkArr();

    partRes *getCurr_CommArr();

    STATUS Sendto_workQ();

    STATUS Sendto_commQ();

    STATUS Tx(INT batchnum);

    STATUS Rx(UINT specs, INT batchnum);

    STATUS ComputeResults();

    STATUS Waitfor_TxData();

    STATUS WaitFor_RxData();

    STATUS SignalExit(BOOL &signal);

    STATUS SignalTx();



};

#endif /* INCLUDE_DSLIM_COMM_H_ */
