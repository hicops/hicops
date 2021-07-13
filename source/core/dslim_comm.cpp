/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#include "dslim_comm.h"

#ifdef USE_MPI

/* Global params */
extern gParams params;

DSLIM_Comm::DSLIM_Comm()
{
    nBatches = RXBUFFERSIZE / (QCHUNK * sizeof(partRes));

    fileArray = new int_t[nBatches];
    sizeArray = new int_t[nBatches];

    nBatches = 0;
    myRXsize = 0;
}

DSLIM_Comm::DSLIM_Comm(int_t tbatches)
{
    auto nodes = params.nodes;
    auto remaining = (tbatches % nodes);

    nBatches = tbatches / nodes;

    if ((remaining > 0) && (remaining > params.myid))
        nBatches += 1;

    if (nBatches > 0)
    {
        fileArray = new int_t[nBatches];
        sizeArray = new int_t[nBatches];
    }

    nBatches = 0;
    myRXsize = 0;
}

DSLIM_Comm::~DSLIM_Comm()
{
    fileArray = NULL;
    sizeArray = NULL;

    nBatches = 0;
    myRXsize = 0;
}

status_t DSLIM_Comm::AddBatch(int_t batchNum, int_t batchSize, int_t fileID)
{
    auto nodes = params.nodes;

    if (batchNum % nodes == params.myid)
    {
        auto position = batchNum/params.nodes;

        fileArray[position] = fileID;
        sizeArray[position] = batchSize;

        nBatches += 1;
        myRXsize += batchSize;
    }

    return SLM_SUCCESS;
}

#endif /* USE_MPI */
