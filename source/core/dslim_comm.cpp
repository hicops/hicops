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
