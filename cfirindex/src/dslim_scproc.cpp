/*
 * This file is part of PCDSFrame software
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
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <unistd.h>
#include "common.h"
#include "dslim.h"
#include "dslim_score.h"
#include "dslim_comm.h"
#include "slmerr.h"
#include "utils.h"
#include "slm_dsts.h"

using namespace std;

/* Global Variables */
DSLIM_Score *ScoreHandle = NULL;
BData       *bdata       = NULL;

STATUS DSLIM_CarryForward(Index *index, DSLIM_Comm *CommHandle, BYICount *Score, hCell *CandidatePSMS, INT cpsmSize)
{
    STATUS status = SLM_SUCCESS;

    bdata = new BData;

    if (bdata != NULL)
    {
        bdata->index     = index;
        bdata->scPtr     = Score;
        bdata->heapArray = CandidatePSMS;
        bdata->fileArray = CommHandle->fileArray;
        bdata->sizeArray = CommHandle->sizeArray;
        bdata->nBatches  = CommHandle->sizeOffset;
        bdata->resPtr    = CommHandle->rxArr;
        bdata->cPSMsize  = cpsmSize;
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }

    return status;
}

INT TestBData()
{
    if (bdata != NULL)
    {
        cout << bdata->index[0].chunksize << endl;
        cout << bdata->scPtr[0].res.bias<< endl;
        cout << bdata->heapArray[0].hyperscore<< endl;
        cout << *bdata->fileArray<< endl;
        cout << *bdata->sizeArray<< endl;
        cout << bdata->nBatches << endl;
        cout << bdata->resPtr[0].N << endl;
        cout << bdata->cPSMsize<< endl;
    }

    return 0;
}

STATUS DSLIM_InitDistScore()
{
    STATUS status = SLM_SUCCESS;

    ScoreHandle = new DSLIM_Score(bdata);

    if (ScoreHandle == NULL)
    {
        status = ERR_INVLD_MEMORY;
    }

    return status;
}

STATUS DSLIM_DeInitDistScore()
{
    if (ScoreHandle != NULL)
    {
        delete ScoreHandle;
        ScoreHandle = NULL;
    }

    return SLM_SUCCESS;
}
