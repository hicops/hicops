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
#include "dslim_fileout.h"
#include "slmerr.h"
#include "utils.h"
#include "slm_dsts.h"

using namespace std;

/* Global Variables */
BOOL         isCarried   = false;

DSLIM_Score *ScoreHandle = NULL;
BData       *bdata       = NULL;

extern gParams           params;

/* Entry function for DSLIM_Score module */
VOID *DSLIM_Score_Thread_Entry(VOID *);

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

        isCarried = true;
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
        cout << bdata->scPtr[0].res.bias << endl;
        cout << bdata->heapArray[0].hyperscore << endl;
        cout << bdata->fileArray[0] << endl;
        cout << bdata->sizeArray[0] << endl;
        cout << bdata->nBatches << endl;
        cout << bdata->resPtr[0].N << endl;
        cout << bdata->cPSMsize << endl;
    }

    return 0;
}

STATUS DSLIM_DistScoreManager()
{
    STATUS status = SLM_SUCCESS;

    /* Check if parameters have been brought */
    if (isCarried == false)
    {
        status = ERR_INVLD_PARAM;
    }

    /* Initialize the ScoreHandle */
    if (status == SLM_SUCCESS)
    {
        ScoreHandle = new DSLIM_Score(bdata);

        if (ScoreHandle == NULL)
        {
            status = ERR_INVLD_MEMORY;
        }
    }

    /* Distributed Scoring Algorithm */
    if (status == SLM_SUCCESS)
    {
        status = ScoreHandle->ComputeDistScores();
    }

    /* Scatter the key-values to all machines */
    if (status == SLM_SUCCESS)
    {
        status = ScoreHandle->ScatterScores();
    }

    /* Display results to files */
    if (status == SLM_SUCCESS)
    {
        status = ScoreHandle->DisplayResults();
    }

    /* Destroy the ScoreHandle */
    if (status == SLM_SUCCESS)
    {
        if (ScoreHandle != NULL)
        {
            delete ScoreHandle;
            ScoreHandle = NULL;
        }
    }

    /* Return the status of execution */
    return status;
}

/* Entry function for the score communicator */
VOID *DSLIM_Score_Thread_Entry(VOID *argv)
{
    /* Get the number of nodes - 1 */
    INT nodes_1   = params.nodes - 1;

    /* MPI pointers */
    MPI_Request *rxRqsts  = NULL;
    MPI_Request *rxRqsts2 = NULL;
    INT         *rxStats  = NULL;
    INT         *rxStats2 = NULL;

    /* 2X per all other machines */
    rxRqsts = new MPI_Request [2 * (nodes_1)];
    rxStats = new INT[2 * (nodes_1)];

    /* Set the pointers to the appropriate offset */
    rxRqsts2 = rxRqsts + (nodes_1);
    rxStats2 = rxStats + (nodes_1);

    /* Fill the rxStats with zeros - not available */
    std::memset(rxStats, 0x0, (sizeof(rxStats[0]) * (2 * (nodes_1))));

    /* Avoid race conditions by waiting for
     * ScoreHandle pointer to initialize */
    while ((VOID *)argv != (VOID *)ScoreHandle);

    /* Check if everything is in order */
    if (rxRqsts != NULL && rxStats != NULL)
    {
        INT cumulate = 0;

        ScoreHandle->RXSizes(rxRqsts);

        /* Wait 500ms between each loop */
        for (cumulate = 0; cumulate < (nodes_1); usleep(500000))
        {
            for (INT ll = 0; ll < nodes_1; ll ++)
            {
                /* Only check if not already received */
                if (!rxStats[ll])
                {
                    MPI_Test(rxRqsts + ll, &rxStats[ll], MPI_STATUS_IGNORE);

                    /* If the RXsize has been received,
                     * start the second RX request.
                     */
                    if (rxStats[ll])
                    {
                        ScoreHandle->RXResults(rxRqsts2 + ll, ll);
                        cumulate++;
                    }
                }
            }
        }

        /* Wait 500ms between each loop */
        for (cumulate = 0; cumulate < (nodes_1); usleep(500000))
        {
            for (INT ll = 0; ll < nodes_1; ll ++)
            {
                /* Only check if not already received */
                if (!rxStats2[ll])
                {
                    MPI_Test(rxRqsts2 + ll, &rxStats2[ll], MPI_STATUS_IGNORE);

                    /* Check if the results
                     * have been received
                     */
                    if (rxStats2[ll])
                    {
                        cumulate++;
                    }
                }
            }
        }
    }

    /* Deallocate the rxRqsts */
    if (rxRqsts != NULL)
    {
        delete [] rxRqsts;
        rxRqsts = NULL;
    }

    /* Deallocate the rxStats */
    if (rxStats != NULL)
    {
        delete [] rxStats;
        rxStats = NULL;
    }

    /* Return the ScoreHandle back */
    return argv;
}
