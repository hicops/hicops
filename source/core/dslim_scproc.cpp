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
#include "dslim_fileout.h"
#include "slmerr.h"
#include "utils.h"
#include "slm_dsts.h"
#include "hicops_instr.hpp"

using namespace std;

/* Global Variables */
BOOL         isCarried   = false;

BData       *bdata       = NULL;

extern gParams           params;

#ifdef USE_MPI
/* Entry function for DSLIM_Score module */
VOID *DSLIM_Score_Thread_Entry(VOID *);
DSLIM_Score *ScoreHandle = NULL;

status_t DSLIM_CarryForward(Index *index, DSLIM_Comm *CommHandle, expeRT *ePtr, hCell *CandidatePSMS, int_t cpsmSize)
{
    status_t status = SLM_SUCCESS;

    bdata = new BData;

    if (bdata != NULL)
    {
        bdata->index     = index;
        bdata->ePtr      = ePtr;
        bdata->heapArray = CandidatePSMS;
        bdata->fileArray = CommHandle->fileArray;
        bdata->sizeArray = CommHandle->sizeArray;
        bdata->nBatches  = CommHandle->nBatches;
        bdata->cPSMsize  = cpsmSize;

        isCarried = true;
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }

    return status;
}
#endif /* USE_MPI */

status_t DSLIM_DistScoreManager()
{
    status_t status = SLM_SUCCESS;

#ifdef USE_MPI
    /* Check if parameters have been brought */
    if (isCarried == false && params.nodes > 1)
    {
        status = ERR_INVLD_PARAM;
    }

    /* No need to do anything if only 1 node */
    if (params.nodes > 1)
    {
        //
        // initialization
        //
        if (status == SLM_SUCCESS)
        {
            ScoreHandle = new DSLIM_Score(bdata);

            if (ScoreHandle == NULL)
            {
                status = ERR_INVLD_MEMORY;
            }
        }

        //
        // combine local results
        //
        if (status == SLM_SUCCESS)
        {
            if (params.myid == 0)
            {
                std::cout << std::endl << "**** Merging Partial Results ****\n" << std::endl;
            }

#if defined (USE_TIMEMORY)
            merge_tuple_t merge_instr("combine");
#endif // USE_TIMEMORY

            status = ScoreHandle->CombineResults();

#if defined (USE_TIMEMORY)
            merge_instr.stop();
#endif // USE_TIMEMORY

            if (params.myid == 0)
            {
                std::cout << std::endl << "Scores Merged with status:\t" << status << std::endl;
            }
        }

        //
        // scatter local scores to relevant nodes
        //

#if defined (USE_TIMEMORY)
        wall_tuple_t sync_penalty("comm_ovhd");
#endif // USE_TIMEMORY

        if (status == SLM_SUCCESS)
        {
            status = ScoreHandle->ScatterScores();

            if (params.myid == 0)
            {
                std::cout << "Scatter Scores with status:\t" << status << std::endl;
            }
        }

        if (status == SLM_SUCCESS)
        {
            status = ScoreHandle->Wait4RX();
        }

#if defined (USE_TIMEMORY)
        sync_penalty.stop();
#endif // USE_TIMEMORY

        //
        // write results to files
        //
        if (status == SLM_SUCCESS)
        {
#if defined (USE_TIMEMORY)
            comm_tuple_t tsv_instr("result_io");
#endif // USE_TIMEMORY

            status = ScoreHandle->DisplayResults();

#if defined (USE_TIMEMORY)
            tsv_instr.stop();
#endif // USE_TIMEMORY

            if (params.myid == 0)
            {
                std::cout << "Writing Results with status:\t" << status << std::endl;
            }
        }

        //
        // synchronization
        //
        if (status == SLM_SUCCESS)
        {
#if defined (USE_TIMEMORY)
            wall_tuple_t sync_penalty("sync_penalty");
            sync_penalty.start();

            // wait for synchronization
            //tim::mpi::barrier(MPI_COMM_WORLD);
            status = MPI_Barrier(MPI_COMM_WORLD);

            sync_penalty.stop();
#else
            MARK_START(sync);

            status = MPI_Barrier(MPI_COMM_WORLD);

            MARK_END(sync)

            if (params.myid == 0)
                std::cout << "Superstep Sync Penalty: " << ELAPSED_SECONDS(sync) << "s" << std::endl<< std::endl;
#endif // USE_TIMEMORY
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
    }

#endif /* USE_MPI */

    /* Return the status of execution */
    return status;
}

#ifdef USE_MPI
/* Entry function for the score communicator */
VOID *DSLIM_Score_Thread_Entry(VOID *argv)
{
    /* Get the number of nodes - 1 */
    int_t nodes   = params.nodes;

    /* MPI pointers */
    MPI_Request *rxRqsts  = NULL;
    int_t         *rxStats  = NULL;

    /* 2X per all other machines */
    rxRqsts = new MPI_Request [nodes];
    rxStats = new int_t[nodes];

    /* Fill the rxStats with ones - available */
    for (int_t kk = 0; kk < nodes; kk++)
    {
        rxStats[kk] = 1;
    }

    /* Avoid race conditions by waiting for
     * ScoreHandle pointer to initialize */
    while ((VOID *)argv != (VOID *)ScoreHandle);

    if (rxRqsts != NULL && rxStats != NULL)
    {
        ScoreHandle->RXSizes(rxRqsts, rxStats);
    }

    /* Wait 500ms between each loop */
    for (int_t cumulate = 0; cumulate < nodes; usleep(500000))
    {
        cumulate = 0;

        for (int_t ll = 0; ll < nodes; ll++)
        {
            /* Only check if not already received */
            if (rxStats[ll] == 0 && ll != (int_t) params.myid)
            {
                //std::cout << (void *)((MPI_Request*)(rxRqsts + ll)) << " rxRqsts@:" << params.myid << std::endl;
                MPI_Test(rxRqsts + ll, &rxStats[ll], MPI_STATUS_IGNORE);

                /* Check if the results
                 * have been received
                 */
                if (rxStats[ll])
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

    /* Fill the rxStats with ones - available */
    for (int_t kk = 0; kk < nodes; kk++)
    {
        rxStats[kk] = 1;
    }

    ScoreHandle->RXResults(rxRqsts, rxStats);

    /* Wait 500ms between each loop */
    for (int_t cumulate = 0; cumulate < nodes; usleep(500000))
    {
        cumulate = 0;

        for (int_t ll = 0; ll < nodes; ll++)
        {
            /* Only check if not already received */
            if (rxStats[ll] == 0 && ll != (int_t) params.myid)
            {
                //std::cout << (void *)((MPI_Request*)(rxRqsts + ll)) << " rxRqsts@:" << params.myid << std::endl;
                MPI_Test(rxRqsts + ll, &rxStats[ll], MPI_STATUS_IGNORE);

                /* Check if the results
                 * have been received
                 */
                if (rxStats[ll])
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

    return NULL;
}
#endif /* USE_MPI */

#ifdef DIAGNOSE
int_t DSLIM_TestBData()
{
    std::cout << "\nTesting BData object: " << std::endl;
    if (bdata != NULL)
    {
        std::cout << bdata->index[0].pepCount << " @: " << params.myid << std::endl;
        std::cout << (void *) &bdata->scPtr[0] << " @: " << params.myid << std::endl;
        std::cout << (void *) &bdata->heapArray[0] << " @: " << params.myid << std::endl;
        std::cout << (void *) &bdata->fileArray[0] << " @: " << params.myid << std::endl;
        std::cout << (void *)&bdata->sizeArray[0] << " @: " << params.myid << std::endl;
        std::cout << (void *) &bdata->indxArray[0] << " @: " << params.myid << std::endl;
        std::cout << bdata->nBatches << std::endl;
        std::cout << (void *) &bdata->resPtr[0] << " @: " << params.myid << std::endl;
        std::cout << bdata->cPSMsize << std::endl;
    }

    std::cout << "Testing BData: DONE " << std::endl;

    return 0;
}
#endif /* DIAGNOSE */
