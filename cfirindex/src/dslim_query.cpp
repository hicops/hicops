/*
 * This file is part of Load Balancing Algorithm for DSLIM
 *  Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
 *  Florida International University, Miami, FL
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "dslim.h"
#include "hyperscore.h"

extern gParams params;

#ifdef BENCHMARK
static DOUBLE duration = 0;
extern DOUBLE compute;
extern DOUBLE fileio;
extern DOUBLE memory;
#endif /* BENCHMARK */


#ifdef FUTURE
extern pepEntry    *pepEntries; /* SLM Peptide Index */
extern PepSeqs          seqPep; /* Peptide sequences */
#ifdef VMODS
extern varEntry    *modEntries;
#endif /* VMODS */
#endif /* FUTURE */

/* Global Variables */
FLOAT *hyperscores = NULL;
UCHAR *sCArr = NULL;

/* FUNCTION: DSLIM_QuerySpectrum
 *
 * DESCRIPTION: Query the DSLIM for all query peaks
 *              and count the number of hits per chunk
 *
 * INPUT:
 * @QA     : Query Spectra Array
 * @len    : Number of spectra in the array
 * @Matches: Array to fill in the number of hits per chunk
 * @threads: Number of parallel threads to launch
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_QuerySpectrum(ESpecSeqs &ss, UINT len, Index *index, UINT idxchunk)
{
    STATUS status = SLM_SUCCESS;
    UINT *QAPtr = NULL;
    FLOAT *iPtr = NULL;
    UINT maxz = params.maxz;
    UINT dF = params.dF;
    UINT threads = params.threads;
    UINT scale = params.scale;

#ifdef BENCHMARK
    DOUBLE tcons[threads] = {0};
#endif

    DOUBLE maxmass = params.max_mass;

#ifndef _OPENMP
    LBE_UNUSED_PARAM(threads);
#endif

    status = HS_InitFile();

#ifdef BENCHMARK
    duration = omp_get_wtime();
#endif

    /* Process all the queries in the chunk */
    for (UINT queries = 0; queries < len; queries++)
    {
        /* Pointer to each query spectrum */
        QAPtr = ss.moz + ss.idx[queries];
        iPtr = ss.intensity + ss.idx[queries];
        UINT qspeclen = ss.idx[queries + 1] - ss.idx[queries];

        std::cout << '\r' << "DONE: " << queries+1 << "/" << len;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif /* _OPENMP */
        for (UINT ixx = 0; ixx < idxchunk; ixx++)
        {
#ifdef BENCHMARK
            DOUBLE stime = omp_get_wtime();
#endif
            UINT speclen = (index[ixx].pepIndex.peplen - 1) * maxz * iSERIES;

            for (UINT chno = 0; chno < index[ixx].nChunks; chno++)
            {
                /* Query each chunk in parallel */
                UINT *bAPtr  = index[ixx].ionIndex[chno].bA;
                UINT *iAPtr  = index[ixx].ionIndex[chno].iA;
                UCHAR *bcPtr = index[ixx].ionIndex[chno].sc.bc;
                UCHAR *ycPtr = index[ixx].ionIndex[chno].sc.yc;
                FLOAT *ibcPtr = index[ixx].ionIndex[chno].sc.ibc;
                FLOAT *iycPtr = index[ixx].ionIndex[chno].sc.iyc;

                /* Query all fragments in each spectrum */
                for (UINT k = 0; k < qspeclen; k++)
                {
                    /* Check for any zeros
                     * Zero = Trivial query */
                    if (QAPtr[k] > dF && QAPtr[k] < ((maxmass * scale) - 1 - dF))
                    {
                        /* Locate iAPtr start and end */
                        UINT start = bAPtr[QAPtr[k] - dF];
                        UINT end = bAPtr[QAPtr[k] + 1 + dF];

                        /* Loop through located iAions */
                        for (UINT ion = start; ion < end; ion++)
                        {
                            UINT raw = iAPtr[ion];

                            /* Calculate parent peptide ID */
                            UINT ppid = (raw / speclen);

                            /* Update corresponding scorecard entries */
                            if ((raw % speclen) < speclen / 2)
                            {
                                bcPtr[ppid] += 1;
                                ibcPtr[ppid] += iPtr[k];
                            }
                            else
                            {
                                ycPtr[ppid] += 1;
                                iycPtr[ppid] += iPtr[k];
                            }
                        }
                    }
                }

                index[ixx].ionIndex[chno].sc.especid = queries;

                INT idaa = -1;
                FLOAT maxhv = 0.0;

                for (UINT it = 0; it < index[ixx].ionIndex[chno].sc.size; it++)
                {
                    if (bcPtr[it] + ycPtr[it] > params.min_shp)
                    {
                        FLOAT hyperscore = log(HYPERSCORE_Factorial(ULONGLONG(bcPtr[it])) *
                                               HYPERSCORE_Factorial(ULONGLONG(ycPtr[it])) *
                                               ibcPtr[it] *
                                               iycPtr[it]);

                        if (hyperscore > maxhv)
                        {
                            idaa = DSLIM_GenerateIndex(&index[ixx], it);
                            maxhv = hyperscore;
                        }
                    }

                    bcPtr[it] = 0;
                    ycPtr[it] = 0;
                    ibcPtr[it] = 0;
                    iycPtr[it] = 0;
                }

                /* Print the highest hyperscore per chunk */
#pragma omp critical
                {
                    /* Printing the hyperscore in OpenMP mode */
                    status = HYPERSCORE_Calculate(queries, idaa, maxhv);
                }
            }
#ifdef BENCHMARK
            tcons[omp_get_thread_num()] += omp_get_wtime() - stime;
#endif
        }
    }

    std::cout << '\n';

#ifdef BENCHMARK
    compute += omp_get_wtime() - duration;

    for (unsigned int thd = 0; thd < params.threads; thd++)
    {
        std:: cout << "Thread #: " << thd << "\t" << tcons[thd] << std::endl;
    }
#endif

    return status;
}

STATUS DSLIM_DeallocateSC(VOID)
{
    if (sCArr != NULL)
    {
        delete[] sCArr;
        sCArr = NULL;
    }

    return SLM_SUCCESS;
}