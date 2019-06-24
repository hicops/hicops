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

extern gParams params;


#ifdef FUTURE
extern pepEntry    *pepEntries; /* SLM Peptide Index */
extern PepSeqs          seqPep; /* Peptide sequences */
#ifdef VMODS
extern varEntry    *modEntries;
#endif /* VMODS */
#endif /* FUTURE */

/* Global Variables */
#define SCSIZE              10000000
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
STATUS DSLIM_QuerySpectrum(UINT *QA, UINT len, Index *index, UINT idxchunk)
{
    STATUS status = SLM_SUCCESS;
    UINT *QAPtr = NULL;
    UINT maxz = params.maxz;
    UINT dF = params.dF;
    UINT threads = params.threads;

    DOUBLE maxmass = params.max_mass;

#ifdef _OPENMP
    if (sCArr == NULL)
    {
        sCArr = new UCHAR[SCSIZE * threads];
        std::memset(sCArr, 0x0, sizeof(UCHAR) * SCSIZE * threads);
    }
#else
    if (sCArr == NULL)
    {
        sCArr = new UCHAR[SCSIZE];
        std::memset(sCArr, 0x0, sizeof(UCHAR) * SCSIZE);
    }

#endif /* _OPENMP */

    if (sCArr == NULL)
    {
        return ERR_BAD_MEM_ALLOC;
    }

#ifdef _OPENMP
#pragma omp parallel  num_threads(threads)
    {
        UINT threadMatches = 0;

#pragma omp for schedule(dynamic, 10) private(threadMatches)
#endif
        /* Process all the queries in the chunk */
        for (UINT queries = 0; queries < len; queries++)
        {
            /* Pointer to each query spectrum */
            QAPtr = QA + (queries * QALEN);

#ifdef _OPENMP
            UCHAR *SCPtr = sCArr + (omp_get_thread_num() * SCSIZE);
#else
            UCHAR *SCPtr = sCArr;
#endif

            for (UINT ixx = 0; ixx < idxchunk; ixx++)
            {
                UINT speclen = (index[ixx].pepIndex.peplen - 1) * maxz * iSERIES;

                for (UINT chno = 0; chno < index[ixx].nChunks; chno++)
                {
                    /* Query each chunk in parallel */
                    UINT *bAPtr = index[ixx].ionIndex[chno].bA;
                    UINT *iAPtr = index[ixx].ionIndex[chno].iA;

                    /* Check if this chunk is the last chunk */
                    UINT size = index[ixx].chunksize;

                    /* Query all fragments in each spectrum */
                    for (UINT k = 0; k < QALEN; k++)
                    {
                        /* Check for any zeros
                         * Zero = Trivial query */
                        if (QAPtr[k] < dF || QAPtr[k] > ((maxmass * scale) - 1 - dF))
                        {
                            continue;
                        }

                        /* Locate iAPtr start and end */
                        UINT start = bAPtr[QAPtr[k] - dF];
                        UINT end = bAPtr[QAPtr[k] + 1 + dF];

                        /* Loop through located iAions */
                        for (UINT ion = start; ion < end; ion++)
                        {
                            /* Calculate parent peptide ID */
                            UINT ppid = (iAPtr[ion] / speclen);

                            /* Update corresponding SC entry */
                            SCPtr[ppid] += 1;

                        }
                    }

                    /* Count the number of candidate peptides
                     * from each chunk
                     */
                    UINT localMatches = 0;

                    for (UINT cntr = 0; cntr < size; cntr++)
                    {
                        if (SCPtr[cntr] >= params.min_shp)
                        {
                            localMatches++;
                        }
                    }

                    /* Avoid too many updates to the Matches */
                    threadMatches += localMatches;

                    /* bitmask not active,
                     * reset the SC instead */
                    std::memset(SCPtr, 0x0, size);
                }
            }
        }

    }

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
