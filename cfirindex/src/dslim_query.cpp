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

extern gParams   params;
extern BYICount  *Score;

#ifdef BENCHMARK
static DOUBLE duration = 0;
extern DOUBLE compute;
extern DOUBLE fileio;
extern DOUBLE memory;
#endif /* BENCHMARK */

/* Global variables */
INT qspecnum = -1;
INT chnum    = -1;


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

static VOID DSLIM_BinarySearch(Index *, FLOAT, INT&, INT&);
static INT DSLIM_BinFindMin(pepEntry *entries, FLOAT pmass1, INT min, INT max);
static INT DSLIM_BinFindMax(pepEntry *entries, FLOAT pmass2, INT min, INT max);

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
    UINT maxz = params.maxz;
    UINT dF = params.dF;
    UINT threads = params.threads;
    UINT scale = params.scale;
    DOUBLE maxmass = params.max_mass;

    if (Score == NULL)
    {
        status = ERR_INVLD_MEMORY;
    }

    if (status == SLM_SUCCESS)
    {
        status = HS_InitFile();
    }

#ifdef BENCHMARK
    DOUBLE tcons[threads];
    std::memset(tcons, 0x0, threads * sizeof(DOUBLE));
#endif /* BENCHMARK */

#ifndef _OPENMP
    LBE_UNUSED_PARAM(threads);
#endif /* _OPENMP */

#ifdef BENCHMARK
    duration = omp_get_wtime();
#endif /* BENCHMARK */

    if (status == SLM_SUCCESS)
    {
        /* Process all the queries in the chunk */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif /* _OPENMP */
         for (UINT queries = 0; queries < len; queries++)
         {
#ifdef BENCHMARK
            DOUBLE stime = omp_get_wtime();
#endif
             /* Pointer to each query spectrum */
             UINT *QAPtr = ss.moz + ss.idx[queries];
             UINT *iPtr = ss.intensity + ss.idx[queries];
             UINT qspeclen = ss.idx[queries + 1] - ss.idx[queries];
             UINT thno = omp_get_thread_num();

             UCHAR *bcPtr  = Score[thno].bc;
             UCHAR *ycPtr  = Score[thno].yc;
             UINT  *ibcPtr = Score[thno].ibc;
             UINT  *iycPtr = Score[thno].iyc;

            if (thno == 0)
            {
                std::cout << "\rDONE: " << (queries * 100) /len << "%";
            }

            for (UINT ixx = 0; ixx < idxchunk; ixx++)
            {
                UINT speclen = (index[ixx].pepIndex.peplen - 1) * maxz * iSERIES;

                for (UINT chno = 0; chno < index[ixx].nChunks; chno++)
                {
                    /* Query each chunk in parallel */
                    UINT *bAPtr = index[ixx].ionIndex[chno].bA;
                    UINT *iAPtr = index[ixx].ionIndex[chno].iA;

                    INT minlimit = 0;
                    INT maxlimit = 0;

                    DSLIM_BinarySearch(index + ixx, ss.precurse[queries], minlimit, maxlimit);

                    /* Spectrum violates limits */
                    if (!(maxlimit - minlimit))
                    {
                        continue;
                    }

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
                                INT ppid = (raw / speclen);

                                if (ppid >= minlimit && ppid < maxlimit)
                                {
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
                    }

                    Score[thno].especid = queries;

                    INT idaa = -1;
                    FLOAT maxhv = 0.0;

/*                    UINT csize = ((chno == index[ixx].nChunks - 1) && (index[ixx].nChunks > 1)) ?
                                    index[ixx].lastchunksize : index[ixx].chunksize;
*/
                    UINT csize = maxlimit - minlimit;

                    for (INT it = minlimit; it < maxlimit; it++)
                    {
                        if (bcPtr[it] + ycPtr[it] > params.min_shp)
                        {
                            FLOAT hyperscore = log(0.001 + HYPERSCORE_Factorial(ULONGLONG(bcPtr[it])) *
                                                   HYPERSCORE_Factorial(ULONGLONG(ycPtr[it])) *
                                                   ibcPtr[it] *
                                                   iycPtr[it]) - 6;

                            if (hyperscore > maxhv)
                            {
                                idaa = it; //DSLIM_GenerateIndex(&index[ixx], (csize * chno) + it);
                                maxhv = hyperscore;
                            }
                        }

                    }

                    /* Print the highest hyperscore per chunk */
#pragma omp critical
                    {
                        /* Printing the hyperscore in OpenMP mode */
                        status = HYPERSCORE_Calculate(queries, idaa, maxhv);
                    }

                    /* Clear the scorecard */
                    std::memset(bcPtr+ minlimit, 0x0, sizeof(UCHAR) * csize);
                    std::memset(ycPtr+ minlimit, 0x0, sizeof(UCHAR) * csize);
                    std::memset(ibcPtr+ minlimit, 0x0, sizeof(UINT) * csize);
                    std::memset(iycPtr+ minlimit, 0x0, sizeof(UINT) * csize);

                }
            }

#ifdef BENCHMARK
            tcons[thno] += omp_get_wtime() - stime;
#endif
        }
    }

    std::cout << '\n';

#ifdef BENCHMARK
    compute += omp_get_wtime() - duration;

    for (unsigned int thd = 0; thd < params.threads; thd++)
    {
        std::cout << "Thread #: " << thd << "\t" << tcons[thd] << std::endl;
    }
#endif

    return status;
}

STATUS DSLIM_DeallocateSC(VOID)
{
    /* Free the Scorecard memory */
    if (Score != NULL)
    {
        for (UINT thd = 0; thd < params.threads; thd++)
        {
            delete[] Score[thd].bc;
            delete[] Score[thd].ibc;
            delete[] Score[thd].yc;
            delete[] Score[thd].iyc;

            Score[thd].bc = NULL;
            Score[thd].ibc = NULL;
            Score[thd].yc = NULL;
            Score[thd].iyc = NULL;

            Score[thd].especid = 0;
            Score[thd].size = 0;
        }

        delete[] Score;
        Score = NULL;

    }

    return SLM_SUCCESS;
}

/*
 * The Binary Search Algorithm
 */
static VOID DSLIM_BinarySearch(Index *index, FLOAT precmass, INT &minlimit, INT &maxlimit)
{
    /* Get the FLOAT precursor mass */
    FLOAT pmass1 = precmass - params.dM;
    FLOAT pmass2 = precmass + params.dM;
    pepEntry *entries = index->pepEntries;

    UINT min = 0;
    UINT max = index->lcltotCnt - 1;

    if (params.dM < 0.0)
    {
        minlimit = min;
        maxlimit = max;

        return;
    }

    /* Check for base case */
    if (pmass1 < entries[min].Mass)
    {
        minlimit = min;
    }
    else if (pmass1 > entries[max].Mass)
    {
        minlimit = max;
        maxlimit = max;
        return;
    }
    else
    {
        /* Find the minlimit here */
        minlimit = DSLIM_BinFindMin(entries, pmass1, min, max);
    }

    min = 0;
    max = index->lcltotCnt - 1;


    /* Check for base case */
    if (pmass2 > entries[max].Mass)
    {
        maxlimit = max;
    }
    else if (pmass2 < entries[min].Mass)
    {
        minlimit = min;
        maxlimit = min;
        return;
    }
    else
    {
        /* Find the maxlimit here */
        maxlimit = DSLIM_BinFindMax(entries, pmass2, min, max);
    }

}

static INT DSLIM_BinFindMin(pepEntry *entries, FLOAT pmass1, INT min, INT max)
{
    INT half = (min + max)/2;

    if (max - min < 500)
    {
        INT current = min;

        while (entries[current++].Mass < pmass1);
        return current;
    }

    if (pmass1 > entries[half].Mass)
    {
        min = half;
        return DSLIM_BinFindMin(entries, pmass1, min, max);
    }
    else if (pmass1 < entries[half].Mass)
    {
        max = half;
        return DSLIM_BinFindMin(entries, pmass1, min, max);
    }

    while (pmass1 >= entries[half--].Mass);

    return (half + 1);

}

static INT DSLIM_BinFindMax(pepEntry *entries, FLOAT pmass2, INT min, INT max)
{
    INT half = (min + max)/2;

    if (max - min < 500)
    {
        INT current = max;

        while (entries[current--].Mass > pmass2);
        return current + 1;
    }

    if (pmass2 > entries[half].Mass)
    {
        min = half;
        return DSLIM_BinFindMax(entries, pmass2, min, max);
    }
    else if (pmass2 < entries[half].Mass)
    {
        max = half;
        return DSLIM_BinFindMax(entries, pmass2, min, max);
    }

    while (pmass2 <= entries[half++].Mass);

    return (half - 1);

}
