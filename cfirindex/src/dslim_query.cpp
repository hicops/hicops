/*
 * This file is part of PCDSFrame software
 * Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
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

#include "dslim_fileout.h"
#include "dslim.h"

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
STATUS DSLIM_QuerySpectrum(Queries &ss, UINT len, Index *index, UINT idxchunk)
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
        status = DFile_InitFiles();
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
             FLOAT pmass = ss.precurse[queries];
             UINT *iPtr = ss.intensity + ss.idx[queries];
             UINT qspeclen = ss.idx[queries + 1] - ss.idx[queries];
             UINT thno = omp_get_thread_num();

             BYC   *bycPtr = Score[thno].byc;
             iBYC *ibycPtr = Score[thno].ibyc;
             Results *resPtr = &Score[thno].res;

            if (thno == 0 && params.myid == 0)
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
                    if ((maxlimit - minlimit) < 1)
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

                                if (ppid >= minlimit && ppid <= maxlimit)
                                {
                                    /* Update corresponding scorecard entries */
                                    if ((raw % speclen) < speclen / 2)
                                    {
                                        bycPtr[ppid].bc += 1;
                                        ibycPtr[ppid].ibc += iPtr[k];
                                    }
                                    else
                                    {
                                        bycPtr[ppid].yc += 1;
                                        ibycPtr[ppid].iyc += iPtr[k];
                                    }
                                }

                            }
                        }
                    }

                    INT csize = maxlimit - minlimit;

                    for (INT it = minlimit; it < maxlimit; it++)
                    {
                        if (bycPtr[it].bc + bycPtr[it].yc >= params.min_shp)
                        {
                            /* Create a heap cell */
                            hCell cell;

                            /* Fill in the information */
                            cell.hyperscore = log10(0.001 +
                                                    DFile_Factorial(ULONGLONG(bycPtr[it].bc)) *
                                                    DFile_Factorial(ULONGLONG(bycPtr[it].yc)) *
                                                    ibycPtr[it].ibc *
                                                    ibycPtr[it].iyc);

                            /* Hyperscore < 0 means either no b- or y- ions were matches */
                            if (cell.hyperscore > 0)
                            {
                                cell.idxoffset = ixx;
                                cell.psid = it;
                                cell.sharedions = bycPtr[it].bc + bycPtr[it].yc;
                                cell.totalions = speclen;

                                /* Insert the cell in the heap dst */
                                resPtr->topK.insert(cell);
                                /* Increase the N */
                                resPtr->cpsms += 1;
                                /* Update the histogram */
                                resPtr->survival[(INT) (cell.hyperscore * 10 + 0.5)] += 1;
                            }
                        }
                    }

                    /* Clear the scorecard */
                    std::memset(bycPtr + minlimit, 0x0, sizeof(BYC) * csize);
                    std::memset(ibycPtr + minlimit, 0x0, sizeof(iBYC) * csize);
                }
            }

            /* Compute expect score if there are any candidate PSMs */
            if (resPtr->cpsms > params.min_cpsm)
            {
                status = DSLIM_ModelSurvivalFunction(resPtr);

                hCell psm = resPtr->topK.getMax();

                /* Estimate the log (s(x)); x = log(hyperscore) */
                DOUBLE lgs_x = resPtr->weight * (psm.hyperscore * 10 + 0.5) + resPtr->bias;

                /* Compute the s(x) */
                DOUBLE s_x = pow(10, lgs_x);

                /* e(x) = n * s(x) */
                DOUBLE e_x = std::min(100.0, resPtr->cpsms * s_x);

                /* Printing the scores in OpenMP mode */
                status = DFile_PrintScore(index, queries, pmass, &psm, e_x, resPtr->cpsms);
            }

            /* Reset the results */
            resPtr->reset();

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

/*
 * DSLIM_DeallocateSC
 */
STATUS DSLIM_DeallocateSC()
{
    /* Free the Scorecard memory */
    if (Score != NULL)
    {
        for (UINT thd = 0; thd < params.threads; thd++)
        {
            delete[] Score[thd].byc;
            delete[] Score[thd].ibyc;
            delete[] Score[thd].res.survival;
            delete[] Score[thd].res.xaxis;

            Score[thd].byc = NULL;
            Score[thd].ibyc = NULL;
            Score[thd].res.survival = NULL;
            Score[thd].res.xaxis = NULL;
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

        while (entries[current].Mass < pmass1)
        {
            current++;
        }

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

    if (pmass1 == entries[half].Mass)
    {
        while (pmass1 == entries[half].Mass)
        {
            half--;
        }

        half++;
    }

    return half;

}

static INT DSLIM_BinFindMax(pepEntry *entries, FLOAT pmass2, INT min, INT max)
{
    INT half = (min + max)/2;

    if (max - min < 500)
    {
        INT current = max;

        while (entries[current].Mass > pmass2)
        {
            current--;
        }

        return current;
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

    if (pmass2 == entries[half].Mass)
    {
        half++;

        while (pmass2 == entries[half].Mass)
        {
            half++;
        }

        half--;
    }

    return half;

}

STATUS DSLIM_ModelSurvivalFunction(Results *resPtr)
{
    STATUS status = SLM_SUCCESS;

    /* Total size and tailends */
    UINT N = resPtr->cpsms;

    /* Choosing the kneePt and endPt to be
     * at 70.7% and 99% respectively */
    UINT kneePt = N - (UINT)((float)N * 0.707);
    UINT endPt = (UINT)((float)N * 0.995);

    /* Copy the slope and bias into local variables */
    DOUBLE slope = resPtr->weight;
    DOUBLE bias = resPtr->bias;

    /* Histogram pointer */
    DOUBLE *histogram = resPtr->survival;
    const INT histosize = 2 + (MAX_HYPERSCORE * 10);

    /* The extracted tail */
    DOUBLE *tail = NULL;
    DOUBLE *axis = NULL;
    INT tailsize = 0;

    /* Loop indexing variable */
    INT ii = 0;

    /* Construct the model */

    /* Initialize the max hyperscore */
    for (ii = histosize - 1; ii > 0; ii--)
    {
        if (histogram[ii] > 0)
        {
            resPtr->maxhypscore = ii;
            ii--;
            break;
        }
    }

    /* Initialize nexthypscore */
    for (; ii > 0; ii--)
    {
        if (histogram[ii] > 0)
        {
            resPtr->nexthypscore = ii;
            break;
        }
    }

    /* Initialize minhypscore */
    for (ii = 0; ii < resPtr->nexthypscore; ii++)
    {
        if (histogram[ii] > 0)
        {
            resPtr->minhypscore = ii;
            break;
        }
    }

    /* Set the minhypscore at kneepoint: ~ 70.7% */
    UINT cumulative = 0;

    for (ii = resPtr->minhypscore; ii < resPtr->nexthypscore; ii++)
    {
        cumulative += histogram[ii];

        /* Mark the lower end of the tail */
        if (cumulative >= kneePt)
        {
            if (ii > resPtr->minhypscore)
            {
                resPtr->minhypscore = ii;
            }

            break;
        }
    }

    /* Set the nexthypscore at endPt: ~99.5% */
    cumulative = N;

    for (ii = resPtr->maxhypscore; ii >= resPtr->minhypscore; ii--)
    {
        cumulative -= histogram[ii];

        /* Mark the upper end of tail */
        if (cumulative <= endPt)
        {
            if (ii < resPtr->nexthypscore)
            {
                resPtr->nexthypscore = ii;
            }

            break;
        }
    }

    /* If both ends at the same point,
     * set the upper end to maxhypscore
     */
    if (resPtr->nexthypscore <= resPtr->minhypscore)
    {
        resPtr->nexthypscore = resPtr->maxhypscore;
    }

    /* Construct s(x) = 1 - CDF in [minhyp, nexthyp] */
    UINT count = histogram[resPtr->nexthypscore];

    for (ii = resPtr->nexthypscore - 1; ii >= resPtr->minhypscore; ii--)
    {
        UINT tmpcount = histogram[ii];

        cumulative -= tmpcount;
        histogram[ii] = (count + histogram[ii + 1]);
        count = tmpcount;
    }

    /* Construct log_10(s(x)) */
    for (ii = resPtr->minhypscore; ii <= resPtr->nexthypscore; ii++)
    {
        histogram[ii] = log10(histogram[ii] / N);
    }

    /* Set the tailPtr, scoreaxis and tailsize */
    tail     = histogram + resPtr->minhypscore;
    axis     = resPtr->xaxis + resPtr->minhypscore;
    tailsize = resPtr->nexthypscore - resPtr->minhypscore + 1;

    /* Perform linear regression (least sq. error) on
     * tail curve and find slope (m) and bias (b)
     */
    (VOID) UTILS_LinearRegression(tailsize, axis, tail, slope, bias);

    /* Assign back from local variables */
    resPtr->weight = slope;
    resPtr->bias   =  bias;


    /* Return the status */
    return status;

}
