/*
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

#include "dslim.h"
using namespace std;

/* Global Variables */
uint_t          *SpecArr = NULL; /* Spectra Array     */
BYICount      *Score   = NULL;
uint_t reduce = 0;

extern gParams params;

#ifdef GENFEATS
BOOL          *features;
#endif


/* FUNCTION: DSLIM_Construct
 *
 * DESCRIPTION: Construct DSLIM chunks
 *
 * INPUT:
 * @threads: Number of parallel threads to launch
 * @modInfo: DSLIM mods Information
 *
 * OUTPUT:
 * @status: status of execution
 */
status_t DSLIM_Construct(Index *index)
{
    status_t status = SLM_SUCCESS;
    uint_t threads = params.threads;
    uint_t maxz = params.maxz;
    uint_t peplen_1 = index->pepIndex.peplen - 1;

    double_t maxmass = params.max_mass;
    uint_t scale = params.scale;

    if (status == SLM_SUCCESS && SpecArr == NULL)
    {
        /* Spectra Array (SA) */
        SpecArr = new uint_t[MAX_IONS];

        /* Check if Spectra Array has been allocated */
        if (SpecArr == NULL)
            status = ERR_BAD_MEM_ALLOC;

    }

    if (status == SLM_SUCCESS)
    {
        /* Allocate memory for all chunks */
        status = DSLIM_AllocateMemory(index);

        /* Construct DSLIM.iA */
        if (status == SLM_SUCCESS)
        {
            /* Distributed SLM Ions Array construction */
            for (uint_t chno = 0; chno < index->nChunks && status == SLM_SUCCESS; chno++)
            {
                /* Construct each DSLIM chunk in Parallel */
                status = DSLIM_ConstructChunk(threads, index, chno);

                /* Apply SLM-Transform on the chunk */
                if (status == SLM_SUCCESS)
                    status = DSLIM_SLMTransform(threads, index, chno);

            }
        }
    }

    if (status == SLM_SUCCESS)
    {
        const uint_t speclen = peplen_1 * maxz * iSERIES;

        /* Construct DSLIM.bA */
#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif /* USE_OMP */
        for (uint_t chunk_number = 0; chunk_number < index->nChunks; chunk_number++)
        {
            uint_t *bAPtr = index->ionIndex[chunk_number].bA;

            /* Get size of the chunk */
            uint_t csize = ((chunk_number == index->nChunks - 1) && (index->nChunks > 1)) ?
                           index->lastchunksize                                 :
                           index->chunksize;

            uint_t count = bAPtr[0];

            /* Initialize first and last bA entries */
            bAPtr[0] = 0;

            for (uint_t li = 1; li <= (uint_t)(maxmass * scale); li++)
            {
                uint_t tmpcount = bAPtr[li];
                bAPtr[li] = bAPtr[li - 1] + count;
                count = tmpcount;
            }

            /* Check if all correctly done */
            if (bAPtr[(int_t)(maxmass * scale)] != (csize * speclen))
            {	
                status = ERR_INVLD_SIZE;	
            }

            assert (bAPtr[(int_t)(maxmass * scale)] == (csize * speclen));
        }
    }


    /* Optimize the CFIR index chunks */
    if (status == SLM_SUCCESS)
    {
        for (uint_t chunk_number = 0; chunk_number < index->nChunks; chunk_number++)
            status = DSLIM_Optimize(index, chunk_number);

    }

    return status;
}

/*
 * FUNCTION: DSLIM_AllocateMemory
 *
 * DESCRIPTION: Allocate memory for DSLIM chunks
 *
 * INPUT:
 * @chsize: Chunk size
 * @Chunks: Number of chunks
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_AllocateMemory(Index *index)
{
    status_t status = SLM_SUCCESS;
    uint_t chsize = index->chunksize;
    uint_t Chunks = index->nChunks;

    uint_t maxmass = params.max_mass;
    uint_t maxz = params.maxz;
    uint_t scale = params.scale;

    uint_t speclen = ((index->pepIndex.peplen-1) * iSERIES * maxz);

    /* Initialize DSLIM pepChunks */
    index->ionIndex = new SLMchunk[Chunks];

    if (index->ionIndex != NULL)
    {

        /* Counter for Chunk size construction */
        int_t totalpeps = (int_t) index->lcltotCnt;

        /* At every loop, check for remaining peps and status */
        for (uint_t i = 0; i < Chunks && totalpeps > 0 && status == SLM_SUCCESS; i++)
        {
            /* Initialize direct hashing bA */
            index->ionIndex[i].bA = new uint_t[(maxmass * scale) + 1];

            if (index->ionIndex[i].bA != NULL)
            {
                /* Calculate the iA chunk size */
                int_t size = ((int_t)(totalpeps - chsize)) > 0 ? chsize : totalpeps;
                totalpeps -= size; // Update the counter

                /* Total Number of Ions = peps * #ion series * ions/ion series */
                index->ionIndex[i].iA = new uint_t[(size * speclen)];

                if (index->ionIndex[i].iA == NULL)
                    status = ERR_INVLD_MEMORY;

            }
            else
            {
                status = ERR_INVLD_MEMORY;
            }
        }
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }

    return status;
}

/*
 * FUNCTION: DSLIM_ConstructChunk
 *
 * INPUT:
 * @threads:      Number of parallel threads
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_ConstructChunk(uint_t threads, Index *index, uint_t chunk_number)
{
    status_t status = SLM_SUCCESS;
    const uint_t peplen_1 = index->pepIndex.peplen - 1;
    const uint_t peplen   = index->pepIndex.peplen;
    const uint_t speclen  = params.maxz * iSERIES * peplen_1;

    const double_t minmass = params.min_mass;
    const double_t maxmass = params.max_mass;
    const uint_t scale = params.scale;

    /* Check if this chunk is the last chunk */
    BOOL lastChunk = (chunk_number == (index->nChunks - 1))? true: false;

    uint_t *BAPtrs[threads];
    uint_t *bA = new uint_t[(int_t)(threads * scale * maxmass)];

#ifdef USE_OMP

    for (uint_t i = 0; i < threads; i++)
        BAPtrs[i] = bA + (int_t)(i * scale * maxmass);

        /* Clear the bAPtrs and SAPtrs */
#pragma omp parallel for num_threads(threads) schedule (static)
    for (uint_t i = 0; i < threads; i++)
        std::memset(BAPtrs[i], 0x0, (scale * maxmass * sizeof(uint_t)));
#else

    std::memset(index->ionIndex[chunk_number].bA, 0x0, (((params.scale * params.max_mass) + 1) * sizeof(uint_t)));

#endif /* USE_OMP */

    if (status == SLM_SUCCESS)
    {
        uint_t start_idx = chunk_number * index->chunksize;
        uint_t interval = index->chunksize;

        /* Check for last chunk */
        if (lastChunk == true && index->nChunks > 1)
            interval = index->lastchunksize;

#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif /* USE_OMP */
        for (uint_t k = start_idx; k < (start_idx + interval); k++)
        {
            /* Filling point */
            uint_t nfilled = (k - start_idx) * speclen;

            /* Temporary Array needed for Theoretical Spectra */
            uint_t* Spectrum = new uint_t[speclen];

#ifdef USE_OMP
            uint_t *bAPtr  = BAPtrs[omp_get_thread_num()];
#else
            uint_t *bAPtr = index->ionIndex[chunk_number].bA;
#endif /* OPENMP */

            /* Extract peptide Information */
            float_t pepMass = 0.0;
            char_t *seq = NULL;
            uint_t pepID = k;

            /* Extract from pepEntries */
            pepEntry *entry = index->pepEntries + pepID;

            seq = &index->pepIndex.seqs[entry->seqID * peplen];

            /* Check if pepID belongs to peps or mods */
            if (entry->sites.modNum == 0)
            {
                /* Generate the Theoretical Spectrum */
                pepMass = UTILS_GenerateSpectrum(seq, peplen, Spectrum);
            }
            else
            {
                /* Generate the Mod. Theoretical Spectrum */
                pepMass = UTILS_GenerateModSpectrum(seq, (uint_t) peplen, Spectrum, entry->sites);
            }

            /* If a legal peptide */
            if (pepMass >= minmass && pepMass <= maxmass)
            {
                /* Fill the ions */
                for (uint_t ion = 0; ion < speclen; ion++)
                {
                    /* Check if legal ion */
                    if (Spectrum[ion] >= (maxmass * scale))
                    {
                        Spectrum[ion] = (maxmass * scale) - 1;
                        //assert (Spectrum[ion] < (maxmass * scale));
                    }

                    SpecArr[nfilled + ion] = Spectrum[ion]; // Fill in the ion
                    bAPtr[Spectrum[ion]]++; // Update the bA counter
                }
            }

            /* Illegal peptide, fill in the container with zeros */
            else
            {
                /* Fill zeros for illegal peptides
                 * FIXME: Should not be filled into the chunk
                 *  and be removed from peptide index as well
                 */
                std::memset(&SpecArr[nfilled], 0x0, sizeof(uint_t) * speclen);
                bAPtr[0] += speclen; // Update the BA counter
            }
        }

#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule (static)
        for (uint_t i = 0; i < (uint_t)(maxmass * scale); i++)
        {
            /* Initialize to zero */
            index->ionIndex[chunk_number].bA[i] = 0;

            for (uint_t j = 0; j < threads; j++)
                index->ionIndex[chunk_number].bA[i] += BAPtrs[j][i];

        }
#endif /* USE_OMP */

    }

#ifdef USE_OMP
    delete[] bA;
    bA = NULL;
#endif /* OPENMP */

    return status;
}

/*
 * FUNCTION: DSLIM_SLMTransform
 *
 * DESCRIPTION: Constructs SLIM Transform
 *
 * INPUT:
 * @threads     : Number of parallel threads
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_SLMTransform(uint_t threads, Index *index, uint_t chunk_number)
{
    status_t status = SLM_SUCCESS;

    /* Check if this chunk is the last chunk */
    uint_t size = ((chunk_number == index->nChunks - 1) && (index->nChunks > 1))?
                   index->lastchunksize : index->chunksize;

    uint_t speclen = (index->pepIndex.peplen - 1) * params.maxz * iSERIES;
    uint_t *iAPtr = index->ionIndex[chunk_number].iA;
    uint_t iAsize = size * speclen;

    /* Construct DSLIM.iA */
#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* USE_OMP */
    for (uint_t k = 0; k < iAsize; k++)
        iAPtr[k] = k;

#ifdef USE_OMP
    /* Parallel Key Value Sort */
    KeyVal_Parallel<uint_t>(SpecArr, iAPtr, iAsize, threads);
#else
    KeyVal_Serial<uint_t>(SpecArr, iAPtr, iAsize);
#endif /* USE_OMP */

    return status;
}

status_t DSLIM_InitializeScorecard(Index *index, uint_t idxs)
{
    status_t status = SLM_SUCCESS;

    /* Get size of the chunk */
    uint_t sAize = params.spadmem / (BYISIZE * params.threads);
    uint_t sz2 = 0;

    for (uint_t ii = 0; ii < idxs; ii++)
    {
        if (index[ii].chunksize > sz2)
            sz2 = index[ii].chunksize;

    }

    sAize = std::min(sz2, sAize);
    Score = new BYICount[params.threads];

    if (Score != NULL)
    {
#ifdef USE_OMP
#pragma omp parallel for schedule(static, 1) num_threads(params.threads)
#endif /* USE_OMP */
        for (uint_t thd = 0; thd < params.threads; thd++)
        {
            Score[thd].byc = new BYC[sAize];
            memset(Score[thd].byc, 0x0, sizeof(BYC) * sAize);

            /* Initialize the histogram */
            Score[thd].res.survival = new double_t[1 + (MAX_HYPERSCORE * 10) + 1]; // +2 for accumulation

            /* Evaluate to the nearest power of 2 */
            uint_t num = (params.topmatches == 0)? 1: params.topmatches;
            num = (uint_t)(floor(log2(num) + 0.999));

            std::memset(Score[thd].res.survival, 0x0, sizeof (double_t) * (2 + MAX_HYPERSCORE * 10));

            /* Initialize the heap with size = pow(2, num) - 1 */
            Score[thd].res.topK.init((1 << num) - 1);
        }
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }

    return status;
}

/*
 * FUNCTION: DSLIM_Optimize
 *
 * DESCRIPTION: Optimization of the index
 *
 * INPUT:
 * @threads     : Number of parallel threads
 * @index       : The SLM Index
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_Optimize(Index *index, uint_t chunk_number)
{
    status_t status = SLM_SUCCESS;

    uint_t *iAPtr = index->ionIndex[chunk_number].iA;
    uint_t *bAPtr = index->ionIndex[chunk_number].bA;

#ifdef USE_OMP
    uint_t interval = 50 *params.scale;
    uint_t threads = params.threads;
#endif /* USE_OMP */

#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, interval)
#endif /* USE_OMP */
    /* Stablize the entries in iAPtr */
    for (uint_t kk = 0; kk < (params.max_mass * params.scale); kk++)
    {
        uint_t offset = bAPtr[kk];
        uint_t size = bAPtr[kk + 1] - offset;

        if (size > 1)
        {
            /* Stablize the KeyVal Sort */
            UTILS_Sort<uint_t>((iAPtr + offset), size, false);
        }

    }

    return status;
}

/*
 * FUNCTION: DSLIM_Analyze
 *
 * DESCRIPTION: Analyze the DSLIM distribution
 *
 * INPUT:
 * @threads: Number of parallel threads
 * @avg    : Pointer to DSLIM mean load
 * @std    : Pointer to DSLIM load distribution
 *
 * OUTPUT:
 * @status: Status of execution
 */
#if 0
status_t DSLIM_Analyze(uint_t threads, double_t &avg, double_t &std)
{
    status_t status = SLM_SUCCESS;
    double_t sum_std = 0;
    ull_t truecount = 1;

#ifdef USE_OMP
    double_t stds[threads] = {};
    ull_t count[threads] = {};
#endif /* USE_OMP */

    if (nchunks <= 1)
        status = ERR_INVLD_SIZE;

    if (status == SLM_SUCCESS)
    {
        double_t *arr = new double_t[(MAX_MASS * SCALE)+1];

        /* Trivial to include zeros or last value */
        arr[(MAX_MASS * SCALE)] = 0x0;
        arr[0] = 0x0;

        /* Compute the means */
#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* USE_OMP */
        for (uint_t i = 1; i < (MAX_MASS * SCALE) -1; i++)
        {
            arr[i] = 0x0;

            for (uint_t j = 0; j < nchunks; j++)
                arr[i] += (dslim.pepChunks[j].bA[i+1] - dslim.pepChunks[j].bA[i]);

            arr[i] /= nchunks;

            if (arr[i] > 0)
            {
#ifdef USE_OMP
                count[omp_get_thread_num()]++;
#else
                truecount++;
#endif /* USE_OMP */
            }
        }

        /* Gather the counts to truecount */
#ifdef USE_OMP
        for (uint_t kk = 0; kk < threads; kk++)
            truecount += count[kk];

#endif /* OPENMP */

        /* Compute the Standard Deviations */
#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* USE_OMP */
        for (uint_t i = 1; i < (MAX_MASS * SCALE) - 1; i++)
        {
            double_t mean = arr[i];
            arr[i] = 0x0;

            if (mean > 0)
            {
                for (uint_t j = 0; j < nchunks; j++)
                {
                    arr[i] += (((double_t) dslim.pepChunks[j].bA[i+1]
                            - (double_t) dslim.pepChunks[j].bA[i] - mean)
                            * ((double_t) dslim.pepChunks[j].bA[i+1]
                            - (double_t) dslim.pepChunks[j].bA[i] - mean));
                }

                arr[i] = sqrt(arr[i] / nchunks);
            }
        }

        /* Compute mean of stdevs */
#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* USE_OMP */
        for (uint_t i = 0; i < MAX_MASS * SCALE; i++)
        {
#ifdef USE_OMP
            stds[omp_get_thread_num()] += arr[i];
#else
            sum_std += arr[i];
#endif /* USE_OMP */
        }

/* Gather the counts to truecount */
#ifdef USE_OMP
        for (uint_t kk = 0; kk < threads; kk++)
        {
            sum_std += stds[kk];
            stds[kk] = 0;
        }
#endif /* OPENMP */

        sum_std /= truecount;
        avg = sum_std;

        sum_std = 0;

        /* Compute stdev of stdevs */
#ifdef USE_OMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* USE_OMP */
        for (uint_t i = 0; i < MAX_MASS * SCALE; i++)
        {
#ifdef USE_OMP
            stds[omp_get_thread_num()] += ((arr[i] - avg) * (arr[i] - avg));
#else
            sum_std += ((arr[i] - avg) * (arr[i] - avg));
#endif /* USE_OMP */
        }

        /* Gather the counts to truecount */
#ifdef USE_OMP
        for (uint_t kk = 0; kk < threads; kk++)
            sum_std += stds[kk];

#endif /* OPENMP */

        sum_std /= truecount;
        sum_std = std::sqrt(sum_std);

        std = sum_std;

        delete[] arr;
        arr = NULL;
    }

    return status;
}

#endif

/*
 * FUNCTION: DSLIM_Deinitialize
 *
 * DESCRIPTION: Deallocate DSLIM memory
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */
status_t DSLIM_Deinitialize(Index *index)
{
    status_t status = SLM_SUCCESS;

    /* Deallocate the Ion Index */
    status = DSLIM_DeallocateIonIndex(index);

    /* Deallocate the peptide index */
    status = DSLIM_DeallocatePepIndex(index);

    return status;
}

status_t DSLIM_DeallocateIonIndex(Index *index)
{
    /* Deallocate all the DSLIM chunks */
    for (uint_t chno = 0; chno < index->nChunks; chno++)
    {
        SLMchunk curr_chunk = index->ionIndex[chno];

        if (curr_chunk.bA != NULL)
        {
            delete[] curr_chunk.bA;
            curr_chunk.bA = NULL;
        }

        if (curr_chunk.iA != NULL)
        {
            delete[] curr_chunk.iA;
            curr_chunk.iA = NULL;
        }
    }

    if (index->ionIndex != NULL)
    {
        delete[] index->ionIndex;
        index->ionIndex = NULL;
    }

    /* Set the chunking to zero */
    index->nChunks = 0;
    index->chunksize = 0;
    index->lastchunksize = 0;

    /* Nothing really to return
     * so success */
    return SLM_SUCCESS;
}

status_t DSLIM_DeallocatePepIndex(Index *index)
{
    if (index->pepEntries != NULL)
    {
        delete[] index->pepEntries;
        index->pepEntries = NULL;
    }

    if (index->pepIndex.seqs != NULL)
    {
        delete[] index->pepIndex.seqs;
        index->pepIndex.seqs = NULL;
    }

    /* Reset peptide Index Variables */
    index->pepCount = 0;
    index->modCount = 0;
    index->totalCount = 0;

    index->pepIndex.AAs = 0;
    index->pepIndex.peplen = 0;

    /* Nothing really to return
     * so success */
    return SLM_SUCCESS;
}

status_t DSLIM_DeallocateSpecArr()
{
    if (SpecArr != NULL)
    {
        delete[] SpecArr;
        SpecArr= NULL;
    }

    /* Nothing really to return
     * so success */
    return SLM_SUCCESS;
}

/*
 * FUNCTION: DSLIM_DeallocateSC
 *
 * DESCRIPTION:
 *
 * INPUT:
 * none
 *
 * OUTPUT:
 * @status: status of execution
 */
status_t DSLIM_DeallocateSC()
{
    /* Free the Scorecard memory */
    if (Score != NULL)
    {
        for (uint_t thd = 0; thd < params.threads; thd++)
        {
            delete[] Score[thd].byc;
            delete[] Score[thd].res.survival;

            Score[thd].byc = NULL;
            Score[thd].res.survival = NULL;
        }

        delete[] Score;
        Score = NULL;

    }

    return SLM_SUCCESS;
}

int_t DSLIM_GenerateIndex(Index *index, uint_t key)
{
    int_t value = -1;

    DistPolicy policy = params.policy;

    uint_t csize = index->lclpepCnt;

    if (policy == _chunk)
        value = (params.myid * csize) + key;
    else if (policy == _cyclic)
        value = (key * params.nodes) + params.myid;
    else if (policy == _zigzag)
        value = false;

    return value;
}
