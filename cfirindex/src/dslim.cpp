/*
 * This file is part of SLM-Transform
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
using namespace std;

/* Global Variables */
UINT          *SpecArr = NULL; /* Spectra Array     */
UINT reduce = 0;


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
STATUS DSLIM_Construct(UINT threads, Index *index)
{
    STATUS status = SLM_SUCCESS;
    UINT maxz = params.maxz;
    UINT peplen_1 = index->pepIndex.peplen - 1;

    DOUBLE maxmass = params.max_mass;
    UINT scale = params.scale;

#ifdef VMODS
#else
    LBE_UNUSED_PARAM(modInfo);
#endif /* VMODS */

#ifndef _OPENMP
    threads = 1;
#endif /* _OPENMP */

    if (status == SLM_SUCCESS && SpecArr == NULL)
    {
        /* Spectra Array (SA) */
        SpecArr = new UINT[MAX_IONS];

        /* Check if Spectra Array has been allocated */
        if (SpecArr == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
    }

    if (status == SLM_SUCCESS)
    {
        /* Allocate memory for SLMChunks and SPI*/
        status = DSLIM_AllocateMemory(index);

        /* Construct DSLIM.iA */
        if (status == SLM_SUCCESS)
        {
            /* Distributed SLM Ions Array construction */
            for (UINT chno = 0; chno < index->nChunks && status == SLM_SUCCESS; chno++)
            {
                /* Construct each DSLIM chunk in Parallel */
                status = DSLIM_ConstructChunk(threads, index, chno);

                /* Apply SLM-Transform on the chunk */
                if (status == SLM_SUCCESS)
                {
                    status = DSLIM_SLMTransform(threads, index, chno);
                }
            }
        }
    }

    if (status == SLM_SUCCESS)
    {
        UINT speclen = peplen_1 * maxz * iSERIES;

        /* Construct DSLIM.bA */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif /* _OPENMP */
        for (UINT chunk_number = 0; chunk_number < index->nChunks; chunk_number++)
        {
            UINT *bAPtr = index->ionIndex[chunk_number].bA;

            /* Get size of the chunk */
            UINT csize = ((chunk_number == index->nChunks - 1) && (index->nChunks > 1)) ?
                           index->lastchunksize                                 :
                           index->chunksize;

            UINT count = bAPtr[0];

            /* Initialize first and last bA entries */
            bAPtr[0] = 0;
            bAPtr[(INT)(maxmass * scale)] = (csize * speclen);

            for (UINT li = 1; li <= (UINT)(maxmass * scale); li++)
            {
                UINT tmpcount = bAPtr[li];
                bAPtr[li] = bAPtr[li - 1] + count;
                count = tmpcount;
#ifdef VALIDATE_SLM
                if (bAPtr[li] < bAPtr[li - 1])
                {
                    std::cout << chunk_number << " " << li << " " << bAPtr[li - 1] << " " << count << " " << bAPtr[li] << std::endl;
                    status = ERR_INVLD_SIZE;

                    while(true);
                }
#endif /* VALIDATE_SLM */
            }

            /* Check if all correctly done */
            if (bAPtr[(INT)(maxmass * scale)] != (csize * speclen))
            {
                status = ERR_INVLD_SIZE;
            }
        }
    }

    /* Optimize DSLIM chunks */
    if (status == SLM_SUCCESS)
    {
        for (UINT chno = 0; chno < index->nChunks; chno++)
        {
            status = DSLIM_Optimize(threads, index, chno);
        }
    }

    /* Remove the temporary SpecArray (SA) */
//    delete[] SpecArr;
//    SpecArr = NULL;

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
STATUS DSLIM_AllocateMemory(Index *index)
{
    STATUS status = SLM_SUCCESS;
    UINT chsize = index->chunksize;
    UINT Chunks = index->nChunks;

    UINT maxmass = params.max_mass;
    UINT maxz = params.maxz;
    UINT scale = params.scale;

    UINT speclen = ((index->pepIndex.peplen-1) * iSERIES * maxz);

    /* Initialize DSLIM pepChunks */
    index->ionIndex = new SLMchunk[Chunks];

    if (index->ionIndex != NULL)
    {

        /* Counter for Chunk size construction */
        INT totalpeps = (INT) index->totalCount;

        /* At every loop, check for remaining peps and status */
        for (UINT i = 0; i < Chunks && totalpeps > 0 && status == SLM_SUCCESS; i++)
        {
            /* Initialize direct hashing bA */
            index->ionIndex[i].bA = new UINT[(maxmass * scale) + 1];

            if (index->ionIndex[i].bA != NULL)
            {
                /* Calculate the iA chunk size */
                INT size = ((INT)(totalpeps - chsize)) > 0 ? chsize : totalpeps;
                totalpeps -= size; // Update the counter

                /* Total Number of Ions = peps * #ion series * ions/ion series */
                index->ionIndex[i].iA = new UINT[(size * speclen)];

                if (index->ionIndex[i].iA == NULL)
                {
                    status = ERR_INVLD_MEMORY;
                }
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

    /* Allocate memory for SPI and Spectra Array */
    if (status == SLM_SUCCESS)
    {
        index->pepEntries = new pepEntry[index->pepCount];

        if (index->pepEntries == NULL)
        {
            status = ERR_INVLD_MEMORY;
        }
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
STATUS DSLIM_ConstructChunk(UINT threads, Index *index, UINT chunk_number)
{
    STATUS status = SLM_SUCCESS;
    UINT peplen_1 = index->pepIndex.peplen - 1;
    UINT peplen   = index->pepIndex.peplen;
    UINT speclen  = params.maxz * iSERIES * peplen_1;

    DOUBLE minmass = params.min_mass;
    DOUBLE maxmass = params.max_mass;
    UINT scale = params.scale;

    /* Check if this chunk is the last chunk */
    BOOL lastChunk = (chunk_number == (index->nChunks - 1))? true: false;

#ifdef _OPENMP
    /* SA ptr for each thread */
    UINT *SAPtrs[threads] = {NULL};

    /* Temporary Array needed for Theoretical Spectra */
    UINT *Spectra = new UINT[threads * speclen];

    UINT *BAPtrs[threads] = {NULL};
    UINT *bA = new UINT[(INT)(threads * scale * maxmass)];

    /* Initialize SAPtrs for each thread */
    if (Spectra != NULL && bA != NULL)
    {
        for (UINT i = 0; i < threads; i++)
        {
            SAPtrs[i] = Spectra + (i * speclen);
            BAPtrs[i] = bA + (INT)(i * scale * maxmass);
        }

        /* Clear the bAPtrs and SAPtrs */
#pragma omp parallel for num_threads(threads) schedule (static)
        for (UINT i = 0; i < threads; i++)
        {
            std::memset(BAPtrs[i], 0x0, (scale * maxmass * sizeof(UINT)));
            std::memset(SAPtrs[i], 0x0, (speclen * sizeof(UINT)));
        }
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }
#else

    UINT SAPtr[speclen] = {};
    std::memset(index->ionIndex[chunk_number].bA, 0x0, (((SCALE * MAX_MASS) + 1) * sizeof(UINT)));

#endif /* _OPENMP */

    if (status == SLM_SUCCESS)
    {
        UINT start_idx = chunk_number * index->chunksize;
        UINT interval = index->chunksize;

        /* Check for last chunk */
        if (lastChunk == true && index->nChunks > 1)
        {
            interval = index->lastchunksize;
        }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1) /* schedule(dynamic, 1) */
#endif /* _OPENMP */
        for (UINT k = start_idx; k < (start_idx + interval); k++)
        {
            /* Filling point */
            UINT nfilled = (k - start_idx) * speclen;

#ifdef _OPENMP
            UINT *Spec = SAPtrs[omp_get_thread_num()];
            UINT *bAPtr  = BAPtrs[omp_get_thread_num()];
#else
            UINT *Spec = SAPtr;
            UINT *bAPtr = index->ionIndex[chunk_number].bA;
#endif /* OPENMP */

            /* Extract peptide Information */
            FLOAT pepMass = 0.0;
            CHAR *seq = NULL;
            UINT pepID = LBE_RevDist(k);

#ifdef VMODS
            /* Check if pepID belongs to peps or mods */
            if (pepID >= index->pepCount)
            {
                /* Extract from Mods */
                varEntry *entry = index->modEntries + (pepID - index->pepCount);
                seq = &index->pepIndex.seqs[entry->seqID * peplen];

                /* Generate the Mod. Theoretical Spectrum */
                pepMass = UTILS_GenerateModSpectrum(seq, (UINT)peplen, Spec, entry->sites);

                /* Fill in the pepMass */
                entry->Mass = pepMass;
            }
            else
#endif /* VMODS */
            {
                /* Extract from Peps */
                pepEntry *entry = index->pepEntries + pepID;
                seq = &index->pepIndex.seqs[peplen * pepID];

                /* Generate the Theoretical Spectrum */
                pepMass = UTILS_GenerateSpectrum(seq, peplen, Spec);

                /* Fill in the pepMass */
                entry->Mass = pepMass;
            }

            /* If a legal peptide */
            if (pepMass >= minmass && pepMass <= maxmass)
            {
                /* Sort by ion Series and Mass */
//                UTILS_Sort<UINT>(Spec, half_len, false);
//                UTILS_Sort<UINT>((Spec + half_len), half_len, false);

                /* Fill the ions */
                for (UINT ion = 0; ion < speclen; ion++)
                {
                    /* Check if legal b-ion */
                    if (Spec[ion] >= (maxmass * scale))
                    {
                        Spec[ion] = (maxmass * scale) - 1;
                    }

                    SpecArr[nfilled + ion] = Spec[ion]; // Fill in the b-ion
                    bAPtr[SpecArr[nfilled + ion]]++; // Update the BA counter
                }
            }

            /* Illegal peptide, fill in the container with zeros */
            else
            {
                /* Fill zeros for illegal peptides
                 * FIXME: Should not be filled into the chunk
                 *        and be removed from SPI as well
                 */
                std::memset(&SpecArr[nfilled], 0x0, sizeof(UINT) * speclen);
                bAPtr[0] += speclen; // Update the BA counter
            }
        }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule (static)
        for (UINT i = 0; i < (UINT)(maxmass * scale); i++)
        {
            /* Initialize to zero */
            index->ionIndex[chunk_number].bA[i] = 0;

            for (UINT j = 0; j < threads; j++)
            {
                index->ionIndex[chunk_number].bA[i] += BAPtrs[j][i];
            }
        }
#endif /* _OPENMP */

    }

#ifdef _OPENMP
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
STATUS DSLIM_SLMTransform(UINT threads, Index *index, UINT chunk_number)
{
    STATUS status = SLM_SUCCESS;

    /* Check if this chunk is the last chunk */
    UINT size = ((chunk_number == index->nChunks - 1) && (index->nChunks > 1))?
                 index->lastchunksize                                  :
                 index->chunksize;

    UINT speclen = (index->pepIndex.peplen - 1) * params.maxz * iSERIES;
    UINT *iAPtr = index->ionIndex[chunk_number].iA;
    UINT iAsize = size * speclen;

#ifdef VALIDATE_SLM
    UINT *integ = new UINT[size * speclen];
    std::memcpy(integ, SpecArr, sizeof(UINT) * iAsize);
#endif /* VALIDATE_SLM */

    /* Construct DSLIM.iA */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
    for (UINT k = 0; k < iAsize; k++)
    {
        iAPtr[k] = k;
    }

#ifdef _OPENMP
    /* Parallel Key Value Sort */
    KeyVal_Parallel<UINT>(SpecArr, iAPtr, iAsize, threads);
#else
    KeyVal_Serial<UINT>(SpecArr, iAPtr, iAsize);
#endif /* _OPENMP */

    /* Check integrity of SLM-Transform */
#ifdef VALIDATE_SLM
    BOOL integrity = true;

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
    for (UINT k = 1; k < iAsize; k++)
    {
        if (integ[iAPtr[k]] < integ[iAPtr[k - 1]])
        {
            integrity = false;
        }
    }

    if (integrity == false)
    {
        while (!integrity);
        status = ERR_INVLD_SORT;
    }
#endif /* VALIDATE_SLM */

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
STATUS DSLIM_Optimize(UINT threads, Index *index, UINT chunk_number)
{
    STATUS status = SLM_SUCCESS;

    UINT *bAPtr = index->ionIndex[chunk_number].bA;
    UINT *iAPtr = index->ionIndex[chunk_number].iA;

#ifdef OPTIMIZE
    UINT *chPtr = new UINT[(MAX_MASS * SCALE)+1];
    std::memset(chPtr, 0x0, (sizeof(UINT) * (MAX_MASS * SCALE)+1));
    UINT speclen = MAXz * (index->pepIndex.peplen - 1) * iSERIES;
#endif /* OPTIMIZE */

    /* Stablize the entries in iAPtr */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 50 * params.scale)
#endif /* _OPENMP */
    for (UINT kk = 0; kk < (UINT)(params.max_mass * params.scale); kk++)
    {
        UINT offset = bAPtr[kk];
        UINT size = bAPtr[kk + 1] - offset;

        if (size > 1)
        {
            /* Stablize the KeyVal Sort */
            UTILS_Sort<UINT>((iAPtr + offset), size, false);
        }

    }

#ifdef OPTIMIZE
    /* For each ionEntries, do the following, check for groups >= 10,
     * Group = ion1, ion1 + peplen, ion1 + 2peplen, ion1+ 3peplen so on
     * Replace group by CELLBLOCK, cell header, size
     * Update the bAPtr size
     */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 50 * SCALE)
#endif /* _OPENMP */
    for (UINT i = 0; i < (MAX_MASS * SCALE); i++)
    {
        UINT start = bAPtr[i];
        UINT end =   bAPtr[i+1];

        UINT nions = end - start;

        /* If less than 500 occurrences, then continue */
        if (nions < 500)
        {
            continue;
        }

        /* Otherwise*/
        UINT fion = start;
        UINT group = 1;

        for (UINT kk = start + 1; kk < end; kk++)
        {

            if (iAPtr[kk] < (iAPtr[fion] + (group * speclen)))
            {
                iAPtr[kk] = -1;
                reduce += 1;
            }
            else if (iAPtr[kk] == (iAPtr[fion] + (group * speclen)))
            {
                iAPtr[kk] = -1;
                group++;
            }
            else
            {
                if (group >= 10)
                {
                    iAPtr[fion + 1] = iAPtr[fion];
                    iAPtr[fion+2] = group;
                    iAPtr[fion] = CELLBLOCK;

                    chPtr[i+1] += group - 3;
                    reduce += group -2;
                }

                fion = kk;
                group = 1;

            }
        }
    }

    delete[] chPtr;
    chPtr = NULL;
    std::cout << "Compression Achieved: " << reduce * 4.0/(1024*1024) << "MiB" << std::endl;
#endif /* OPTIMIZE */

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
STATUS DSLIM_Analyze(UINT threads, DOUBLE &avg, DOUBLE &std)
{
    STATUS status = SLM_SUCCESS;
    DOUBLE sum_std = 0;
    ULONGLONG truecount = 1;

#ifdef _OPENMP
    DOUBLE stds[threads] = {};
    ULONGLONG count[threads] = {};
#endif /* _OPENMP */

    if (nchunks <= 1)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        DOUBLE *arr = new DOUBLE[(MAX_MASS * SCALE)+1];

        /* Trivial to include zeros or last value */
        arr[(MAX_MASS * SCALE)] = 0x0;
        arr[0] = 0x0;

        /* Compute the means */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT i = 1; i < (MAX_MASS * SCALE) -1; i++)
        {
            arr[i] = 0x0;

            for (UINT j = 0; j < nchunks; j++)
            {
                arr[i] += (dslim.pepChunks[j].bA[i+1] - dslim.pepChunks[j].bA[i]);
            }

            arr[i] /= nchunks;

            if (arr[i] > 0)
            {
#ifdef _OPENMP
                count[omp_get_thread_num()]++;
#else
                truecount++;
#endif /* _OPENMP */
            }
        }

        /* Gather the counts to truecount */
#ifdef _OPENMP
        for (UINT kk = 0; kk < threads; kk++)
        {
            truecount += count[kk];
        }
#endif /* OPENMP */

        /* Compute the Standard Deviations */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT i = 1; i < (MAX_MASS * SCALE) - 1; i++)
        {
            DOUBLE mean = arr[i];
            arr[i] = 0x0;

            if (mean > 0)
            {
                for (UINT j = 0; j < nchunks; j++)
                {
                    arr[i] += (((DOUBLE) dslim.pepChunks[j].bA[i+1]
                            - (DOUBLE) dslim.pepChunks[j].bA[i] - mean)
                            * ((DOUBLE) dslim.pepChunks[j].bA[i+1]
                            - (DOUBLE) dslim.pepChunks[j].bA[i] - mean));
                }

                arr[i] = sqrt(arr[i] / nchunks);
            }
        }

        /* Compute mean of stdevs */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT i = 0; i < MAX_MASS * SCALE; i++)
        {
#ifdef _OPENMP
            stds[omp_get_thread_num()] += arr[i];
#else
            sum_std += arr[i];
#endif /* _OPENMP */
        }

/* Gather the counts to truecount */
#ifdef _OPENMP
        for (UINT kk = 0; kk < threads; kk++)
        {
            sum_std += stds[kk];
            stds[kk] = 0;
        }
#endif /* OPENMP */

        sum_std /= truecount;
        avg = sum_std;

        sum_std = 0;

        /* Compute stdev of stdevs */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT i = 0; i < MAX_MASS * SCALE; i++)
        {
#ifdef _OPENMP
            stds[omp_get_thread_num()] += ((arr[i] - avg) * (arr[i] - avg));
#else
            sum_std += ((arr[i] - avg) * (arr[i] - avg));
#endif /* _OPENMP */
        }

        /* Gather the counts to truecount */
#ifdef _OPENMP
        for (UINT kk = 0; kk < threads; kk++)
        {
            sum_std += stds[kk];
        }
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
STATUS DSLIM_Deinitialize(Index *index)
{
    STATUS status = SLM_SUCCESS;

    status = DSLIM_DeallocateSC();

    /* Deallocate all the DSLIM chunks */
    for (UINT chno = 0; chno < index->nChunks; chno++)
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

    if (index->pepEntries != NULL)
    {
        delete[] index->pepEntries;
        index->pepEntries = NULL;
    }

#ifdef VMODS
    if (index->modEntries != NULL)
    {
        delete[] index->modEntries;
        index->modEntries = NULL;
    }
#endif /* VMODS */

    if (index->pepIndex.seqs != NULL)
    {
        delete[] index->pepIndex.seqs;
        index->pepIndex.seqs = NULL;
    }

    /* Reset Index Variables */
    index->nChunks = 0;
    index->pepCount = 0;
    index->modCount = 0;
    index->totalCount = 0;
    index->chunksize = 0;
    index->lastchunksize = 0;

    index->pepIndex.AAs = 0;
    index->pepIndex.peplen = 0;

    return status;
}

STATUS DSLIM_DeallocateSpecArr(VOID)
{
    if (SpecArr != NULL)
    {
        delete[] SpecArr;
        SpecArr= NULL;
    }

    return SLM_SUCCESS;
}
#if 0
/* FUNCTION: DSLIM_WriteLIBSVM
 *
 * DESCRIPTION: Write the MS/MS spectra data in libsvm format
 *
 * INPUT:
 * @path   : Path for output data
 * @chno   : Chunk number
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_WriteLIBSVM(STRING path, UINT peplen, UINT chno)
{
    STATUS status = SLM_SUCCESS;

#ifdef GENDATA
    ofstream file;
#endif

#ifdef GENFEATS
    ofstream ffile;
#endif

#ifdef GENBA
    ofstream bAfile;
#endif

    UINT peplen_1 = peplen - 1;

#ifdef GENDATA
    /* Create the filename string */
    STRING filename = path + std::to_string(peplen) + "_" + std::to_string(chno) + ".libsvm";
#endif

#ifdef GENBA
    STRING bAfname  = path + "bA" + std::to_string(peplen) + "_" + std::to_string(chno) + ".csv";
#endif

#ifdef GENFEATS
    STRING ffname = path + "fA" + std::to_string(peplen) + "_" + std::to_string(chno) + ".csv";
#endif

    /* Check if this chunk is the last chunk */
    BOOL lastChunk = (chno == (index->nChunks - 1))? true: false;

#ifdef GENDATA
    /* Do not apply the append flag in here */
    file.open((const CHAR *) filename.c_str(), std::ofstream::out);
#endif

#ifdef GENBA
    bAfile.open((const CHAR *) bAfname.c_str(), std::ofstream::out);
#endif

#ifdef GENFEATS
    ffile.open((const CHAR *) ffname.c_str(), std::ofstream::out);
#endif

    if (status == SLM_SUCCESS)
    {
        UINT start = chno * chunksize;
        UINT interval = chunksize;

        /* Check for last chunk */
        if (lastChunk == true && nchunks > 1)
        {
            interval = lastchunksize;
        }

#ifdef GENDATA
        for (UINT k = start; k < (start + interval); k++)
        {
            /* Filling point */
            UINT nfilled = (k - start) * (iSERIES * peplen_1 * MAXz);
            UINT fsize = (iSERIES * peplen_1 * MAXz);
            UINT bl = 0, yl = 0;
            UINT feat = 0;

            STRING line = std::to_string(k);

            for (UINT z = 0; z < MAXz; z++)
            {
                for (bl = 1; bl <= B_FEATURES; bl++)
                {
                    line += " " + std::to_string(feat) + ":"
                            + std::to_string(SpecArr[chno][nfilled + bl - 1]);
                    feat++;
                }
            }

            for (UINT z = 0; z < MAXz; z++)
            {
                for (yl = B_FEATURES + 1; yl <= (B_FEATURES + Y_FEATURES); yl++)
                {
                    line += " " + std::to_string(feat) + ":"
                         + std::to_string(SpecArr[chno][nfilled + (int) (fsize / 2)
                                          + yl - B_FEATURES - 1]);
                    feat++;
                }
            }

            line += "\n";

            file << line;
        }
#endif

#ifdef GENBA
        STRING line = std::to_string(0);

        line += "," + std::to_string(0);
        line += "\n";

        bAfile << line;

        for (UINT k = 1; k <= (MAX_MASS * SCALE) + 1; k++)
        {
            STRING line = std::to_string(k);

            line += "," + std::to_string(bAPtr[k]);
            line += "\n";

            bAfile << line;
        }
#endif

#ifdef GENFEATS
        for (UINT k = 1; k <= (iSERIES *peplen_1 * MAXz); k++)
        {
            UINT sum = 0;

            ffile << std::to_string(k) + ",";

            for (UINT ll = 0; ll < (MAX_MASS * SCALE); ll++)
            {
                sum += (CHAR)features[(k-1) * (MAX_MASS * SCALE) + ll];
            }

            ffile << std::to_string(sum) + "\n";
        }
#endif

    }



    return status;
}

/* FUNCTION: DSLIM_WriteLIBSVM
 *
 * DESCRIPTION: Write the MS/MS spectra data in libsvm format
 *
 * INPUT:
 * @path   : Path for output data
 * @chno   : Chunk number
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_WriteCSV(STRING path, UINT peplen, UINT chno)
{
    STATUS status = SLM_SUCCESS;

#ifdef GENDATA
    ofstream file;
#endif

#ifdef GENFEATS
    ofstream ffile;
#endif

#ifdef GENBA
    ofstream bAfile;
#endif

    UINT peplen_1 = peplen - 1;
    UINT *bAPtr = dslim.pepChunks[chno].bA;

#ifdef GENDATA
    /* Create the filename string */
    STRING filename = path + std::to_string(peplen) + "_" + std::to_string(chno) + ".csv";
#endif

#ifdef GENBA
    STRING bAfname  = path + "bA" + std::to_string(peplen) + "_" + std::to_string(chno) + ".csv";
#endif

#ifdef GENFEATS
    STRING ffname = path + "fA" + std::to_string(peplen) + "_" + std::to_string(chno) + ".csv";
#endif
    /* Check if this chunk is the last chunk */
    BOOL lastChunk = (chno == (nchunks - 1))? true: false;

#ifdef GENDATA
    /* Do not apply the append flag in here */
    file.open((const CHAR *) filename.c_str(), std::ofstream::out);
#endif

#ifdef GENBA
    bAfile.open((const CHAR *) bAfname.c_str(), std::ofstream::out);
#endif

#ifdef GENFEATS
    ffile.open((const CHAR *) ffname.c_str(), std::ofstream::out);
#endif

    if (status == SLM_SUCCESS)
    {
        UINT start = chno * chunksize;
        UINT interval = chunksize;

        /* Check for last chunk */
        if (lastChunk == true && nchunks > 1)
        {
            interval = lastchunksize;
        }

#ifdef GENDATA
        for (UINT k = start; k < (start + interval); k++)
        {
            /* Filling point */
            UINT nfilled = (k - start) * (iSERIES * peplen_1 * MAXz);
            UINT fsize = (iSERIES * peplen_1 * MAXz);
            UINT bl = 0, yl = 0;

            STRING line = std::to_string(SpecArr[chno][nfilled]);

            for (UINT z = 0; z < MAXz; z++)
            {
                for (bl = 1; bl < B_FEATURES; bl++)
                {
                    line += "," + std::to_string(SpecArr[chno][nfilled + bl]);
                }
            }

            for (UINT z = 0; z < MAXz; z++)
            {
                for (yl = B_FEATURES; yl < (B_FEATURES + Y_FEATURES); yl++)
                {
                    line += "," + std::to_string(SpecArr[chno][nfilled + (int) (fsize / 2)
                                                               + yl - B_FEATURES]);
                }
            }

            line += "\n";

            file << line;
        }
#endif

#ifdef GENBA
        STRING line = std::to_string(0);

        line += "," + std::to_string(0);
        line += "\n";

        bAfile << line;

        for (UINT k = 1; k <= (MAX_MASS * SCALE) + 1; k++)
        {
            STRING line = std::to_string(k);

            line += "," + std::to_string(bAPtr[k]);
            line += "\n";

            bAfile << line;
        }
#endif

#ifdef GENFEATS
        for (UINT k = 1; k <= (iSERIES *peplen_1 * MAXz); k++)
        {
            UINT sum = 0;

            ffile << std::to_string(k) + ",";

            for (UINT ll = 0; ll < (MAX_MASS * SCALE); ll++)
            {
                sum += (CHAR)features[(k-1)* (MAX_MASS * SCALE) + ll];
            }

            ffile << std::to_string(sum) + "\n";
        }
#endif
    }

    if (status == SLM_SUCCESS)
    {
#ifdef GENDATA
        file.close();
#endif

#ifdef GENBA
        bAfile.close();
#endif

#ifdef GENFEATS
        ffile.close();
#endif
    }

    return status;
}
#endif /* 0 */
