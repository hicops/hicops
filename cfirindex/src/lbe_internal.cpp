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

#include "lbe.h"
using namespace std;

vector<STRING> Seqs;
UINT cumusize = 0;
UINT *varCount = NULL;
ifstream file;

extern gParams params;

#ifdef BENCHMARK
static DOUBLE duration = 0;
extern DOUBLE compute;
extern DOUBLE fileio;
extern DOUBLE memory;
#endif /* BENCHMARK */

/* Static function Prototypes */
static STATUS LBE_AllocateMem(Index *index);
static BOOL CmpPepEntries(pepEntry e1, pepEntry e2);
/*
 * FUNCTION: LBE_AllocateMem
 *
 * DESCRIPTION: Allocate memory for Data Structures
 *
 * INPUT:
 * @N: Number of peptides
 * @M: Number of mods
 *
 * OUTPUT:
 * @status: Status of execution
 */
static STATUS LBE_AllocateMem(Index *index)
{
    STATUS status = SLM_SUCCESS;
    UINT M = index->lcltotCnt;

    index->pepEntries = NULL;

    /* Allocate Memory for seqPep */
    index->pepIndex.seqs = new AA[index->pepIndex.AAs];

    if (index->pepIndex.seqs == NULL)
    {
        status = ERR_BAD_MEM_ALLOC;
    }

#ifdef VMODS
    /* Allocate the seqMod */
    if (status == SLM_SUCCESS)
    {
        index->pepEntries = new pepEntry[M];

        if (index->pepEntries == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
    }
#endif /* VMODS */


    return status;
}

BOOL LBE_ApplyPolicy(Index *index,  BOOL pepmod, UINT key)
{
    BOOL value = false;

    DistPolicy policy = params.policy;

    UINT csize = index->lclmodCnt;

    if (pepmod == false)
    {
        csize = index->lclpepCnt;
    }

    if (policy == _chunk)
    {
        value = (key / csize) == (params.myid);
    }
    else if (policy == _cyclic)
    {
        value = key % (params.nodes) == params.myid;
    }
    else
    {
        cout << "This policy is not implemented yet\n";
        value = false;
    }

    return value;
}

/*
 * FUNCTION: LBE_Initialize
 *
 * DESCRIPTION: Initialize internal peptides
 *              database from FASTA file
 *
 * INPUT:
 * @threads      : Number of parallel threads
 * @filename     : Path to FASTA file
 * @modconditions: String with mod conditions
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_Initialize(Index *index)
{
    STATUS status = SLM_SUCCESS;
    UINT iCount = 1;
    STRING seq;
    STRING modconditions = params.modconditions;

#ifdef _OPENMP
    UINT threads = params.threads;
#endif /* _OPENMP */

    /* Check if ">" entries are > 0 */
    if (index->lclpepCnt > 0)
    {
#ifdef BENCHMARK
        duration = omp_get_wtime();
#endif
        status = LBE_AllocateMem(index);

#ifdef BENCHMARK
        memory += omp_get_wtime() - duration;
#endif
    }
    else
    {
        status = ERR_INVLD_PARAM;
    }

    /* If Seqs was successfully filled */
    if (Seqs.size() != 0 && status == SLM_SUCCESS)
    {
        UINT seqlen = Seqs.at(0).length();

#ifdef BENCHMARK
        duration = omp_get_wtime();
#endif

#ifdef DEBUG
        cout << seq << endl;
#endif /* DEBUG */

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule (static) reduction(+: iCount)
#endif
        for (UINT i = 0; i < Seqs.size(); i++)
        {
            /* Extract Sequences */
            STRING seq = Seqs.at(i);

            /* Copy into the seqPep.seqs array */
            memcpy((void *) &index->pepIndex.seqs[i * (seqlen)], (const void *) seq.c_str(), seqlen);

            /* Increment the counters */
            iCount += 2;

#ifdef DEBUG
            cout << seq << endl;
#endif /* DEBUG */
        }

        /* Get the peptide count */
        iCount /= 2;

#ifdef BENCHMARK
        memory += omp_get_wtime() - duration;
#endif
    }
    else
    {
        status = ERR_FILE_NOT_FOUND;
    }

    /* Fill in the peps entry */
    if (status == SLM_SUCCESS)
    {
        status = LBE_GeneratePeps(index);
    }

#ifdef VMODS

    if (index->lclmodCnt > 0)
    {
#ifdef BENCHMARK
        duration = omp_get_wtime();
#endif
        /* Fill in the mods entry */
        status = MODS_GenerateMods(index);
#ifdef BENCHMARK
        compute += omp_get_wtime() - duration;
#endif
    }
#endif /* VMODS */

    /* Make sure Seqs is clear anyway */
    Seqs.clear();

    /* Sort the peptide index based on peptide precursor mass */
    if (status == SLM_SUCCESS && params.dM > 0.0)
    {
        std::sort(index->pepEntries, index->pepEntries + index->lcltotCnt, CmpPepEntries);
    }

    /* Deinitialize varCount */
    if (varCount != NULL)
    {
        delete[] varCount;
        varCount = NULL;
    }

    if (status != SLM_SUCCESS)
    {
        (VOID) LBE_Deinitialize(index);
    }

    return status;
}

STATUS LBE_GeneratePeps(Index *index)
{
    STATUS status = SLM_SUCCESS;
    UINT interval = index->lclpepCnt;
    pepEntry *entries = index->pepEntries;
    UINT seqlen = Seqs.at(0).length();

#ifdef _OPENMP
    UINT threads = params.threads;
#endif /* _OPENMP */

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
    for (UINT fill = 0; fill < interval; fill++)
    {
        UINT idd = DSLIM_GenerateIndex(index, fill);

        entries[fill].Mass = UTILS_CalculatePepMass((AA *)Seqs.at(idd).c_str(), seqlen);
        entries[fill].seqID = idd;
        entries[fill].sites.sites = 0x0;
        entries[fill].sites.modNum = 0x0;
    }

    return status;
}
/*
 * FUNCTION: LBE_Deinitialize
 *
 * DESCRIPTION: Deallocate all memory and
 *              reset variables
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_Deinitialize(Index *index)
{
    return DSLIM_Deinitialize(index);
}

/*
 * FUNCTION: LBE_Distribute
 *
 * DESCRIPTION: Apply Distribution Policy on peptides
 *
 * INPUT:
 * @threads   : Number of parallel threads
 * @policy    : Policy ID (enum)
 * @slm_chunks: Number of distribution chunks
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_Distribute(Index *index)
{
    STATUS status = 0;
    UINT N = index->lcltotCnt;
    UINT speclen = (index->pepIndex.peplen-1) * params.maxz * iSERIES;
    UINT maxchunksize = (MAX_IONS / speclen);
    UINT maxchunksize2 = params.spadmem / (BYISIZE * params.threads);
    UINT nchunks = 0;
    UINT chunksize = 0;
    UINT lastchunksize = 0;

    /* Calculate the chunksize */
    chunksize = std::min(N, maxchunksize);
    chunksize = std::min(chunksize, maxchunksize2);

    /* Set the number of chunks to p */
    nchunks = (N / chunksize);

    if ((N % chunksize) != 0)
    {
        nchunks += 1;
    }

    /* Calculate the size of last chunk */
    UINT factor = N / chunksize;

    lastchunksize = ((N % chunksize) == 0)?
                     chunksize            :
                     N - (chunksize * factor);

    if (status == SLM_SUCCESS)
    {
        /* Return the number of chunks created */
        index->nChunks = nchunks;
        index->chunksize = chunksize;
        index->lastchunksize = lastchunksize;
    }

    return status;
}

/*
 * FUNCTION: LBE_CreatePartitions
 *
 * DESCRIPTION: Creates the partition size for the current node
 *
 * INPUT:
 * @index : (Distributed) Virtual Peptide ID.
 *
 * OUTPUT:
 * @status: Actual SPI peptide ID
 */
STATUS LBE_CreatePartitions(Index *index)
{
    STATUS status = SLM_SUCCESS;

    UINT N = index->pepCount;
    UINT p = params.nodes;
    UINT myid = params.myid;

    UINT chunksize = 0;

    /* More than one nodes in the system ? */
    if (p > 1)
    {
        /* Partition the pepCount */
        chunksize = N / p;

        if (((N % p) > 0) && ((N % p) > myid))
        {
            chunksize += 1;
        }

        index->lclpepCnt = chunksize;

        /* Partition the modCount */
        N = index->modCount;

        chunksize = N / p;

        if (((N % p) > 0) && ((N % p) > myid))
        {
            chunksize += 1;
        }

        index->lclmodCnt = chunksize;

        index->lcltotCnt = index->lclpepCnt + index->lclmodCnt;
    }
    else
    {
        index->lcltotCnt = index->totalCount;
        index->lclpepCnt = index->pepCount;
        index->lclmodCnt = index->modCount;
    }

    return status;
}

/*
 * FUNCTION: LBE_CountPeps
 *
 * DESCRIPTION: Count peptides in FASTA and the
 *              number of mods that will be generated
 *
 * INPUT:
 * @threads      : Number of parallel threads
 * @filename     : Path to FASTA file
 * @modconditions: Mod generation conditions
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_CountPeps(CHAR *filename, Index *index, UINT explen)
{
    STATUS status = SLM_SUCCESS;
    STRING line;
    FLOAT pepmass = 0.0;
    STRING modconditions = params.modconditions;
    UINT maxmass= params.max_mass;
    UINT minmass= params.min_mass;

    /* Initialize Index parameters */
    index->pepIndex.AAs = 0;
    index->pepCount = 0;
    index->modCount = 0;

#ifndef VMODS
    LBE_UNUSED_PARAM(modconditions);
#endif /* VMODS */
    
#ifdef BENCHMARK
    duration = omp_get_wtime();
#endif

    /* Open file */
    file.open(filename);

    if (file.is_open())
    {
        while (getline(file, line))
        {
            if (line.at(0) != '>')
            {
                /* Linux has a weird \r at end of each line */
                if (line.at(line.length() - 1) == '\r')
                {
                    line = line.substr(0, line.size() - 1);
                }

                /* Transform to all upper case letters */
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);

                /* Check the integrity of the peptide
                 * sequence being read */
                if (line.length() != explen)
                {
                    status = ERR_INVLD_SIZE;
                    break;
                }

                /* Calculate mass of peptide */
                pepmass = UTILS_CalculatePepMass((AA *)line.c_str(), line.length());

                /* Check if the peptide mass is legal */
                if (pepmass >= minmass && pepmass <= maxmass)
                {
                    index->pepCount++;
                    index->pepIndex.AAs += line.length();
                    Seqs.push_back(line);
                }
            }
        }
    }
    else
    {
        cout << endl << "FATAL: Could not read FASTA file" << endl;
        status = ERR_INVLD_PARAM;
    }

#ifdef BENCHMARK
    fileio += omp_get_wtime() - duration;
#endif

    /* Allocate teh varCount array */
    if (status == SLM_SUCCESS)
    {
        varCount = new UINT[index->pepCount + 1];
    }

#ifdef VMODS
    /* Count the number of variable mods given
     * modification information */
    if (status == SLM_SUCCESS)
    {
#ifdef BENCHMARK
        duration = omp_get_wtime();
#endif

        index->modCount = MODS_ModCounter(index);

#ifdef BENCHMARK
        compute += omp_get_wtime() - duration;
#endif
    }

#endif /* VMODS */

    /* Check if any errors occurred in MODS_ModCounter */
    if (index->modCount == (UINT)(-1) || index->pepIndex.AAs != index->pepCount * explen)
    {
        status = ERR_INVLD_SIZE;
    }

    /* Print if everything is okay */
    if (status == SLM_SUCCESS)
    {
        /* Return the total count */
        index->totalCount = index->pepCount + index->modCount;
        cumusize += index->totalCount;

        if (params.myid == 0)
        {
            cout << "Number of Peptides    =\t\t" << index->pepCount << endl;
            cout << "Number of Variants    =\t\t" << index->modCount << endl;
            cout << "Total Index Size      =\t\t" << index->totalCount << endl;
            cout << "Cumulative Index Size =\t\t" << cumusize << endl << endl;
        }

        /* Close the file once done */
        file.close();
    }

    return status;
}

/*
 * FUNCTION: LBE_PrintHeader
 *
 * DESCRIPTION: Prints the LBE header
 *
 * INPUT : none
 * OUTPUT: none
 */
VOID LBE_PrintHeader()
{
    cout << "\n"
            "***************************************\n"
            "* HiCOPS: HPC Database Peptide Search *\n"
            "* School of Computing & Info Sciences *\n"
            "*   Florida International University  *\n"
            "*         Miami, Florida, USA         *\n"
            "***************************************\n"
          << endl << endl;

    return;
}

static BOOL CmpPepEntries(pepEntry e1, pepEntry e2)
{
    return e1 < e2;
}
