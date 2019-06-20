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

#include "lbe.h"
using namespace std;

ULONGLONG    pepCount;
ULONGLONG    modCount;
ULONGLONG    totalCount;
PepSeqs      seqPep;
ifstream     dbfile;
/* The output file */
ofstream outfile;
STRING       line;
vector<STRING> Seqs;


/* External Variables */
#ifdef VMODS
varEntry *modEntries;
#endif /* VMODS */

/* Static function Prototypes */
static STATUS LBE_AllocateMem(UINT N, UINT M);

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
static STATUS LBE_AllocateMem(UINT N, UINT M)
{
    STATUS status = SLM_SUCCESS;
#ifdef VMODS
    modEntries = NULL;
#endif /* VMODS */

    /* Allocate Memory for seqPep */
    seqPep.seqs = new AA[seqPep.AAs];
    seqPep.idx = new UINT[N+1];

    /* Set the seqPep.idx[N] to AAs for consistency*/
    seqPep.idx[N] = seqPep.AAs;

    if (seqPep.seqs == NULL || seqPep.idx == NULL)
    {
        status = ERR_BAD_MEM_ALLOC;
    }

#ifdef VMODS
    /* Allocate the seqMod */
    if (status == SLM_SUCCESS)
    {
        modEntries = new varEntry[M];

        if (modEntries == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
    }
#endif /* VMODS */


    return status;
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
STATUS LBE_Initialize(UINT threads, STRING modconditions, SLM_vMods* varMod)
{
    STATUS status = SLM_SUCCESS;
    UINT iCount = 1;
    STRING seq;

    /* Check if ">" entries are > 0 */
    if (pepCount > 0)
    {
        status = LBE_AllocateMem(pepCount, modCount);
    }
    else
    {
        status = ERR_INVLD_PARAM;
    }

    /* If Seqs was successfully filled */
    if (Seqs.size() != 0 && status == SLM_SUCCESS)
    {
        UINT idx = 0;

        /* Extract First Sequence manually */
        seq = Seqs.at(0);

        /* Increment counter */
        iCount += 2;

        memcpy((void *)&seqPep.seqs[idx], (const void *)seq.c_str(), seq.length());
        seqPep.idx[0] = idx;
        idx+=seq.length();

#ifdef DEBUG
        cout << seq << endl;
#endif /* DEBUG */

        for (UINT i = 1; i < Seqs.size(); i++)
        {
            /* Extract Sequences */
            STRING seq = Seqs.at(i);
            /* Copy into the seqPep.seqs array */
            memcpy((void *) &seqPep.seqs[idx], (const void *) seq.c_str(), seq.length());

            /* Increment the counters */
            iCount += 2;
            seqPep.idx[iCount / 2 - 1] = idx;
            idx += seq.length();

#ifdef DEBUG
            cout << seq << endl;
#endif /* DEBUG */
        }

        /* Get the peptide count */
        iCount /= 2;
    }
    else
    {
        status = ERR_FILE_NOT_FOUND;
        cout << endl << "ABORT: File could not be opened" << endl;
    }

    /* Make sure that the '>' == PEPTIDES */
    if (iCount != pepCount)
    {
        cout << endl << "pepCount != iCount - Please check the FASTA file";
        status = ERR_INVLD_SIZE;
    }

#ifdef VMODS

    /* Initialize gModInfo */
    status = UTILS_InitializeModInfo(varMod);

    if (modCount > 0)
    {
        /* Fill in the mods entry */
        status = MODS_GenerateMods(threads, modCount, modconditions);
    }
#endif /* VMODS */

    /* Make sure Seqs is clear anyway */
    Seqs.clear();

    if (status != SLM_SUCCESS)
    {
        (VOID) LBE_Deinitialize();
    }

    return status;
}

/*
 * FUNCTION: LBE_WritePEPREC
 *
 * DESCRIPTION: Prints the LBE header
 *
 * INPUT : threads
 * OUTPUT: status of execution
 */
STATUS LBE_WritePEPREC(UINT threads, CHAR *outpath, UINT peplen, UINT charge)
{
    STATUS status = SLM_SUCCESS;

    line = "spec_id modifications peptide charge\n";

    /* Do not apply the append flag in here */
	std::cout << outpath << std::endl;
    outfile.open((const CHAR *) outpath, std::ofstream::out);

    outfile << line;

    /* For each peptide entry, write a PEPREC entry */
    for (UINT ii = 0; ii < pepCount; ii++)
    {
        line = "pep" + std::to_string(ii+1) + " - ";

        for (UINT jj = 0; jj < peplen; jj++)
        {
            line += seqPep.seqs[(peplen * ii) + jj];
        }

        line += " " + std::to_string(charge) + "\n";

        outfile << line;
    }

    for (UINT ii = 0; ii < modCount; ii++)
    {
        UINT modsSeen = 0;
        line = "var" + std::to_string(ii+1) + " ";

        /* Write the mod and position info */
        for (UINT mm = 0; mm < peplen; mm++)
        {
            if (ISBITSET(modEntries[ii].sites.sites, mm))
            {
                int modNum = ((modEntries[ii].sites.modNum & (0x0F << (4 * modsSeen))) >> (4 * modsSeen));

                if (modsSeen != 0)
                {
                    line += "|";
                }

                modsSeen++;

                line += std::to_string(mm+1) + "|mod" + std::to_string(modNum) + seqPep.seqs[(peplen * modEntries[ii].seqID) + mm];
            }
        }

        line += " ";

        /* Write the peptide character */
        for (UINT jj = 0; jj < peplen; jj++)
        {
            line += seqPep.seqs[(peplen * modEntries[ii].seqID) + jj];
        }

        /* Add the max charge info */
        line += " " + std::to_string(charge) + "\n";

        outfile << line;
    }

	outfile.close();

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
STATUS LBE_Deinitialize(VOID)
{
   //(VOID) DSLIM_Deinitialize();

    /* Reset all counters */
    seqPep.AAs = 0;
    pepCount = 0;
    modCount = 0;
    totalCount = 0;

    /* Deallocate memory */
    if (seqPep.idx != NULL)
    {
        delete[] seqPep.idx;
        seqPep.idx = NULL;
    }

    if (seqPep.seqs != NULL)
    {
        delete[] seqPep.seqs;
        seqPep.seqs = NULL;
    }

#ifdef VMODS
    if (modEntries != NULL)
    {
        delete[] modEntries;
        modEntries = NULL;
    }
#endif /* VMODS */

    return SLM_SUCCESS;
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
STATUS LBE_CountPeps(UINT threads, CHAR *filename, STRING modconditions)
{
    STATUS status = SLM_SUCCESS;
    FLOAT pepmass = 0.0;

    pepCount = 0;
    modCount = 0;
    seqPep.AAs = 0;

#ifndef VMODS
    LBE_UNUSED_PARAM(modconditions);
#endif /* VMODS */

    cout << filename << endl;

    /* Open file */
    dbfile.open(filename);

    if (dbfile.is_open())
    {
        while (getline(dbfile, line))
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

                /* Calculate mass of peptide */
                pepmass = UTILS_CalculatePepMass((AA *)line.c_str(), line.length());

                /* Check if the peptide mass is legal */
                if (pepmass >= MIN_MASS && pepmass <= MAX_MASS)
                {
                    pepCount++;
                    seqPep.AAs += line.length();
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

#ifdef VMODS
    /* Count the number of variable mods given
     * modification information */
    if (status == SLM_SUCCESS)
    {
        modCount = MODS_ModCounter(threads, modconditions);
    }

#endif /* VMODS */

    /* Print if everything is okay */
    if (status == SLM_SUCCESS)
    {
        /* Return the total count */
        totalCount = pepCount + modCount;

        cout << "Number of Peptides   = \t\t" << pepCount << endl;
        cout << "Number of Variants   = \t\t" << modCount << endl;
        cout << "Total Index Size     = \t\t" << totalCount << endl << endl;

        /* Close the file once done */
        dbfile.close();
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
VOID LBE_PrintHeader(VOID)
{
    cout << "\n"
            "********************************\n"
            "      Peptide Data Generator    \n"
            "Florida International University\n"
            "        Miami, FL, USA          \n"
            "********************************"
          << endl << endl;

    return;
}
