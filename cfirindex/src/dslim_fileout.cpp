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
 * GNU General Public License for more detailSpectrum.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <vector>
#include "common.h"
#include "dslim_fileout.h"

/* Global parameters */
extern gParams params;

BOOL FilesInit = false;

/* Data structures for the output file */
std::ofstream *tsvs = NULL; /* The output files */

/*
 * FUNCTION: DFile_InitFile
 *
 * DESCRIPTION:
 *
 * INPUT:
 * none
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DFile_InitFiles()
{
    STATUS status = SLM_SUCCESS;

    if (FilesInit == false)
    {
        status = ERR_BAD_MEM_ALLOC;

        tsvs = new std::ofstream[params.threads];

        STRING common = params.workspace + '/' + DFile_Datetime();

        if (tsvs != NULL)
        {
            status = SLM_SUCCESS;

            FilesInit = true;

            for (UINT f = 0; f < params.threads; f++)
            {
                STRING filename = common + "_" + std::to_string(f) + ".tsv";
                tsvs[f].open(filename);

#ifndef ANALYSIS
                tsvs[f] << "scan_num\t" << "prec_mass\t" << "peptide\t"
                        << "matched_ions\t" << "total_ions\t"
                        << "calc_pep_mass\t" << "mass_diff\t" << "mod_info\t"
                        << "hyperscore\t" << "expectscore\t" << "num_hits" << std::endl;
#else
                tsvs[f] << "scan_num\t" << "cpsm\t" << "weight\t"
                        << "bias\t" << "min\t"
                        << "max" << std::endl;
#endif /* ANALYSIS */
            }
        }
    }

    return status;
}

STATUS DFile_DeinitFiles()
{
    for (UINT i = 0; i < params.threads; i++)
    {
        tsvs[i].close();
    }

    delete[] tsvs;

    tsvs = NULL;

    FilesInit = false;

    return SLM_SUCCESS;
}

STATUS  DFile_PrintPartials(UINT specid, Results *resPtr)
{
    STATUS status = SLM_SUCCESS;
    UINT thno = omp_get_thread_num();

    tsvs[thno]         << std::to_string(specid + 1);
    tsvs[thno] << '\t' << std::to_string(resPtr->cpsms);
    tsvs[thno] << '\t' << std::to_string(resPtr->weight);
    tsvs[thno] << '\t' << std::to_string(resPtr->bias);
    tsvs[thno] << '\t' << std::to_string(resPtr->minhypscore);
    tsvs[thno] << '\t' << std::to_string(resPtr->nexthypscore);

    tsvs[thno] << std::endl;

    return status;

}

STATUS DFile_PrintScore(Index *index, UINT specid, FLOAT pmass, hCell* psm, DOUBLE e_x, UINT npsms)
{
    UINT thno = omp_get_thread_num();

    Index * lclindex = index + psm->idxoffset;
    INT peplen = lclindex->pepIndex.peplen;
    INT pepid = psm->psid;

    /* The size is peplen + 1 to add the \0 character at the end */
    CHAR pepseq[peplen + 1];

    /* Adds the string followed by the \0 character into the buffer */
    snprintf((CHAR *)&(pepseq[0]), peplen + 1, lclindex->pepIndex.seqs +
             (lclindex->pepEntries[psm->psid].seqID * peplen));

    /* Make a string from the char [] */
    STRING pep = pepseq;

    /* Print the PSM info to the file */
    tsvs[thno]         << std::to_string(specid + 1);
    tsvs[thno] << '\t' << std::to_string(pmass);
    tsvs[thno] << '\t' << pep;
    tsvs[thno] << '\t' << std::to_string(psm->sharedions);
    tsvs[thno] << '\t' << std::to_string(psm->totalions);
    tsvs[thno] << '\t' << std::to_string(lclindex->pepEntries[pepid].Mass);
    tsvs[thno] << '\t' << std::to_string((pmass - lclindex->pepEntries[pepid].Mass));
    tsvs[thno] << '\t'; // TODO: mod_info
    tsvs[thno] << '\t' << std::to_string(psm->hyperscore);
    tsvs[thno] << '\t' << std::to_string(e_x);
    tsvs[thno] << '\t' << std::to_string(npsms);

    tsvs[thno] << std::endl;

    return SLM_SUCCESS;
}

/*
 * FUNCTION: DFile_Factorial
 *
 * DESCRIPTION: Calculate the factorial of a number
 *
 * INPUT:
 * @n : input value for which to calculate factorial
 *
 * OUTPUT:
 * @factorial : the factorial of the input number n
 */
ULONGLONG DFile_Factorial(ULONGLONG n)
{
    return (n == 1 || n == 0) ? 1 : DFile_Factorial(n - 1) * n;
}

/*
 * FUNCTION: DFile_Datetime
 *
 * DESCRIPTION: Get the date time in a readable format.
 *
 * INPUT:
 *
 * OUTPUT:
 * @datetime : date & time in string format
 */
STRING DFile_Datetime()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%m.%d.%Y_%H.%M.%S_", timeinfo);

    std::string b(buffer);
    return b + std::to_string(params.myid);
}
