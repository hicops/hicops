/*
 *  Copyright (C) 2019  Muhammad Usman Tariq
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
 *  GNU General Public License for more detailSpectrum.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "hyperscore.h"

/* Data structures related to the current experimental spectrum */
extern FLOAT *hyperscores;          /* Array to store the scores for the current especid */
UINT size;                          /* Size of the above array */
UINT especid;                       /* Current especid (experimental spectrum id) */

/* Data structures for the output file */
BOOL FileInitiated = false;      /* Flag to indicate that the output file has been initialized */
BOOL HeadersDone = false;        /* Flag that the headers have bee written to the output file */
std::ofstream myfile;            /* The output file */

/*
 * FUNCTION: HYPERSCORE_Calculate
 *
 * DESCRIPTION:
 *
 * INPUT:
 * @counts : structure to hold b and y counts and sum of intensities
 * @threshold : cutoff value for the number of peak matches to considered a candidate
 *
 * OUTPUT:
 * @status: Status of execution
 */

STATUS HS_InitFile()
{
    if (!FileInitiated)
    {
        std::stringstream sst;
        sst << HYPERSCORE_Datetime() << ".tsv";
        myfile.open(sst.str());
        FileInitiated = true;
    }

    return SLM_SUCCESS;
}

STATUS HS_DeinitFile()
{
    if (FileInitiated)
    {
        myfile.close();
        FileInitiated = false;
    }

    return SLM_SUCCESS;
}

STATUS HYPERSCORE_Calculate(UINT specid, INT psid, FLOAT maxhv)
{
    if (FileInitiated)
    {
        myfile << "Spectrum ID: " << specid + 1;

        myfile << '\t' << std::to_string(psid);

        myfile << '\t' << std::to_string(maxhv);

        myfile << std::endl;
    }

    return SLM_SUCCESS;
}

/*
 * FUNCTION: HYPERSCORE_WriteHeaders
 *
 * DESCRIPTION: Write the headers (column names) to the output file
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */

/*
 * FUNCTION: HYPERSCORE_WriteToFile
 *
 * DESCRIPTION: Write the hyperscores stored in the hypscores array to the output file
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS HYPERSCORE_WriteToFile ()
{
    myfile << "Spectrum ID: " << especid+1;

    myfile << ',' << std::to_string(hyperscores[0]);

    for (UINT i = 1; i < 10; i++)
    {
        myfile << "," << hyperscores[i];
    }

    myfile << std::endl;

    return SLM_SUCCESS;
}

/*
 * FUNCTION: HYPERSCORE_Factorial
 *
 * DESCRIPTION: Calculate the factorial of a number
 *
 * INPUT:
 * @n : input value for which to calculate factorial
 *
 * OUTPUT:
 * @factorial : the factorial of the input number n
 */
ULONGLONG HYPERSCORE_Factorial(ULONGLONG n)
{
    return (n == 1 || n == 0) ? 1 : HYPERSCORE_Factorial(n - 1) * n;
}


/*
 * FUNCTION: HYPERSCORE_Datetime
 *
 * DESCRIPTION: Get the date time in a readable format.
 *
 * INPUT:
 *
 * OUTPUT:
 * @datetime : date & time in string format
 */
std::string HYPERSCORE_Datetime()
{

    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,80,"%d-%m-%Y %H-%M-%S",timeinfo);
    return std::string(buffer);
}
