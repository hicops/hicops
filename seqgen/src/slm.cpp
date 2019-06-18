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

#ifdef _PROFILE
 #include <gperftools/profiler.h>
#endif /* _PROFILE */

using namespace std;

/* Global large strings */
STRING dbpath;
STRING params;
STRING dbfile;
STRING outfile;

/* FUNCTION: SLM_Main (main)
 *
 * DESCRIPTION: Driver Application
 *
 * INPUT: none
 *
 * OUTPUT
 * @status: Status of execution
 */
STATUS SLM_Main(INT argc, CHAR* argv[])
{
    STATUS status = SLM_SUCCESS;

    /* Print start time */
    auto start_tim = chrono::system_clock::now();
    time_t start_time = chrono::system_clock::to_time_t(start_tim);

    SLM_vMods vModInfo;

    /* Benchmarking */
    auto start = chrono::system_clock::now();
    auto end   = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    UINT maxz = MAXz;
    STRING modconditions;
    STRING line;

    /* Print Header */
    LBE_PrintHeader();

    cout << endl << "Start Time: " << ctime(&start_time) << endl;

    if (argc < 5)
    {
        cout << "Missing arguments\n";
        cout << "Format: ./pseq.exe /<dbpath>/ /<mods.txt> <minlen> <maxlen> <opt:charge> \n";
        status = ERR_INVLD_PARAM;
        exit (status);
    }

    /* Database and Dataset files */
    dbpath = argv[1];
    STRING extension = ".peps";

    STRING params = argv[2];
    UINT minlen = atoi(argv[3]);
    UINT maxlen = atoi(argv[4]);

    /* Max Charge info provided */
    if (argc == 6)
    {
        maxz = atoi(argv[5]);
    }

    /* Check for a dangling / character */
    if (dbpath.at(dbpath.length()- 1) == '/')
    {
        dbpath = dbpath.substr(0, dbpath.size() - 1);
    }

    /* Open the params file and parse for mods */
    ifstream pfile(params);

    /* Check if mods file is open */
    if (pfile.is_open())
    {
        /* Get number of mods */
        getline(pfile, line);
        vModInfo.num_vars = std::atoi((const char *) line.c_str());

        /* If no mods then init to 0 M 0 */
        if (vModInfo.num_vars == 0)
        {
            modconditions = "0 M 0";
        }
        else
        {
            /* Get max vmods per peptide sequence */
            getline(pfile, line);
            vModInfo.vmods_per_pep = std::atoi((const char *) line.c_str());
            modconditions = std::to_string(vModInfo.vmods_per_pep);

            /* Fill in information for each vmod */
            for (USHORT md = 0; md < vModInfo.num_vars; md++)
            {
                /* Get and set the modAAs */
                getline(pfile, line);
                modconditions += " " + line;

                std::strncpy((char *) vModInfo.vmods[md].residues, (const char *) line.c_str(),
                        std::min(4, (const int) line.length()));

                /* Get and set the modmass */
                getline(pfile, line);
                vModInfo.vmods[md].modMass = (UINT) (std::atof((const char *) line.c_str()) * SCALE);

                /* Get and set the modAAs_per_peptide */
                getline(pfile, line);
                modconditions += " " + line;

                vModInfo.vmods[md].aa_per_peptide = std::atoi((const char *) line.c_str());
            }
        }

        pfile.close();
    }
    else
    {
//        status = ERR_FILE_NOT_FOUND;
        cout << "WARNING: mods file not found, generating seqs without mods";
        modconditions = "0 M 0";
    }

#ifdef _PROFILE
    ProfilerStart("C:/work/lbe.prof");
#endif /* _PROFILE */

    /* Get thread count */
    UINT threads = UTILS_GetNumProcs();

    cout << "Num Threads: " << threads << endl;
    cout << modconditions << endl;

#ifndef _OPENMP
    threads = 1;
#endif /* _OPENMP */

    for (UINT peplen = minlen; peplen <= maxlen; peplen++)
    {
        dbfile = dbpath + "/" + std::to_string(peplen) + extension;
        outfile = dbpath + "/" + std::to_string(peplen) + ".peprec";

        /* Count the number of ">" entries in FASTA */
        if (status == SLM_SUCCESS)
        {
            status = LBE_CountPeps(threads, (CHAR *)dbfile.c_str(), modconditions);
        }

        /* Initialize internal structures */
        if (status == SLM_SUCCESS)
        {
            start = chrono::system_clock::now();

            /* Initialize the LBE */
            status = LBE_Initialize(threads, modconditions, &vModInfo);

            end = chrono::system_clock::now();

            /* Compute Duration */
            elapsed_seconds = end - start;
            cout << "Initialized with status:\t" << status << endl;
            cout << "Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;
        }

        if (status == SLM_SUCCESS)
        {
            start = chrono::system_clock::now();
            status = LBE_WritePEPREC(threads, (CHAR *) outfile.c_str(), peplen, maxz);
            end = chrono::system_clock::now();

            /* Compute Duration */
            elapsed_seconds = end - start;
            cout << "Writing PEPREC file with status:\t" << status << endl;
            cout << "Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;
        }

        /* Deinitialize DSLIM and LBE */
        if (status == SLM_SUCCESS)
        {
            status = LBE_Deinitialize();
            cout << endl << "Deinitialized with status:\t" << status << endl;
        }
    }
    /* Print final program status */
    cout << "\n\nEnded with status: \t\t" << status << endl;

    /* Print end time */
    auto end_tim = chrono::system_clock::now();
    time_t end_time = chrono::system_clock::to_time_t(end_tim);
    cout << endl << "End Time: " << ctime(&end_time) << endl;
    elapsed_seconds = end_tim - start_tim;
    cout << "Total Elapsed Time: " << elapsed_seconds.count() << "s" <<endl;

#ifdef _PROFILE
    ProfilerFlush();
    ProfilerStop();
#endif /* _PROFILE */

    /* Make sure stdout is empty at the end */
    fflush(stdout);

    return status;
}
