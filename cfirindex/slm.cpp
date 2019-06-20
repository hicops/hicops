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

Index *slm_index = NULL;
STRING dbpath;
STRING querypath;
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
    ULONGLONG Matches = 0;

    SLM_vMods vModInfo;

    /* Benchmarking */
    auto start = chrono::system_clock::now();
    auto end   = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    chrono::duration<double> qtime = end - start;
    UINT maxz = MAXz;
    STRING modconditions;
    STRING line;

    /* Print Header */
    LBE_PrintHeader();

    cout << endl << "Start Time: " << ctime(&start_time) << endl;

    if (argc < 5)
    {
        cout << "ERROR: Missing Params\n"
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

    /* Check for a dangling / character */
    if (querypath.at(querypath.length()- 1) == '/')
    {
        querypath = querypath.substr(0, querypath.size() - 1);
    }

    /* Open all the query files */
    dir = opendir(querypath.c_str());

    /* Check if opened */
    if (dir != NULL)
    {
        while ((pdir = readdir(dir)) != NULL)
        {
            string cfile(pdir->d_name);

            /* Only add if there is a matching file */
            if (cfile.find(patt) != std::string::npos)
            {
                queryfiles.push_back(querypath + '/' + pdir->d_name);
            }
        }
    }
    else
    {
        status = ERR_FILE_NOT_FOUND;
    }

#ifndef _OPENMP
    threads = 1;
#endif /* _OPENMP */

    for (UINT peplen = minlen; peplen <= maxlen; peplen++)
    {
        dbfile = dbpath + "/" + std::to_string(peplen) + extension;

        /* Set the peptide length in the pepIndex */
        slm_index[peplen-minlen].pepIndex.peplen = peplen;

        /* Count the number of ">" entries in FASTA */
        if (status == SLM_SUCCESS)
        {
            status = LBE_CountPeps(threads, (CHAR *) dbfile.c_str(), modconditions, (slm_index + peplen-minlen));
        }

        /* Initialize internal structures */
        if (status == SLM_SUCCESS)
        {
            start = chrono::system_clock::now();

            /* Initialize the LBE */
            status = LBE_Initialize(threads, modconditions, (slm_index + peplen-minlen));

            end = chrono::system_clock::now();

            /* Compute Duration */
            elapsed_seconds = end - start;
            cout << "Initialized with status:\t" << status << endl;
            cout << "Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;
        }

        /* Distribution Algorithm */
        if (status == SLM_SUCCESS)
        {
            start = chrono::system_clock::now();

            /* Distribute peptides among cores */
            status = LBE_Distribute(threads, _chunk, (slm_index + peplen - minlen));

        }

        /* DSLIM-Transform */
        if (status == SLM_SUCCESS)
        {
            start = chrono::system_clock::now();

            /* Construct DSLIM by SLM Transformation */
            status = DSLIM_Construct(threads, &vModInfo, dbpath, (slm_index + peplen - minlen));

            end = chrono::system_clock::now();

            /* Compute Duration */
            elapsed_seconds = end - start;
            cout << "SLM-Transform with status:\t" << status << endl;
            cout << "Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;
        }
    }

    if (status == SLM_SUCCESS)
    {
        status = DSLIM_DeallocateSpecArr();
    }

    /* Query the index */
    if (status == SLM_SUCCESS)
    {
        /* Allocate the Query Array */
        UINT *QA = new UINT[QCHUNK * QALEN];

        /* Initialize and process Query Spectra */
        for (UINT qf = 0; qf < queryfiles.size(); qf++)
        {
            /* Initialize Query MS/MS file */
            status = MSQuery_InitializeQueryFile((CHAR *) queryfiles[qf].c_str());
            cout << "Query File: " << queryfiles[qf] << endl;

            UINT spectra = 0;

            /* DSLIM Query Algorithm */
            if (status == SLM_SUCCESS)
            {
                /* Extract a chunk of MS/MS spectra and
                 * query against DSLIM Index */
                for (; QA != NULL;)
                {
                    /* Extract a chunk and return the chunksize */
                    UINT ms2specs = MSQuery_ExtractQueryChunk(QA, threads);

                    /* If the chunksize is zero, all done */
                    if (ms2specs <= 0)
                    {
                        break;
                    }
                    spectra += ms2specs;

                    start = chrono::system_clock::now();

                    /* Query the chunk */
                    status = DSLIM_QuerySpectrum(QA, ms2specs, Matches, threads, slm_index, (maxlen-minlen+1));
                    end = chrono::system_clock::now();

                    /* Compute Duration */
                    qtime += end - start;

                }
            }

            /* Compute Duration */
            cout << "Queried Spectra:\t\t" << spectra << endl;
            cout << "Query Time: " << qtime.count() << "s" << endl;
            cout << "Queried with status:\t\t" << status << endl << endl;
            end = chrono::system_clock::now();
        }

        if (QA != NULL)
        {
            delete[] QA;
            QA = NULL;
        }
    }

    /* Print end time */
    auto end_tim = chrono::system_clock::now();
    time_t end_time = chrono::system_clock::to_time_t(end_tim);
    cout << endl << "End Time: " << ctime(&end_time) << endl;
    elapsed_seconds = end_tim - start_tim;
    cout << "Total Elapsed Time: " << elapsed_seconds.count() << "s" <<endl;

    for (UINT peplen = minlen; peplen <= maxlen; peplen++)
    {
        status = LBE_Deinitialize(slm_index + peplen - minlen);
    }

    /* Print final program status */
    cout << "\n\nEnded with status: \t\t" << status << endl;

#ifdef _PROFILE
    ProfilerFlush();
	ProfilerStop();
#endif /* _PROFILE */

    /* Make sure stdout is empty at the end */
    fflush(stdout);

    return status;
}

