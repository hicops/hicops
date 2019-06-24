/*
 *  This file is part of SLM-Transform
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
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "lbe.h"

#ifdef _PROFILE
 #include <gperftools/profiler.h>
#endif /* _PROFILE */

using namespace std;

/* Global Variables */
Index *slm_index = NULL;
DIR*    dir;
dirent* pdir;
vector<STRING> queryfiles;
STRING dbfile;

gParams params;

static STATUS ParseParams(CHAR* paramfile);

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

    /* Benchmarking */
    auto start = chrono::system_clock::now();
    auto end   = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    chrono::duration<double> qtime = end - start;

    STRING patt[3] = {".ms2", ".mzML", "mzXML"};
    CHAR extension[] = ".peps";

    /* Print Header */
    LBE_PrintHeader();

    cout << endl << "Start Time: " << ctime(&start_time) << endl;

    if (argc < 2)
    {
        cout << "ERROR: Missing arguments\n";
        cout << "Format: ./cfir.exe <uparams.txt>\n";
        status = ERR_INVLD_PARAM;
        exit (status);
    }


    /* Parse the parameters */
    status = ParseParams(argv[1]);

#ifdef _PROFILE
    ProfilerStart("C:/work/lbe.prof");
#endif /* _PROFILE */


    /* Add all the query files to the vector */
    dir = opendir(params.datapath.c_str());

    /* Check if opened */
    if (dir != NULL)
    {
        while ((pdir = readdir(dir)) != NULL)
        {
            STRING cfile(pdir->d_name);

            /* Add the matching files */
            if (cfile.find(patt[0]) != std::string::npos)
            {
                queryfiles.push_back(params.datapath + '/' + pdir->d_name);
            }
            if (cfile.find(patt[1]) != std::string::npos)
            {
                queryfiles.push_back(params.datapath + '/' + pdir->d_name);
            }
            if (cfile.find(patt[2]) != std::string::npos)
            {
                queryfiles.push_back(params.datapath + '/' + pdir->d_name);
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

    /* Create local variables to avoid trouble */
    UINT minlen = params.min_len;
    UINT maxlen = params.max_len;
    UINT threads = params.threads;

    slm_index = new Index[maxlen - minlen + 1];

    if (slm_index == NULL)
    {
        status = ERR_INVLD_MEMORY;
    }

    for (UINT peplen = minlen; peplen <= maxlen; peplen++)
    {
        dbfile = params.dbpath + "/" + std::to_string(peplen) + extension;

        /* Set the peptide length in the pepIndex */
        slm_index[peplen-minlen].pepIndex.peplen = peplen;

        /* Count the number of ">" entries in FASTA */
        if (status == SLM_SUCCESS)
        {
            status = LBE_CountPeps(threads, (CHAR *) dbfile.c_str(), (slm_index + peplen-minlen));
        }

        /* Initialize internal structures */
        if (status == SLM_SUCCESS)
        {
            start = chrono::system_clock::now();

            /* Initialize the LBE */
            status = LBE_Initialize(threads, (slm_index + peplen-minlen));

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
            status = DSLIM_Construct(threads, (slm_index + peplen - minlen));

            end = chrono::system_clock::now();

            /* Compute Duration */
            elapsed_seconds = end - start;
            cout << "SLM-Transform with status:\t" << status << endl;
            cout << "Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;
        }
    }

    /* We don't need the original data anymore - Deallocate */
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
                    status = DSLIM_QuerySpectrum(QA, ms2specs, slm_index, (maxlen-minlen+1));
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

/* FUNCTION: ParseParams
 *
 * DESCRIPTION: Parse the input params and initialize
 *
 * INPUT: none
 *
 * OUTPUT
 * @status: Status of execution
 */
static STATUS ParseParams(CHAR* paramfile)
{
    STATUS status = SLM_SUCCESS;

    STRING line;

    ifstream pfile(paramfile);

    /* Check if mods file is open */
    if (pfile.is_open())
    {
        /* Get path to DBparts */
        getline(pfile, line);

        /* Check for a dangling / character */
        if (line.at(line.length()- 1) == '/')
        {
            line = line.substr(0, line.size() - 1);
        }

        params.dbpath = line;

        /* Get path to MS2 data */
        getline(pfile, line);

        /* Check for a dangling / character */
        if (line.at(line.length()- 1) == '/')
        {
            line = line.substr(0, line.size() - 1);
        }

        params.datapath = line;

        /* Get the max threads to use */
        getline(pfile, line);
        params.threads = std::atoi(line.c_str());

        /* Get the min peptide length */
        getline(pfile, line);
        params.min_len = std::atoi(line.c_str());

        /* Get the max peptide length */
        getline(pfile, line);
        params.max_len = std::atoi(line.c_str());

        /* Get the max fragment charge */
        getline(pfile, line);
        params.maxz = std::atoi(line.c_str());

        /* Get the fragment mass tolerance */
        getline(pfile, line);
        params.dF = (UINT)(std::atof(line.c_str()) * params.scale);

        /* Get the precursor mass tolerance */
        getline(pfile, line);
        params.dM = std::atof(line.c_str());

        /* Get the m/z axis resolution */
        getline(pfile, line);
        params.res = std::atof(line.c_str());

        /* Get the scaling factor */
        getline(pfile, line);
        params.scale = std::atoi(line.c_str());

        /* Get the min mass */
        getline(pfile, line);
        params.min_mass = std::atoi(line.c_str());

        /* Get the max mass */
        getline(pfile, line);
        params.max_mass = std::atoi(line.c_str());

        /* Get the top matches to report */
        getline(pfile, line);
        params.topmatches = std::atoi(line.c_str());

        /* Get the shp threshold */
        getline(pfile, line);
        params.min_shp = std::atoi(line.c_str());

        /* Get number of mods */
        getline(pfile, line);
        params.vModInfo.num_vars = std::atoi((const CHAR *) line.c_str());

        /* If no mods then init to 0 M 0 */
        if (params.vModInfo.num_vars == 0)
        {
            params.modconditions = "0 M 0";
        }
        else
        {
            /* Get max vmods per peptide sequence */
            getline(pfile, line);
            params.vModInfo.vmods_per_pep = std::atoi((const CHAR *) line.c_str());
            params.modconditions = std::to_string(params.vModInfo.vmods_per_pep);

            /* Fill in information for each vmod */
            for (USHORT md = 0; md < params.vModInfo.num_vars; md++)
            {
                /* Get and set the modAAs */
                getline(pfile, line);
                params.modconditions += " " + line;

                std::strncpy((char *) params.vModInfo.vmods[md].residues, (const char *) line.c_str(),
                        std::min(4, (const int) line.length()));

                /* get and set the modmass */
                getline(pfile, line);
                params.vModInfo.vmods[md].modMass = (UINT) (std::atof((const char *) line.c_str()) * params.scale);

                /* Get and set the modAAs_per_peptide */
                getline(pfile, line);
                params.modconditions += " " + line;

                params.vModInfo.vmods[md].aa_per_peptide = std::atoi((const char *) line.c_str());
            }
        }

        pfile.close();
    }
    else
    {
        status = ERR_FILE_NOT_FOUND;
    }

    return status;
}



