/*
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

#include "counter.hpp"

using namespace std;

/* Global Variables */
string_t dbfile;

extern ull_t cumusize;
extern ull_t ions;

gParams params;

static status_t ParseParams(char_t* paramfile);

/* FUNCTION: SLM_Main (main)
 *
 * DESCRIPTION: Driver Application
 *
 * INPUT: none
 *
 * OUTPUT
 * @status: Status of execution
 */
status_t main(int_t argc, char_t* argv[])
{
    status_t status = SLM_SUCCESS;

    char_t extension[] = ".peps";

    if (argc < 2)
    {
        std::cout << "ERROR: Missing arguments\n";
        std::cout << "Format: ./counter.exe <uparams.txt>\n";
        status = ERR_INVLD_PARAM;
        exit (status);
    }

    /* Parse the parameters */
    status = ParseParams(argv[1]);

    /* Create local variables to avoid trouble */
    uint_t minlen = params.min_len;
    uint_t maxlen = params.max_len;

    for (uint_t peplen = minlen; peplen <= maxlen; peplen++)
    {
        dbfile = params.dbpath + "/" + std::to_string(peplen) + extension;

        /* Count the number of ">" entries in FASTA */
        status = DBCounter((char_t *) dbfile.c_str());

    }

    /* The only output should be the cumulative size of the index */
    std::cout << "spectra:" << cumusize << std::endl;
	std::cout << "ions:" << ions << std::endl;

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
static status_t ParseParams(char_t* paramfile)
{
    status_t status = SLM_SUCCESS;

    string_t line;

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
#ifdef USE_OMP
        params.threads = std::atoi(line.c_str());
#else
        params.threads = 1;
#endif /* USE_OMP */

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
        params.dF = (uint_t)(std::atof(line.c_str()) * params.scale);

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

        getline(pfile, line);
        params.spadmem = std::atoi(line.c_str());

        params.spadmem *= 1024 * 1024; // Convert to MBs (max scratch space for score card)

        /* Get the distribution policy */
        getline(pfile, line);

        if (!line.compare("cyclic"))
        {
            params.policy = _cyclic;
        }
        else if (!line.compare("chunk"))
        {
            params.policy = _chunk;
        }
        else if (!line.compare("zigzag"))
        {
            params.policy = _zigzag;
        }

        /* Get number of mods */
        getline(pfile, line);
        params.vModInfo.num_vars = std::atoi((const char_t *) line.c_str());

        /* If no mods then init to 0 M 0 */
        if (params.vModInfo.num_vars == 0)
        {
            params.modconditions = "0 M 0";
        }
        else
        {
            /* Get max vmods per peptide sequence */
            getline(pfile, line);
            params.vModInfo.vmods_per_pep = std::atoi((const char_t *) line.c_str());
            params.modconditions = std::to_string(params.vModInfo.vmods_per_pep);

            /* Fill in information for each vmod */
            for (ushort_t md = 0; md < params.vModInfo.num_vars; md++)
            {
                /* Get and set the modAAs */
                getline(pfile, line);
                params.modconditions += " " + line;

                std::strncpy((char *) params.vModInfo.vmods[md].residues, (const char *) line.c_str(),
                        std::min(4, (const int) line.length()));

                /* get and set the modmass */
                getline(pfile, line);
                params.vModInfo.vmods[md].modMass = (uint_t) (std::atof((const char *) line.c_str()) * params.scale);

                /* Get and set the modAAs_per_peptide */
                getline(pfile, line);
                params.modconditions += " " + line;

                params.vModInfo.vmods[md].aa_per_peptide = std::atoi((const char *) line.c_str());
            }
        }

        params.myid = 0;
        params.nodes = 1;

        pfile.close();
    }
    else
    {
        status = ERR_FILE_NOT_FOUND;
    }

    return status;
}
