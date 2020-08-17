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
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "hicops.hpp"

using namespace std;

/* Global Variables */
Index *slm_index = NULL;
DIR*    dir;
dirent* pdir;
vector<string_t> queryfiles;
string_t dbfile;

gParams params;

static status_t ParseParams(char_t* paramfile);

/* FUNCTION: main
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

    // file extensions to find
    const string_t patt = {".ms2"};
    const char_t extension[] = ".peps";

    /* Benchmarking */
    double_t elapsed_seconds;

    /* Print start time */
    MARK_START(start_tim);
    const time_t start_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

#if defined (USE_TIMEMORY)
    // reset any previous configuration
    bundle_t::reset();

    // get extra tools from the environment variable
    auto env_tool = tim::get_env<std::string>("HICOPS_INST_COMPONENTS", "");
    auto env_enum = tim::enumerate_components(tim::delimit(env_tool));
    env_enum.erase(std::remove_if(env_enum.begin(), env_enum.end(),
                                  [](int c) { return c == WALL_CLOCK || 
                                                     c == CPU_UTIL; }),
                                  env_enum.end());

    // configure PAPI events
    const std::string def_papi_evts = "PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L3_TCM, PAPI_L2_TCA, PAPI_L3_TCA, PAPI_MEM_WCY, PAPI_RES_STL, PAPI_STL_CCY, PAPI_BR_CN, PAPI_BR_PRC, PAPI_FUL_ICY";
    std::string papi_evts = tim::get_env<std::string>("HICOPS_PAPI_EVENTS", def_papi_evts);
    tim::settings::papi_events() = papi_evts;

    // configure the bundle
    tim::configure<bundle_t>(env_enum);

    // start a timer
    time_tuple_t total("total_time");
#endif

    if (argc < 2)
    {
        std::cout << "ABORT: Missing parameters\n";
        std::cout << "Format: mpirun -np <k> hicops <uparams.txt>\n";
        status = ERR_INVLD_PARAM;
        exit (status);
    }

#ifdef USE_MPI
#    if defined (USE_TIMEMORY)
    // timemory MPI settings
    tim::settings::mpi_thread() = true;
    tim::settings::mpi_thread_type() = "multiple";

    // init MPI via timemory
    tim::mpi::initialize(argc, argv);

#   else
    // init MPI
    int provided = -1;
    status = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    
    // Check if desired MPI level available 
    if (provided != MPI_THREAD_MULTIPLE)
    {
        std::cout << "************************ Warning: ************************\n";
                "Your MPI distribution does not support MPI_THREAD_MULTIPLE\n";
                "HiCOPS may not work properly\n";
                "**********************************************************\n";
    }
#   endif // USE_TIMEMORY
#endif // USE_MPI

    // --------------------------------------------------------------------------------------------- //

    //
    // Initialization
    //

#if defined (USE_TIMEMORY)
    tim::timemory_init(argc, argv);
#endif // USE_TIMEMORY

    // parse parameters
    if (status == SLM_SUCCESS)
    {
        status = ParseParams(argv[1]);
    }

    // Print HiCOPS header after the ranks have been assigned
    if (params.myid == 0)
    {
        LBE_PrintHeader();
        std::cout << std::endl << "Start Time: " << ctime(&start_time) << std::endl;
    }

    /* Add all the query files to the vector */
    dir = opendir(params.datapath.c_str());

    /* Check if opened */
    if (dir != NULL)
    {
        while ((pdir = readdir(dir)) != NULL)
        {
            string_t cfile(pdir->d_name);

            /* Add the matching files */
            if (cfile.find(patt/* patt[0] */) != std::string::npos)
            {
                queryfiles.push_back(params.datapath + '/' + pdir->d_name);
            }

#ifdef FUTURE
            if (cfile.find(patt[1]) != std::string::npos)
            {
                queryfiles.push_back(params.datapath + '/' + pdir->d_name);
            }
            if (cfile.find(patt[2]) != std::string::npos)
            {
                queryfiles.push_back(params.datapath + '/' + pdir->d_name);
            }
#endif /* FUTURE */
        }
    }

    /* No file to query - Abort */
    if (queryfiles.size() < 1)
    {
        status = ERR_FILE_NOT_FOUND;
    }

    /* Create local variables to avoid trouble */
    uint_t minlen = params.min_len;
    uint_t maxlen = params.max_len;

    /* Create (max - min + 1) instances of SLM_Index */
    if (status == SLM_SUCCESS)
    {
        slm_index = new Index[maxlen - minlen + 1];

        /* Check if successful memory allocation */
        if (slm_index == NULL)
        {
            status = ERR_INVLD_MEMORY;
        }
    }

    // Initialize the mod information
    if (status == SLM_SUCCESS)
    {
        status = UTILS_InitializeModInfo(&params.vModInfo);
    }

    // Initialize the ModGen Engine
    if (status == SLM_SUCCESS)
    {
        status = MODS_Initialize();
    }

    // --------------------------------------------------------------------------------------------- //

    //
    // Indexing
    //

#if defined (USE_TIMEMORY)
    time_tuple_t index_inst("indexing");
#endif

    for (uint_t peplen = minlen; peplen <= maxlen; peplen++)
    {
        dbfile = params.dbpath + "/" + std::to_string(peplen) + extension;

        /* Set the peptide length in the pepIndex */
        slm_index[peplen-minlen].pepIndex.peplen = peplen;

        /* Count the number of ">" entries in FASTA */
        if (status == SLM_SUCCESS)
        {
            MARK_START(lbe_cnt);

            status = LBE_CountPeps((char_t *) dbfile.c_str(), (slm_index + peplen-minlen), peplen);

            MARK_END(lbe_cnt);

            /* Compute Duration */
            elapsed_seconds = ELAPSED_SECONDS(lbe_cnt);

            if (params.myid == 0)
            {
                std::cout << "DONE: Peptide Counting:\tstatus: " << status << std::endl << std::endl;
                PRINT_ELAPSED(elapsed_seconds);
            }
        }

        if (status == SLM_SUCCESS)
        {
            MARK_START(parts);

            status  = LBE_CreatePartitions((slm_index + peplen-minlen));

            MARK_END(parts);

            /* Compute Duration */
            elapsed_seconds = ELAPSED_SECONDS(parts);

            if (params.myid == 0)
            {
                std::cout << "DONE: Partition:\tstatus: " << status << std::endl << std::endl;
                PRINT_ELAPSED(elapsed_seconds);
            }
        }

        /* Initialize internal structures */
        if (status == SLM_SUCCESS)
        {
            MARK_START(lbe_init);

            /* Initialize the LBE */
            status = LBE_Initialize((slm_index + peplen - minlen));

            MARK_END(lbe_init);

            /* Compute Duration */
            elapsed_seconds = ELAPSED_SECONDS(lbe_init);

            if (params.myid == 0)
            {
                std::cout << "DONE: Initialize:\tstatus: " << status << std::endl << std::endl;
                PRINT_ELAPSED(elapsed_seconds);
            }
        }

        /* Distribution Algorithm */
        if (status == SLM_SUCCESS)
        {
            MARK_START(lbe_dist);

            /* Distribute peptides among cores */
            status = LBE_Distribute((slm_index + peplen - minlen));

            MARK_END(lbe_dist);

            /* Compute Duration */
            elapsed_seconds = ELAPSED_SECONDS(lbe_dist);

            if (params.myid == 0)
            {
                std::cout << "DONE: Internal Partition:\tstatus: " << status << std::endl << std::endl;
                PRINT_ELAPSED(elapsed_seconds);
            }
        }

        /* DSLIM-Transform */
        if (status == SLM_SUCCESS)
        {
            MARK_START(dslim);

            /* Construct DSLIM by SLM Transformation */
            status = DSLIM_Construct((slm_index + peplen - minlen));

            MARK_END(dslim);

            /* Compute Duration */
            elapsed_seconds = ELAPSED_SECONDS(dslim);

            if (params.myid == 0)
            {
                std::cout << "DONE: Indexing:\tstatus: " << status << std::endl << std::endl;
                PRINT_ELAPSED(elapsed_seconds);
            }
        }
    }

        /* We don't need the original data anymore */
    if (status == SLM_SUCCESS)
    {
        status = DSLIM_DeallocateSpecArr();
    }

    /* Initialize the Scorecard */
    if (status == SLM_SUCCESS)
    {
        status = DSLIM_InitializeScorecard(slm_index, (maxlen - minlen + 1));
    }

#if defined (USE_TIMEMORY)
    // stop measurements for indexing
    index_inst.stop();
#endif

    // print indexing time
    if (status == SLM_SUCCESS && params.myid == 0)
    {
#if defined (USE_TIMEMORY)
        auto wc = index_inst.get<wall_clock>();
        std::cout << "Indexing Time: " << wc->get() << "s" << std::endl;
#else
        std::cout << "Indexing Time: " << ELAPSED_SECONDS(start_tim) << "s" << std::endl;
#endif // USE_TIMEMORY
    }

    // --------------------------------------------------------------------------------------------- //

    //
    // Distributed Search
    //
#if defined (USE_TIMEMORY)
    time_tuple_t search_inst("search");
    mem_tuple_t  search_mem_inst("search_mem");
#endif // USE_TIMEMORY

    // Perform the distributed database search 
    if (status == SLM_SUCCESS)
    {
        MARK_START(dslim_search);
        status = DSLIM_SearchManager(slm_index);
        MARK_END(dslim_search);
        elapsed_seconds = ELAPSED_SECONDS(dslim_search);

        if (params.myid == 0)
        {
            std::cout << "DONE: Search:\tstatus: " << status << std::endl;
            PRINT_ELAPSED(elapsed_seconds);
        }
    }

    /* Deinitialize the scorecard */
    if (status == SLM_SUCCESS)
    {
        /* Deallocate the scorecard */
        status = DSLIM_DeallocateSC();
    }

    /* De-initialize the ion index */
    for (uint_t peplen = minlen; peplen <= maxlen; peplen++)
    {
        status = DSLIM_DeallocateIonIndex(slm_index + peplen - minlen);
    }

#if defined (USE_TIMEMORY)
    // stop instrumentation
    search_inst.stop();
    search_mem_inst.stop();
#endif

    // --------------------------------------------------------------------------------------------- //
    
    //
    // Merging
    //
#if defined (USE_TIMEMORY)
    time_tuple_t merge_inst("merge");
    mem_tuple_t merge_mem_inst("merge_mem");
#endif

    /* Compute the distributed scores */
    if (status == SLM_SUCCESS)
    {
        MARK_START(dist_score);
        status = DSLIM_DistScoreManager();
        MARK_END(dist_score);

        elapsed_seconds = ELAPSED_SECONDS(dist_score);

        if (params.myid == 0)
        {
            std::cout << "\nDONE: Merge:\tstatus: " << status << std::endl;
            PRINT_ELAPSED(elapsed_seconds);
        }
    }
#if defined (USE_TIMEMORY)
    // stop instrumentation
    merge_inst.stop();
    merge_mem_inst.stop();
#endif

    // --------------------------------------------------------------------------------------------- //

    //
    // Deinitialize
    //

    /* De-initialize the remaining index */
    if (status == SLM_SUCCESS)
    {
        for (uint_t peplen = minlen; peplen <= maxlen; peplen++)
        {
            status = DSLIM_DeallocatePepIndex(slm_index + peplen - minlen);
        }
    }

#ifdef DIAGNOSE2
        std::cout << "SCProc DONE@ " << params.myid << std::endl;
#endif /* DIAGNOSE2 */

    /* Delete the index Handle */
    if (status == SLM_SUCCESS)
    {
        delete [] slm_index;
        slm_index = NULL;
    }

#ifdef USE_MPI
    if (status == SLM_SUCCESS && params.nodes > 1)
    {
#   if defined (USE_TIMEMORY)
        tim::mpi::barrier();
#   else
        /* Wait for everyone to synchronize */
        status = MPI_Barrier(MPI_COMM_WORLD);
#endif // USE_TIMEMORY
    }
#endif /* USE_MPI */

    /* Print end time */
    const auto end_tim = chrono::system_clock::now();
    const time_t end_time = chrono::system_clock::to_time_t(end_tim);

    if (params.myid == 0)
    {
        std::cout << std::endl << "End Time: " << ctime(&end_time) << std::endl;
#if defined (USE_TIMEMORY)

    total.stop();
    auto tt = total.get<wall_clock>();
    auto es = tt->get();
#else
        auto es = ELAPSED_SECONDS(start_tim);
#endif
        std::cout << "Total Elapsed Time: " << es << "s" << std::endl;

        /* Print final program status */
        std::cout << "\n\nEnded with status: \t\t" << status << std::endl << std::endl;
    }

    /* Make sure stdout is empty at the end */
    fflush(stdout);

#if defined (USE_TIMEMORY)
    tim::timemory_finalize();
#endif

#ifdef USE_MPI
#   if defined (USE_TIMEMORY)
    tim::mpi::finalize();
#   else
    status = MPI_Finalize();
#   endif // USE_TIMEMORY
#endif // USE_MPI

    /* Return the status of execution */
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

        /* Get path to workspace path */
        getline(pfile, line);

        /* Check for a dangling / character */
        if (line.at(line.length()- 1) == '/')
        {
            line = line.substr(0, line.size() - 1);
        }

        params.workspace = line;

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

        /* Get the max expect score to report */
        getline(pfile, line);
        params.expect_max = std::atof(line.c_str());

        /* Get the shp threshold */
        getline(pfile, line);
        params.min_shp = std::atoi(line.c_str());

        /* Get the minhits threshold */
        getline(pfile, line);
        params.min_cpsm = std::atoi(line.c_str());

        /* Base Intensity */
        getline(pfile, line);
        params.base_int = std::atoi(line.c_str());

        /* Cutoff intensity ratio */
        getline(pfile, line);
        params.min_int = (int_t)(std::atof(line.c_str()) * (double_t)(params.base_int));

        /* Get the scorecard + scratch memory */
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

                /* Convert to all uppercases if not already */
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);

                params.modconditions += " " + line;

                std::strncpy((char *) params.vModInfo.vmods[md].residues, (const char *) line.c_str(),
                        std::min(4, static_cast<const int>(line.length())));

                /* get and set the modmass */
                getline(pfile, line);
                params.vModInfo.vmods[md].modMass = (uint_t) (std::atof((const char *) line.c_str()) * params.scale);

                /* Get and set the modAAs_per_peptide */
                getline(pfile, line);
                params.modconditions += " " + line;

                params.vModInfo.vmods[md].aa_per_peptide = std::atoi((const char *) line.c_str());
            }
        }

        params.perf = new double_t[params.nodes];

        /* Initialize to 1.0 for now */
        for (uint_t nn = 0; nn < params.nodes; nn++)
        {
            params.perf[nn] = 1.0;
        }

#ifdef USE_MPI
        status = MPI_Comm_rank(MPI_COMM_WORLD, (int_t *)&params.myid);
        status = MPI_Comm_size(MPI_COMM_WORLD, (int_t *)&params.nodes);
#else
        params.myid = 0;
        params.nodes = 1;

#endif /* USE_MPI */

        pfile.close();
    }
    else
    {
        status = ERR_FILE_NOT_FOUND;
    }

    return status;
}
