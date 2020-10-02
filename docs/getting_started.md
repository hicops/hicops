---
title: Getting Started
---

# Getting Started
For this document, we will be assuming that the HiCOPS has been installed at : `$HICOPS_INSTALL`

## Setup Database


## Setup MS/MS dataset
HiCOPS currently only supports the `MS2` format for experimental MS/MS data. Please convert all experimental MS/MS data files into this format using the `raw2ms2` command line tool available with HiCOPS. Read more about the usage of `raw2ms2` tool [here]({{ site.baseurl }}/tools/ms2prep/raw2ms2).

## Setup instrumentation
Optional: If HiCOPS instrumentation was enabled during build, it can be configured and modified using the following environment variables. See how to enable HiCOPS instrumentation in the [Installation]({{ site.baseurl }}/installations) document:

| Variable                 | Description                                                                                                          |
|--------------------------|----------------------------------------------------------------------------------------------------------------------|
| `TIMEMORY_ENABLED`       | Enable/disable Timemory instrumentation interface. Set to : ON (default), OFF                                        |
| `HICOPS_MPIP_INSTR`      | Enable MPI data communication instrumentation. Set to: ON (default), OFF                                             |
| `HICOPS_INST_COMPONENTS` | Append instrumentation components. Set to: HICOPS_INST_COMPONENTS="c1,c2,.." where ci is a Timemory component  |
| `HICOPS_PAPI_EVENTS`     | Modify the hardware counters. Set to: HICOPS_PAPI_EVENTS="h1,h2,.." where hi is a PAPI counter                 |

To list all available timemory components [here](https://timemory.readthedocs.io/en/develop/tools/timemory-avail/README.html?highlight=user_bundle#available-components). By default, the following hardware counters are inserted into the `HICOPS_PAPI_EVENTS`.

```bash
HICOPS_PAPI_EVENTS="PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L3_TCM, \\
PAPI_L2_TCA, PAPI_L3_TCA, PAPI_MEM_WCY, PAPI_RES_STL, \\
PAPI_STL_CCY, PAPI_BR_CN, PAPI_BR_PRC, PAPI_FUL_ICY"
``` 

To see which hardware counters are available on your system and their description, use the `papi_avail` or `timemory-avail` tool. Refer to the PAPI documentation [here](https://icl.utk.edu/papi/) for more information. 

**NOTE:** If a PAPI counter is not available on the system but is added to the `HICOPS_PAPI_EVENTS` anyway, the profiler will not instrument any of the counters in the list regardless of their availability.

## Run HiCOPS on XSEDE Comet
Follow the instructions in [XSEDE]({{ site.baseurl }}/getting_started/xsede) document.

## Run HiCOPS
Follow the steps mentioned below. We will use the `hicops_config` tool to generate input parameters file, called `uparams.txt` for HiCOPS. Information on `hicops_comet` tool can be found [here]({{ site.baseurl }}/tools/runtime/hicops_config).

### Parameter Generation
1. Ensure that the hicops-core library path has been added to `LD_LIBRARY_PATH`.      

```bash
# append hicops-core lib path to LD_LIBRARY_PATH.
$ export LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH
```

2. Generate HiCOPS template runtime parameters file using the `hicops_config` tool located at `$HICOPS_INSTALL/bin`.     

```bash
$ $HICOPS_INSTALL/bin/hicops_config -g
Generated: ./sampleparams.txt

SUCCESS
```

3. Edit the generated sampleparams.txt file and setup HiCOPS' runtime parameters, database and data paths.     

4. Generate the HiCOPS runtime parameters file (`uparams.txt`) using `hicops_config` as:     

```bash
$ $HICOPS_INSTALL/bin/hicops_config ./sampleparams.txt
# uparams.txt will be generated
Generated: uparams.txt
```

**Note:** Repeat Steps 3-4 when you modify parameters in the `sampleparams.txt`.       

### Local Computer
1. Generate `uparams.txt` file using the steps [above](###Parameter-Generation).         

2. Run HiCOPS with `uparams.txt` as input argument with or without MPI depending on HiCOPS install [options]({{ site.baseurl }}/installation/##CMake-Options).       

```bash
# without using MPI
$ $HICOPS_INSTALL/bin/hicops $HICOPS_INSTALL/bin/uparams.txt

# using MPI
$ mpirun -np [P] [OPTIONS] $HICOPS_INSTALL/bin/hicops \\
  $HICOPS_INSTALL/bin/uparams.txt
```

**Note:**  Configure the mpirun options as follows: set binding level to `socket` and binding policy to `scatter`.

3. After HiCOPS execution is complete, run the `psm2excel` tool with `workspace` output directory (set in the sampleparams.txt file) as arguments.       

```bash
$ $HICOPS_INSTALL/tools/psm2excel [/path/to/hicops/workspace/output]
```

### SLURM
1. Generate `uparams.txt` file using the steps [above](###Parameter-Generation).         

2. Run HiCOPS with `uparams.txt` as input argument using SLURM with or without MPI depending on HiCOPS install [options]({{ site.baseurl }}/installation/##CMake-Options).        

```bash
# without using MPI
$ srun [OPTIONS] $HICOPS_INSTALL/bin/hicops $HICOPS_INSTALL/bin/uparams.txt

# using MPI
$ srun [OPTIONS] mpirun -np [P] [OPTIONS] $HICOPS_INSTALL/bin/hicops \\
  $HICOPS_INSTALL/bin/uparams.txt
```

**Note:** We highly recommend running HiCOPS through batch job submission via `sbatch` instead of `srun`. Make sure to follow the Hybrid MPI/OpenMP batch submission template when doing so.

**Note:**  Configure the mpirun options as follows: set binding level to `socket` and binding policy to `scatter`.

3. After HiCOPS execution, run the `psm2excel` tool using SLURM with `workspace` output directory (set in the sampleparams.txt file) as arguments.

```bash
$ srun [OPTIONS] --nodes=1 $HICOPS_INSTALL/tools/psm2excel -i \\
  [/path/to/hicops/workspace/output]
```