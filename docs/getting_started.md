---
title: Getting Started
---

# Getting Started

Follow the below steps to get started with HiCOPS:

<!-- TOC -->

- [Setup](#setup)
  - [Setup Database](#setup-database)
  - [Setup MS/MS Dataset](#setup-msms-dataset)
  - [Setup Instrumentation](#setup-instrumentation)
- [Run HiCOPS](#run-hicops)
  - [XSEDE Comet](#xsede-comet)
  - [Any Other System](#any-other-system)
    - [Generate Params](#generate-params)
    - [Execute](#execute)
- [Precautions](#precautions)

<!-- /TOC -->

## Setup
Setup the peptide database, experimental MS/MS dataset and HiCOPS instrumentation using the below instructions.

### Setup Database
Get the desired protein sequence database from UniProt/Swissprot. Digest the protein sequence database into a peptide sequence database using Digestor tool available with [OpenMS](https://www.openms.de/) or using [ProteoWizard](http://proteowizard.sourceforge.net/). Make sure that the generated peptide sequence database is in FASTA format. Use the `db_prep` tool to separate coarse-grained peptide sequence clusters. This tool will generate many files in `./parts/len.pep` directory. Read more about the usage of `db_prep` tool [here]({{ site.baseurl }}/tools/dbtools/dbprep).

### Setup MS/MS Dataset
HiCOPS currently only supports the `MS2` format for experimental MS/MS data. Please convert all experimental MS/MS data files into this format using the `raw2ms2` command line tool available with HiCOPS. Read more about the usage of `raw2ms2` tool [here]({{ site.baseurl }}/tools/ms2prep/raw2ms2).

### Setup Instrumentation
Optional: If HiCOPS instrumentation was enabled during build, it can be configured and modified using the following environment variables. See how to enable HiCOPS instrumentation in the [Install]({{ site.baseurl }}/installation) document.

| Variable                 | Description                                                                                  |
|--------------------------|----------------------------------------------------------------------------------------------|
| `TIMEMORY_ENABLED`       | Enable/disable Timemory instrumentation interface. Set to : ON (default), OFF                |
| `HICOPS_MPIP_INSTR`      | Enable MPI data communication instrumentation. Set to: ON (default), OFF                     |
| `HICOPS_INST_COMPONENTS` | Add more inst components. Set as: `HICOPS_INST_COMPONENTS="ci,.."` where `ci`= Timemory comp |
| `HICOPS_PAPI_EVENTS`     | Modify hardware counters. Set as: `HICOPS_PAPI_EVENTS="hi,.."` where `hi`= PAPI counter      |

To list all available timemory components [here](https://timemory.readthedocs.io/en/develop/tools/timemory-avail/README.html?highlight=user_bundle#available-components). By default, the following hardware counters are inserted into the `HICOPS_PAPI_EVENTS`.

```bash
HICOPS_PAPI_EVENTS="PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L3_TCM, \\
PAPI_L2_TCA, PAPI_L3_TCA, PAPI_MEM_WCY, PAPI_RES_STL, \\
PAPI_STL_CCY, PAPI_BR_CN, PAPI_BR_PRC, PAPI_FUL_ICY"
``` 

To see which hardware counters are available on your system and their description, use the `papi_avail` or `timemory-avail` tool. Refer to the PAPI documentation [here](https://icl.utk.edu/papi/) for more information. 

**NOTE:** If a PAPI counter is not available on the system but is added to the `HICOPS_PAPI_EVENTS` anyway, the profiler will not instrument any of the counters in the list regardless of their availability.

## Run HiCOPS
Follow the instructions relevant to your compute environment to seamlessly run HiCOPS. We categorize the compute environments into two categories

### XSEDE Comet
If you are running on XSEDE Comet environment, skip the rest of this document and follow the instructions [here]({{ site.baseurl }}/getting_started/xsede).

### Any Other System
Follow the below instructions if you are running on any system but XSEDE Comet.

#### Generate Params
**i.** Ensure that the hicops-core library path is added to `LD_LIBRARY_PATH`.      

```bash
# append hicops-core lib path to LD_LIBRARY_PATH.
$ export LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH
```

**ii.** Generate HiCOPS template runtime parameters file using the `hicops_config` tool located at `$HICOPS_INSTALL/bin`.     

```bash
# run hicops_comet with -g
$ $HICOPS_INSTALL/bin/hicops_config -g
# generated: ./sampleparams.txt
```

**iii.** Edit the generated sampleparams.txt file and setup HiCOPS' runtime parameters, database and data paths.     

**iv.** Generate the HiCOPS runtime parameters file (`uparams.txt`) using `hicops_config` as:     

```bash
# run hicops_comet with sampleparams.txt
$ $HICOPS_INSTALL/bin/hicops_config ./sampleparams.txt

# generated: uparams.txt
```

#### Execute
**v.** Run HiCOPS with `uparams.txt` as input argument with or without MPI depending on HiCOPS install [options]({{ site.baseurl }}/installation#cmake-options). Use the resource manager (SLURM, LSH) if working on a managed cluster system.       

**Note:**  Configure the mpirun options as follows: set binding level to `socket` and binding policy to `scatter`.     

**Note:** We highly recommend running HiCOPS through batch job submission `sbatch` instead of `srun`. Make sure to follow the relevant Hybrid MPI/OpenMP batch submission template when doing so.     

```bash
# without MPI
$ $HICOPS_INSTALL/bin/hicops $HICOPS_INSTALL/bin/uparams.txt

# SLURM without MPI
$ srun [OPTIONS] $HICOPS_INSTALL/bin/hicops $HICOPS_INSTALL/bin/uparams.txt

# with MPI
$ mpirun -np [N] [OPTIONS] $HICOPS_INSTALL/bin/hicops \\
  $HICOPS_INSTALL/bin/uparams.txt

# SLURM with MPI
$ srun [OPTIONS] mpirun -np [N] [OPTIONS] $HICOPS_INSTALL/bin/hicops \\
  $HICOPS_INSTALL/bin/uparams.txt
```    

**vi.** After HiCOPS execution is complete, run the `psm2excel` tool with `workspace` output directory (set in the sampleparams.txt file) as arguments.       

```bash
# psm2excel
$ $HICOPS_INSTALL/tools/psm2excel [/path/to/hicops/workspace/output]

# psm2excel with SLURM
$ srun [OPTIONS] --nodes=1 $HICOPS_INSTALL/tools/psm2excel -i \\
  [/path/to/hicops/workspace/output]
```    

**vii.** Repeat Steps **iii.** to **vi.** when you modify parameters in the `sampleparams.txt`.

## Precautions
Please read and follow the following precautions to avoid any errors.

* Always use a unique workspace directory for each experiment, specially for the simultaneously running HiCOPS instances to avoid overwriting intermediate results and other errors.      
* Always convert the TSV results into Excel file using `psm2excel` (last step) before using the same workspace folder for another experiment. It is always better to use a unique workspace folder for each experimental run.     
* Do not run too many simultaneous HiCOPS instances with large number of nodes allocated to each instance to avoid I/O bandwidth contention and thus, performance degradation.      
* Do not modify the generated files such as `uparams.txt` manually and instead re-generate using the relevant runtime tool.    
* Avoid using relative paths in the sampleparams.txt file to avoid any errors.    