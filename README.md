![build](https://github.com/hicops/hicops/workflows/build/badge.svg) [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/hicops/hicops/blob/develop/README.md#contributing) [![pages yes](https://img.shields.io/badge/pages-yes-blue.svg)](https://hicops.github.io) [![GitHub forks](https://img.shields.io/github/forks/hicops/hicops.svg?style=social&label=Fork&maxAge=2592000)](https://GitHub.com/hicops/hicops/network/) [![GitHub stars](https://img.shields.io/github/stars/hicops/hicops.svg?style=social&label=Star&maxAge=2592000)](https://GitHub.com/hicops/hicops/stargazers/) [![GitHub contributors](https://img.shields.io/github/contributors/hicops/hicops.svg)](https://GitHub.com/hicops/hicops/graphs/contributors/) [![GitHub issues](https://img.shields.io/github/issues/hicops/hicops.svg)](https://GitHub.com/hicops/hicops/issues/) [![Github all releases](https://img.shields.io/github/downloads/hicops/hicops/total.svg)](https://GitHub.com/hicops/hicops/releases/)

# HiCOPS
*HiCOPS*: A computational framework for accelerated peptide identification from LC-MS/MS data on HPC systems.

This document contains a quick starter guide for HiCOPS. For full documentation, please refer to the HiCOPS [website](https://hicops.github.io).

# Installation

## Recommended Environment
* GCC compiler version 7.2.0 or later supporting C++14, OpenMP and threads. You may use Intel or LLVM compilers but make sure to follow through the requirements & installation steps accordingly. 
* Linux OS (Ubuntu v16.04, v18.04 and CentOS-7)
* MPI with threads support.
* Python 3.7+

## Install required packages
Install and load the following packages preferably using [Spack](https://spack.readthedocs.io). Read more about how to install Spack, and how to install and load packages using Spack [here](https://spack.readthedocs.io/en/latest/getting_started.html).

```bash
Required Packages:
boost@1.73.0       cmake@3.18.1        python@3.7.8              py-numpy@1.19.1         py-setuptools-scm@4.1.2     py-kiwisolver@1.1.0    py-python-dateutil@2.8.0
pkgconf@1.7.3      py-numexpr@2.7.0    py-setuptools@49.2.0      py-et-xmlfile@1.0.1     py-pillow@7.2.0             py-bottleneck@1.2.1    papi@6.0.0.1
py-jdcal@1.3       py-pyparsing@2.4.2  py-cython@0.29.21         py-pandas@1.1.0         py-subprocess32@3.5.4       py-cycler@0.10.0       py-openpyxl@3.0.3    
py-six@1.14.0      py-argparse@1.4.0   py-matplotlib@3.3.0       py-pytz@2019.3
```

**NOTE**: The package versions listed in the above list are not compulsory. You may install the latest versions of the packages. Make sure that you installed all the above listed packages using the same compiler that you will use to install HiCOPS as well. 


## Install Timemory (Instrumentation)
Optionally install timemory using CMake or Spack using the instructions [here](https://timemory.readthedocs.io/en/develop/installation.html) to enable HiCOPS instrumentation. After installation, make sure that the path to timemory installation has been appended to the enviornment variable `CMAKE_PREFIX_PATH`. Also make sure that the Timemory's python package (Pytimemory) is in your `PYTHONPATH`. Also make sure that the environment variable `MPLCONFIGDIR` is set to a writeable directory. e.g. `$HOME/mplconfigdir`

```bash
$ export CMAKE_PREFIX_PATH=$TIMEMORY_INSTALL:$CMAKE_PREFIX_PATH
$ export PYTHONPATH=$PYTIMEMORYPATH:$PYTHONPATH
$ export MPLCONFIGDIR=$HOME/mplconfigdir
```

### Timemory Install Example
If using Spack, you can install and load timemory and its dependencies using:

**Using Spack**
```bash
$ spack install timemory%gcc@version +ompt +tools +ompt_library ~dyninst +gotcha +python +papi ~caliper +mpi +mpip_library
$ spack load -r timemory
```
Add the `--only dependencies` flag right after `spack install` in above command if you want to only install the dependencies of Timemory.

**Using CMake**
```bash
$ spack --only dependencies install timemory%gcc@version +ompt +tools +ompt_library ~dyninst +gotcha +python +papi ~caliper +mpi +mpip_library
$ git clone https://github.com/NERSC/timemory.git && cd timemory && mkdir build && cd build
$ cmake .. -DTIMEMORY_USE_MPI=ON -DTIMEMORY_BUILD_MPIP_LIBRARY=ON -DTIMEMORY_USE_OMPT=ON -DTIMEMORY_USE_GOTCHA=ON -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_CXX_STANDARD=14 -DTIMEMORY_USE_PYTHON=ON -DTIMEMORY_BUILD_TOOLS=ON -DUSE_MPI=ON -DUSE_OPENMP=ON -DTIMEMORY_USE_PAPI=ON
$ make install
$ export CMAKE_PREFIX_PATH=$PWD/../install:$CMAKE_PREFIX_PATH
$ export PYTHONPATH=$PWD:$PYTHONPATH
```

**NOTE:** Timemory and its dependencies can take upto 2-3 hours to install depending on the system so please be patient.

## Install HiCOPS
Install HiCOPS using Git & CMake using the following steps:

```bash
$ git clone https://github.com/hicops/hicops
$ cd hicops
$ mkdir build && cd build
$ CC=$(which gcc) CXX=$(which g++) cmake .. [CMAKE_OPTIONS] -G [BUILD_SYSTEM] [HICOPS_OPTIONS]
$ make install -j [JOBS]
```

Available HiCOPS options:

```bash
USE_MPI                 Enable MPI support. Set to: ON (default), OFF
USE_TIMEMORY            Enable timemory interface. Set to: ON, OFF (default) => Requires timemory installation.
USE_MPIP_LIBRARY        Enables the MPIP data_tracker via Timemory. Set to: ON, OFF (default) => Requires timemory installation. 
TAILFIT                 Use the tailfit method instead of Gumbelfit for e-value computation. Set to: ON (default), OFF
PROGRESS                Display progress marks. Set to: ON (default), OFF
MAX_SEQ_LEN             Allowed maximum peptide sequence length. Set to: 7 to 60. Default: 60
QALEN                   Maximum number of top K peaks to keep when spectrum preprocess. Default: 100
QCHUNK                  Max size of each batch extracted from the dataset. Default: 10000
MAX_HYPERSCORE          Maximum allowed hyperscore computed. Default: 100
```

Available CMake options:

```bash
CMAKE_INSTALL_PREFIX    Set the installation path for HiCOPS, must be a writable directory without sudo
CMAKE_BUILD_TYPE        Build type. Set to: Release, Debug, MinSizeRel, RelWithDebInfo (default)
CMAKE_CXX_STANDARD      C++ standard. Set to: 11, 14 (default), 17
```

**NOTE:** Compiling HiCOPS with Timemory enabled may take some time (~ 3-5 minutes)

# Run HiCOPS
We will be assuming that the HiCOPS was installed at : `$HICOPS_INSTALL`. Further, make sure that the `HICOPS_INSTALL/lib` has been added to the environment variable `LD_LIBRARY_PATH`.

```bash
$ export LD_LIBRARY_PATH=$HICOPS_INSTALL/lib:$LD_LIBRARY_PATH
```

## Setup Experimental Data, Database and Instrumentation
Please follow through the instructions detailed [here](https://hicops.github.io/getting_started#setup) to prepare experimental data and database, and setup instrumentation.

## Running HiCOPS (Anywhere)
1. Generate HiCOPS sample runtime parameters file using the `hicops_config` located at `$HICOPS_INSTALL/bin`.

```bash
$ $HICOPS_INSTALL/bin/hicops_config -g
Generated: ./sampleparams.txt

SUCCESS
```

2. Edit the generated sampleparams.txt file and add/modify HiCOPS' runtime parameters.
3. Generate the HiCOPS runtime parameters file (called `uparams.txt`) as:

```bash
$ $HICOPS_INSTALL/bin/hicops_config ./sampleparams.txt
Generated: uparams.txt
```

4. Run HiCOPS with `uparams.txt` as arguments and optionally with MPI.

```bash
$ mpirun -np 4 [OPTIONS] $HICOPS_INSTALL/bin/hicops $HICOPS_INSTALL/bin/uparams.txt
```

5. Run the `psm2excel` tool and pass the HiCOPS workspace output directory (that was set in the sampleparams.txt file) as parameter.

```bash
$ $HICOPS_INSTALL/tools/psm2excel [/path/to/hicops/workspace/output]
```

**NOTE:** Repeat Steps # 2 to 5 if you modify the `sampleparams.txt`.

## Running HiCOPS (on XSEDE Comet)
1. Generate HiCOPS sample runtime parameters file using the `hicops_comet` tool located at `$HICOPS_INSTALL/bin/tools`.

```bash
$ $HICOPS_INSTALL/bin/tools/hicops_comet -g
Generated: ./sampleparams.txt

SUCCESS
```

2. Edit the generated sampleparams.txt file and add/modify HiCOPS' runtime parameters.

3. Run HiCOPS using the same tool i.e. `hicops_comet`, however, this time providing the updated sampleparams.txt as parameter to the tool.

```bash
$ $HICOPS_INSTALL/bin/tools/hicops_comet sampleparams.txt
```

4. Run the `psm2excel` tool using SLURM and pass the HiCOPS workspace output directory (that was set in the sampleparams.txt file) as parameters.

```bash
$ srun --partition=compute  --nodes=1 --ntasks-per-node=1 -t 00:15:00 --export=ALL $HICOPS_INSTALL/tools/psm2excel -i [/path/to/hicops/workspace/output]
```

**NOTE:** Repeat Steps # 2 to 4 if you modify the `sampleparams.txt`.

# Precautions
Please read and follow the following precautions to avoid any errors.

* Always use a unique workspace directory for each experiment, specially for the simultaneously running HiCOPS instances to avoid overwriting intermediate results and other errors.
* Always convert the TSV results into Excel file using `psm2excel` (last step) before using the same workspace folder for another experiment. It is always better to use a unique workspace folder for each experimental run.
* Do not run too many simultaneous HiCOPS instances with large number of nodes allocated to each instance to avoid I/O bandwidth contention and thus, performance degradation.
* Do not modify the generated files such as uparams.txt manually and instead regenerate using the relevant tool.
* Avoid using relative paths in the `sampleparams.txt` file to avoid any errors.

# Integrating Code with HiCOPS
Please refer to the documentation [here](https://hicops.github.io/getting_started/integrate).

# About this repository

## Credits
1. [Muhammad Haseeb](https://sites.google.com/fiu.edu/muhammadhaseeb)
2. [Fahad Saeed](http://www.saeedfahad.com)

## Code Maintainers
1. [Muhammad Haseeb](https://github.com/mhaseeb123)

## Issues and Feature Requests
Open an issue [here](https://github.com/hicops/hicops/issues).

Please include any logs, screenshots and/or helpful observations. Also, do not forget to describe the dataset(s), database, steps performed etc. so that the issue can be reproduced.

## Contributing
All contributions are welcome including new features, documentation and bug fixes through standard GitHub pull request method. Generic guidelines:

1. Fork this repository to your local GitHub, checkout a new branch from the `develop` branch.
2. Make your changes/updates.
3. Make sure that you pull the latest changes from `hicops:develop` into your branch and merge before committing your changes.
4. Commit your changes and push your branch to `origin`. i.e. `your_hicops_fork`.
5. Open a pull request (PR) from `your_hicops_fork:new_branch` to `hicops:develop`.
