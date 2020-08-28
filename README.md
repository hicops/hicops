![HiCOPS C/C++ CI](https://github.com/mhaseeb123/hicops/workflows/HiCOPS%20C/C++%20CI/badge.svg)

# HiCOPS
*HiCOPS*: A computational framework for accelerated peptide identification from LC-MS/MS data on HPC systems.

# Installation

## Recommended Compiler
GCC compiler version 7.2.0 or later supporting C++14. You may use Intel or LLVM compilers but make sure to follow through the installation steps accordingly. We have tested the HiCOPS on Linux OS (Ubuntu v16.04, v18.04 and CentOS-7) with GCC v7.2.0, v8.4.0 and v9.3.0 running on Haswell, Broadwell, Kabylake and Skylake processors.

## Install and Load the required packages
Install and load the following packages preferably using [Spack](https://spack.readthedocs.io). Read more about how to install Spack, and how to install and load packages using Spack [here](https://spack.readthedocs.io/en/latest/getting_started.html).

```bash
Required Packages:
boost@1.73.0       cmake@3.18.1        python@3.7.8              py-numpy@1.19.1         py-setuptools-scm@4.1.2     py-kiwisolver@1.1.0    py-python-dateutil@2.8.0
pkgconf@1.7.3      py-numexpr@2.7.0    py-setuptools@49.2.0      py-et-xmlfile@1.0.1     py-pillow@7.2.0             py-bottleneck@1.2.1    papi@6.0.0.1
py-jdcal@1.3       py-pyparsing@2.4.2  py-cython@0.29.21         py-pandas@1.1.0         py-subprocess32@3.5.4       py-cycler@0.10.0       py-openpyxl@3.0.3    
py-six@1.14.0      py-argparse@1.4.0   py-matplotlib@3.3.0       py-pytz@2019.3
```

**NOTE**: The package versions listed in the above list are not compulsory. You may install the latest versions of the packages. 

Make sure that you installed all the above listed packages using the same compiler that you will use to install HiCOPS as well. See [Recommended Compiler](##Recommended-Compiler). The packages currently installed via Spack can be checked using `spack find`. For example:

```bash
$ spack find
==> 63 installed packages
-- linux-centos7-haswell / gcc@7.2.0 ----------------------------
autoconf@2.69        elfutils@0.180   intel-tbb@2020.3     libsigsegv@2.12  papi@6.0.0.1         py-kiwisolver@1.1.0  py-python-dateutil@2.8.0  readline@8.0
automake@1.16.2      expat@2.2.9      libbsd@0.10.0        libtool@2.4.6    perl@5.30.3          py-matplotlib@3.3.0  py-pytz@2019.3            sqlite@3.31.1
berkeley-db@18.1.40  findutils@4.6.0  libffi@3.3           libxml2@2.9.10   pkgconf@1.7.3        py-numexpr@2.7.0     py-setuptools@49.2.0      tar@1.32
boost@1.73.0         freetype@2.10.1  libiberty@2.33.1     m4@1.4.18        py-bottleneck@1.2.1  py-numpy@1.19.1      py-setuptools-scm@4.1.2   texinfo@6.5
bzip2@1.0.8          gdbm@1.18.1      libiconv@1.16        nasm@2.14.02     py-cycler@0.10.0     py-openpyxl@3.0.3    py-six@1.14.0             util-macros@1.19.1
cmake@3.18.1         gettext@0.20.2   libjpeg-turbo@2.0.4  ncurses@6.2      py-cython@0.29.21    py-pandas@1.1.0      py-subprocess32@3.5.4     xz@5.2.5
diffutils@3.7        gotcha@1.0.3     libpciaccess@0.13.5  openblas@0.3.10  py-et-xmlfile@1.0.1  py-pillow@7.2.0      python@3.7.8              zlib@1.2.11
dyninst@10.2.0       hwloc@2.2.0      libpng@1.6.37        openssl@1.1.1g   py-jdcal@1.3         py-pyparsing@2.4.2   qhull@2020.1
```

### On a regular computer (skip if using XSEDE Comet)
Install the `mpich` package using `spack install mpich%gcc@version`

### On XSEDE Comet
Load the MPI and GNU modules
```bash
$ module purge
$ module load gnu/7.2.0
$ module load openmpi_ib
```

## Install timemory for HiCOPS instrumentation/profiling - Optional
Install timemory using CMake or Spack using the instructions [here](https://timemory.readthedocs.io/en/develop/installation.html). After installation, make sure that the path to timemory installation has been appended to the enviornment variable `CMAKE_PREFIX_PATH`. Also make sure that the Timemory's python package (Pytimemory) is in your `PYTHONPATH`. Also make sure that the environment variable `MPLCONFIGDIR` is set to a writeable directory. e.g. `$HOME/mplconfigdir`

```bash
$ export CMAKE_PREFIX_PATH=$TIMEMORY_INSTALL:$CMAKE_PREFIX_PATH
$ export PYTHONPATH=$PYTIMEMORYPATH:$PYTHONPATH
$ export MPLCONFIGDIR=$HOME/mplconfigdir
```

### Timemory Install Example
All Timemory installation examples can be found [here](https://github.com/NERSC/timemory/wiki/Installation-Examples). Below we demonstrate a couple of them.

#### Using Spack
If using Spack, you can install and load timemory and its dependencies using:

```bash
$ spack install timemory%gcc@version +ompt +tools +ompt_library ~dyninst +gotcha +python +papi ~caliper +mpi +mpip_library
$ spack load -r timemory
```
Add the `--only dependencies` flag right after `spack install` in above command if you want to only install the dependencies of Timemory.

#### Using CMake
If using CMake, then assuming all that all required Timemory dependency packages have been installed and loaded:

```bash
$ git clone https://github.com/NERSC/timemory.git && cd timemory && mkdir build && cd build
$ cmake .. -DTIMEMORY_USE_MPI=ON -DTIMEMORY_BUILD_MPIP_LIBRARY=ON -DTIMEMORY_USE_OMPT=ON -DTIMEMORY_USE_GOTCHA=ON -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_CXX_STANDARD=14 -DTIMEMORY_USE_PYTHON=ON -DTIMEMORY_BUILD_TOOLS=ON -DUSE_MPI=ON -DUSE_OPENMP=ON -DTIMEMORY_USE_PAPI=ON
$ make install
$ export CMAKE_PREFIX_PATH=$PWD/../install:$CMAKE_PREFIX_PATH
$ export PYTHONPATH=$PWD:$PYTHONPATH
```

**NOTE:** Timemory and its dependencies can take upto 2-3 hours to install depending on the system so please be patient.

## Install HiCOPS
Install HiCOPS using Git & CMake using the following steps:

### Configure

```bash
$ git clone https://github.com/pcdslab/hicops
$ cd hicops
$ mkdir build && cd build
$ CC=$(which gcc) CXX=$(which g++) cmake .. [CMAKE_OPTIONS] -G [BUILD_SYSTEM] [HICOPS_OPTIONS]
```

Available HiCOPS options:

```bash
USE_OMP                 Enable the use of OpenMP multithreading. Set to: ON(default), OFF
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

You can check the value of each setting in your build using `ccmake`. For example: 

```bash
$ ccmake ..

 Boost_NO_BOOST_CMAKE             ON
 CMAKE_BUILD_TYPE                 Release
 CMAKE_ENABLE_EXPORTS             ON
 CMAKE_INSTALL_PREFIX             /home/mhaseeb/repos/hicops/install
 CMAKE_INSTALL_RPATH_USE_LINK_P   ON
 CMAKE_POSITION_INDEPENDENT_COD   ON
 USE_MPI                          ON
 USE_OMP                          ON
 USE_TIMEMORY                     ON
 timemory_DIR                     /home/mhaseeb/repos/timemory/install/share/cmake/timemory
 timemory_ONETIME_MESSAGE_DELIV   ON
```

## Build and Install

Depending on the build system generated, build and install HiCOPS. For example in case of Makefile

```bash
$ make install -j [JOBS]
```

**NOTE:** Compiling HiCOPS with Timemory enabled may take some time (~ 3-5 minutes)

# Run HiCOPS
For the rest of the document, we will be assuming that the HiCOPS was installed at : `$HICOPS_INSTALL`

## Update LD_LIBRARY_PATH
Make sure that the `HICOPS_INSTALL/lib` has been added to the environment variable `LD_LIBRARY_PATH`.

```bash
$ export LD_LIBRARY_PATH=$HICOPS_INSTALL/lib:$LD_LIBRARY_PATH
```

## Setup Instrumentation (with Timemory) - Optional
If the `USE_TIMEMORY=ON` option was set in [Configure](###Configure) step, the HiCOPS instrumentation can be configured and updated via the following environment variables:

```bash
HICOPS_MPIP_INSTR           Enable MPI data communication instrumentation. Set to: ON (default), OFF
HICOPS_INST_COMPONENTS      Append to the list of Timemory components (metrics) used for instrumenting the HiCOPS parallel search algorithm. 
                            Set to: HICOPS_INST_COMPONENTS="<c1>,<c2>,.." where each <ci> is a Timemory component.
HICOPS_PAPI_EVENTS          Modify (not append) the vector of PAPI hardware counters used for instrumenting the HiCOPS parallel search algorithm.
                            Set to: HICOPS_PAPI_EVENTS="<hw1>, <hw2>,.." where each <hwi> is a PAPI hardware counter.
```

See more about how to list available timemory components [here](https://timemory.readthedocs.io/en/develop/tools/timemory-avail/README.html?highlight=user_bundle#available-components). 

To see which hardware counters are available on your system and their description, run the `papi_avail` tool or refer to the documentation [here](https://icl.utk.edu/papi/). By default, the following hardware counters are inserted into the `HICOPS_PAPI_EVENTS`.

```bash
HICOPS_PAPI_EVENTS="PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L3_TCM, PAPI_L2_TCA, PAPI_L3_TCA, PAPI_MEM_WCY, PAPI_RES_STL, PAPI_STL_CCY, PAPI_BR_CN, PAPI_BR_PRC, PAPI_FUL_ICY"
```

**NOTE:** If a PAPI counter is not available on the system but is added to the `HICOPS_PAPI_EVENTS` anyway, the profiler will not instrument any of the counters in the list regardless of their availability.

## On a regular computer (skip if using XSEDE Comet)
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

Run HiCOPS with `uparams.txt` as arguments and optionally MPI if `USE_MPI=ON` was set in [Configure](###Configure) stage.

```bash
$ mpirun -np 4 [OPTIONS] $HICOPS_INSTALL/bin/hicops $HICOPS_INSTALL/bin/uparams.txt
```

**NOTE:** Repeat Steps # 2 and 3 if you modify parameters in the `sampleparams.txt`.

## On XSEDE Comet
1. Generate HiCOPS sample runtime parameters file using the `hicops_comet` wrapper script located at `$HICOPS_INSTALL/bin/wrappers`.

```bash
$ $HICOPS_INSTALL/bin/wrappers/hicops_comet -g
Generated: ./sampleparams.txt

SUCCESS
```

2. Edit the generated sampleparams.txt file and add/modify HiCOPS' runtime parameters.

3. Run HiCOPS using the same wrapper script i.e. `hicops_comet`, however, this time providing the updated sampleparams.txt as parameter to the wrapper.

```bash
$ $HICOPS_INSTALL/bin/wrappers/hicops_comet sampleparams.txt
```

**NOTE:** Repeat Steps # 2 and 3 if you modify parameters in the `sampleparams.txt`.

# Post-processing HiCOPS output
HiCOPS generates PSM data in partial TSV files that can be merged using the `psm2excel` tool located at: `$HICOPS_INSTALL/wrappers`. The tool generates a combined Excel file called `Concat.xlsx` containing the final PSM data (no-FDR).

## On a regular computer (skip if using XSEDE Comet)
Run the `psm2excel` tool and pass the HiCOPS workspace output directory (that was set in the sampleparams.txt file) as parameter.

```bash
$ $HICOPS_INSTALL/wrappers/psm2excel [/path/to/hicops/workspace/output]
```

## On XSEDE Comet
Run the `psm2excel` tool using SLURM and pass the HiCOPS workspace output directory (that was set in the sampleparams.txt file) as parameters.

```bash
$ srun --partition=shared  --nodes=1 --ntasks-per-node=1 -t 00:10:00 --export=ALL $HICOPS_INSTALL/wrappers/psm2excel [/path/to/hicops/workspace/output]
```

# About this repository

## Authors
1. [Muhammad Haseeb](https://sites.google.com/fiu.edu/muhammadhaseeb)
2. [Fahad Saeed](http://www.saeedfahad.com)

## Code Maintainers
1. [Muhammad Haseeb](https://github.com/mhaseeb123)

## Issue Reporting
Open an issue [here](https://github.com/pcdslab/hicops/issues) (preferred) or email: mhase003@fiu.edu or fsaeed@fiu.edu. Use the template subject line: \[BUG:HiCOPS]: *One line summary of the issue*.  

Please include any logs, screenshots and/or helpful observations. Also, do not forget to describe the dataset(s), database, steps performed etc. so that the issue can be reproduced.

## Contributing

All contributions are welcome including new features, documentation updates and bug fixes.

1. Fork this repository to your local GitHub, checkout a new branch from the `develop` branch.
2. Make your changes/updates.
3. Make sure that you pull the latest changes from `hicops:develop` into your branch and merge before committing your changes.
4. Commit your changes and push your branch to `origin`. i.e. `your_hicops_fork`.
5. Open a pull request (PR) from `your_hicops_fork:new_branch` to `hicops:develop`.
