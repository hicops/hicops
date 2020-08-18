![HiCOPS C/C++ CI](https://github.com/mhaseeb123/hicops/workflows/HiCOPS%20C/C++%20CI/badge.svg)

# HiCOPS
HiCOPS: Software framework for accelerated peptide identification from LC-MS/MS data on HPC systems.

# Installation

## Recommended Compiler
Preferred: GCC compiler version 7.2.0 or later supporting C++14.

## Install the required packages
Install the following packages preferably using [Spack](https://spack.readthedocs.io). Read more about how to install Spack and how to install packages using Spack [here](https://spack.readthedocs.io/en/latest/getting_started.html). Make sure that you install all following packages using the same `GCC` compiler that you will use to install HiCOPS as well. See [Recommended Compiler](##Recommended_Compiler). The packages currently installed via Spack can be checked by:

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

`Note`: The package versions listed in the above list are not compulsory. You may install the latest versions of the packages.

### On a regular computer (skip if using XSEDE Comet)
Install the `mpich` package using `spack install mpich%gcc@version`

### On XSEDE Comet
Load the MPI and GNU modules
```bash
$ module purge
$ module load gnu/7.2.0
$ module load openmpi_ib
```

## Install timemory for HiCOPS instrumentation/profiling (optional)
Install timemory using CMake or Spack using the instructions [here](https://timemory.readthedocs.io/en/develop/installation.html). After installation, make sure that the path to timemory installation has been appended to the enviornment variable `CMAKE_PREFIX_PATH`.

```bash
$ export CMAKE_PREFIX_PATH=$timemory_DIR:$CMAKE_PREFIX_PATH
```

## Install HiCOPS
Install HiCOPS using Git & CMake using the following steps:

### Get HiCOPS

```bash
$ git clone https://github.com/pcdslab/hicops
$ cd hicops
```

### Configure

```bash
$ mkdir build & cd build
$ CC=$(which gcc) CXX=$(which g++) cmake .. [CMAKE_OPTIONS] -G [BUILD_SYSTEM] -D<VARIABLE>=<VALUE>
```

The available variables are as follows:

```bash
CMAKE_INSTALL_PREFIX: Set the installation path for HiCOPS, must be a writable directory and should not require sudo
CMAKE_BUILD_TYPE: Build type. Set to: Release, Debug, MinRelSize, RelWithDebInfo (default)
CMAKE_CXX_STANDARD: C++ standard. Set to: 11, 14 (default), 17
USE_OMP: Enable the use of OpenMP multithreading. Set to: ON(default)/OFF
USE_MPI: Enable MPI= support. Set to: ON(default)/OFF
USE_TIMEMORY: Enable timemory interface. Set to: ON/OFF(default). Requires timemory installation
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

# Run HiCOPS
For the rest of the document, we will be assuming that the HiCOPS was installed at : `$HICOPS_INSTALL`

## Update LD_LIBRARY_PATH
Append the `HICOPS_INSTALL/lib` to the environment variable `LD_LIBRARY_PATH`.

```bash
$ export LD_LIBRARY_PATH=$HICOPS_INSTALL/lib:$LD_LIBRARY_PATH
```

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

`Note:` Repeat Steps # 2 and 3 if you modify parameters in the `sampleparams.txt`.

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

`Note:` Repeat Steps # 2 and 3 if you modify parameters in the `sampleparams.txt`.

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

# Issue Reporting 
Open an issue on HiCOPS [GitHub](https://github.com/pcdslab/hicops/issues) (preferred) or email: `{mhase003, fsaeed}@fiu.edu`

# About this repository

## Authors
1. [Muhammad Haseeb](https://sites.google.com/fiu.edu/muhammadhaseeb)
2. [Fahad Saeed](https://saeedfahad.com)

## Contributors
1. [Muhammad Haseeb](https://github.com/mhaseeb123)
