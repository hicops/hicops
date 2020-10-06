---
title: Install
---

# Install

HiCOPS uses a standard CMake installation. Set the environment variables `CC` and `CXX` before executing CMake to specify the compilers. By default, HiCOPS enables `OpenMP` and `MPI` and disables instrumentation support. HiCOPS will search the `CMAKE_PREFIX_PATH` for Timemory installation.

## Environment
* GCC 7.2+ compiler with C++14, OpenMP and threading
* MPI
* Python 3.7+
* CMake 3.11+
* Linux based OS such as Ubuntu 16.04+, CentOS7+.

**Important**: Please make sure that your MPI distribution supports multiple threads (MPI thread multiple) option.

## Required Packages
Install all required packages by following the instructions [here]({{ site.baseurl }}/installation/packages#required-packages).

## Instrumentation Support
Optional: Install Timemory to enable instrumentation by following the instructions [here]({{ site.baseurl }}/installation/packages#optional-packages).

## Install HiCOPS
Install HiCOPS using Git & CMake using the following steps:

```bash
# clone
$ git clone https://github.com/mhaseeb123/hicops
$ cd hicops
$ mkdir build && cd build

# configure with CMake
$ CC=$(which gcc) CXX=$(which g++) cmake .. [CMAKE_OPTIONS] -G [BUILD_SYSTEM] [HICOPS_OPTIONS]

# build and install
$ make install -j [JOBS]

# append  hicops-core lib path to LD_LIBRARY_PATH.
$ export LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH
```

**Note:** Compiling HiCOPS with Timemory enabled may take some time (~ 3-5 minutes)

## CMake Options

The available CMake options for HiCOPS include:

| Option             | Description                                                                                             |
|--------------------|---------------------------------------------------------------------------------------------------------|
| `USE_MPI`          | Enable MPI support. Set to: ON (default), OFF                                                           |
| `USE_TIMEMORY`     | Enable timemory interface. Set to: ON, OFF (default) => Requires timemory installation.                 |
| `USE_MPIP_LIBRARY` | Enables the MPIP data_tracker via Timemory. Set to: ON, OFF (default) => Requires timemory installation |
| `TAILFIT`          | Use the tailfit method instead of Gumbelfit for e-value computation. Set to: ON (default), OFF          |
| `PROGRESS`         | Display  HiCOPS progress marks. Set to: ON (default), OFF                                               |
| `MAX_SEQ_LEN`      | Allowed maximum peptide sequence length. Set to: 7 to 60. Default: 60                                   |
| `QALEN`            | Maximum number of top K peaks to keep when spectrum preprocess. Default: 100                            |
| `QCHUNK`           | Max size of each batch extracted from the dataset. Default: 10000                                       |
| `MAX_HYPERSCORE`   | Maximum allowed hyperscore computed. Default: 100                                                       |

Available CMake build control options include:

| Option                 | Description                                                                     |
|------------------------|---------------------------------------------------------------------------------|
| `CMAKE_INSTALL_PREFIX` | Set the installation path for HiCOPS, must be a writable directory without sudo |
| `CMAKE_BUILD_TYPE`     | Build type. Set to: Release, Debug, MinSizeRel, RelWithDebInfo (default)        |
| `CMAKE_CXX_STANDARD`   | C++ standard. Set to: 14 (default), 17                                          |