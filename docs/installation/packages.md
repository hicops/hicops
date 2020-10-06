---
title: Required Packages
---

# Required Packages

Install and load the latest versions of the following packages, preferably through [Spack](https://spack.readthedocs.io). Read more about how to install Spack, and how to install and load packages using Spack [here](https://spack.readthedocs.io/en/latest/getting_started.html).

| Required Packages  |                   |               |
|--------------------|-------------------|---------------|
| boost              | cmake             | python        |
| py-numpy           | py-setuptools-scm | py-kiwisolver |
| py-python-dateutil | pkgconf           | py-numexpr    |
| py-setuptools      | py-et-xmlfile     | py-pillow     |
| py-bottleneck      | papi              | py-jdcal      |
| py-pyparsing       | py-cython         | py-pandas     |
| py-subprocess32    | py-cycler         | py-openpyxl   |
| py-six             | py-argparse       | py-matplotlib |
| py-pytz            | mpich             | python        |

**Note:** You may skip MPI installation if you are running on a cluster system with built-in MPI distribution and simply load that into your environment.

## Optional Packages
Install the optional packages, preferably using Spack.

| Optional Packages |          |
|-------------------|----------|
| papi              | timemory |

### Timemory Install
Install timemory with Python, MPI, MPIP, OMPT and PAPI support using the following instructions:

```bash
# install timemory dependencies via Spack
$ spack install --only dependencies timemory%gcc@<version> +ompt +tools \\
  +ompt_library ~dyninst +gotcha +python +papi \\
  ~caliper +mpi +mpip_library

# load dependencies
$ spack load --only dependencies timemory

# clone timemory
$ git clone https://github.com/NERSC/timemory.git && cd timemory && mkdir build && cd build

# configure using CMake
$ cmake .. -DTIMEMORY_USE_MPI=ON -DTIMEMORY_BUILD_MPIP_LIBRARY=ON -DTIMEMORY_USE_OMPT=ON -DTIMEMORY_USE_GOTCHA=ON -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_CXX_STANDARD=14 -DTIMEMORY_USE_PYTHON=ON -DTIMEMORY_BUILD_TOOLS=ON -DUSE_MPI=ON -DUSE_OPENMP=ON -DTIMEMORY_USE_PAPI=ON

# build timemory
$ make install

# export required environment variables
$ export CMAKE_PREFIX_PATH=$PWD/../install:$CMAKE_PREFIX_PATH
$ export PYTHONPATH=$PWD:$PYTHONPATH
$ export MPLCONFIGDIR=$HOME/mplconfigdir
```
**NOTE:** Timemory and its dependencies can take upto 2-3 hours to install depending on the system so please be patient.

Ensure that the timemory installation path has correctly been appended to `CMAKE_PREFIX_PATH`. Also make sure that the Timemory's python package (PyTimemory) is in `PYTHONPATH`. Set `MPLCONFIGDIR` to a writeable directory. e.g. `$HOME/mplconfigdir`

#### More Information
See Timemory installation instructions [here](https://timemory.readthedocs.io/en/develop/installation.html). More Timemory installation examples can be found [here](https://github.com/NERSC/timemory/wiki/Installation-Examples).
