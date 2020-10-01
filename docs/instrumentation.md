# Instrumentation
The instrumentation and performance analysis is provided via Timemory toolkit. HiCOPS incorporates instrumentation using C++14 `extern templates` which allow custom instantiation of required instrumentation metrics. Please refer to [Timemory](https://timemory.readthedocs.io/en/develop/about.html) documentation for more information.

## Provided Interfaces
HiCOPS provides instrumentation interfaces using Timemory bundles which are C++ tuples containing instrumentation metrics. Time-only measurements are performed via `MACROS` which direct instrumentation to HiCOPS' native instrumentation or Timemory depending on the configuration.

## Install Timemory
Install timemory using CMake or Spack using the instructions [here](https://timemory.readthedocs.io/en/develop/installation.html). After installation, make sure that the path to timemory installation has been appended to the enviornment variable `CMAKE_PREFIX_PATH`. Also make sure that the Timemory's python package (Pytimemory) is in your `PYTHONPATH`. Also make sure that the environment variable `MPLCONFIGDIR` is set to a writeable directory. e.g. `$HOME/mplconfigdir`

```bash
$ export CMAKE_PREFIX_PATH=$TIMEMORY_INSTALL:$CMAKE_PREFIX_PATH
$ export PYTHONPATH=$PYTIMEMORYPATH:$PYTHONPATH
$ export MPLCONFIGDIR=$HOME/mplconfigdir
```

### Timemory install example
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

## Setup Instrumentation
If the `USE_TIMEMORY=ON` option is ena, the HiCOPS instrumentation can be configured and updated via the following environment variables:

```bash
HICOPS_MPIP_INSTR           Enable MPI data communication instrumentation. Set to: ON (default), OFF
HICOPS_INST_COMPONENTS      Append to the list of Timemory components (metrics) used for instrumenting the HiCOPS parallel search algorithm. 
                            Set to: HICOPS_INST_COMPONENTS="<c1>,<c2>,.." where each <ci> is a Timemory component.
HICOPS_PAPI_EVENTS          Modify (not append) the vector of PAPI hardware counters used for instrumenting the HiCOPS parallel search algorithm.
                            Set to: HICOPS_PAPI_EVENTS="<hw1>, <hw2>,.." where each <hwi> is a PAPI hardware counter.
TIMEMORY_ENABLED            Enable/disable Timemory instrumentation interface. Set to : ON (default), OFF
```

See more about how to list available timemory components [here](https://timemory.readthedocs.io/en/develop/tools/timemory-avail/README.html?highlight=user_bundle#available-components). 

To see which hardware counters are available on your system and their description, use the `papi_avail` or `timemory-avail` tool. Refer to the PAPI documentation [here](https://icl.utk.edu/papi/) for more information. By default, the following hardware counters are inserted into the `HICOPS_PAPI_EVENTS`.

```bash
HICOPS_PAPI_EVENTS="PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L3_TCM, PAPI_L2_TCA, PAPI_L3_TCA, PAPI_MEM_WCY, PAPI_RES_STL, PAPI_STL_CCY, PAPI_BR_CN, PAPI_BR_PRC, PAPI_FUL_ICY"
```

**NOTE:** If a PAPI counter is not available on the system but is added to the `HICOPS_PAPI_EVENTS` anyway, the profiler will not instrument any of the counters in the list regardless of their availability.
