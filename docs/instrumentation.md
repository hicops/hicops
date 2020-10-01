# Instrumentation
The instrumentation and performance analysis is provided via Timemory toolkit. HiCOPS incorporates instrumentation using C++14 `extern templates` which allow custom instantiation of required instrumentation metrics. Please refer to [Timemory](https://timemory.readthedocs.io/en/develop/about.html) documentation for more information.

## Provided Interfaces
HiCOPS provides instrumentation interfaces using Timemory bundles which are C++ tuples containing instrumentation metrics. Time-only measurements are performed via `MACROS` which direct instrumentation to HiCOPS' native instrumentation or Timemory depending on the configuration.

## Enable Instrumentation

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
