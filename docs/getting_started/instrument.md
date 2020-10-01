---
title: Setup Instrumentation
---

## Setup Instrumentation
If HiCOPS instrumentation was enabled during build, it can be configured and updated by setting the following environment variables. See how to enable HiCOPS instrumentation in the [Installation]({{ site.baseurl }}/installations) document:

| Variable               | Description                                                                                                          |
|------------------------|----------------------------------------------------------------------------------------------------------------------|
| TIMEMORY_ENABLED       | Enable/disable Timemory instrumentation interface. Set to : ON (default), OFF                                        |
| HICOPS_MPIP_INSTR      | Enable MPI data communication instrumentation. Set to: ON (default), OFF                                             |
| HICOPS_INST_COMPONENTS | Append instrumentation components. Set to: HICOPS_INST_COMPONENTS="c1,c2,.." where ci is a Timemory component  |
| HICOPS_PAPI_EVENTS     | Modify the hardware counters. Set to: HICOPS_PAPI_EVENTS="h1,h2,.." where hi is a PAPI counter                 |

To list all available timemory components [here](https://timemory.readthedocs.io/en/develop/tools/timemory-avail/README.html?highlight=user_bundle#available-components). By default, the following hardware counters are inserted into the `HICOPS_PAPI_EVENTS`.

```bash
HICOPS_PAPI_EVENTS="PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L3_TCM, PAPI_L2_TCA, PAPI_L3_TCA, PAPI_MEM_WCY, PAPI_RES_STL, PAPI_STL_CCY, PAPI_BR_CN, PAPI_BR_PRC, PAPI_FUL_ICY"
``` 

To see which hardware counters are available on your system and their description, use the `papi_avail` or `timemory-avail` tool. Refer to the PAPI documentation [here](https://icl.utk.edu/papi/) for more information. 

**NOTE:** If a PAPI counter is not available on the system but is added to the `HICOPS_PAPI_EVENTS` anyway, the profiler will not instrument any of the counters in the list regardless of their availability.
