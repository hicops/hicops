# Getting Started

For this document, we will be assuming that the HiCOPS was installed at : `$HICOPS_INSTALL`

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
TIMEMORY_ENABLED            Enable/disable Timemory instrumentation interface. Set to : ON (default), OFF
```

See more about how to list available timemory components [here](https://timemory.readthedocs.io/en/develop/tools/timemory-avail/README.html?highlight=user_bundle#available-components). 

To see which hardware counters are available on your system and their description, use the `papi_avail` or `timemory-avail` tool. Refer to the PAPI documentation [here](https://icl.utk.edu/papi/) for more information. By default, the following hardware counters are inserted into the `HICOPS_PAPI_EVENTS`.

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

**NOTE:** Repeat Steps # 2 and 3 if you modify parameters in the `sampleparams.txt`.

# Post-processing HiCOPS output
HiCOPS generates PSM data in partial TSV files that can be merged using the `psm2excel` tool located at: `$HICOPS_INSTALL/tool`. The tool generates a combined Excel file called `Concat.xlsx` containing the final PSM data (no-FDR).

## On a regular computer (skip if using XSEDE Comet)
Run the `psm2excel` tool and pass the HiCOPS workspace output directory (that was set in the sampleparams.txt file) as parameter.

```bash
$ $HICOPS_INSTALL/tools/psm2excel [/path/to/hicops/workspace/output]
```

## On XSEDE Comet
Run the `psm2excel` tool using SLURM and pass the HiCOPS workspace output directory (that was set in the sampleparams.txt file) as parameters.

```bash
$ srun --partition=compute  --nodes=1 --ntasks-per-node=1 -t 00:15:00 --export=ALL $HICOPS_INSTALL/tools/psm2excel -i [/path/to/hicops/workspace/output]
```