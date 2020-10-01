---
title: Getting Started
---

# Getting Started
For this document, we will be assuming that the HiCOPS was installed at : `$HICOPS_INSTALL`

## Set environment variables
Append  `HICOPS_INSTALL/lib` to the `LD_LIBRARY_PATH`.

```bash
$ export LD_LIBRARY_PATH=$HICOPS_INSTALL/lib:$LD_LIBRARY_PATH
```

## Setup Database

## Setup MS dataset

## Setup instrumentation
Optionally setup instrumentation by following the instructions in [Instrumentation]({{ site.baseurl }}/instrumentation##Setup-Instrumentation) document.


## Run HiCOPS
Generate HiCOPS sample runtime parameters file using the `hicops_config` located at `$HICOPS_INSTALL/bin`.

```bash
$ $HICOPS_INSTALL/bin/hicops_config -g
Generated: ./sampleparams.txt

SUCCESS
```

Edit the generated sampleparams.txt file and add/modify HiCOPS' runtime parameters.
Generate the HiCOPS runtime parameters file (called `uparams.txt`) as:

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
Generate HiCOPS sample runtime parameters file using the `hicops_comet` tool located at `$HICOPS_INSTALL/bin/tools`.

```bash
$ $HICOPS_INSTALL/bin/tools/hicops_comet -g
Generated: ./sampleparams.txt

SUCCESS
```

Edit the generated sampleparams.txt file and add/modify HiCOPS' runtime parameters.

Run HiCOPS using the same tool i.e. `hicops_comet`, however, this time providing the updated sampleparams.txt as parameter to the tool.

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