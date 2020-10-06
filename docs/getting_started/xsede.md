---
title: Run on XSEDE Comet
---

# Run on XSEDE Comet
Follow the below steps to run HiCOPS on XSEDE Comet:

<!-- TOC -->

- [Setup](#setup)
- [Run HiCOPS](#run-hicops)
- [Precautions](#precautions)

<!-- /TOC -->

**Important:** The instructions in this document are ***ONLY*** for the XSEDE Comet system. For any other system, please use the instructions detailed [here]({{ site.baseurl }}/getting_started).

## Setup
Complete the setup steps by following through the instructions in [Getting Started]({{ site.baseurl }}/getting_started#setup) document.

## Run HiCOPS
Assuming that the HiCOPS has been installed at : `$HICOPS_INSTALL`.

**o.** Are you working on XSEDE Comet system? If no, then follow this [link]({{ site.baseurl }}/getting_started) instead.

**i.** Ensure that the hicops-core library path has been added to `LD_LIBRARY_PATH`.      

```bash
# append hicops-core lib path to LD_LIBRARY_PATH.
$ export LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH
```

**ii.** Generate HiCOPS template runtime parameters file using the `hicops_comet` tool located at `$HICOPS_INSTALL/bin/tools`.       

```bash
# run hicops_comet with -g
$ $HICOPS_INSTALL/bin/tools/hicops_comet -g
# generates: ./sampleparams.txt
```

**iii.** Edit the generated sampleparams.txt file and setup HiCOPS' runtime parameters, database and data paths.     

**iv.** Run HiCOPS through `hicops_comet`, however, this time providing the sampleparams.txt as argument to the tool.        

```bash
# run hicops_comet with sampleparams.txt
$ $HICOPS_INSTALL/bin/tools/hicops_comet sampleparams.txt
```

**v.** Run HiCOPS with `uparams.txt` as input argument using SLURM with or without MPI depending on HiCOPS install [options]({{ site.baseurl }}/installation/#cmake-options).        

```bash
# run psm2excel using SLURM
$ srun --partition=compute  --nodes=1 --ntasks-per-node=1 -t 00:15:00 \\
  --export=ALL $HICOPS_INSTALL/tools/psm2excel -i \\
  [/path/to/hicops/workspace/output]
```

**vi.** Repeat Steps **iii** to **v** when you modify parameters in the `sampleparams.txt`.

## Precautions
Please read and understand the precautions mentioned [here]({{ site.baseurl }}/getting_started#precautions) before running HiCOPS.