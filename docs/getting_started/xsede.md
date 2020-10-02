---
title: Run on XSEDE
---
# Run HiCOPS on XSEDE Comet

**Important:** The instructions contained in this document must be used *ONLY* for the XSEDE Comet system. For any other system, please use the instructions [here]({{ site.baseurl }}/getting_started###Any-Other-System).

Assuming that the HiCOPS has been installed at : `$HICOPS_INSTALL` and the pre-requisite [Setup]({{ site.baseurl }}/getting_started##Setup) steps have been performed.

## Steps
1. Ensure that the hicops-core library path has been added to `LD_LIBRARY_PATH`.      

```bash
# append hicops-core lib path to LD_LIBRARY_PATH.
$ export LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH
```

2. Generate HiCOPS template runtime parameters file using the `hicops_comet` tool located at `$HICOPS_INSTALL/bin/tools`.       

```bash
# run hicops_comet with -g
$ $HICOPS_INSTALL/bin/tools/hicops_comet -g
# generates: ./sampleparams.txt
```

3. Edit the generated sampleparams.txt file and setup HiCOPS' runtime parameters, database and data paths.     

4. Run HiCOPS through `hicops_comet`, however, this time providing the sampleparams.txt as argument to the tool.        

```bash
# run hicops_comet with sampleparams.txt
$ $HICOPS_INSTALL/bin/tools/hicops_comet sampleparams.txt
```

5. Run HiCOPS with `uparams.txt` as input argument using SLURM with or without MPI depending on HiCOPS install [options]({{ site.baseurl }}/installation/##CMake-Options).        

```bash
# run psm2excel using SLURM
$ srun --partition=compute  --nodes=1 --ntasks-per-node=1 -t 00:15:00 \\
  --export=ALL $HICOPS_INSTALL/tools/psm2excel -i \\
  [/path/to/hicops/workspace/output]
```

**Note:** Repeat Steps 3 to 5 if you modify parameters in the `sampleparams.txt`.

## Precautions
Please read the precautions mentioned [here]({{ site.baseurl }}/getting_started##Precautions) to ensure correct and efficient HiCOPS experimental runs.