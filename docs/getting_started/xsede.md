---
title: Run on XSEDE
---
# Run HiCOPS on XSEDE Comet
For this document, we will be assuming that the HiCOPS has been installed at : `$HICOPS_INSTALL`

We will use the `hicops_comet` runtime tool to generate and parse user parameters, generate `uparams.txt`, check enough memory in resource pool, and create and submit HiCOPS SLURM jobs. For more information about the `hicops_comet` tool, refer to [tools]({{ site.baseurl }}/tools/runtime/hicops_comet).

## Steps
1. Ensure that the hicops-core library path has been added to `LD_LIBRARY_PATH`.      

```bash
# append hicops-core lib path to LD_LIBRARY_PATH.
$ export LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH
```

2. Generate HiCOPS template runtime parameters file using the `hicops_comet` tool located at `$HICOPS_INSTALL/bin/tools`.       

```bash
$ $HICOPS_INSTALL/bin/tools/hicops_comet -g
Generated: ./sampleparams.txt
```

3. Edit the generated sampleparams.txt file and setup HiCOPS' runtime parameters, database and data paths.     

4. Run HiCOPS through `hicops_comet`, however, this time providing the sampleparams.txt as argument to the tool.        

```bash
$ $HICOPS_INSTALL/bin/tools/hicops_comet sampleparams.txt
```

5. Run HiCOPS with `uparams.txt` as input argument using SLURM with or without MPI depending on HiCOPS install [options]({{ site.baseurl }}/installation/##CMake-Options).        

```bash
$ srun --partition=compute  --nodes=1 --ntasks-per-node=1 -t 00:15:00 \\
  --export=ALL $HICOPS_INSTALL/tools/psm2excel -i \\
  [/path/to/hicops/workspace/output]
```

**NOTE:** Repeat Steps 3-5 if you modify parameters in the `sampleparams.txt`.