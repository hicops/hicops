---
title: hicops_comet
---

# hicops_comet

Command-line tool that parses user parameters, and creates and runs the HiCOPS parallel tasks on the XSEDE Comet cluster. This tool also generates a template user paramter file which can be modified as needed and then used used with `hicops_comet` to run HiCOPS.

**Note:** This tool is specialized only when running HiCOPS on XSEDE Comet. For any other system, please follow the instructions [here]({{ site.baseurl }}/getting_started#any-other-system) and use the [hicops_config]({{ site.baseurl }}/tools/runtime/hicops_config) tool instead.

## General Syntax

### Generate a template user parameters file

```bash
$ $HICOPS_INSTALL/bin/tools/hicops_comet -g
Generated: ./sampleparams.txt

SUCCESS
```

**Parameters:**
* `-g`: flag that specifies to generate a template user paramters file.     

### Run HiCOPS on XSEDE Comet

```
$ $HICOPS_INSTALL/bin/tools/hicops_comet <params.txt>
```

**Parameters:**
* `params`: user parameters file, template generated via the `-g` flag and then edited for experiment.     

## Example Usage

```bash
$ $HICOPS_INSTALL/bin/tools/hicops_comet sampleparams.txt

***************************
*  HiCOPS for XSEDE Comet *
*   PCDS Lab, SCIS, FIU   *
***************************

Provided Parameters:
Workspace   = /oasis/scratch/comet/mhaseeb/temp_project/workspaces/settest/workspace1
Job time = 01:00:00
Using nodes = 16
Using cores/node  = 24
MPI binding policy = scatter
MPI binding level = socket
Optimizations = 1
Email events = FAIL
Processed DB Parts = /home/mhaseeb/database/swiss/parts
MS/MS dataset = /oasis/scratch/comet/mhaseeb/temp_project/PXD015890
Max mods/pep  = 0
Adding mod   = M 15.99 3
Adding mod   = Y 361.146 3
Min pep len  = 7
Max pep len  = 46
Min precursor mass = 500
Max precursor mass = 5000
Max frag chg = 3
Required min PSM hits = 4
Base Normalized Intensity = 100000
Intensity Cutoff Ratio = 0.01
Scratch Memory (MB) = 2048
Resolution   = 0.01
Peptide mass tolerance (dM) = 500.0
Fragment mass tolerance (dF) = 0.01
Top matches = 10
Max expect value to report = 20.0

Initializing Workspace at:  /oasis/scratch/comet/mhaseeb/temp_project/workspaces/settest/workspace1


****** Autotuning parameters *******

Submitted batch job 36204166

Waiting for job scheduler

****** Autotuning parameters *******


Extracted System Settings

Available cores/machine  = 24
Available cores/socket  = 12
Available sockets/machine  = 2
Available NUMA nodes/machine = 2
Available max NUMA memory (- 1GB ) = 58984


Submitted batch job 36204175


Estimating Index Size

Estimated Index Size (x 1E6 Spectra) = 28


Optimized HiCOPS settings...

Setting Threads/MPI = 12
Setting Max Prep threads/MPI = 4
Setting MPI/machine = 2
Setting MPI Policy  = scatter
Setting MPI Binding = socket
Setting Index / MPI = 894414

SUCCESS

Submitted batch job 36204191

HiCOPS is running now

You can check the job progress by:

$ squeue -u $USER

The output will be written at: /oasis/scratch/comet/mhaseeb/temp_project/workspaces/settest/workspace1/output

SUCCESS

After job completion, run:

$ srun --partition=compute --nodes=1 --ntasks-per-node=1 -t 00:25:00 --export=ALL /home/mhaseeb/repos/hicops/install/bin/tools/../tools/psm2excel -i /oasis/scratch/comet/mhaseeb/temp_project/workspaces/settest/workspace1/output


#----------------------------------------------------------------------------------------------------#
     Read More: https://github.com/pcdslab/hicops/blob/develop/README.md
#----------------------------------------------------------------------------------------------------#

```