---
title: Features
toc: false
---
# Features
This various features packed by the HiCOPS include:
* TOC
{:toc}

## Parallel Core Library
The parallel core workflow based on Bulk Synchronous Parallel supersteps model. Parallel versions of many core algorithms along with necessary boiler-plate code has been implemented and available for faster development of new tools.

## Parallel Performance
Maximum usage of underlying resources with minimal overheads enable serious super-linear speedups in peptide identification speed.

## Seamless Integration
The generic workflow design allows any existing database peptide search algorithms such as data indexing algorithms, search algorithms, scoring algorithms, post-processing, data analytics algorithms to be easily integrated with HiCOPS. 

## Command-line Tools
HiCOPS provides several command-line tools for auxiliary purposes including database and experimental data preparation and extraction, post-processing results and runtime support. Please refer to [Tools]({{ site.baseurl }}/tools) for detailed information on tools. The `hicops_config` tool generates the parameter files used to run HiCOPS.

```bash
$ hicops_comet sampleparamters.txt

*****************************
*  HiCOPS Configuration Gen *
*    PCDS Lab, SCIS, FIU    *
*****************************

Provided Parameters:
Workspace   = /lclhome/mhase003/workspace
Using cores/node  = 22
Processed DB Parts = /lclhome/mhase003/Data/Database/Digested/parts
MS/MS dataset = /lclhome/mhase003/Data/Dataset/PXD009072
Max mods/pep  = 3
Adding mod   = M 15.99 2
Min pep len  = 7
Max pep len  = 40
Min precursor mass = 500
Max precursor mass = 5000
Max frag chg = 3
Required min PSM hits = 4
Base Normalized Intensity = 100000
Intensity Cutoff Ratio = 0.01
Scratch Memory (MB) = 2048
Resolution   = 0.01
Peptide mass tolerance (dM) = 10.0
Fragment mass tolerance (dF) = 0.01
Top matches = 10
Max expect value to report = 20.0

Initializing Workspace at:  /lclhome/mhase003/workspace

SUCCESS

Written: /disk/raptor/lclhome/mhase003/repos/hicops/install-no-inst/bin/uparams.txt

You can now run HiCOPS as:

$ mpirun -np [N] [OPTIONS] /disk/raptor/lclhome/mhase003/repos/hicops/install-no-inst/bin/hicops /disk/raptor/lclhome/mhase003/repos/hicops/install-no-inst/bi
n/uparams.txt

After HiCOPS, run post-processing as:

$ /disk/raptor/lclhome/mhase003/repos/hicops/install-no-inst/bin/tools/psm2excel -i /lclhome/mhase003/workspace/output

#----------------------------------------------------------------------------------------------------#
     Read More: https://github.com/pcdslab/hicops/blob/develop/README.md
#----------------------------------------------------------------------------------------------------#

```

## Performance Analysis Support
Performance analysis and instrumentation support has been provided at multiple levels via Timemory toolkit using C++14 `extern templates` allowing custom instantiation of only the required instrumentation components resulting in reduced compile times. Please refer to [Timemory%20Integration](https://github.com/NERSC/timemory#c-template-interface) document on more information on instrumentation integration.