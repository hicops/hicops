---
title: hicops_config
---

# hicops_config

Command line tool that takes in the user paramters and generates a file, called `uparams.txt`, which is used as arguments to the `hicops` binary. This tool also generates a template user paramter file which can be modified as needed and then used used with `hicops_config` to generate the `uparams.txt` file. 

**Note:** Regenerate a `uparams.txt` everytime you modify the user parameters.

**Note:** If you are working on XSEDE Comet, you can follow an easier workflow by following instructions [here]({{ site.baseurl }}/getting_started/xsede) and using the [hicops_comet]({{ site.baseurl }}/tools/runtime/hicops_comet) tool.

## General Syntax

### Generate a template user parameters file

```bash
$ $HICOPS_INSTALL/bin/hicops_config -g
Generated: ./sampleparams.txt

SUCCESS
```

**Parameters:**
* `-g`: flag that specifies to generate a template user paramters file.     

### Generate `uparams.txt` from user params

```bash
$ $HICOPS_INSTALL/bin/hicops_config <params.txt>
```

**Parameters:**
* `params`: user parameters file, template generated via the `-g` flag and then edited for experiment.     

## Example Usage

```bash
$ $HICOPS_INSTALL/bin/hicops_config userparams.txt

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

$ mpirun -np [N] [OPTIONS] /disk/raptor/lclhome/mhase003/repos/hicops/install-no-inst/bin/hicops /disk/raptor/lclhome/mhase003/repos/hicops/install-no-inst/bin/uparams.txt

After HiCOPS, run post-processing as:

$ /disk/raptor/lclhome/mhase003/repos/hicops/install-no-inst/bin/tools/psm2excel -i /lclhome/mhase003/workspace/output

#----------------------------------------------------------------------------------------------------#
     Read More: https://github.com/pcdslab/hicops/blob/develop/README.md
#----------------------------------------------------------------------------------------------------#
```