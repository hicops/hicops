---
title: raw2ms2
---

# raw2ms2

Command-line tool that converts the raw experimental MS/MS data to the `MS2` format. This tool is a Linux wrapper for the ProteoWizard's [msconvert](http://proteowizard.sourceforge.net/tools.shtml) tool and requires `docker` to be installed to run. Please see [this](https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses) for additional information. HiCOPS only supports the `MS2` format as of now thus, all data must be converted to MS2 before using with HiCOPS.

**Note:** If you are working on a Windows machine, please use the `msconvert` tool directly and not use this `docker` based tool for better performance.

## General Syntax

```bash
$ sudo HICOPS_INSTALL/tools/ms2prep/raw2ms2 [fmt] [absolute/path/to/data*.fmt|.FMT]
```
**Paramaters:**
* `fmt`: The file format of the input experimental data. For example: `raw`, `mzML`. Default: `raw`.
* `abs_path`: Absolute path to the folder where the data.fmt or data.FMT are located.       

## Example Usage

```bash
$ sudo HICOPS_INSTALL/tools/ms2prep/raw2ms2 mgf /disk/raptor-2/mhase003/datasets/test/
Using default tag: latest
latest: Pulling from chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
Digest: sha256:f0f70674df9f180f45607ecd8500ca9ca0d091b251fca9ab1cd8710de8fa5132
Status: Image is up to date for chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:latest
format: MS2
outputPath: /disk/raptor-2/mhase003/datasets/test//converted
extension: .ms2
contactFilename:
runIndexSet:

spectrum list filters:

chromatogram list filters:

filenames:
  /disk/raptor-2/mhase003/datasets/test//14Sep18_Olson_S01.mgf

processing file: /disk/raptor-2/mhase003/datasets/test//14Sep18_Olson_S01.mgf
calculating source file checksums
writing output file: /disk/raptor-2/mhase003/datasets/test//converted\14Sep18_Olson_S01.ms2
```