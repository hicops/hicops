---
title: Tools
---

# Tools
<a id="markdown-tools" name="tools"></a>
This section covers the executable tools that are distributed as part of HiCOPS. The provided tools aid in peptide database processing, MS/MS data conversions, and runtime support. The subsequent sections provide an overview of the listed tools.

<!-- TOC -->

- [Runtime](#runtime)
  - [hicops_config](#hicops_config-sitebaseurl-toolsruntimehicops_config)
  - [hicops_comet](#hicops_comet-sitebaseurl-toolsruntimehicops_comet)
- [Database Preprocessing](#database-preprocessing)
  - [db_prep](#db_prep-sitebaseurl-toolsdbtoolsdbprep)
- [MS Data Processing](#ms-data-processing)
  - [raw2ms2](#raw2ms2-sitebaseurl-toolsms2prepraw2ms2)
  - [extractms2](#extractms2-sitebaseurl-toolsms2prepextractms2)

<!-- /TOC -->

## Runtime
<a id="markdown-runtime" name="runtime"></a>
The following tools help generate input parameters and run HiCOPS.

### [hicops_config]({{ site.baseurl }}/tools/runtime/hicops_config)
<a id="markdown-hicops_config" name="hicops_config"></a>
This tool generates a compact params file, called `uparams.txt` from a user-friendly parameters/setting file. The generated `uparams.txt` is used with HiCOPS as arguments. A template of user-friendly parameters file can also be generated. More details [here]({{ site.baseurl }}/tools/runtime/hicops_config).

### [hicops_comet]({{ site.baseurl }}/tools/runtime/hicops_comet)
<a id="markdown-hicops_comet" name="hicops_comet"></a>
This tool serves as the main utility for running HiCOPS on XSEDE Comet using a user-friendly parameters file. A template parameters file can also be generated. More details [here]({{ site.baseurl }}/tools/runtime/hicops_comet).

## Database Preprocessing
<a id="markdown-database-preprocessing" name="database-preprocessing"></a>
The following tools help prepare the peptide sequence database for HiCOPS.

### [db_prep]({{ site.baseurl }}/tools/dbtools/dbprep)
<a id="markdown-db_prep" name="db_prep"></a>
This tool prepares coarse-grained peptide sequence clusters from a peptide sequence database (in FASTA format). The generated files are used as input database to HiCOPS. More details [here]({{ site.baseurl }}/tools/database/db_prep).

## MS Data Processing
<a id="markdown-ms-data-processing" name="ms-data-processing"></a>
The following tools help prepare the experimental MS/MS data for HiCOPS.

### [raw2ms2]({{ site.baseurl }}/tools/ms2prep/raw2ms2)
<a id="markdown-raw2ms2" name="raw2ms2"></a>
This tool is a wrapper for the ProteoWizard's msconvert tool on Linux systems. This tools requires `docker` to be installed on the system. Alternatively, the users may directly use msconvert tool for MS/MS data conversion to MS2 format. More details [here]({{ site.baseurl }}/tools/ms2prep/raw2ms2).

### [extractms2]({{ site.baseurl }}/tools/ms2prep/extractms2)
<a id="markdown-extractms2" name="extractms2"></a>
This tool extracts a range of MS/MS spectra from an `MS2` file and writes them to a new file. This tools may be used in debugging purposes or data sampling. More details [here]({{ site.baseurl }}/tools/ms2prep/extractms2).