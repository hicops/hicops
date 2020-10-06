---
title: extractms2
---

# extractms2

Command-line tool that extracts the specified (by index numbers beginning with 1) MS/MS spectra from a `MS2` file and writes them to a new `MS2` file. This tool may be used for sampling dataset, debugging, testing and so on.

## General Syntax

```bash
$ $HICOPS_INSTALL/tools/ms2prep/extractms2 -h
USAGE: extractms2 [-h] -i IFILE -d N [N ...]

Extract spectra from MS2 files

optional arguments:
  -h, --help            show this help message and exit
  -i IFILE, --infile IFILE
                        The MS2 file to extract spectra from
  -d N [N ...], --id N [N ...]
                        Spectra numbers in the file
```


## Example Usage

```bash
$ $HICOPS_INSTALL/tools/ms2prep/extractms2 -i /disk/raptor-2/mhase003/datasets/test/converted/14Sep18_Olson_S01.ms2 -d $(seq 1 10)

Writing extracted data at: /disk/raptor-2/mhase003/datasets/test/converted/14Sep18_Olson_S01_extracted.ms2

$ cat /disk/raptor-2/mhase003/datasets/test/converted/14Sep18_Olson_S01_extracted.ms2 | head -n 20
H       CreationDate Mon Oct  5 21:20:06 2020
H       Extractor       ProteoWizard
H       Extractor version       pwiz_3.0.20239
H       Source file     14Sep18_Olson_S01.mgf
S       0       0       0
I       NativeID        index=0
I       RTime   0.03216476
I       BPI     125204.6
I       BPM     357.0676
I       TIC     360397.3
Z       4       1425.175
120.3386 271.9986
125.4088 291.5909
145.5829 296.1437
149.231 322.5773
149.2636 418.1497
149.2701 539.7873
149.2766 525.7592
149.2828 993.4144
149.2897 1670.596
```