---
title: db_prep
---

# db_prep

Command line tool that separates coarse-grained peptide sequence clusters of the specified range of peptide sequence lengths from a peptide sequence database (FASTA format). The output `.pep` files are directly used by the HiCOPS for database constructions and processing.

## General Syntax

```bash
$ $HICOPS_INSTALL/tools/dbtools/db_prep [peptide_seqs.fasta] [min_len] [max_len] [Optional: output_path]
```

**Paramaters:**
* `peptide_seqs`: path to peptide sequence database.        
* `min_len`: minimum length peptide sequences to extract.      
* `max_len`: maximum length peptide sequences to extract.      
* `output_path`: path where the .pep files will be written.       

## Example Usage

```bash
$ $HICOPS_INSTALL/tools/dbtools/db_prep ./human_peptides.fasta 6 20 ./output
Processing peptides of length: 6
Processing peptides of length: 7
Processing peptides of length: 8
Processing peptides of length: 9
Processing peptides of length: 10
Processing peptides of length: 11
Processing peptides of length: 12
Processing peptides of length: 13
Processing peptides of length: 14
Processing peptides of length: 15
Processing peptides of length: 16
Processing peptides of length: 17
Processing peptides of length: 18
Processing peptides of length: 19
Processing peptides of length: 20

$ ls output
10.peps  11.peps  12.peps  13.peps  14.peps  15.peps  16.peps  17.peps  18.peps  19.peps  20.peps  6.peps  7.peps  8.peps  9.peps
```