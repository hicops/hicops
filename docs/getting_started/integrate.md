---
title: Integration
---

# Integration
HiCOPS provides a function interface for integrating database peptide search algorithms its parallel core library. Currently, all HiCOPS data structures, boiler-plate code, parallel algorithms are exposed publicly to allow full control. 

## Methods
The shared memory versions of algorithms are implemented inside relevant functions using OpenMP or `C++::threads` depending on the function's specification. The following sub-sections will provide an overview of integrating new or existing shared memory database peptide search algorithms with HiCOPS parallel library.

### Database Indexing

### MS/MS Data Preprocessing

### Search Scoring

### Statistical Scoring

### False Discovery Rate


## Active Development (Future)
We are actively developing a new C++ template meta-programming based interface for HiCOPS allowing simple modular application following the build your own workflow philosophy. The new interface will be similar to LLVM and Timemory allowing users to easily build and customize their workflows without having to parallelize (even with OpenMP or threads) or directly interact with the core library.

