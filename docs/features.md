---
title: Features
---

# Features

This various features provided by the HiCOPS include:

<!-- TOC -->

- [Parallel Core Library](#parallel-core-library)
- [Performance](#performance)
- [Performance Analysis Support](#performance-analysis-support)
- [Command-line Tools](#command-line-tools)
- [Seamless Integration](#seamless-integration)

<!-- /TOC -->

## Parallel Core Library
The HiCOPS parallel core library is based on the Bulk Synchronous Parallel (BSP) model and divides the database peptide search workflow in ***four*** supersteps. The algorithmic work in each superstep is executed in asynchronous parallel fashion and synchronization is performed between supersteps if required.

 The first superstep partitions the database in a load balanced fashion constructing partial databases. The second superstep preprocesses the experimental MS data in parallel. The third superstep performs the partial database peptide search producing intermediate results which are stored at the shared file system. The final superstep assembles the intermediate results yielding complete results which are written back to the shared file system. Several optimization techniques including buffering, sampling, overlapped computations and communications and task scheduling are implemented to aid in optimal parallel performance. 
 
 Further, the core library packs parallel versions of many core algorithms and the essential boiler-plate code for quicker development of new tools based on the discussed parallel design. Read more about the HiCOPS parallel design, its mathematical analysis and optimizations [here]().

## Performance
The proposed parallel design along with appropriate optimizations allow maximum utilization of underlying resources with minimal overheads enabling ultra-fast peptide search speeds. The highlights of performance evaluation results include:

* **Speedup vs Other tools:** More than $100\times$ against all tools and up to $420\times$ for large workloads compared to MSFragger.      
* **Parallel Efficiency:** 75-85% when the workload size is sufficiently large.     
* **Load Imabalance:** Less than 10%.    
* **Roofline:** Super-linear speedups with increasing nodes when the problem size is large.    
* **Comm/Compute Ratio:** Less than 5%.      

Read more about HiCOPS results, scalability and performance [here]().

## Performance Analysis Support
Performance analysis and instrumentation support has been provided at multiple levels via Timemory toolkit using C++14 `extern templates` allowing custom instantiation of only the required instrumentation components resulting in reduced compile times. Please refer to [Timemory%20Integration](https://github.com/NERSC/timemory#c-template-interface) document on more information on instrumentation integration.

## Command-line Tools
HiCOPS provides several command-line tools for auxiliary purposes including database and experimental data preparation and extraction, post-processing results and runtime support. Please refer to [Tools]({{ site.baseurl }}/tools) for detailed information on tools.

## Seamless Integration
The generic workflow design allows any existing database peptide search algorithms such as data indexing algorithms, search algorithms, scoring algorithms, post-processing, data analytics algorithms to be easily integrated with HiCOPS. 