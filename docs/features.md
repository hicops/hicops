---
title: Features
---

# Features
<a id="markdown-features" name="features"></a>

This various features packed by the HiCOPS include:

<!-- TOC -->

- [Parallel Core Library](#parallel-core-library)
- [Parallel Performance](#parallel-performance)
- [Seamless Integration](#seamless-integration)
- [Command-line Tools](#command-line-tools)
- [Performance Analysis Support](#performance-analysis-support)

<!-- /TOC -->

## Parallel Core Library
<a id="markdown-parallel-core-library" name="parallel-core-library"></a>
The parallel core workflow based on Bulk Synchronous Parallel supersteps model. Parallel versions of many core algorithms along with necessary boiler-plate code has been implemented and available for faster development of new tools.

## Parallel Performance
<a id="markdown-parallel-performance" name="parallel-performance"></a>
Maximum usage of underlying resources with minimal overheads enable serious super-linear speedups in peptide identification speed.

## Seamless Integration
<a id="markdown-seamless-integration" name="seamless-integration"></a>
The generic workflow design allows any existing database peptide search algorithms such as data indexing algorithms, search algorithms, scoring algorithms, post-processing, data analytics algorithms to be easily integrated with HiCOPS. 

## Command-line Tools
<a id="markdown-command-line-tools" name="command-line-tools"></a>
HiCOPS provides several command-line tools for auxiliary purposes including database and experimental data preparation and extraction, post-processing results and runtime support. Please refer to [Tools]({{ site.baseurl }}/tools) for detailed information on tools.

## Performance Analysis Support
<a id="markdown-performance-analysis-support" name="performance-analysis-support"></a>
Performance analysis and instrumentation support has been provided at multiple levels via Timemory toolkit using C++14 `extern templates` allowing custom instantiation of only the required instrumentation components resulting in reduced compile times. Please refer to [Timemory%20Integration](https://github.com/NERSC/timemory#c-template-interface) document on more information on instrumentation integration.