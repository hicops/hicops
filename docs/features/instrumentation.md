---
title: Instrumentation
---
# Instrumentation
The instrumentation and performance analysis is provided via Timemory toolkit. HiCOPS incorporates instrumentation using C++14 `extern templates` which allow custom instantiation of required instrumentation metrics. Please refer to [Timemory](https://timemory.readthedocs.io/en/develop/about.html) documentation for more information.

## Provided Interfaces
HiCOPS provides instrumentation interfaces using Timemory bundles which are C++ tuples containing instrumentation metrics. Time-only measurements are performed via MACROS which direct instrumentation to HiCOPS' native instrumentation or Timemory depending on the configuration.