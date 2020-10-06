---
title: Welcome
toc: true
---

![build](https://github.com/mhaseeb123/hicops/workflows/build/badge.svg) [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)]({{ site.baseurl }}/contributing#Contribute) [![pages yes](https://img.shields.io/badge/pages-yes-blue.svg)](https://mhaseeb123.github.io/hicops) [![GitHub forks](https://img.shields.io/github/forks/mhaseeb123/hicops.svg?style=social&label=Fork&maxAge=2592000)](https://GitHub.com/mhaseeb123/hicops/network/) [![GitHub stars](https://img.shields.io/github/stars/mhaseeb123/hicops.svg?style=social&label=Star&maxAge=2592000)](https://GitHub.com/mhaseeb123/hicops/stargazers/) [![GitHub contributors](https://img.shields.io/github/contributors/mhaseeb123/hicops.svg)](https://GitHub.com/mhaseeb123/hicops/graphs/contributors/) [![GitHub issues](https://img.shields.io/github/issues/mhaseeb123/hicops.svg)](https://GitHub.com/mhaseeb123/hicops/issues/) [![Github all releases](https://img.shields.io/github/downloads/mhaseeb123/hicops/total.svg)](https://GitHub.com/mhaseeb123/hicops/releases/)

# About

HiCOPS is a software framework designed to accelerate database peptide search workflows on HPC environments for large-scale peptide identification from the mass spectrometry (MS/MS) data. The HiCOPS parallel implements many optimization techniques, algorithms for data indexing, processing, analysis and regression, and novel data structures that can be used to develop new database peptide search algorithms or parallelize existing ones.

HiCOPS has been implemented using C++14, Python 3.7 and Bash. Detailed information about HiCOPS development, upcoming changes can be found in [README](https://github.com/mhaseeb123/hicops) on GitHub.

## Research
Please read the following research paper for extensive introduction to the research work, parallel design and technical concepts behind HiCOPS. A preprint of the paper is available [here]().

*Placeholder for the paper*

Please cite us if you use our work. Thank you.

## Application
Computational Proteomics researchers and developers can utilize HiCOPS core library to accelerate their worklflows. Integration is as simple as implementing shared-memory versions of database indexing, filtering, peptide-to-spectrum scoring, post-processing etc. algorithms within HiCOPS.

## Supported Environments
HiCOPS can seamlessly run on any symmetric multicore multinode (the most common) HPC environment. Sufficient amount of memory resources must be allocated to HiCOPS for database indexing.

## Credits
HiCOPS is being developed by the Parallel Computing and Data Science Laboratory at Florida International University.

| Name               |                                        Affiliation                                        |                    GitHub                     |
| ------------------ | :---------------------------------------------------------------------------------------: | :-------------------------------------------: |
| Muhammad Haseeb       |       [FIU](https://tinyurl.com/mhaseeb22)       | [mhaseeb123](https://github.com/mhaseeb123) |

## Contributions
We welcome contributions via GitHub pull requests. For more information, please refer to [Contributing]({{ site.baseurl }}/contributing) document.

## License
GPL v3.0 for all academic users. Commercial users may acquire a license from the [FIU Technology Transfer Office](http://research.fiu.edu/ored/)