---
title: About
---
# About
HiCOPS is a software framework designed to accelerate database peptide search workflows on HPC environments for large-scale peptide identification from the mass spectrometry (MS/MS) data. The HiCOPS parallel implements many optimization techniques, algorithms for data indexing, processing, analysis and regression, and novel data structures that can be used to develop new database peptide search algorithms or parallelize existing ones.

HiCOPS has been implemented using C++14, Python 3.7 and Bash.

## Application
Computational Proteomics researchers and developers may integrate their algorithms with HiCOPS by simply implementing their serial indexing/serach/post-processing algorithms in HiCOPS.

## Environments
HiCOPS can seamlessly run in symmetric multicore multinode HPC system (most common HPC systems) environment. Serious speedups compared to all existing (shared and distributed parallel) database peptide search tools. The only constraint is that there should be enough memory in the pool of nodes allocated to HiCOPS to contain the indexed database.

## Credits
Timemory is actively developed by Parallel Computing and Data Science Lab at the Florida International University.

| Name               |                                        Affiliation                                        |                    GitHub                     |
| ------------------ | :---------------------------------------------------------------------------------------: | :-------------------------------------------: |
| Muhammad Haseeb       |       [FIU](https://tinyurl.com/mhaseeb22)       | [mhaseeb123](https://github.com/mhaseeb123) |

## License
GPL v3.0 for all academic users. Commercial users may acquire a license online from the Florida International University Technology Transfer Office.