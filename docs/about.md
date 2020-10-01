# About
HiCOPS is a software framework designed to accelerate database peptide search workflows on HPC environments for large-scale peptide identification from the mass spectrometry (MS/MS) data. The HiCOPS parallel implements many optimization techniques, algorithms for data indexing, processing, analysis and regression, and novel data structures that can be used to develop new database peptide search algorithms or parallelize existing ones.

# Credits
Timemory is actively developed by Parallel Computing and Data Science Lab at the Florida International University.

# Applications
Computational Proteomics researchers and developers may integrate their algorithms with HiCOPS by simply implementing their serial indexing/serach/post-processing algorithms in HiCOPS.

# Features
Maximum usage of underlying resources with minimal overheads enable serious speedups in peptide identification speed. The generic workflow design allows any existing database peptide search algorithms such as data indexing algorithms, search algorithms, scoring algorithms, post-processing, data analytics algorithms to be easily integrated with HiCOPS. Parallel versions of many core algorithms along with necessary boiler-plate code has been implemented and available for faster development of new tools. Performance instrumentation and analysis (both macro and micro levels) has never been easier using simple template interfaces several metrics provided via [Timemory](https://github.com/NERSC/timemory.git).

# Usage
HiCOPS can seamlessly run in symmetric multicore multinode HPC system (most common HPC systems) environment. Serious speedups compared to all existing (shared and distributed parallel) database peptide search tools. The only constraint is that there should be enough memory in the pool of nodes allocated to HiCOPS to contain the indexed database.

# Limitations (for now)
Scale seamlessly in asymmetric multinode HPC systems. Offload algorithmic workload to GPUs, accelerators etc. Multi-language support

# Contributions
HiCOPS encourages contributions via GitHub pull requests. Please follow the following method for contributing new algorithms, bug fixes and so on:

1. Fork this repository to your local GitHub account
2. Create a new branch from the `develop` branch.
3. Make your changes/updates.
4. Make sure that you have the latest changes from `hicops:develop` into your branch
5. Commit and push your changes to your fork.
6. Open a pull request (PR) from `your_hicops_fork:new_branch` to `hicops:develop`.
7. We will review the changes and merge them into `hicops:develop`.

# Issue Reporting
Open an [issue](https://github.com/pcdslab/hicops/issues). Please include any logs, screenshots and/or helpful observations. Also, do not forget to describe the dataset(s), database, steps performed etc. so that the issue can be reproduced.

# License
GPL v3.0 for all academic users. Commercial users may acquire a license online from the Florida International University Technology Transfer Office.