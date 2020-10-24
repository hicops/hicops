![build](https://github.com/hicops/hicops/workflows/build/badge.svg) [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/hicops/hicops/blob/develop/README.md#contributing) [![pages yes](https://img.shields.io/badge/pages-yes-blue.svg)](https://hicops.github.io) [![GitHub forks](https://img.shields.io/github/forks/hicops/hicops.svg?style=social&label=Fork&maxAge=2592000)](https://GitHub.com/hicops/hicops/network/) [![GitHub stars](https://img.shields.io/github/stars/hicops/hicops.svg?style=social&label=Star&maxAge=2592000)](https://GitHub.com/hicops/hicops/stargazers/) [![GitHub contributors](https://img.shields.io/github/contributors/hicops/hicops.svg)](https://GitHub.com/hicops/hicops/graphs/contributors/) [![GitHub issues](https://img.shields.io/github/issues/hicops/hicops.svg)](https://GitHub.com/hicops/hicops/issues/) [![Github all releases](https://img.shields.io/github/downloads/hicops/hicops/total.svg)](https://GitHub.com/hicops/hicops/releases/)

# HiCOPS
*HiCOPS*: A computational framework for accelerated peptide identification from LC-MS/MS data on HPC systems.

For full documentation [here](https://hicops.github.io).

## Install
Please follow the documentation [here](https://hicops.github.io/installation).

## Getting Started with HiCOPS
Please follow the documentation [here](https://hicops.github.io/getting_started).

## Precautions
Please read and follow the following precautions to avoid any errors.

* Always use a unique workspace directory for each experiment, specially for the simultaneously running HiCOPS instances to avoid overwriting intermediate results and other errors.
* Always convert the TSV results into Excel file using `psm2excel` (last step) before using the same workspace folder for another experiment. It is always better to use a unique workspace folder for each experimental run.
* Do not run too many simultaneous HiCOPS instances with large number of nodes allocated to each instance to avoid I/O bandwidth contention and thus, performance degradation.
* Do not modify the generated files such as uparams.txt manually and instead regenerate using the relevant tool.
* Avoid using relative paths in the `sampleparams.txt` file to avoid any errors.

## Integrating Code with HiCOPS
Please refer to the documentation [here](https://hicops.github.io/getting_started/integrate).

# Published Research
*Placeholder for the published research paper*

# Credits
1. [Muhammad Haseeb](https://sites.google.com/fiu.edu/muhammadhaseeb)
2. [Fahad Saeed](http://www.saeedfahad.com)

# Issues and Feature Requests
Open an issue [here](https://github.com/hicops/hicops/issues). Please include any logs, screenshots and/or helpful observations. Also, do not forget to describe the dataset(s), database, steps performed etc. so that the issue can be reproduced.

# Contributing
All contributions are welcome including new features, documentation and bug fixes through standard GitHub pull request method. Generic guidelines:

1. Fork this repository to your local GitHub, checkout a new branch from the `develop` branch.
2. Make your changes/updates.
3. Make sure that you pull the latest changes from `hicops:develop` into your branch and merge before committing your changes.
4. Commit your changes and push your branch to `origin`. i.e. `your_hicops_fork`.
5. Open a pull request (PR) from `your_hicops_fork:new_branch` to `hicops:develop`.
