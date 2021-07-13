#!/bin/bash
#
# SLURM duplicator for MSFragger experiments
# 
# Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
# Florida International University, Miami, FL# 
# This program is licensed under the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
# 

# print usage
function usage() {
    echo "USAGE: gen_slurms <slurm_file> <nodes1> <nodes2>"
}

# if no input params
if [ -z "$1" ]; then
    usage
    exit 0
fi

# if no input params
if [ "$1" = -h ]; then
    usage
    exit 0
fi

# if everything is alright then
echo "Duplicating slurm files"

# make a temporary copy
echo "Source = ${1}"

# loop in the range of dst workspaces
for i in $(seq ${2} ${3})
do
    cp ${1} ${i}_${1} -r
    sync
    echo "Setting data partition = $i"
    sed -i "s/mzML4Fragger/mzML4Fragger\/part_${i}/g" ${i}_${1}
done

# echo "Submitting jobs"

#for i in $(seq ${2} ${3})
#do
#    cd workspace${i}
#    sbatch autogen/hicops
#    cd ..
#done

echo "SUCCESS"