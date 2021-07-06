#!/bin/bash
#
# SLURM duplicator for MSFragger experiments
# Copyrights(C) 2019 PCDS Laboratory
# Muhammad Haseeb, and Fahad Saeed
# School of Computing and Information Sciences
# Florida International University (FIU), Miami, FL
# Email: {mhaseeb, fsaeed}@fiu.edu
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