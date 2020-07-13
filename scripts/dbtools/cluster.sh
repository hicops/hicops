# 
# Peptide Sequence Clusterer
# Copyrights(C) 2019 PCDS Laboratory
# Muhammad Haseeb, and Fahad Saeed
# School of Computing and Information Sciences
# Florida International University (FIU), Miami, FL
# Email: {mhaseeb, fsaeed}@fiu.edu
#

#!/bin/bash

# Path to input peptide database
PEPDB=$1
# Min length to filter
MIN=$2
# Max length to filter
MAX=$3

# Perform C-term K/R data wrangling?
WRNG=$4

# Output Directory
OUTPUT=$5/parts

# Temp File
TEMPR=$5/PEPDB.tmp

# Remove the previously existing ./parts if present
rm -rf ${OUTPUT}
wait
sync
mkdir -p ${OUTPUT}
sync

# Remove the FASTA header files and create a temp file
sed -e '/>.*/d' ${PEPDB} > ${TEMPR}

# for loop to dump peptide sequences of increasing lengths
for LOL in $(seq $MIN $MAX); do

    sed -nr '/^.{'$LOL','$LOL'}$/p' $TEMPR | sed 's/K/*/g' | sed 's/L/K/g'| sed 's/*/L/g' | sort | uniq | sed 's/L/*/g' | sed 's/K/L/g'| sed 's/*/K/g' > $OUTPUT/$LOL.peps

    if [[ $WRNG == "y" ]]; then
        grep 'K$' $OUTPUT/$LOL.peps > $OUTPUT/$LOL.K.peps
        grep 'R$' $OUTPUT/$LOL.peps > $OUTPUT/$LOL.R.peps

        sed -i '/K$/d' $OUTPUT/$LOL.peps
        sed -i '/R$/d' $OUTPUT/$LOL.peps

        cat $OUTPUT/$LOL.R.peps $OUTPUT/$LOL.K.peps > tmp
        cat $OUTPUT/$LOL.peps  >> tmp
        mv tmp $OUTPUT/$LOL.peps
        rm -rf $OUTPUT/$LOL.K.peps $OUTPUT/$LOL.R.peps
    fi
done

# Remove the temporary files
wait
sync
rm -rf ${TEMPR}
sync
