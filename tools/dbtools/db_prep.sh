#!@BASH_EXECUTABLE@
# 
# Peptide Sequence Database Preprocessor
# 
# Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
# Florida International University, Miami, FL# 
# This program is licensed under the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
# 

# Path to input peptide database
PEPDB=$1
# Min length to filter
MIN=$2
# Max length to filter
MAX=$3

# print usage
function usage() {
    echo "USAGE: db_prep [peptide_seqs.fasta] [min_len] [max_len] [Optional: output_path]"
}

# Output directory
OUTPUT=$4

# if no input database then error
if [ -z "$1" ]; then 
    usage
    exit 0
fi

# if no min length, then error
if [ -z "$2" ]; then 
    usage
    exit 0
fi

# if no max length then error
if [ -z "$3" ]; then 
    usage
    exit 0
fi

# if not then use current directory
if [ -z "$4" ]; then 
    OUTPUT=$PWD
fi

# Temp File
TEMPR=$OUTPUT/PEPDB.tmp

#
# Proceed
#

# Remove the previously existing ./parts if present
rm -rf ${OUTPUT}
wait
sync
mkdir -p ${OUTPUT}
sync

# Remove the FASTA headers and leading/trailing whitespaces, and dump to a file
sed -e '/>.*/d' ${PEPDB} | sed 's/^[ \t]*//;s/[ \t]*$//' > ${TEMPR}

# for loop to dump peptide sequences of increasing lengths
for LOL in $(seq $MIN $MAX); do
    # Print status
    echo "Processing peptides of length: $LOL"

    # separate by length, interchange K and L, lex sort, remove duplicates, interchange K and L again
    sed -nr '/^.{'$LOL','$LOL'}$/p' $TEMPR | sed 's/K/*/g' | sed 's/L/K/g'| sed 's/*/L/g' | sort | uniq | sed 's/L/*/g' | sed 's/K/L/g'| sed 's/*/K/g' > $OUTPUT/$LOL.peps

    # do we want to separate peptides based on K/R termini? No for now
#    if [[ $WRNG == "y" ]]; then
#        grep 'K$' $OUTPUT/$LOL.peps > $OUTPUT/$LOL.K.peps
#        grep 'R$' $OUTPUT/$LOL.peps > $OUTPUT/$LOL.R.peps

#        sed -i '/K$/d' $OUTPUT/$LOL.peps
#        sed -i '/R$/d' $OUTPUT/$LOL.peps

#        cat $OUTPUT/$LOL.R.peps $OUTPUT/$LOL.K.peps > tmp
#        cat $OUTPUT/$LOL.peps  >> tmp
#        mv tmp $OUTPUT/$LOL.peps
#        rm -rf $OUTPUT/$LOL.K.peps $OUTPUT/$LOL.R.peps
#    fi
done

# Remove the temporary files
wait
sync
rm -rf ${TEMPR}
sync
