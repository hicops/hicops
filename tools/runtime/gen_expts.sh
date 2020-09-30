#!@BASH_EXECUTABLE@
# 
# Workspace duplicator to generate experiment sets
# Copyrights(C) 2019 PCDS Laboratory
# Muhammad Haseeb, and Fahad Saeed
# School of Computing and Information Sciences
# Florida International University (FIU), Miami, FL
# Email: {mhaseeb, fsaeed}@fiu.edu
#

# print usage
function usage() {
    echo "USAGE: gen_expts <src> <dst1> <dst2> <dM>"
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

# if everything is fine
echo "Duplicating workspaces..."

# make a temporary copy
echo "Source = workspace${1}"

cp workspace${1} workspacetmp -r
sync

# remove garbage data
rm workspacetmp/output/* -rf
rm workspacetmp/timemory* -rf

# loop in the range of dst workspaces
for i in $(seq ${2} ${3})
do
    cp workspacetmp workspace${i} -r
    sync
    n=$((2**$i))
    m=$((2**${1}))
    echo "Setting nodes = $n"
    sed -i "s/\/workspace${1}/\/workspace${i}/g" workspace${i}/autogen/hicops
    sed -i "s/nodes=${m}/nodes=${n}/g" workspace${i}/autogen/hicops
    sed -i "s/\/workspace${1}/\/workspace${i}/g" workspace${i}/autogen/uparams.txt

    echo "Setting dM = ${4}"
    sed -i "10s/.*/${4}/g" workspace${i}/autogen/uparams.txt
done

# remove temporary data
rm -rf workspacetmp

echo "Submitting jobs"

for i in $(seq ${2} ${3})
do
    cd workspace${i}
    sbatch autogen/hicops
    cd ..
done

echo "SUCCESS"
