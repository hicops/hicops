#!@BASH_EXECUTABLE@
# 
# Workspace duplicator to generate experiment sets
# Copyrights(C) 2019 PCDS Laboratory
# Muhammad Haseeb, and Fahad Saeed
# School of Computing and Information Sciences
# Florida International University (FIU), Miami, FL
# Email: {mhaseeb, fsaeed}@fiu.edu
#

echo "Creating duplicate workspaces"

# make a copy
cp workspace1 workspace2 -r

# remove garbage data
rm workspace2/output/*
rm workspace2/timemory* -rf

# make more copies
cp workspace2 workspace3 -r
cp workspace2 workspace4 -r

echo "Making path modifications"

# make appropriate edits
for i in $(seq 2 4) 
do
    sed -i "s/\/workspace1/\/workspace$i/g" workspace$i/autogen/hicops

    sed -i "s/\/workspace1/\/workspace$i/g" workspace$i/autogen/uparams.txt
done

echo "Making edits for dM"


# edits for dM
sed -i 's/10.0/100.0/g' workspace2/autogen/uparams.txt

sed -i 's/10.0/200.0/g' workspace3/autogen/uparams.txt

sed -i 's/10.0/500.0/g' workspace4/autogen/uparams.txt

echo "Submitting jobs"

for i in $(seq 1 4)
do
    cd workspace$i
    sbatch autogen/hicops
    cd ..
done

echo "SUCCESS"
