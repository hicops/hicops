#!@BASH_EXECUTABLE@

# Convert RAW files to MS2 for HiCOPS on Linux

# Copyrights(C) 2019 PCDS Laboratory
# Muhammad Haseeb and Fahad Saeed
# School of Computing and Information Sciences (SCIS)
# Florida International University (FIU), Miami, FL
# Email: {mhaseeb, fsaeed}@fiu.edu

# ------------------------------------------------------------- #
# NOTE: Always run this script with sudo
# USAGE: sudo msconvert absolute/path/to/data*[.raw|.RAW]
# DEPENDENCY: docker installed on the system
# ------------------------------------------------------------- #

# Set working directory from command line
WDIR=$1

# if not then use current directory
if [ -z "$1" ]; then 
    WDIR=$PWD
fi

# make a list of all the .raw and .RAW files
ls -rt -d -1 "$PWD"/{*,.*} | grep -E "\.raw|\.RAW" > $WDIR/list

# pull msconvert docker image
docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses

# run msconvert
docker run -it --rm -v $WDIR:$WDIR chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert -f $WDIR/list --ms2 -o $WDIR/converted

# make the output directory user accessible
chown $USER:user $WDIR/converted -R

# wait for everything to sync up
wait
sync

# remove the list
rm -rf $WDIR/list