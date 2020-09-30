#!@BASH_EXECUTABLE@

# Convert RAW files to MS2 for HiCOPS on Linux

# Copyrights(C) 2019 PCDS Laboratory
# Muhammad Haseeb and Fahad Saeed
# School of Computing and Information Sciences (SCIS)
# Florida International University (FIU), Miami, FL
# Email: {mhaseeb, fsaeed}@fiu.edu

# --------------------------------------------------------------------------------- #
# NOTE: Always run this script with sudo
# USAGE: sudo raw2ms2 absolute/path/to/data*[.raw|.RAW]
# DEPENDENCY: docker installed on the system
# HELP: https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
# --------------------------------------------------------------------------------- #

# Set working directory from command line
WDIR=$1

# if not then use current directory
if [ -z "$1" ]; then 
    WDIR=$PWD
fi

# make a list of all the .raw and .RAW files
ls -rt -d -1 "$WDIR"/{*,.*} | grep -E "\.raw|\.RAW" > $WDIR/list

# pull msconvert docker image
docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses

# run msconvert
docker run -it --rm -v $WDIR:$WDIR chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert -f $WDIR/list --ms2 -o $WDIR/converted

# make the output directory user accessible
chown ${SUDO_USER:-${USER}}:user $WDIR/converted -R

# wait for everything to sync up
wait
sync

# remove the list
rm -rf $WDIR/list
