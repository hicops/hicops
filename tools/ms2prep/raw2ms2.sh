#!@BASH_EXECUTABLE@

# Convert RAW files to MS2 for HiCOPS on Linux

# Copyrights(C) 2019 PCDS Laboratory
# Muhammad Haseeb and Fahad Saeed
# School of Computing and Information Sciences (SCIS)
# Florida International University (FIU), Miami, FL
# Email: {mhaseeb, fsaeed}@fiu.edu

# --------------------------------------------------------------------------------- #
# NOTE: Always run this script with sudo
# USAGE: sudo raw2ms2 [fmt] absolute/path/to/data*[.fmt|.FMT]
# DEPENDENCY: docker installed on the system
# HELP: https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
# --------------------------------------------------------------------------------- #

# print usage
function usage() {
    echo "USAGE: sudo raw2ms2 [fmt] [absolute/path/to/data*.fmt|.FMT]"
    echo "DEPENDENCY: docker installed on the system"
    echo "HELP: https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses"
}

# make sure of sudo acess
if [ "$(id -u)" != "0" ] ; then
   echo "ERROR: Run this script as sudo"
   usage
   exit 0
fi

# if not then use current directory
if [ -z "$1" ]; then
    usage
    exit 0
fi

# if not then use current directory
if [ -z "$2" ]; then 
    usage
    exit 0
fi

# --------------------------------------------------------------------------------- #

# input format
FMT=$1

# UPPERCASE the input format
FMT2=${FMT^^}

# get data directory
WDIR=$2

# make a list of all the .fmt and .FMT files
ls -rt -d -1 "$WDIR"/{*,.*} | grep -E "\.${FMT}|\.${FMT2}" > $WDIR/list

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
