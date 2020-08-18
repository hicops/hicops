#!@PYTHON_EXECUTABLE@
#  This file is a part of HiCOPS software
#
#  Copyright (C) 2019 Parallel Computing and Data Science (PCDS) Laboratory
#                       School of Computing and Information Sciences
#                         Florida International University   (FIU)
#                          Authors: Muhammad Haseeb, Fahad Saeed
#                            Email: {mhaseeb, fsaeed}@fiu.edu
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#


# Required Imports

import os
import sys
import math
import glob
import time
import os.path
import datetime
import filecmp
from shutil import copyfile

# The main function
if __name__ == '__main__':

    # Check the Python version
    pyversion = float(str(sys.version_info[0]) + '.' + str(sys.version_info[1]))
    if (pyversion < 3.5):
        print ('\nERROR: This software requries Python v3.5 or later')
        print (' Your Python version is: ' + str(pyversion) + '\n')
        exit(-1)

    if len(sys.argv) > 1:
        paramfile = sys.argv[1]
    else:
        print ("ERROR: Missing arguments to the HiCOPS_COMET")
        print ("USAGE: python3.5+ hicops_config <params>")
        print ("GENERATE: sampleparams.txt using: \n")
        print (" $ hicops_config -g")
        sys.exit(-1)

    # Generate the sampleparams.txt file
    if (paramfile == '-g'):
        sample = open("./sampleparams.txt","w+")

        sample.write('# \n')
        sample.write('# HiCOPS for Comet\n')
        sample.write('# Copyrights(c) 2020 PCDS Laboratory\n')
        sample.write('# Muhammad Haseeb, and Fahad Saeed\n')
        sample.write('# School of Computing and Information Sciences\n')
        sample.write('# Florida International University (FIU), Miami, FL\n')
        sample.write('# Email: {mhaseeb, fsaeed}@fiu.edu\n')
        sample.write('# \n')
        sample.write('# Auto generated sampleparams.txt\n')
        sample.write('# Sample parameters generated for HiCOPS (XSEDE COMET) version\n')
        sample.write('# For more information: https://portal.xsede.org/sdsc-comet\n')
        sample.write('# \n')
        sample.write('# Generated on: ' + (datetime.datetime.now()).strftime("%Y-%m-%d %H:%M %Z") + '\n')
        sample.write('# \n')
        sample.write('# IMPORTANT: DO NOT put any spaces between variable=value\n')
        sample.write('# \n\n')
        
        sample.write('# Path (absolute or relative) to Workspace directory\n')
        sample.write('workspace=/path/to/workspace\n\n')

        sample.write('# OpenMP Multithreads per MPI process\n')
        sample.write('cores=24\n\n')

        sample.write('# ABSOLUTE path to processed proteome database parts\n')
        sample.write('# NOTE: Run the dbprocessor.py to process a FASTA \n')
        sample.write('#       proteome database and generate dbparts\n')
        sample.write('dbparts=/path/to/processed/database/parts\n\n')

        sample.write('# ABSOLUTE path to MS/MS dataset\n')
        sample.write('ms2data=/path/to/msms/dataset\n\n')

        sample.write('# Mods to include per peptide sequence\n')
        sample.write('nmods=3\n\n')

        sample.write('# Mods Information: AA(max 4) mod_mass mods_per_pep\n')
        sample.write('mod1=M 15.99 2\n')
        sample.write('mod2=X 0 0\n')
        sample.write('mod3=X 0 0\n')
        sample.write('mod4=X 0 0\n')
        sample.write('mod5=X 0 0\n')
        sample.write('mod6=X 0 0\n')
        sample.write('mod7=X 0 0\n')
        sample.write('mod8=X 0 0\n')
        sample.write('mod9=X 0 0\n')
        sample.write('mod10=X 0 0\n')
        sample.write('mod11=X 0 0\n')
        sample.write('mod12=X 0 0\n')
        sample.write('mod13=X 0 0\n')
        sample.write('mod14=X 0 0\n')
        sample.write('mod15=X 0 0\n\n')

        sample.write('# Min peptide length\n')
        sample.write('min_length=7\n\n')

        sample.write('# Max peptide length\n')
        sample.write('max_length=40\n\n')

        sample.write('# Min precursor mass (Da)\n')
        sample.write('min_prec_mass=500\n\n')

        sample.write('# Max precursor mass (Da)\n')
        sample.write('max_prec_mass=5000\n\n')

        sample.write('# Max fragment charge\n')
        sample.write('maxz=3\n\n')

        sample.write('# Min shared peaks\n')
        sample.write('shp=4\n\n')

        sample.write('# Required min shared peaks\n')
        sample.write('min_hits=4\n\n')

        sample.write('# Base normalized Intensity for MS/MS data \n')
        sample.write('base_int=100000\n\n')

        sample.write('# Cutoff ratio wrt base intensity \n')
        sample.write('cutoff_ratio=0.01\n\n')

        sample.write('# Scratch pad memory for scorecard in MBs (min: 2048MB)\n')
        sample.write('spadmem=2048\n\n')

        sample.write('# Resolution (Da)\n')
        sample.write('res=0.01\n\n')

        sample.write('# Precursor Mass Tolerance (+-Da): -1 means infinity \n')
        sample.write('dM=500\n\n')

        sample.write('# Fragment Mass Tolerance (+-Da)\n')
        sample.write('dF=0.01\n\n')

        sample.write('# Top Matches to report\n')
        sample.write('top_matches=10\n\n')

        sample.write('# Max expect value to report\n')
        sample.write('expect_max=20.0\n')


        print('Generated: ./sampleparams.txt')
        print ("\nSUCCESS")
        
        sys.exit(0)

#
# ----------------------------------------------------------------------------------------------------
#

    # Initialize the parameters
    cores = 24
    expect_max = 20.0
    dbparts = ''
    ms2data = ''
    nmods = 0
    madded = 0
    mods = []
    min_length = 7
    max_length = 40
    maxz       = 3
    dF = 0.01
    dM = 500
    res = 0.01
    scale = int(1/res)
    min_prec_mass = 500
    max_prec_mass = 5000
    top_matches = 10
    shp_cnt = 4
    min_hits = 4
    base_int = 100000
    cutoff_ratio = 0.01
    workspace = './workspace'
    policy = 'cyclic'
    spadmem = 2048
    hicopspath = os.path.dirname(os.path.realpath(__file__))
    newparams = False
    uparams = ''

#
# ----------------------------------------------------------------------------------------------------
#

    print ('\n*****************************')
    print   ('*  HiCOPS Configuration Gen *')
    print   ('*    PCDS Lab, SCIS, FIU    *')
    print   ('*****************************\n')

    print ('Provided Parameters:')

    # Parse the params file
    with open(paramfile) as params:
        for line in params:

            # Ignore the empty or comment lines
            if (line[0] == '\r' or line[0] == '#' or line[0] == '\n'):
                continue

            # Split line into param and value
            param, val = line.split('=', 1)

            # Set database file 
            if (param == 'dbparts'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                dbparts = os.path.abspath(val)
                if (os.path.exists(dbparts) == False):
                    print ("ERROR: Enter valid path to processed proteome database parts directory")
                    sys.exit(-2)
                else:
                    print ('Processed DB Parts =', dbparts)

            # Set the ms2data path
            elif (param == 'ms2data'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                ms2data = os.path.abspath(val)
                if (os.path.exists(ms2data) == False):
                    print ("ERROR: Enter valid path to MS2 dataset")
                    sys.exit(-3)
                else:
                    print ('MS/MS dataset =', ms2data)

            # Cores per node
            elif (param == 'cores'):
                cores = int(val)
                if (cores <= 0 or cores > 44):
                    cores = 44
                print ('Using cores/node  =', cores)

            # Set max mods
            elif (param == 'nmods'):
                nmods = int (val)
                if (nmods < 0):
                    nmods = 0 
                if (nmods > 8):
                    nmods = 8 
                print ('Max mods/pep  =', nmods)

            # There is a \n at the end of each string
            elif (param[:-1] == 'mod'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                if (val != 'X 0 0'):
                    mods.append(val)
                    print ('Adding mod   =', str(val))
                    madded += 1

            # Set the min digestion length
            elif (param == 'min_length'):
                min_length = int(val)

                if (min_length < 4):
                    min_length = 4 
                if (min_length > 60):
                    min_length = 60 
                print ('Min pep len  =', min_length)

            # Set the max digestion length
            elif (param == 'max_length'):
                max_length = int(val)
                if (max_length < 4):
                    max_length = 4 
                if (max_length > 60):
                    max_length = 60
                print ('Max pep len  =', max_length)

            # Set the max fragment ion charge
            elif (param == 'maxz'):
                maxz = int(val)
                if (maxz < 1):
                    maxz = 1 
                if (maxz > 5):
                    maxz = 5
                print ('Max frag chg =', maxz)

            # Fragment mass tolerance
            elif (param == 'dF'):
                if (dF < 0.0):
                    dF = 0.0
                if (dF > 0.02):
                    dF = 0.02
                dF = float(val)
                print ('dF           =', dF)

            # Peptide precursor mass tolerance
            elif (param == 'dM'):
                dM = float(val)
                if (dM < 0.001):
                    dM = 0.001 
                print ('dM           =', dM)

            # m/z axis resolution
            elif (param == 'res'):
                res = float(val)
                if (res <= 0):
                    res = 0.01 
                if (res > 5.0):
                    res = 5.0
                print ('resolution   =', res)

            # Minimum precursor mass
            elif (param == 'min_prec_mass'):
                min_prec_mass = int(val)
                if (min_prec_mass < 0):
                    min_prec_mass = 0 
                if (min_prec_mass > 10000):
                    min_prec_mass = 10000
                print ('min_prec_mass =', min_prec_mass)

            # Maximum precursor mass
            elif (param == 'max_prec_mass'):
                max_prec_mass = int(val)
                if (max_prec_mass < 0):
                    max_prec_mass = 0 
                if (max_prec_mass > 10000):
                    max_prec_mass = 10000
                print ('max_prec_mass =', max_prec_mass)

            # Minimum Shared Peaks
            elif (param == 'shp_cnt'):
                shp_cnt = int(val)
                if (shp_cnt < 1):
                    shp_cnt = 1 
                if (shp_cnt > 20):
                    shp_cnt = 20
                print ('Min Shared Peaks =', shp_cnt)
                
            # Minimum required PSM hits
            elif (param == 'min_hits'):
                min_hits = int(val)
                if (min_hits < 4):
                    min_hits = 4
                print ('Required min PSM hits =', min_hits)

            # Base Normalized Intensity
            elif (param == 'base_int'):
                base_int = int(val)
                if (base_int < 1000):
                    base_int = 1000

                print ('Base Normalized Intensity =', base_int)
                
            # Intensity Cutoff Ratio
            elif (param == 'cutoff_ratio'):
                cutoff_ratio = float(val)
                if (cutoff_ratio >= 0.10):
                    cutoff_ratio = 0.10
                if (cutoff_ratio <= 0):
                    cutoff_ratio = 0.01
                print ('Intensity Cutoff Ratio =', cutoff_ratio)

            # Scorecard memory
            elif (param == 'spadmem'):
                spadmem = int(val)
                if (spadmem < 2048):
                    spadmem = 2048
                print ('Scratch Memory (MB) =', spadmem)

            # Workspace Path
            elif (param == 'workspace'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                if (val[-1] == '/'):
                    val = val[:-1]

                workspace = os.path.abspath(str(val))
                print ('workspace   =', workspace)

            # Maximum precursor mass
            elif (param == 'top_matches'):
                top_matches = int(val)
                if (top_matches < 1):
                    top_matches = 1
                print ('Top matches =', top_matches)

            elif (param == 'expect_max'):
                expect_max = float(val)
                if (expect_max < 0):
                    expect_max = 0
                if (expect_max > 100):
                    expect_max = 100
                print ('Max expect value to report =', expect_max)

#    print ('Mods Added', mods)

    if (len(mods) == 0):
        mods.append("X 0 0")
        nmods = 0

    # Close the params file
    params.close()

#
# ----------------------------------------------------------------------------------------------------
#

    # Create a workspace directory
    print ('\nInitializing Workspace at: ', workspace)

    if (os.path.exists(workspace) == False):
        os.mkdir(workspace)

    # Create the output directory for results
    if (os.path.exists(workspace + '/output') == False):
        os.mkdir(workspace + '/output')

    # Check if the params have been changed from the last run
    if (os.path.isfile(workspace + '/runtime_settings.txt') == False or filecmp.cmp(workspace + '/runtime_settings.txt', paramfile) == False):
        newparams = True
        copyfile(paramfile, workspace + '/runtime_settings.txt')

    # Sanity check
    if (min_length > max_length):
        temp = min_length
        min_length = max_length
        max_length = temp
        print('WARNING: min_length > max_length. Swapping them\n')

    # Check if all database parts are available
    for k in range(min_length, max_length + 1):
        if (os.path.isfile(dbparts + '/' + str(k) + '.peps') == False):
            print ('ABORT: Database part(s) are missing\n')
            exit(-3)

#
# ----------------------------------------------------------------------------------------------------
#

    # Prepare the uparams.txt file for PCDSFrame
    uparams = hicopspath + '/uparams.txt'
    modfile = open(uparams, "w+")

    # Write params for the CFIR index
    modfile.write(dbparts + '\n')
    modfile.write(ms2data + '\n')
    modfile.write(workspace + '/output\n')
    modfile.write(str(cores) + '\n')
    modfile.write(str(min_length) + '\n')
    modfile.write(str(max_length) + '\n')
    modfile.write(str(maxz) + '\n')
    modfile.write(str(dF) + '\n')
    modfile.write(str(dM) + '\n')
    modfile.write(str(res) + '\n')
    modfile.write(str(scale) + '\n')
    modfile.write(str(min_prec_mass) + '\n')
    modfile.write(str(max_prec_mass) + '\n')
    modfile.write(str(top_matches) + '\n')
    modfile.write(str(expect_max) + '\n')
    modfile.write(str(shp_cnt) + '\n')
    modfile.write(str(min_hits) + '\n')
    modfile.write(str(base_int) + '\n')
    modfile.write(str(cutoff_ratio) + '\n')
    modfile.write(str(spadmem) + '\n')
    modfile.write(str(policy) + '\n')

    modfile.write(str(len(mods)) + '\n')
    modfile.write(str(nmods) + '\n')
    for info in mods:
        aa,ms,num = info.split(' ', 2)
        modfile.write(aa + '\n')
        modfile.write(str(ms) + '\n')
        modfile.write(str(num) + '\n')

    modfile.close()

    print ('\nSUCCESS\n')
    print('Written:', hicopspath + '/uparams.txt')
    print('You can now run HiCOPS as: mpirun -np [N] [OPTIONS] '+ hicopspath + '/hicops ' + hicopspath + '/uparams.txt\n')

    print ('Issue Reporting: https://github.com/pcdslab/hicops/issues\n')

    print (' # ---------------------------------------------------------------------------------------------------- #\n\n')
