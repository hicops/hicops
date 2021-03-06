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
import filecmp
import datetime
import subprocess
from subprocess import call
from shutil import copyfile
from functools import reduce

#
# ------------------------------ Helper Functions -------------------------------------------
#

# Computes factors of a number in descending order
def factors(n):    
    return list(set(reduce(list.__add__, ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))))

# Checks if any jobs are running
def checkRunningJobs(username):
    squeue = subprocess.run('squeue -u ' + username + ' | wc -l', stdout=subprocess.PIPE, shell=True)
    proc = int(squeue.stdout.decode('utf-8'))
    if (proc == 1):
        return False
    else:
        return True

# Generates a normal unicore job script
def genNormalScript(wkspc, jobname, outname, partition, nds, ntask_per_node, minust, comd, mail, user, evts):
    script = open(wkspc + '/autogen/' + jobname, 'w+')
    script.write('#!/bin/bash\n')
    script.write('\n')
    script.write('#SBATCH --job-name="' + jobname +'"\n')
    script.write('#SBATCH --output="' + wkspc + '/autogen/' + outname + '.out"\n')
    script.write('#SBATCH --partition=' + partition + '\n')
    script.write('#SBATCH --nodes=' + nds + '\n')
    script.write('#SBATCH --ntasks-per-node=' + ntask_per_node + '\n')
    script.write('#SBATCH --export=ALL\n')
    script.write('#SBATCH -t ' + minust + '\n')

    if (mail):
        script.write('#SBATCH --mail-type=' + evts + '\n')
        script.write('#SBATCH --mail-type=' + user + '\n')

    script.write('\n')
    script.write(comd + '\n')

    return

# Generates a multithreaded OpenMP job script
def genOpenMPScript(wkspc, jobname, outname, partition, nds, ntask_per_node, minust, ompthrds, command, args, mail, user, evts):
    script = open(wkspc + '/autogen/' + jobname, 'w+')
    script.write('#!/bin/bash\n')
    script.write('\n')
    script.write('#SBATCH --job-name="' + jobname +'"\n')
    script.write('#SBATCH --output="' + wkspc + '/autogen/' + outname + '.out"\n')
    script.write('#SBATCH --partition=' + partition + '\n')
    script.write('#SBATCH --nodes=' + nds + '\n')
    script.write('#SBATCH --ntasks-per-node=' + ntask_per_node + '\n')
    script.write('#SBATCH --export=ALL\n')
    script.write('#SBATCH -t ' + minust + '\n')
    if (mail):
        script.write('#SBATCH --mail-type=' + evts + '\n')
        script.write('#SBATCH --mail-type=' + user + '\n')

    script.write('\n')
    script.write ('export OMP_NUM_THREADS      ' + ompthrds + '\n')
    script.write('\n')
    script.write(command + ' ' + args)

    return

# Generates a Hybrid MPI/OpenMP job script
def genMPI_OpenMPScript(wkspc, jobname, outname, partition, nds, ntask_per_node, minust, ompthrds, command, npernode, blevel, bpolicy, args, mail, user, evts):
    script = open(wkspc + '/autogen/' + jobname, 'w+')
    script.write('#!/bin/bash\n')
    script.write('\n')
    script.write('#SBATCH --job-name="' + jobname +'"\n')
    script.write('#SBATCH --output="' + wkspc + '/output/' + outname + '.%j.%N.out"\n')
    script.write('#SBATCH --partition=' + partition + '\n')
    script.write('#SBATCH --nodes=' + nds + '\n')
    script.write('#SBATCH --ntasks-per-node=' + ntask_per_node + '\n')
    script.write('#SBATCH --export=ALL\n')
    script.write('#SBATCH -t ' + minust + '\n')
    if (mail):
        script.write('#SBATCH --mail-type=' + evts + '\n')
        script.write('#SBATCH --mail-type=' + user + '\n')
    script.write('\n')
    script.write ('export OMP_NUM_THREADS      ' + ompthrds + '\n')
    script.write('\n')
    script.write('ibrun --npernode ' + npernode + ' -bl ' + blevel + ' -bp ' + bpolicy + ' ' + command + ' ' + args)

    return

#
# ------------------------------ Main Function -------------------------------------------
#

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
        print ("USAGE: python3.5+ hicops_comet.py <params>")
        print ("GENERATE: A sample params.txt using python3.5+ hicops_comet.py -g")
        sys.exit(-1)

    # Generate the sampleparams.txt file
    if (paramfile == '-g'):
        # get user name to make workspace path
        username = os.environ['USER']

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
        sample.write('# \n\n\n')
        sample.write('# IMPORTANT: DO NOT put any spaces between variable=value\n')
        sample.write('# \n\n')

        sample.write('# XSEDE (Comet) Username\n')
        sample.write('username='+ username + '\n\n')

        sample.write('# Get emails for job? 1/0? \n')
        sample.write('mail=0\n')

        sample.write('# Events to get emails? Options: BEGIN, END, FAIL or ALL\n')
        sample.write('mailtype=FAIL\n')

        sample.write('# ABSOLUTE Path to Workspace directory\n')
        sample.write('workspace=/oasis/scratch/comet/'+ username + '/temp_project/hicops_workspace\n\n')

        sample.write('# Job Time: hh:mm:ss (max: 48:00:00)\n')
        sample.write('jobtime=00:45:00\n\n')

        sample.write('# Nodes available\n')
        sample.write('nodes=2\n\n')

        sample.write('# Cores per machine\n')
        sample.write('cores=24\n\n')

        sample.write('# Cores per MPI process\n')
        sample.write('cores_per_mpi=12\n\n')

        sample.write('# MPI binding policy: scatter, compact \n')
        sample.write('bp=scatter\n\n')

        sample.write('# MPI binding level: core, socket, numanode\n')
        sample.write('bl=socket\n\n')

        sample.write('# Optimize settings? 1/0? Highly Recommended \n')
        sample.write('optimize=1\n\n')

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
# ------------------------------ Initialization -------------------------------------------
#

    # Initialize the parameters
    cores = 24
    sockets = 2
    numa = 2
    nodes = 2
    expect_max = 20.0
    numamem = math.inf
    mpi_per_node = sockets
    cores_per_socket = int(cores/sockets)
    cores_per_numa = int (cores/numa)
    threads = cores_per_socket
    prep_threads = int(threads/4)
    bp = 'scatter'
    bl = 'socket'
    optimize = 1
    dbparts = ''
    ms2data = ''
    nmods = 0
    madded = 0
    mods = []
    min_length = 6
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
    workspace = '/oasis/scratch/comet/$USER/temp_project/hicops_workspace'
    policy = 'cyclic'
    spadmem = 2048
    indexsize = 0
    nions     = 0
    size_mb   = 0
    mb_per_numa = 0
    mb_per_mpi  = 0
    nparts      = 0
    jobtime='00:45:00'
    hicopspath = os.path.dirname(os.path.realpath(__file__)) + '/..'
    newparams = False
    username = os.environ['USER']
    uparams = ''
    pparams = ''
    mail = False
    evts = 'FAIL'
    possible_evts = ['BEGIN', 'END', 'ALL', 'FAIL']

#
# ------------------------------ Parse Parameters -------------------------------------------
#

    print ('\n***************************')
    print   ('*  HiCOPS for XSEDE Comet *')
    print   ('*   PCDS Lab, SCIS, FIU   *')
    print   ('***************************\n')

    print ('Provided Parameters:')

    # Parse the params file
    with open(paramfile) as params:
        for line in params:

            # Ignore the empty or comment lines
            if (line[0] == '\r' or line[0] == '#' or line[0] == '\n'):
                continue

            # Split line into param and value
            param, val = line.split('=', 1)

            # Set XSEDE username 
            if (param == 'username'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                username = val

            # Set up emails
            elif (param == 'mail'):
                mail = int(val)
                if (mail <= 0):
                    mail = 0
                if (mail > 0):
                    mail = 1
                print ('Emails =', mail)

            # Email events
            if (param == 'dbparts'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                evts = val
                if not evts in possible_evts:
                    evts = 'FAIL'
                print ('Email events =', evts)

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

            # Set the job time
            elif (param == 'jobtime'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                hh,mm,ss = map(int, val.split(':',2))
                if (hh == 0 and mm == 0 and ss == 0):
                    val = '00:45:00'
                else:
                    jobtime = val
                print ('Job time =', jobtime)

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

            # Set max nodes in the system [1,72] on COMET
            elif (param == 'nodes'):
                nodes = int(val)
                if (nodes <= 0):
                    nodes = 1
                if (nodes > 72):
                    nodes = 72
                print ('Using nodes =', nodes)

            # Cores per node
            elif (param == 'cores'):
                cores = int(val)
                if (cores <= 0 or cores > 24):
                    cores = 24
                print ('Using cores/node  =', cores)

            # Optimize number of cores and MPI processes to run?
            elif (param == 'optimize'):
                optimize = int(val)
                if (optimize <= 0):
                    optimize = 0
                if (optimize > 0):
                    optimize = 1
                print ('Optimizations =', optimize)

            # Set the MPI binding level
            elif (param == 'bl'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                if (val == 'socket' or val == 'numanode' or val == 'core'):
                    bl = val
                print ('MPI binding level =', bl)

                # Set the MPI binding policy
            elif (param == 'bp'):
                if (val[-1] == '\n'):
                    val = val[:-1]
                if (val[-1] == '\r'):
                    val = val[:-1]

                if (val == 'scatter' or val == 'compact'):
                    bp = val
                print ('MPI binding policy =', bp)

            # Set OMP cores per MPI
            elif (param == 'cores_per_mpi'):
                threads = int(val)
                if (threads <= 0 or threads > cores):
                    threads = int(cores)

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
                print ('Fragment mass tolerance (dF) =', dF)

            # Peptide precursor mass tolerance
            elif (param == 'dM'):
                dM = float(val)
                if (dM < 0.001):
                    dM = 0.001 
                print ('Peptide mass tolerance (dM) =', dM)

            # m/z axis resolution
            elif (param == 'res'):
                res = float(val)
                if (res <= 0):
                    res = 0.01 
                if (res > 5.0):
                    res = 5.0
                print ('Resolution   =', res)

            # Minimum precursor mass
            elif (param == 'min_prec_mass'):
                min_prec_mass = int(val)
                if (min_prec_mass < 0):
                    min_prec_mass = 0 
                if (min_prec_mass > 10000):
                    min_prec_mass = 10000
                print ('Min precursor mass =', min_prec_mass)

            # Maximum precursor mass
            elif (param == 'max_prec_mass'):
                max_prec_mass = int(val)
                if (max_prec_mass < 0):
                    max_prec_mass = 0 
                if (max_prec_mass > 10000):
                    max_prec_mass = 10000
                print ('Max precursor mass =', max_prec_mass)

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
                print ('Workspace   =', workspace)

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
# ------------------------------ Workspace Creation -------------------------------------------
#

    # Create a workspace directory
    print ('\nInitializing Workspace at: ', workspace)

    os.makedirs(workspace, exist_ok=True)

    # Create the output directory for results
    os.makedirs(workspace + '/output', exist_ok=True)

    # Create directory where autogen stuff will be placed
    os.makedirs(workspace + '/autogen', exist_ok=True)

    # Check if the params have been changed from the last run
    if (os.path.isfile(workspace + '/autogen/settings.txt') == False or filecmp.cmp(workspace + '/autogen/settings.txt', paramfile) == False):
        newparams = True
        copyfile(paramfile, workspace + '/autogen/settings.txt')

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
# ------------------------------ Optimization -------------------------------------------
#

    # Optimizer
    if (optimize == 1):
        print ("\n\n****** Optimizing parameters *******\n")

        # Call the lsinfo to gather CPU information
        if (os.path.isfile(workspace + '/autogen/info.out') == False):
            genNormalScript(workspace, 'info', 'info', 'compute', '1','1', '00:00:10', 'lscpu | tr -d " \\r" && numactl --hardware | tr -d " \\r"', False, username, evts)

            optimize = call('sbatch ' + workspace + '/autogen/info', shell=True)
            print ('\nWaiting for job scheduler\n')

        # Wait for the lscpu process to complete 
        while (os.path.isfile(workspace + '/autogen/info.out') == False or checkRunningJobs(username) == True):
            time.sleep(0.5)

        print ('\nExtracted System Settings\n')

        # Parse the machine info file
        with open(workspace + '/autogen/info.out') as minfo:
            for line in minfo:

                # Ignore the empty or comment lines
                if (line[0] == '\r' or line[0] == '#' or line[0] == '\n'):
                    continue

                # Split line into param and value
                param, val = line.split(':', 1)        

                # Set the sockets per node
                if (param == 'Socket(s)'):
                    sockets = int(val)
                    print ('Available sockets/machine  =', sockets)

                elif (param == 'CPU(s)'):
                    cores = int(val)
                    print ('Available cores/machine  =', cores)

                elif (param == 'Core(s)persocket'):
                    cores_per_socket = int(val)
                    print ('Available cores/socket  =', cores_per_socket)

                elif (param == 'NUMAnode(s)'):
                    numa = int(val)
                    print ('Available NUMA nodes/machine =', numa)

                # Get the available NUMA memory
                elif (param[:4] == 'node' and param[-4:] == 'free'):
                    mem = int(val[:-3]) - 1024
                    if (mem < numamem):
                        numamem = mem

                elif (param == 'nodedistances'):
                    break

        print ('Available max NUMA memory (- 1GB ) =', numamem)

        cores_per_numa = int(cores/numa)

        minfo.close()

        # Check if params file was modified
        if (newparams == True or os.path.isfile(workspace + '/autogen/counter.out') == False):

            # Prepare the pparams.txt file for seq generator
            pparams = workspace + '/autogen/pparams.txt'

            modfile = open(pparams, "w+")

            # Write params for the CFIR index
            modfile.write(dbparts + '\n')
            modfile.write(ms2data + '\n')
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
            modfile.write(str(shp_cnt) + '\n')
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

            print ('\n')

            # Remove the previous counter file
            if (os.path.isfile(workspace + '/autogen/counter.out')):
                os.remove(workspace + '/autogen/counter.out')

            # Make counter.exe if missing
            if (os.path.isfile(hicopspath + '/counter') == False):
                cleancntr = call("make -C counter allclean", shell=True)
                makecntr = call("make -C counter", shell=True)
            
            genOpenMPScript(workspace, 'counter', 'counter', 'compute', '1', str(cores), '00:30:00', str(cores), hicopspath + '/counter', pparams, False, username, evts)

            # Call the counter process            
            optimize3 = call('sbatch ' + workspace + '/autogen/counter', shell=True)

        # Wait for the counter process to complete
        while (os.path.isfile(workspace + '/autogen/counter.out') == False or checkRunningJobs(username) == True):
            time.sleep(0.5)

        print ('\nEstimating Index Size\n')

        # Parse the index size file
        with open(workspace + '/autogen/counter.out') as minfo:
            for line in minfo:

                # Ignore the empty or comment lines
                if (line[0] == '\r' or line[0] == '#' or line[0] == '\n' or line[0] == '/'):
                    continue

                # Split line into param and value
                param, val = line.split(':', 1)    

                if (param == 'spectra'):
                    indexsize = int(val)
                    print ('Estimated Index Size (x 1E6 Spectra) =', int(indexsize/(1000 * 1000)))

                elif (param == 'ions'):
                    nions = int(val)

        minfo.close()

        # Remove the temporary pparams.txt
        if (os.path.isfile(pparams)):
            os.remove(pparams)
        
        if (nions == 0 or indexsize == 0):
            print ('\nFATAL: counter.exe failed. Please check the ./cnterr.txt\n')

            if (os.path.isfile(workspace + '/autogen/counter.out') == True):
                copyfile(workspace + '/autogen/counter.out', hicopspath + '/cnterr.txt')
                os.remove(workspace + '/autogen/counter.out')

            # Exit abnormally
            exit(-3)

        print ('\n')

#
# ------------------------------ Apply Optimizations -------------------------------------------
#

        # Apply the optimizations 
            
        # Case 1: Sockets >= NUMA nodes (one or multiple sockets/NUMA)
        if (sockets >= numa):

            # Set the BL to socket, BP to scatter, mpi_per_node to sockets, and threads_per_mpi to cores_per_socket
            threads = cores_per_socket
            mpi_per_node = sockets
        

            # Estimate size of index in MBs
            size_mb = ((nions * 4 + (mpi_per_node * nodes * max_prec_mass * scale * 4)) / (1024 * 1024))  + (spadmem * mpi_per_node * nodes)

            size_per_mpi = size_mb/mpi_per_node
            
            bl = 'socket'
            bp = 'scatter'

        # Case 2: Socket mapped to multiple NUMA nodes
        elif (sockets < numa):
            threads = int(cores_per_socket/numa)
            mpi_per_node = int(sockets * numa)

            # Estimate size of index in MBs
            size_mb = ((nions * 4 + (mpi_per_node * nodes * max_prec_mass * scale * 4)) / (1024 * 1024))  + (spadmem * mpi_per_node * nodes)

            bl = 'numanode'
            bp = 'scatter'
        
        # Optimize based on the index size (in spectra) per MPI
        # If partition size > 25 million, then increase number of partitions
        min_threads = 6
        max_mpi_per_node = cores / min_threads

        if (indexsize/(mpi_per_node * nodes) > 25E6):
            
            # Get set of factors
            possible = factors(threads)
            possible = sorted(possible, reverse=True)

            for cc in possible:

                if (indexsize/(mpi_per_node * nodes) <= 48E6 or cc < min_threads or mpi_per_node > max_mpi_per_node):
                    break
                else:
                    threads = cc
                    mpi_per_node = int(cores/cc)

        # Get MBs per NUMA node
        mbs_per_numa = size_mb/(numa * nodes)

        # We hope this never happens :)
        if (mbs_per_numa > numamem):
            print ('WARNING: Memory required per NUMA node = ' + str(mbs_per_numa) + ' >' + str(numamem) + ' = available NUMA mem\n')
            print ('         Either increase the number of nodes or expect performance degradation due to NUMA access and page faults\n')

        # if very small index then the preprocessing threads may be increased to 50%
        if (indexsize/(mpi_per_node * nodes) < 10E6 or dM < 50):
            prep_threads = int(threads/3)

        print('Optimized HiCOPS settings...\n')
        print('Setting Threads/MPI =', threads)
        print('Setting Max Prep threads/MPI =', prep_threads)
        print('Setting MPI/machine =', mpi_per_node)
        print('Setting MPI Policy  =', bp)
        print('Setting MPI Binding =', bl)
        print('Setting Index / MPI =', int(indexsize/(mpi_per_node * nodes)))

        print('\nSUCCESS\n')

#
# ------------------------------ Write uparams.txt -------------------------------------------
#

    # Prepare the uparams.txt file for PCDSFrame
    uparams = workspace + '/autogen/uparams.txt'
    modfile = open(uparams, "w+")

    # Write params for the CFIR index
    modfile.write(dbparts + '\n')
    modfile.write(ms2data + '\n')
    modfile.write(workspace + '/output\n')
    modfile.write(str(threads) + '\n')
    modfile.write(str(prep_threads) + '\n')
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

    # Generate the job script
    genMPI_OpenMPScript(workspace, 'hicops', 'hicops', 'compute', str(nodes), str(cores), jobtime, str(threads), hicopspath + '/hicops', str(mpi_per_node), bl, bp, uparams, mail, username, evts)

    # Generate the post-processing script
    genNormalScript(workspace, 'postprocess', 'postprocess', 'shared', '1', '1', "00:20:00", 'psm2excel ' + workspace + '/output', mail, username, evts)
#
# ------------------------------ Schedule HiCOPS job -------------------------------------------
#

    # Run HiCOPS
    cfir = call('sbatch ' + workspace + '/autogen/hicops', shell=True)

    print ('\nHiCOPS is running now\n')
    print ('You can check the job progress by: \n')
    print ('$ squeue -u $USER\n')
    print ('The output will be written at: '+ workspace + '/output')

    print ('\nSUCCESS\n')

    print ('After job completion, run:\n')
    print ('$ srun --partition=compute --nodes=1 --ntasks-per-node=1 -t 00:25:00 --export=ALL ' + hicopspath + '/tools/psm2excel -i ' + workspace + '/output')

    print ('\n')

    print ('#---------------------------------------------#')
    print ('     Read More: https://hicops.github.io')
    print ('#---------------------------------------------#\n\n')
