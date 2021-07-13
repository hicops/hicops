#!@PYTHON_EXECUTABLE@

# # HDF5-based I/O module for MS/MS Data
# 
# Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
# Florida International University, Miami, FL# 
# This program is licensed under the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
# 

#
# imports
#
import os
import sys
import glob
import copy
import math
import h5py
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from mpi4py import MPI

# --------------------------------- Helper functions --------------------------------- #

#
# Parse MS2 files
#
def parseMS2(file, meta, spectra, extra):
    # open the MS2 file in read only
    file = open(ms2file, 'r')
    # number of spectra in the file
    nspecs = 0
    # empty list to contain tuples
    peaks = []

    # clear the list if not already
    meta.clear()
    # clear list to contain spectra as np.arrays
    spectra.clear()

    # loop through all lines
    for line in file:
        # test for empty lines
        if not line.strip():
            continue
        # If irrelevant lines
        if line[0] == 'H':
            continue
        elif line[0] == 'S':
            if nspecs > 0:
                # copy the peaks list to spectra                
                spectra.append(np.asarray(peaks))

                # clear the peaks
                peaks.clear()

            # S    1    2    m/z
            info = line.split()
            meta.append([['S', int(info[1]), int(info[2])]])

            # add the number of spectra
            nspecs += 1

        elif line[0] == 'I' or line[0] == 'D':
            # extract extra fields only if needed and present (of course)
            if extra ==True:
                # I/D    Key    Value
                info = line.split()
                meta[-1].append([info[0]+'_'+info[1], info[2:]])

        elif line[0] == 'Z':
            # Z    charge    mass
            info = line.split()
            meta[-1].append(['mz', float(info[2]), int(info[1])])

        # else append data lines to peaks
        else:
            # split by whitespace
            peak = line.split()
            # convert to float and append as a tuple
            peaks.append(tuple(map(lambda x : float(x), peak)))

    # copy the last peaks list to spectra as well
    spectra.append(np.asarray(peaks))

    # clear the peaks
    peaks.clear()

    # close the MS2 file
    file.close()


#
# Custom convert string to bool
# Ref: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
#
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


#
# argument parser
#
def parseArgs(comm):
    # parse user parameters and broadcast
    ms2path = None
    outpath = None
    qchunk = None
    extra = None

    # parse only at rank 0
    if comm.Get_rank() == 0:
        # add description and the data model
        parser = argparse.ArgumentParser(description= 'Convert the MS2 dataset into HDF5 dataset for HiCOPS.\n\n'
                                                      '******* DATA MODEL *******\n'
                                                      'file: 1 (MS2 file)\n'
                                                      '   group: 1 (batch of spectra)\n'
                                                      '        attributes (metadata)\n'
                                                      '\n'
                                                      '        dataset: 1 (spectrum)\n'
                                                      '            attrs\n'
                                                      '        dataset: 2\n'
                                                      '            attrs\n'
                                                      '        |\n'
                                                      '        dataset: [CHUNK]\n'
                                                      '            attrs\n'
                                                      '    group: 2\n'
                                                      '    |\n'
                                                      '    group: j\n'
                                                      '\n'
                                                      'file: 2\n'
                                                      '|\n'
                                                      'file: n\n'
                                                      '******* \\DATA MODEL *******',
                                            formatter_class = RawTextHelpFormatter,
                                            epilog =  'MPI PARALLEL usage: mpirun -n [N] [OPTIONS] %(prog)s [OPTIONS] [ARGS]')

        # make a new argument group
        requiredNamed = parser.add_argument_group('required arguments')

        # argument for input MS2 data
        requiredNamed.add_argument('-i', '--idir', dest='idir', type=str, required=True,
                        help='path to the MS2 dataset')

        # argument for output HDF5 data
        parser.add_argument('-o', '--odir', dest='odir', type=str, required=False,
                        help='path to the output HDF5 dataset, default: input_dir/hdf5')

        # argument for chunk size
        parser.add_argument('-q', '--chunk', '--group', dest='chunk', type=int, required=False, default=10000,
                        help='max spectra group size, default: 10000')

        # argument for chunk size
        parser.add_argument('-e', '--extra', dest='extra', action='store_true',
                        help='extract I and D fields from spectra, default: False')

        # parse arguments
        args = parser.parse_args()

        # path to input MS2 data (required)
        ms2path = os.path.abspath(args.idir.strip())

        # path to output hdf5 data (default to input/hdf5)
        if args.odir is not None:
            outpath = os.path.abspath(args.odir.strip())
        else:
            outpath = ms2path + '/hdf5'

        # dataset chunk size per H5group
        qchunk = args.chunk

        # extract extra fields from spectra
        extra = args.extra

        # make the output dir if does not exist
        print ('Creating output at: {}'.format(outpath), flush=True)
        os.makedirs(outpath, exist_ok=True)

    # broadcast user settings to all processes
    ms2path = comm.bcast(ms2path, root=0)
    outpath = comm.bcast(outpath, root=0)
    qchunk = comm.bcast(qchunk, root=0)
    extra  = comm.bcast(extra, root=0)

    return ms2path, outpath, qchunk, extra


# --------------------------------- The parallel MPI main function --------------------------------- #

#
# MPI4Py and HDF5 ***core*** driver to convert a MS2 dataset 
# into HDF5 dataset in parallel
#
# Parallelization Scheme: Distribute individual MS2 files among 
#                         parallel processes in a round-robin fashion
#
if __name__ == '__main__':

    # MPI settings
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    world = comm.Get_size()

    # parse user parameters and broadcast
    ms2path = None
    outpath = None
    qchunk = None
    extra = None

    # parse input parameters
    ms2path, outpath, qchunk, extra = parseArgs(comm)

    # get all MS2 files in here
    files = []

    # read all MS2 files from ms2path
    os.chdir(ms2path)
    for ms2file in glob.glob("*.ms2"):
        files.append(ms2file)

    # sanity check
    if len(files) < 1:
        print("FATAL: No MS2 file found at:" + ms2path)
        exit(-1)

    # read all .MS2 files in ms2path
    nfiles = len(files)

    # files to process locally
    myfiles = []

    # if world size > 1
    if world > 1:
       # extract my files and add to myfiles
        fno = rank
        while fno < nfiles:
            myfiles.append(files[fno])
            fno += world
    else:
        myfiles = files

    #
    # parse all local MS2 files and create HDF5 output
    #
    for ms2file in myfiles:
        # read the file into spectra array
        spectra = []
        # read the metadata into meta
        meta = []

        print('Rank: {}, Processing {}'.format(rank, ms2file))

        # parse the file
        parseMS2(ms2path + '/' + ms2file, meta, spectra, extra)

        # get the maxmimum length of spectra in file.
        longest_spec = max(len(x) for x in spectra)

        # number of batches (groups) in the spectra file
        ngroups = math.ceil(len(spectra)/float(qchunk))

        # global index of batches (groups)
        gindex = 0

        # global index of spectra
        gspectra = 0

        # 
        # Create the HDF5 configuration using tuple type to store the spectra.
        #
        ofilename = outpath + '/' + ms2file[:-4] + '.h5'
        f = h5py.File(ofilename, mode='w')

        # remaining spectra
        remaining_specs = len(spectra)
        curr_specs = 0
        # total_specs = len(spectra)

        # batch attributes
        # nbatches = file.keys().shape()
        # spectra in batch = specs
        # total spectra = len(spectra) - do we really need that?
        # longest spectrum size = maxlen

        # create groups for batches
        for i in range(ngroups):
            g = f.create_group('g_' + str(i))
            specs = qchunk
            if remaining_specs <= qchunk:
                specs = remaining_specs
                remaining_specs = 0
            else:
                remaining_specs -= qchunk

            # add the attributes to groups
            #g.attrs.create('gindex', gindex + i, dtype='int32')
            #g.attrs.create('gspectra', gspectra, dtype='int32')
            g.attrs.create('nspecs', specs, dtype='int32')

            # TODO: add the following attrs only in the first batch of file?
            # if i == 0:
            # g.attrs.create('total_specs', total_specs, dtype='int32')
            g.attrs.create('maxlen', longest_spec, dtype='int32')


            #
            # Create datasets from MS2 data
            #
            for j in range(len(spectra[curr_specs:curr_specs + specs])):
                # create a dataset with spectral data
                si = g.create_dataset(str(j+curr_specs), data=spectra[j+curr_specs], dtype='float32')

                #
                # add dataset attributes: 
                # scan number: S 1 1
                # mass charge: mz
                # I_ and D_ attributes (--extra)
                #
                for info in meta[j+curr_specs]:
                    si.attrs.create(info[0], info[1:])

            # update curr_specs
            curr_specs += specs

            # update gspectra
            gspectra += specs

        # update the global group index
        gindex += ngroups

        #
        # close the HDF5 file
        #
        f.close()

    # Wait for all parallel processes to complete
    comm.Barrier()

    # flush the stdio
    sys.stdout.flush()

    # print status
    if rank == 0:
        # print help
        print('\nSUCCESS\n')

        # print more information
        print('INFO: Use the: "h5dump [-h] [file.h5]" to read HDF5 files')