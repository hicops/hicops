#
# imports
#

import os
import sys
import glob
import copy
import numpy as np
import argparse

#
# main function
#

if __name__ == '__main__':

    # initialize arg parser
    parser = argparse.ArgumentParser(description='Split/merge the MS/MS data set into/from partitions')

    # data root directory
    parser.add_argument('-i', '--idir', dest='inpath', type=str, required=True,
                        help='Path to data set files')

    # data file extensions
    parser.add_argument('-e', '--ext', dest='extension', type=str, required=True,
                        help='Data set file extension e.g. ms2, mzML, mzXML, mgf')

    # number of nodes/partitions
    parser.add_argument('-N', dest='nodes', type=int, required=True,
                        help='Number of partitions: >= 2 if -m not set, else anything')

    # merge operation?
    parser.add_argument('-m', '--merge', dest='merge', action='store_true',
                        help='Merge all existing partitions: -N <anything is fine>, default: False')

    # parse arguments
    args = parser.parse_args()

    # input path to the MS2 data directory
    inpath = args.inpath.lstrip(' ')
    inpath = os.path.expanduser(inpath)

    # check if directory exists
    if not os.path.isdir(inpath):
        print ('ERROR: Directory does not exist\n')
        sys.exit (-1)

    # get file extension
    extension = args.extension.lstrip(' ')

    # extract extra fields from spectra
    merge = args.merge

    # if merge operation
    if (merge == True):
        
        # get all files in part_*/*.ext
        files = np.array(glob.glob(inpath + '/part_*/*.' + extension))
        
        # move data parts back to the root
        for file in files:
            os.rename(file, inpath + '/' + os.path.split(file)[1])
            
    # scatter operation
    else:
        # get number of nodes
        nodes = args.nodes
    
        # check for number of nodes
        if nodes < 2:
            print ('ERROR: N >= 2 required')
            sys.exit(-1)
    
        # get all files with *.ext
        files = np.array(glob.glob(inpath + '/*.' + extension))

        # check if files found
        if len(files) < 1:
            print ('ERROR: No ' + extension + ' files found\n')
            sys.exit(-1)

        # check number of nodes and files
        if nodes > len(files):
            print ('WARNING: # partitions > # files')
            print ('Some partitions may be empty\n')
    
        # shuffle the files
        np.random.shuffle(files)
    
        # split files in 'nodes' partitions
        splits = np.array_split(files, nodes)
    
        # create new directories
        for i in range(nodes):
            os.makedirs(inpath + '/part_' + str(i+1), exist_ok=True)

            # move data to parts
            for file in splits[i]:
                os.rename(file, inpath + '/part_' +  str(i+1) + '/' + os.path.split(file)[1])

    # print final status
    print ('DONE\n')