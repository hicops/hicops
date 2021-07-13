#!@PYTHON_EXECUTABLE@
# 
# Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
# Florida International University, Miami, FL# 
# This program is licensed under the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
# 


# Import Packages
import os
import sys
import glob
import time
import argparse
import pandas as pd

# Sanity Checking

# The main function
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Merge partial TSV files')

    parser.add_argument('-i', '--idir', dest='data_dir', type=str, required=True,
                        help='Path to partial TSV files')

    parser.add_argument('-o', '--ofile', dest='ofile', type=str, required=False, 
                        help='Path to Output file (default: idir/results.tsv')

    # parse args
    args = parser.parse_args()

    # get the data directory
    data_dir = args.data_dir.lstrip(' ')

    # get the output file name
    output = ''
    
    if args.ofile is not None:
        output = args.ofile.lstrip(' ')
    else:
        output = data_dir + '/results.tsv'

    # Open the TSV files
    data_dir = os.path.expanduser(data_dir)

    # Check if directory exists
    if not os.path.isdir(data_dir):
        print ("Directory does not exist\n")
        sys.exit (-1)
        
    # Get all files with TSV
    tsv_files = glob.glob(data_dir + '/*.tsv')
    #print (tsv_files)


    # Extract TSV data

    # Matrix where the df will be collected
    matrix = []

    # Read all TSVs into data matrix
    if (len(tsv_files) > 1):
        for kk in tsv_files:
            if kk.find('results') != -1:
                t = time.localtime()
                current_time = time.strftime("%H.%M.%S", t)
                output = data_dir + '/results.' + current_time + '.tsv'
            else:
                print ('Loading File: ', kk)
                dat = pd.read_csv(kk, sep='\t', index_col=None, header=0)
                matrix.append(dat)


    # Construct data frame
    print ("Constructing DataFrame...")

    # Concatenate data into a single data frame
    frame = pd.concat(matrix, axis=0, ignore_index=True)

    # Print the new data frame shape
    if (frame.shape[1] == 0):
        print ('ERROR: Empty data frame')
        sys.exit(-2)
    # else:
    #    print(frame.shape)

    # Write to Excel file
    print ('Writing concatenated TSV...')

    # Write the TSV file
    frame.to_csv(output, sep = '\t')

    # Delete all partial TSVs
    if (len(tsv_files) > 1):
        for kk in tsv_files:
            if kk.find('results') != -1:
                continue
            else:
                os.remove(kk)

    # Print the output address
    print ('DONE: ', output)