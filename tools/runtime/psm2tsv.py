#!@PYTHON_EXECUTABLE@
#   Copyright (c) 2020 Muhammad Haseeb, Fahad Saeed
#    School of Computing and Information Sciences
#      Florida International University   (FIU)
#         Email: {mhaseeb, fsaeed} @fiu.edu
# 
#  License
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