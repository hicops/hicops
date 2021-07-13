#!@PYTHON_EXECUTABLE@

# 
# Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
# Florida International University, Miami, FL# 
# This program is licensed under the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
# 

# ## Import Packages

import os
import re
import sys
import math
import time
import glob
import numpy as np
import pandas as pd
import operator
import argparse

#
# main function
#
if __name__ == '__main__':
    # parse user paramters
    parser = argparse.ArgumentParser(description='Extract spectra from MS2 files')
    parser.add_argument('-i', '--infile', dest='ifile', type=str, required=True,
                        help='The MS2 file to extract spectra from')

    parser.add_argument('-d', '--id', metavar='N', dest='specs', required=True, 
                        type=int, nargs='+', help='Spectra numbers in the file')

    args = parser.parse_args()

    # List of spectrum ID to extract
    ms2file = args.ifile.lstrip(' ')
    spectralist = np.array(args.specs)

    spectralist = spectralist[spectralist >= 0]

    # Check if file exists
    if not os.path.isfile(ms2file):
        print ("Error: The file %s does not exist\n", ms2file)
        exit (-1)
    
    # Extract Spectra

    # Empty list to store lines
    lines = []
    
    # Spectrum number currently being extracted
    currSpec = 0
    
    # Current spectrum number
    specno = 0
    
    # Boolean to keep the data or not
    keep = False
    
    # Extract the required data
    with open(ms2file) as file:
        for line in file:
            # Check if we want to keep the line
            if line[:2] == 'H\t':
                lines.append(line)
    
            # If a spectrum is starting
            elif line[:2] == 'S\t':
    
                # If all acquired then break
                if currSpec == len(spectralist):
                    break
    
                # If we want to keep it
                if specno in spectralist:
                    keep = True
                    currSpec += 1
                    lines.append(line)
    
                # We don't want to keep it
                else:
                    keep = False
                
                # Increment the spectrum counter
                specno += 1
    
            # Normal data line
            else:
                if keep == True:
                        lines.append(line)

    # Write data to a new file
    ms2file2 = ms2file[:-4] + '_extracted.ms2'
    print ('Writing extracted data at: ' + ms2file2)

    nfile = open(ms2file2, "w+")
    for line in lines:
        nfile.write(line)
    
    nfile.close()
