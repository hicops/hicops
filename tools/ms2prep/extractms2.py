#!@PYTHON_EXECUTABLE@
# coding: utf-8

# ####   Copyright (c) 2020 Muhammad Haseeb, Fahad Saeed
# ####    School of Computing and Information Sciences
# ####      Florida International University   (FIU)
# ####         Email: {mhaseeb, fsaeed} @fiu.edu
# 
# ### License
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

# ## Import Packages

# In[1]:


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


# ## Sanity Checking

# In[3]:

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract spectra from MS2 files')
    parser.add_argument('-i', '--infile', dest='ifile', type=str, required=True,
                        help='The MS2 file to extract spectra from')

    parser.add_argument('-d', '--id', metavar='N', dest='specs', required=True, 
                        type=int, nargs='+', help='Spectra numbers in the file')

    args = parser.parse_args()

    # List of spectrum ID to extract
    ms2file = args.ifile.lstrip(' ')
    spectralist = np.array(args.specs)

    spectralist[spectralist < 1] = 1

    # Check if file exists
    if not os.path.isfile(ms2file):
        print ("Error: The file %s does not exist\n", ms2file)
        exit (-1)
    
    # ## Extract Spectra
    
    # Empty list to store lines
    lines = []
    
    # Spectrum number currently being extracted
    currSpec = 0
    
    # Current spectrum number
    specno = 1
    
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
                if specno == spectralist[currSpec]:
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
    
    
    # ## Write data to a new file
    
    # In[29]:
    
    
    ms2file2 = ms2file[:-4] + '_extracted.ms2'
    print ('Writing extracted data at: ' + ms2file2)
    
    
    # In[30]:
    
    
    nfile = open(ms2file2, "w+")
    for line in lines:
        nfile.write(line)
    
    nfile.close()
