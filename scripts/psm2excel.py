
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

# In[3]:


import os
import re
import sys
import math
import time
import glob
import numpy as np
import pandas as pd
import operator
import matplotlib.pyplot as plt


# ## Sanity Checking

# In[9]:


# The main function
if __name__ == '__main__':

	if len(sys.argv) > 1:
		data_dir = sys.argv[1]
	else:
		data_dir = './'
		
	# Open the TSV files
	data_dir = os.path.expanduser(data_dir)

	# Check if directory exists
	if not os.path.isdir(data_dir):
		print ("Directory does not exist\n")
		sys.exit (-1)
	
	# Get all files with TSV
	tsv_files = glob.glob(data_dir + '/*.tsv')
	#print (tsv_files)
	
	
	# ## Extract TSV data
	
	# In[11]:
	
	
	# Matrix where the df will be collected
	matrix = []
	
	# Read all TSVs into data matrix
	if (len(tsv_files) > 1):
		for kk in tsv_files:
			dat = pd.read_csv(kk, sep='\t', index_col=None, header=0)
			matrix.append(dat)
			os.remove(kk)
	
	
	# ## Construct data frame
	
	# In[12]:
	
	
	# Concatenate data into a single data frame
	frame = pd.concat(matrix, axis=0, ignore_index=True)
	
	
	# In[19]:
	
	
	# Print the new data frame shape
	if (frame.shape[1] == 0):
		print ('ERROR: Empty data frame')
		sys.exit(-2)
	# else:
	#	print(frame.shape)
	
	
	# ## Write to Excel file
	
	# In[20]:
	
	
	# Write to Excel format
	frame.to_excel(data_dir + '/Concat.xlsx')
	
	# Print the output address
	print ('Writing: ' + data_dir + '/Concat.xlsx')
