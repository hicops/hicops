#!@PYTHON_EXECUTABLE@
#
#  This file is a part of HPC PCDSFrame software
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
from shutil import which

# The main function
if __name__ == '__main__':

	# Read the arguments
	if len(sys.argv) > 1:
		paramfile = sys.argv[1]
	else:
		print ("ERROR: Missing arguments to the DBProcessor")
		print ("USAGE: python3.5+ dbprocessor.py <dbparams>")
		print ("GENERATE: A sample dbparams.txt using python3.5+ dbprocessor.py -g")
		sys.exit(-1)

	# Generate the sampleparams.txt file
	if (paramfile == '-g'):
		sample = open("./dbparams.txt","w+")

		sample.write('# \n')
		sample.write('# Protein Database Processor\n')
		sample.write('# Copyrights(C) 2019 PCDS Laboratory\n')
		sample.write('# Muhammad Haseeb, and Fahad Saeed')
		sample.write('# School of Computing and Information Sciences\n')
		sample.write('# Florida International University (FIU), Miami, FL\n')
		sample.write('# Email: {mhaseeb, fsaeed}@fiu.edu\n')
		sample.write('# \n')
		sample.write('# Auto generated sampledbparams.txt\n')
		sample.write('# Sample parameters generated for Protein Database Processor\n')
		sample.write('# \n')
		sample.write('# REQUIREMENTS: Digestor tool available with OpenMSv2.4.0 (www.openms.de)\n')
		sample.write('# \n')
		sample.write('# Generated on: ' + (datetime.datetime.now()).strftime("%Y-%m-%d %H:%M %Z") + '\n')
		sample.write('# \n')
		sample.write('# IMPORTANT: DO NOT put any spaces between variable=value\n')
		sample.write('# \n\n')

		sample.write('# ABSOLUTE path to protein_database.fasta\n')
		sample.write('database=/path/to/database.fasta\n\n')
		
		sample.write('# Path to the output files\n')
		sample.write('output=/path/to/output\n\n')

		sample.write('# Digestion Enzyme\n')
		sample.write('enzyme=Trypsin\n\n')

		sample.write('# Missed cleavages\n')
		sample.write('missed_cleavages=1\n\n')

		sample.write('# Min peptide length\n')
		sample.write('min_length=6\n\n')

		sample.write('# Max peptide length\n')
		sample.write('max_length=40\n\n')

		print('Generated: ./dbparams.txt')
		print ("\nSUCCESS")
		
		sys.exit(0)

# ##################################################################################

	# Initialize the parameters
	database = ''
	mcleavages = 2
	min_length = 6
	max_length = 40
	enzyme = 'Trypsin'
	output = ''
	cores = subprocess.run("grep ^cpu\\scores /proc/cpuinfo | uniq | awk '{print $4}'", 
								stdout=subprocess.PIPE, shell=True)
	cores = int(cores.stdout.decode('utf-8'))
# ##################################################################################

	print ('\n************************************\n')
	print   ('*    Protein Database Processor    *\n')
	print   ('*  Copyrights PCDS Lab, SCIS, FIU  *\n')
	print   ('************************************\n\n')

	# Parse the params file
	with open(paramfile) as params:
		for line in params:

			# Ignore the empty or comment lines
			if (line[0] == '\r' or line[0] == '#' or line[0] == '\n'):
				continue

			# Split line into param and value
			param, val = line.split('=', 1)

			# Set database file 
			if (param == 'database'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]

				database = os.path.abspath(val)
				if (os.path.isfile(database) == False):
					print ("ERROR: Enter valid path to Protein Database.FASTA")
					sys.exit(-2)
				else:
					print ('Database =', database)
					
			# Set output directory
			if (param == 'output'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]

				output = os.path.abspath(val)
				print ('Output =', output)
			# Set the enzyme for digestion
			elif (param == 'enzyme'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]
				enzyme = val
				print ('Using enzyme  =', enzyme)

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

			# Set the max missed cleavages
			elif (param == 'missed_cleavages'):
				mcleavages = int(val)
				if (mcleavages < 0):
					mcleavages = 0 
				if (mcleavages > 5):
					mcleavages = 5
				print ('Missed Cleavages =', mcleavages)

	# Close the params file
	params.close()

# ##################################################################################

	if ((which('Digestor') is not None) == False):
		print ('FATAL: The Digestor tool is missing\n')
		print ('ACTION: Please install OpenMS toolkit v2.4\n')
		print ('        and make sure that its binutils are\n')
		print ('        are included in the PATH\n')

		exit(-3)

	if (os.path.exists(output) == False):	
		os.mkdir(output)

	# Sanity check
	if (min_length > max_length):
		temp = min_length
		min_length = max_length
		max_length = temp
		print('WARNING: min_length > max_length. Swapping them\n')

# ##################################################################################

	# Run the digestor now
	digesteddb = output + '/digested_db.fasta'

	digestcommand = "Digestor.exe -in " + database + " -out " + digesteddb + " -out_type fasta -threads " + str(cores) + " -missed_cleavages " + str(mcleavages) + " -enzyme " + enzyme +  " -min_length " + str(min_length) + " -max_length " + str(max_length) + " -FASTA:ID number -FASTA:description remove"

	print ('\nRunning: OpenMS\'s Digestor.exe \n')

	# Run the Digester.exe
	digestor = call(digestcommand, shell=True)

	print ("\nSUCCESS\n")

	# Print the next steps
	print ('\nRunning: Peptide Sequence Clusterer')

	# Create the cluster command
	clustercommand = 'cluster.sh ' + digesteddb + ' ' + str(min_length) + ' ' + str(max_length) + ' no ' + output

	# Run the cluster command and pass arguments
	cluster = call(clustercommand, shell=True)

	print ("\nSUCCESS\n")

# ##################################################################################