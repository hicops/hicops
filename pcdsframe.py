# Required Imports

import os
import sys
import glob
import os.path
import datetime
import numpy as np
import scipy as scp
import pandas as pd
import subprocess
from subprocess import call

# The main function
if __name__ == '__main__':

	# Read the argument
	if len(sys.argv) > 1:
		paramfile = sys.argv[1]
	else:
		print ("ERROR: Enter the path to params file")
		print ("USAGE: py preprocess.py ./params.txt")
		print ("Generate a sampleparams.txt: py preprocess.py -g")
		sys.exit(-1)

	# Generate the sampleparams.txt file
	if (paramfile == '-g'):
		sample = open("./sampleparams.txt","w+")

		sample.write('# \n')
		sample.write('# DDA MS/MS Proteomics Pipeline\n')
		sample.write('# Copyrights 2019 PCDS Lab\n') 
		sample.write('# School of Computing and Information Sciences\n')
		sample.write('# Florida International University, Miami, FL\n')	
		sample.write('# \n')
		sample.write('# Auto generated sampleparams.txt\n')
		sample.write('# Generated on: ' + (datetime.datetime.now()).strftime("%Y-%m-%d %H:%M %Z")+'\n')
		sample.write('# \n')
		sample.write('# NOTE: Please do not put any spaces between variable=value\n')
		sample.write('# \n\n')
		
		sample.write('# Maxmimum threads: 0 for max\n')
		sample.write('threads=0\n\n')
		
		sample.write('# Path to proteome database\n')
		sample.write('database=/path/to/database.fasta\n\n')
		
		sample.write('# Mods to include per peptide sequence\n')
		sample.write('nmods=3\n\n')

		sample.write('# Mods Information: AAs mass mods_per_pep\n')
		sample.write('mod1=M 15.99 2\n')
		sample.write('mod2=X 0 0\n')
		sample.write('mod3=X 0 0\n')
		sample.write('mod4=X 0 0\n')
		sample.write('mod5=X 0 0\n')
		sample.write('mod6=X 0 0\n')
		sample.write('mod7=X 0 0\n\n')
		
		sample.write('# Path to MS/MS dataset\n')
		sample.write('ms2data=/path/to/msms/dataset\n\n')

		sample.write('# Missed cleavages\n')
		sample.write('missed_cleavages=1\n\n')
		
		sample.write('# Min peptide length\n')
		sample.write('min_length=6\n\n')
		
		sample.write('# Max peptide length\n')
		sample.write('max_length=40\n\n')

		sample.write('# Max fragment charge\n')
		sample.write('maxz=3\n\n')

		sample.write('# Digestion Enzyme\n')
		sample.write('enzyme=Trypsin\n')
		
		print('Generated: ./sampleparams.txt')
		print ("\nSUCCESS")
		
		sys.exit(0)
		
	# Max threads and cpus
	max_threads = subprocess.run("grep -c ^processor /proc/cpuinfo", 
									stdout=subprocess.PIPE, shell=True)
	max_threads = int(max_threads.stdout.decode('utf-8'))

	max_cpus = subprocess.run("grep ^cpu\\scores /proc/cpuinfo | uniq | awk '{print $4}'", 
								stdout=subprocess.PIPE, shell=True)
	max_cpus = int(max_cpus.stdout.decode('utf-8'))

	# Print stuff
#	print(max_cpus)
#	print(max_threads)

	# Initialize the parameters
	threads = -1
	database = ''
	ms2data = ''
	nmods = 0
	madded = 0
	mods = []
	mcleavages = 2
	min_length = 6
	max_length = 40
	maxz       = 3
	enzyme = 'Trypsin'	
	
	print ('\n************************************\n')
	print   ('*  DDA MS/MS Proteomics Pipeline   *\n')
	print   ('*  Copyrights PCDS Lab, SCIS, FIU  *\n')
	print   ('************************************\n\n')
	# Parse the params file
	with open(paramfile) as params:

		for line in params:

			# Ignore the empty or comment lines
			if (line[0] == 'r' or line[0] == '#' or line == '\n'):
				continue

			# Split line into param and value
			param, val = line.split('=', 1)

			# Set database file 
			if (param == 'database'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]

				database = val
				print ('Proteome DB   =', database)
				if (os.path.isfile(database) == False):
					print ("ERROR: Enter valid path to database.fasta")
					sys.exit(-2)
			
			elif (param == 'ms2data'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]

				ms2data = val
				print ('MS/MS dataset =', ms2data)
				if (os.path.exists(ms2data) == False):
					print ("ERROR: Enter valid path to MS2 dataset")
					#sys.exit(-3)
				
			# Set max threads to use
			elif (param == 'threads'):
				threads = int(val)
				if (threads <= 0 or threads > max_threads):
					threads = max_threads
				print ('Using threads = ' + str(threads))
				
			# Set the enzyme for digestion
			elif (param == 'enzyme'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]
				enzyme = val
				print ('Using enzyme  =', enzyme)
				
			# Set max mods
			elif (param == 'nmods'):
				nmods = int (val)
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

			# Set the max digestion length
			elif (param == 'max_length'):
				max_length = int(val)
				print ('Max pep len  =', max_length)

			# Set the max missed cleavages
			elif (param == 'missed_cleavages'):
				mcleavages = int(val)
				print ('Min pep len  =', min_length)

			# Set the max fragment ion charge
			elif (param == 'maxz'):
				maxz = int(val)
				print ('Max frag chg =', maxz)
			
	print ('Mods Added', mods)
	
	# Close the params file
	params.close()
	
	# Run the digestor now
	digesteddb = database[:-6] + '_digested.fasta'
	digestcommand = "Digestor.exe -in " + database + " -out " + digesteddb + " -out_type fasta -threads " + str(threads) + " -missed_cleavages " + str(mcleavages) + " -enzyme " + enzyme +  " -min_length " + str(min_length) + " -max_length " + str(max_length) + " -FASTA:ID number -FASTA:description remove"
	
	print ('\nRunning: '+ digestcommand + '\n')

	# Run the Digester.exe
	digestor = call(digestcommand, shell=True)
	
	print ("\nSUCCESS\n")
	
	
	# Print the next steps
#	print ('\nRunning: '+ 'Separate by Peptide Length')
#	print ('\nRunning: '+ 'Custom Lexicographical Sort\n')

	# Print the clustering command
	clustercommand = './bash/sep_by_len.sh ' + digesteddb + ' ' + str(min_length) + ' ' + str(max_length)

	print ('Running: ' + clustercommand)

	# Run the cluster command and pass arguments
	cluster = subprocess.run(['./bash/sep_by_len.sh ', digesteddb, str(min_length), str(max_length)], stdout=subprocess.PIPE, shell=True)

	print ("\nSUCCESS\n")

	# Prepare the mods.txt file for seq generator
	modfile = open('./mods.txt', "w+")
	modfile.write(str(len(mods)) + '\n')
	modfile.write(str(nmods) + '\n')

	for info in mods:
		aa,ms,num = info.split(' ', 2)
		modfile.write(aa + '\n')
		modfile.write(str(ms) + '\n')
		modfile.write(str(num) + '\n')

	modfile.close()

	# Generate sequences


	# Construct CFIR index and compute shared peak count