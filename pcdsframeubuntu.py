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

		sample.write('# Workspace directory \n')
		sample.write('workspace=./workspace\n\n')

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
		sample.write('mod7=X 0 0\n')
		sample.write('mod8=X 0 0\n')
		sample.write('mod9=X 0 0\n')
		sample.write('mod10=X 0 0\n')
		sample.write('mod11=X 0 0\n')
		sample.write('mod12=X 0 0\n')
		sample.write('mod13=X 0 0\n')
		sample.write('mod14=X 0 0\n')
		sample.write('mod15=X 0 0\n\n')
		
		sample.write('# Path to MS/MS dataset\n')
		sample.write('ms2data=/path/to/msms/dataset\n\n')

		sample.write('# Missed cleavages\n')
		sample.write('missed_cleavages=1\n\n')
		
		sample.write('# Min peptide length\n')
		sample.write('min_length=6\n\n')
		
		sample.write('# Max peptide length\n')
		sample.write('max_length=40\n\n')

		sample.write('# Min precursor mass (Da)\n')
		sample.write('min_prec_mass=500\n\n')
		
		sample.write('# Max precursor mass (Da)\n')
		sample.write('max_prec_mass=5000\n\n')

		sample.write('# Digestion Enzyme\n')
		sample.write('enzyme=Trypsin\n\n')

		sample.write('# Index Distribution Policy: chunk, cyclic, zigzag\n')
		sample.write('policy=cyclic\n\n')

		sample.write('# Max fragment charge\n')
		sample.write('maxz=3\n\n')

		sample.write('# Min shared peak\n')
		sample.write('shp=4\n\n')
		
		sample.write('# Resolution (Da)\n')
		sample.write('res=0.01\n\n')
		
		sample.write('# Precursor Mass Tolerance (+-Da): -1 means infinity \n')
		sample.write('dM=500\n\n')

		sample.write('# Fragment Mass Tolerance (+-Da)\n')
		sample.write('dF=0.05\n\n')

		sample.write('# Top Matches to report\n')
		sample.write('top_matches=10\n')

		
		print('Generated: ./sampleparams.txt')
		print ("\nSUCCESS")
		
		sys.exit(0)
		
	# Max threads and cpus
	max_threads = subprocess.run("grep -c ^processor /proc/cpuinfo", 
									stdout=subprocess.PIPE, shell=True)
	max_threads = int(max_threads.stdout.decode('utf-8'))/2
	
	max_threads = int(max_threads)

	max_cpus = subprocess.run("grep ^cpu\\scores /proc/cpuinfo | uniq | awk '{print $4}'", 
								stdout=subprocess.PIPE, shell=True)
#	print(max_cpus.stdout.decode('utf-8'))
#	max_cpus = int(max_cpus.stdout.decode('utf-8'))

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
	dF = 0.02
	dM = 500
	res = 0.01
	scale = int(1/res)
	min_prec_mass = 500
	max_prec_mass = 5000
	top_matches = 10
	shp_cnt = 4
	workspace = './workspace'
	policy = 'cyclic'


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
	#				sys.exit(-2)
			
			elif (param == 'ms2data'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]

				ms2data = val
				print ('MS/MS dataset =', ms2data)
				if (os.path.exists(ms2data) == False):
					print ("ERROR: Enter valid path to MS2 dataset")
	#				sys.exit(-3)
				
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

			# Set the distribution policy
			elif (param == 'policy'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]
				policy = val
				print ('Using policy =', policy)

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

			# Set the max missed cleavages
			elif (param == 'missed_cleavages'):
				mcleavages = int(val)
				if (mcleavages < 0):
					mcleavages = 0 
				if (mcleavages > 5):
					mcleavages = 5
				print ('Missed Cleavages =', mcleavages)

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
				if (dF < 0.001):
					dF = 0.001 
				if (dF > 5.0):
					dF = 5.0
				dF = float(val)
				print ('dF           =', dF)

			# Peptide precursor mass tolerance
			elif (param == 'dM'):
				dM = float(val)
				if (dM < 0.001):
					dM = 0.001 
				print ('dM           =', dM)
				
			# m/z axis resolution
			elif (param == 'res'):
				res = float(val)
				if (min_prec_mass <= 0):
					min_prec_mass = 0.01 
				if (min_prec_mass > 5.0):
					min_prec_mass = 5.0
				print ('resolution   =', res)
	
			# Minimum precursor mass
			elif (param == 'min_prec_mass'):
				min_prec_mass = int(val)
				if (min_prec_mass < 0):
					min_prec_mass = 0 
				if (min_prec_mass > 10000):
					min_prec_mass = 10000
				print ('min_prec_mass =', min_prec_mass)

			# Maximum precursor mass
			elif (param == 'max_prec_mass'):
				max_prec_mass = int(val)
				if (max_prec_mass < 0):
					max_prec_mass = 0 
				if (max_prec_mass > 10000):
					max_prec_mass = 10000
				print ('max_prec_mass =', max_prec_mass)				

			# Minimum Shared Peaks
			elif (param == 'shp_cnt'):
				shp_cnt = int(val)
				if (shp_cnt < 1):
					shp_cnt = 1 
				if (shp_cnt > 20):
					shp_cnt = 20
				print ('Min Shared Peaks =', shp_cnt)

			# Workspace Path
			elif (param == 'workspace'):
				if (val[-1] == '\n'):
					val = val[:-1]
				if (val[-1] == '\r'):
					val = val[:-1]

				if (val[-1] == '/'):
					val = val[:-1]

				workspace = str(val)
				print ('workspace   =', workspace)
				
			# Maximum precursor mass
			elif (param == 'top_matches'):
				top_matches = int(val)
				print ('Top matches =', top_matches)

#	print ('Mods Added', mods)

	if (len(mods) == 0):
		mods.append("X 0 0")
		nmods = 0
	
	# Close the params file
	params.close()

	# Create a workspace directory
	print ('\nInitializing Workspace at: ', workspace)
	
	if (os.path.exists(workspace) == False):	
		os.mkdir(workspace)
	
	# Run the digestor now
	digesteddb = workspace + '/digested_db.fasta'
	digestcommand = "Digestor.exe -in " + database + " -out " + digesteddb + " -out_type fasta -threads " + str(threads) + " -missed_cleavages " + str(mcleavages) + " -enzyme " + enzyme +  " -min_length " + str(min_length) + " -max_length " + str(max_length) + " -FASTA:ID number -FASTA:description remove"
	
	print ('\nRunning: '+ digestcommand + '\n')

	# Run the Digester.exe
#	digestor = call(digestcommand, shell=True)
	
	print ("\nSUCCESS\n")
	
	
	# Print the next steps
#	print ('\nRunning: '+ 'Separate by Peptide Length')
#	print ('\nRunning: '+ 'Custom Lexicographical Sort\n')

	# Print the clustering command
	clustercommand = './bash/sep_by_len.sh ' + digesteddb + ' ' + str(min_length) + ' ' + str(max_length)

	print ('Running: ' + clustercommand)

	# Run the cluster command and pass arguments
#	cluster = subprocess.run(['./bash/sep_by_len.sh ', digesteddb, str(min_length), str(max_length)], stdout=subprocess.PIPE, shell=True)

	print ("\nSUCCESS\n")

	# Prepare the uparams.txt file for seq generator
	modfile = open(workspace + '/uparams.txt', "w+")

	modfile.write('/lclhome/mhase003/Data/Database/Digested/parts' + '\n')
	modfile.write(ms2data + '\n')
	modfile.write(str(max_threads) + '\n')
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
	modfile.write(str(policy) + '\n')

	modfile.write(str(len(mods)) + '\n')
	modfile.write(str(nmods) + '\n')
	for info in mods:
		aa,ms,num = info.split(' ', 2)
		modfile.write(aa + '\n')
		modfile.write(str(ms) + '\n')
		modfile.write(str(num) + '\n')

	modfile.close()

	# Generate sequences (All set for future use with MS2PIP)
#	makeseqs = call("make -C seqgen clean", shell=True)
#	makeseqs = call("make -C seqgen", shell=True)
#	genseqs = subprocess.run(['./seqgen/seqgen.exe ', './parts', './mods.txt', str(min_length), str(max_length), str(maxz)], stdout=subprocess.PIPE, shell=True)

	# Install deps for MS2PIP, set up the config.txt and run for all PEPREC files.
	# TODO:

	# Construct CFIR index and compute shared peak count
	uparams = './workspace/uparams.txt\n'
	cleancfir = call("make -C cfirindex clean", shell=True)
	makecfir = call("make -C cfirindex", shell=True)

	cfir = subprocess.run(['./cfirindex/cfir.exe ', uparams], stdout=subprocess.PIPE, shell=True)
	
#	print (cfir.stdout.decode('utf-8'))
