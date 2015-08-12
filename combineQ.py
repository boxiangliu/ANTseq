#!/usr/bin/env python 
import sys
import numpy as np 

##read command line args:
flist = sys.argv[1]

##every line in flist is a filename,
with open(flist, 'r') as f:

	##read filename: 
	for i, fname in enumerate(f):
		fname = fname.strip()
		split_fname = fname.split(".")
		try: 
			##read idx for individual:
			idx = int(split_fname[-3])
		except ValueError:
			print "file name in wrong format; individual index must be third to last field"

		##read Q file into numpy array:
		Q = np.loadtxt(fname) 
		nrow = Q.shape[0]
		ncol = Q.shape[1]

		##ADMIXTURE sometimes switch the order of the columns,
		##store the first file so that later files can be ordered according to it. 
		if i == 0:
			firstQ = Q
			output = Q

		##calculate correlation between columns of firstQ and current Q file:
		else: 
			concat = np.concatenate((firstQ,Q), axis = 1)
			corr = np.corrcoef(concat.T)[0:ncol, -ncol:]
			
			##reorder columns of files according to the first Q file if needed:
			if corr[0,0] < corr[0,1]:
				sys.stder.write("columns flipped in %s"%fname)
				output = Q[:,[1,0]]
			else:
				output = Q

		##output the line from input file corresponding to the index:
		sys.stdout.write("%f %f\n"%(output[idx-1,0],output[idx-1,1]))
