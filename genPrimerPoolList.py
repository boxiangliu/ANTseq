#!/usr/bin/env python 

import pandas as pd
import sys 
##########
#functions
##########
def unique(input):
	'''return unique elements from input without changing order'''
	output = []
	for x in input:
		if x not in output:
			output.append(x)
	return output

######
#main
######
primerPool = pd.read_table(sys.stdin)
primerPoolNum = unique(primerPool['pool'])
for i in primerPoolNum:
	sys.stdout.write("primerPool_%s.txt\n"%i)


