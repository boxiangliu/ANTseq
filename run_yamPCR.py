#!/usr/bin/env python 

##usage: 
##python run_yamPCR.py path/to/aimsFile primerPoolSize workdingDir

##example:
##python run_yamPCR.py "/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs/AFR.EUR/AFR_EUR_500k_500.aims" 12 "/srv/persistent/bliu2/ancestry/AIMS_selection/multiplexPrimers/test/"

######
#setup
######
from subprocess import call 
import sys 
import os 
import json
import timeit
import datetime

############
# constants
############
yamPCR = '/srv/persistent/bliu2/ancestry/multiplex_primer_design/mmPCR_gDNA_BL/yamPCR/yamPCR.pl'
debug = True
##########
#function
##########
def read_aims(fname):
	'''
	read *.aims file into a list of dictionaries. 
	each dictionary contains information about 1 AIMs marker. 
	'''
	aims_list = []
	with open(fname, 'r') as f:
		rank = 0 
		for line in f: 
			if line.strip().startswith('snp'): continue
			rank = rank + 1 
			split_line = line.strip().split()
			aims_dict = {}
			aims_dict['snp'] = split_line[0]
			aims_dict['chr'] = split_line[1]
			aims_dict['pos'] = split_line[2]
			aims_dict['In'] = split_line[5]
			aims_dict['rank'] = rank 
			aims_list.append(aims_dict)
	return aims_list

def hasEnoughAIMs(aims_list, n):
	'''check if aims_list contains more than n AIMs'''
	return len(aims_list) >= n

def retrieveAIMs(aims_list, n):
	'''retrieve n AIMs from the aims_list, simultaneously removing them from aims_list'''
	aims_subset = aims_list[0:n]
	aims_list = aims_list[n:]
	return aims_subset, aims_list

def removeAIMs(aims_list, aimsToRemove):
	'''remove aimsToRemove from aims_list'''
	if len(aimsToRemove) == 1:
		aims_list.remove(aimsToRemove[0])
	elif len(aimsToRemove) > 1:
		for aim in aimsToRemove:
			aims_list.remove(aim)
	else: 
		pass 
	return aims_list

def addAIMs(aims_list, aimsToAdd):
	'''add aimsToAdd to the end of aims_list
	   aimsToAdd should be a list of aims_dict objects'''
	aims_list = aims_list + aimsToAdd
	return aims_list

def gen_yamPCR_input(aims_list, fname):
	'''generate yamPCR targetFile'''
	with open(fname, 'w') as f: 
		for target in aims_list:
			f.write("%s\tchr%s\t%s\t%s\n"%(target["snp"], target["chr"], str(int(target["pos"]) - 1), target["pos"]))

def call_yamPCR(yamPCR, input_file, output_file, error_file, includedPrimers = None, numberCandidatePrimers = 10, maxPrimerTm = 64, minPrimerTm = 57, maxAmpliconSize = 400, minAmpliconSize = 150, execute = True):
	'''function to run yamPCR'''
	if execute == False: return 
	with open(output_file,'w') as out, open(error_file,'w') as err:
		if includedPrimers == None: 
			call([yamPCR,'--targetFile=%s'%input_file, '--numberCandidatePrimers=%i'%numberCandidatePrimers, '--maxPrimerTm=%i'%maxPrimerTm, '--minPrimerTm=%i'%minPrimerTm, '--maxAmpliconSize=%i'%maxAmpliconSize, '--minAmpliconSize=%i'%minAmpliconSize], stdout = out, stderr = err)
		else: 
			call([yamPCR,'--targetFile=%s'%input_file, '--numberCandidatePrimers=%i'%numberCandidatePrimers, '--maxPrimerTm=%i'%maxPrimerTm, '--minPrimerTm=%i'%minPrimerTm, '--maxAmpliconSize=%i'%maxAmpliconSize, '--minAmpliconSize=%i'%minAmpliconSize, '--includedPrimers=%s'%includedPrimers], stdout = out, stderr = err)

def extract_good_primers(yamPCR_output):
	'''function to read yamPCR output'''
	read = False
	primers_list = []
	with open(yamPCR_output, 'r') as f:
		for line in f: 

			if line.strip() == "Searching for a set of multiplex primers ...":
				read = True
			if line.strip() == "No compatible primer set found for all amplicons.": 
				read = False 
			if "candidate amplicons" in line.strip():
				read = True 
			if "Total length of region cover by these amplicons" in line.strip():
				read = False 

			if (read == False) or ("candidate amplicons" in line.strip()): continue
			
			primers = {}
			split_line = line.strip().split()

			try: 
				primers['snp'] = split_line[0]
				primers['forward'] = split_line[1]
				primers['forward_temp'] = split_line[2]
				primers['reverse'] = split_line[3]
				primers['reverse_temp'] = split_line[4]
				primers['range'] = split_line[5]
				primers['length'] = split_line[6]
				if (line.strip() != "Searching for a set of multiplex primers ...") and ("candidate amplicons" not in line.strip()): 
					primers_list.append(primers)
			except IndexError:
				pass
				# print "line %s has less number of fields than required"%line.strip()

	return primers_list

def remove_Primer3FailedSNPs(yamPCR_output, aims_list):
	'''remove AIMs from aims_list for which Primer3 failed to design primers for'''
	bad_snp = []
	with open(yamPCR_output, 'r') as f:
		for line in f:
			if "Primer3 couldn't find any primer" in line: 
				bad_snp.append(line.strip().split('for ')[1].split(',')[0])

	aims_list = [aims for aims in aims_list if aims['snp'] not in bad_snp]
	numPrimer3FailedSNPs = len(bad_snp)
	return aims_list, numPrimer3FailedSNPs, bad_snp

def noInteraction(yamPCR_output):
	'''determine whether there is interaction between candidate primers and included primers'''
	noInteraction = False 
	with open(yamPCR_output, 'r') as f:
		for line in f: 
			if line.strip() == "Number of candidate primers removed due to having interaction with included primers: 0.":
				noInteraction = True 
	return noInteraction

def remove_otherFailedSNPs(yamPCR_output, aims_list):
	'''remove SNPs for which there is no good primer
	   Primer3 is able to design primers for these SNPs,
	   but due to having more than one amplicon or unsuitable
	   amplicon size, no candidate primers can be found for these 
	   SNPs.'''
	bad_snp = []
	with open(yamPCR_output, 'r') as f:
		for line in f: 
			if "No good candidate primers could be found" in line:
				bad_snp.append(line.strip().split('for ')[1].split('! ')[0])
	aims_list = [aims for aims in aims_list if aims['snp'] not in bad_snp]
	numFailedSNPs = len(bad_snp)
	return aims_list, numFailedSNPs, bad_snp

def write_includedPrimers(primers_list, includedPrimersFile):
	'''function to write good primers into a includedPrimers file'''
	with open(includedPrimersFile, 'w') as f:
		for primers in primers_list:
			f.write(primers['forward'] + "\n")
			f.write(primers['reverse'] + "\n")

def removeAIMsWithPrimer(aims_list, primers_list):
	'''remove AIMs with good primers from the target list
	   return the remaining aims'''
	snp_list = [aims['snp'] for aims in aims_list]
	snp_withPrimer = [primers['snp'] for primers in primers_list]
	snp_withoutPrimer = list(set(snp_list) - set(snp_withPrimer))
	remaining_aims = [aims for aims in aims_list if aims['snp'] in snp_withoutPrimer]
	return remaining_aims

def allAimsHasPrimers(aims_list, primers_list):
	'''check if all AIMs has primers'''
	snp_withPrimer = set([primers['snp'] for primers in primers_list])
	snp_list = set([aims['snp'] for aims in aims_list])
	return snp_withPrimer == snp_list

def writePrimerPool(primers_list, primerPoolFile):
	'''write primer pool to file'''
	with open(primerPoolFile, 'w') as f:
		f.write("snp\tforward\tfwd_temp\treverse\trev_temp\trange\tlength\n")
		for primers in primers_list:
			f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(primers['snp'],primers['forward'],primers['forward_temp'],primers['reverse'],primers['reverse_temp'],primers['range'],primers['length']))

def repeatCounter(aims_list, previous_aims_list, timesRepeated):
	'''count the number of times the same set of aims has been repeated within the for loop'''
	if aims_list == previous_aims_list:
		timesRepeated += 1 
	else: 
		##if any change occured to aims_list, reset counter: 
		timesRepeated = 0 
	previous_aims_list = aims_list
	return timesRepeated, previous_aims_list

def reportRemainingAIMs(remaining_aims):
	'''report status on remaining aims''' 
	print "%i AIMs remaining:"%len(remaining_aims)
	print "\n".join([aims['snp'] for aims in remaining_aims])

def designMultiplexPrimers(aims_list, poolSize, poolNum):
	##check if more than 1-poolSize aims left in the list: 
	if len(aims_list) < poolSize: 
		aims_list = []
		primers_list = []
		return primers_list, aims_list

	print "## Iteration 0 ##"
	##retrieve aims from the master list:
	aims_currentPool, aims_list = retrieveAIMs(aims_list, poolSize)

	##generate yamPCR targetFile:
	if not os.path.exists('./tmp'):
		os.makedirs('./tmp')
	targetFile = "./tmp/targetFile.txt" # this will be a temporary file. 
	gen_yamPCR_input(aims_currentPool, targetFile)

	##call yamPCR:
	yamPCR_output = "./tmp/yamPCR_output_pool%i_iter0.txt"%poolNum
	yamPCR_error = "./tmp/yamPCR_output_pool%i_iter0.err"%poolNum
	call_yamPCR(yamPCR, targetFile, yamPCR_output, yamPCR_error)

	##extract compatible primers yamPCR output: 
	primers_list = []
	primers_list = extract_good_primers(yamPCR_output)

	##remove SNPs for which Primer3 cannot find any primer:
	aims_currentPool, numPrimer3FailedSNPs, primer3FailedSNPs = remove_Primer3FailedSNPs(yamPCR_output, aims_currentPool)
	print "Primer 3 cannot find primer for:"
	print "\n".join(primer3FailedSNPs)

	##replenish current aims list.
	##if there are not enough AIMs left, return primers_list and empty aims_list:  
	if hasEnoughAIMs(aims_list, numPrimer3FailedSNPs):
		retrieved_aims, aims_list = retrieveAIMs(aims_list, numPrimer3FailedSNPs)
		aims_currentPool += retrieved_aims
	else: 
		aims_list = []
		return primers_list, aims_list
	assert len(aims_currentPool) == poolSize, "current AIMs list has %i elements!"%len(aims_currentPool)

	##initialize: take out SNPs with good primers from target list, write remaining SNPs to target file: 
	remainingTargetFile = "./tmp/remainingTargetFile.txt"
	remaining_aims = removeAIMsWithPrimer(aims_currentPool, primers_list)
	gen_yamPCR_input(remaining_aims, remainingTargetFile)
	reportRemainingAIMs(remaining_aims)

	##loop until hitting poolSize: 
	n = 0
	timesRepeated = 0 
	previous_remaining_aims = []
	while len(primers_list) < poolSize: 
		n += 1
		print "## Iteration %i ##"%n

		##write includedPrimers
		includedPrimersFile = "./tmp/includedPrimers.txt"
		write_includedPrimers(primers_list, includedPrimersFile)

		##generate primers for remaining aims:
		numberCandidatePrimers = 10*n
		yamPCR_output = "./tmp/yamPCR_output_pool%i_iter%i.txt"%(poolNum,n)
		yamPCR_error = "./tmp/yamPCR_output_pool%i_iter%i.err"%(poolNum,n)
		call_yamPCR(yamPCR, input_file = remainingTargetFile, output_file = yamPCR_output, error_file = yamPCR_error, includedPrimers = includedPrimersFile, numberCandidatePrimers = numberCandidatePrimers)

		##append good primers to primers_list:
		primers_list += extract_good_primers(yamPCR_output)

		##remove SNPs for which Primer3 cannot find any primer:
		aims_currentPool, numPrimer3FailedSNPs, primer3FailedSNPs = remove_Primer3FailedSNPs(yamPCR_output, aims_currentPool)
		
		if numPrimer3FailedSNPs > 0: 
			print "Primer 3 cannot find primer for:"
			print "\n".join(primer3FailedSNPs)

			##replenish current aims list.
			##if there are not enough AIMs left, return primers_list and empty aims_list:  
			if hasEnoughAIMs(aims_list, numPrimer3FailedSNPs):
				retrieved_aims, aims_list = retrieveAIMs(aims_list, numPrimer3FailedSNPs)
				aims_currentPool += retrieved_aims
			else: 
				aims_list = []
				return primers_list, aims_list
			
		##check if remaining_aims primers interact with included primers,
		##if no, replace remaining_aims with new ones from master list:
		if noInteraction(yamPCR_output) == True: 
			aims_currentPool, numFailedSNPs, otherFailedSNPs = remove_otherFailedSNPs(yamPCR_output, aims_currentPool)
			
			##debug:
			if debug: print "within noInteraction conditionl, after removal, aims_currentPool is %s"%", ".join([aims['snp'] for aims in aims_currentPool])  
			
			print "Discarded SNPs:"
			print "\n".join(otherFailedSNPs)

			##add new SNPs from the master list: 
			if hasEnoughAIMs(aims_list, numFailedSNPs):
				retrieved_aims, aims_list = retrieveAIMs(aims_list, numFailedSNPs)
				aims_currentPool += retrieved_aims
			else: 
				aims_list = []
				return primers_list, aims_list

		
		#debug
		if debug: 
			print "previous remainning AIMs:"
			print "\n".join([aims['snp'] for aims in previous_remaining_aims])


		##count how many times the same aims_currentPool has appeared:
		timesRepeated, previous_remaining_aims = repeatCounter(remaining_aims, previous_remaining_aims, timesRepeated)
		
		##report:
		print "Repeated current pool for %i time(s)"%timesRepeated

		##if the same aims_currentPool has appeared more than 3 times
		##replace remaining_aims in aims_currentPool with new ones from master list,
		##This happens when remaining_aims candidate primer has interaction with includedPrimers. 
		if timesRepeated > 3: 
			##remove remaining_aims from aims_currentPool: 
			aims_currentPool = removeAIMs(aims_currentPool,remaining_aims)

			##debug:
			if debug: 
				print "after removeal, aims_currentPool is %s"%", ".join([aims['snp'] for aims in aims_currentPool])  
			
			##replace with new aims from master list: 
			if hasEnoughAIMs(aims_list, len(remaining_aims)):
				retrieved_aims, aims_list = retrieveAIMs(aims_list,len(remaining_aims))
				aims_currentPool += retrieved_aims

				##report: 
				print "Replaced %s with %s"%(", ".join([aims['snp'] for aims in remaining_aims]), ", ".join([aims['snp'] for aims in retrieved_aims]))

			else:
				aims_list = []
				return primers_list, aims_list

			##put remaining_aims back to master list.
			##these aims will be reused later: 
			aims_list = addAIMs(aims_list, remaining_aims)

			##reset times repeated:
			timesRepeated = 0

		##take out SNPs with good primers from target list, write remaining SNPs to target file: 
		remainingTargetFile = "./tmp/remainingTargetFile.txt"
		remaining_aims = removeAIMsWithPrimer(aims_currentPool, primers_list)
		gen_yamPCR_input(remaining_aims, remainingTargetFile)
		reportRemainingAIMs(remaining_aims)


		##sanity check whether aims_currentPool has 1-poolSize SNPs: 
		assert len(aims_currentPool) == poolSize, "current AIMs list has %i elements!"%len(aims_currentPool)

		##debug: 
		# with open('./debug_%i.txt'%n, 'w') as f:
		# 	json.dump(aims_currentPool, f)
		# 	json.dump(primers_list, f)

	return primers_list, aims_list

######
#main
######
##start timer: 
tic = timeit.default_timer()

##read command line arguments:
aims_fname = sys.argv[1]
poolSize = int(sys.argv[2])
wd = sys.argv[3]

##change to working directory: 
os.chdir(wd)

##read AIMs sites into a master list:
# aims_fname = '/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs/AFR.EUR/AFR_EUR_500k_500.aims'
aims_list = read_aims(aims_fname)

##design multiplex primers: 
poolNum = 0
while len(aims_list) > 0:
	ticc = timeit.default_timer()
	poolNum += 1
	print "Starting primer pool %i"%poolNum
	primers_list, aims_list = designMultiplexPrimers(aims_list, poolSize, poolNum)
	primerPoolFile = "./primerPool_%i.txt"%poolNum
	writePrimerPool(primers_list, primerPoolFile)
	tocc = timeit.default_timer()
	elapsed = str(datetime.timedelta(seconds = tocc - ticc))
	print "Primer pool %i done!"%poolNum
	print "Time elapsed: %s"%elapsed

##stop timer:
toc = timeit.default_timer()
total_elaped = str(datetime.timedelta(seconds = toc - tic))

##report done:
print "Designed %i primers pools"%poolNum
print "Total time elapsed: %s"%total_elaped