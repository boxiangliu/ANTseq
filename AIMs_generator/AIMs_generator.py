'''
AIMs generator
Joshua Galanter and Chris Gignoux
December 16, 2010

Modified by Michelle Daya to support more than three ancestral populations
April 19, 2012

This program will take genotype files and allele frequency files from plink as well
as a pairwise LD file and will generate sets of AIMs based on informativeness.
'''

import sys
from sys import stdout
import os
from math import log
from numpy import *
from scipy import stats
import operator
from pdb import set_trace as t



def calc_sigma(af1, af2):
	return af1 + af2



def calc_delta(af1, af2):
	return abs(af1 - af2)



def calc_Fst(af1,af2):
	'''Calculates pairwise Fst given allele frequencies'''
	if (af1 == 0 and af2 == 0) or (af1 == 1 and af2 == 1):
		return 0
	return calc_delta(af1, af2) ** 2 / (calc_sigma(af1, af2) * (2 - calc_sigma(af1, af2)))



def calc_In(af1, af2):
	'''calculate Rosenberg et al.'s In statistic for pairs of populations'''
	if (af1 == 0 and af2 == 0) or (af1 == 1 and af2 == 1):
		return 0
	if af1 == 1:
		af1 = 0
		af2 = 1 - af2
	elif af2 == 1:
		af1 = 1 - af1
		af2 = 0
	sigma = calc_sigma(af1, af2)
	delta = calc_delta(af1, af2)
	a_exp = -0.5 * log(sigma ** sigma * (2 - sigma) ** (2 - sigma))
	b_exp = 0.25 * log( (sigma + delta) ** (sigma + delta) * (2 - sigma - delta) ** (2 - sigma - delta) * (sigma - delta) ** (sigma - delta) * (2 - sigma + delta) ** (2 - sigma + delta) )
	# return sum
	return a_exp + b_exp



def calc_lsbl(ab, ac, bc):
	'''Calculates locus specific branch length given three pairwise statistics 
	(either In or Fst)'''
	a = (ab + ac - bc)/2
	b = (ab + bc - ac)/2
	c = (ac + bc - ab)/2
	return (a, b, c)



def calculate_ld(ldfile, rsq_threshold):
	'''
	Read in LD and position calculated in PLINK for all snps in LD, given a 
	window sizepassed as a parameter.
	'''
	lddict = {}
	while True:
		try:
			# count the number of lines in ldfile
			with open(ldfile) as f: 
				for total, line in enumerate(f):
					pass
			with open(ldfile) as f: 
				for i, line in enumerate(f):
					if i > 0:
						print_progress('Importing LD file', 50, i, total)
						line = line.strip().split()
						if float(line[-1]) > rsq_threshold:
							snpA = line[2]
							snpB = line[5] # if PLINK r2--with-freqs is specified show change snpB = line[5] -> snpB = line[6]; Bosh 20150709 
							if snpA in lddict:
								lddict[snpA].append(snpB)
							else:
								lddict[snpA] = [snpB]
							if snpB in lddict:
								lddict[snpB].append(snpA)
							else:
								lddict[snpB] = [snpA]
			print
			return lddict
		except:
			sys.exit('File error; LD file is not in the correct format.')



def get_bim(filename):
	'''
	This function will read in a plink "bim" file to get the chromosomal positions.
	It will output a SNP - (chromosome, position) dictionary, as well as a 
	dictionary of alleles (minor, major).
	'''
	try:
		print 'Reading in position data from %s...'%filename
		with open(filename) as f: 
			posdict = dict([[line.strip().split()[1], \
				(line.strip().split()[0], line.strip().split()[3])] \
				for line in f])
		with open(filename) as f: 
			alleledict = dict([[line.strip().split()[1], \
				(line.strip().split()[4], line.strip().split()[5])] 
				for line in f])
		return (posdict, alleledict)
	except:
		sys.exit('File error; bim file is not in the correct format.')



def get_freq(population, nsnps, filename):
	'''
	This function will read in a plink "frq" file to get minor allele frequencies.
	'''
	try:
		print 'Reading allele frequency data from %s...'%filename
		freq = open(filename)
		freq.readline() # get rid of header line. 
		if nsnps == -1:
			return [line.strip().split() for line in freq.readlines()]
		else:
			snps = []
			for i, line in enumerate(freq.readlines()):
				if i > nsnps:
					sys.exit('Error, more SNPs in frequency file than in bim file.')
				print_progress('Importing frequency file', 50, i, nsnps)
				snp = line.strip().split()
				if len(snp) == 6:
					snps.append(snp)
				else:
					sys.exit('File error; incorrect allele frequency file format.')
			return snps	
		freq.close() 
	except:
		sys.exit('File error; incorrect allele frequency file format.')



def correct_freq(freq, mafdict):
	n = 0
	for i, snp in enumerate(freq):
		if snp[2] == '0':
			if snp[3] == mafdict[snp[1]][1]:
				freq[i][2] = mafdict[snp[1]][0]
			else:
				freq[i][4] = '1'
				freq[i][3] = mafdict[snp[1]][1]
				freq[i][2] = mafdict[snp[1]][0]
		elif snp[2] == mafdict[snp[1]][1]:
			n += 1
			freq[i][4] = str(1-float(snp[4]))
			freq[i][2] = mafdict[snp[1]][0]
			freq[i][3] = mafdict[snp[1]][1]
	print '\nFixed %s out of %s minor allele frequency flips' %(n, i)
	return freq


def sort_snps(snps, positions):
	print 'Sorting snps...'
	return sorted(snps, key = lambda order: (positions[order[1]][0], positions[order[1]][1]))


def calc_pairwise_aims(positions, pop1_frq, pop2_frq):
	'''
	Calculates the pairwise aims statistics for two populations given their allele 
	frequencies.
	Returns a list of aims, that include SNP, chromosome, position,
	allele frequency in population 1, allele frequency in population 2,
	sigma (the sum of allele frequencies), delta (the difference in allele 
	frequencies), pairwise Fst, and pairwise In
	'''
	aims = []
	nsnps = len(pop1_frq)
	for i, snp in enumerate(pop1_frq):
		if int(snp[0]) > 22:
			# ignore non-autosomal SNPs
			continue
		else:
			print_progress('Calculating pairwise AIMs statistics', 61, i, nsnps)
			af1 = float(pop1_frq[i][4])
			af2 = float(pop2_frq[i][4])
			aims.append([snp[1], positions[snp[1]][0], positions[snp[1]][1], af1, af2, calc_sigma(af1, af2), calc_delta(af1, af2), calc_Fst(af1,af2), calc_In(af1,af2)])
	return aims



def sort_pairwise_aims(unsorted_aims, stat = 'In'):
	'''
	Sorts the pairwise aims statistics given a set of aims.  It can sort on the basis
	of In or Fst.  Defaults to In.
	'''
	if stat == 'In':
		print '\nsorting output...'
		aims = sorted(unsorted_aims, key = lambda order: order[8])
		aims.reverse()
		return aims
	elif stat == 'Fst':
		print '\nsorting output...'
		aims = sorted(unsorted_aims, key = lambda order: order[7])
		aims.reverse()
		return aims
	else:
		print 'Invalid sort command, returning unsorted aims.'
		return aims



def output_pairwise_aims(aims, pop1, pop2, outfilename, n = -1):
	'''
	Writes the AIMs statistics to a file, given the filename.  n is an optional
	parameter that specifies how many AIMs to output.  It defaults to -1 (all)
	'''
	print 'Writing pairwise AIMs statistics for %s/%s populations to file %s' %(pop1, pop2, outfilename)
	outfile = open(outfilename, 'w')
	outfile.write('snp\tchr\tposition\t%s_allele_freq\t%s_allele_freq\tsigma\tdelta\tFst\tIn\n' %(pop1, pop2))
	lines = 0
	for line in aims:
		if lines == n:
			return
		else:
			outfile.write('%s\n' % ('\t'.join([str(val) for val in line])))
			lines += 1
	outfile.close()
	return



def calc_lsbl_pop_aims(positions, pop1_frq, pop2_frq, pop3_frq, AB_pairwise_aims, AC_pairwise_aims, BC_pairwise_aims, populations, outstem='LACE'):
	'''
	Calculates the locus specific brach length statistics for three populations, 
	given the pairwise AIMs statistics of three populations.  It returns a tuple 
	containing three sorted lists (sorted by best In statistic), one for each of the 
	three populations. The lists contain variables for:
	rsID, chromosome, position, allele frequency, locus specific branch length 
	(LSBL) by Fst, and LSBL by In
	'''
	pop1_aims = []
	pop2_aims = []
	pop3_aims = []
	n = 0
	nsnps = len(pop1_frq)
	for i, snp in enumerate(pop1_frq):
		if int(snp[0]) > 22:
			# ignore non-autosomal SNPs
			continue
		else:
			print_progress('Calculating branch length AIM statistics', 66, i, nsnps)
			threeway_Fst = calc_lsbl(AB_pairwise_aims[n][7], AC_pairwise_aims[n][7], BC_pairwise_aims[n][7])
			threeway_In = calc_lsbl(AB_pairwise_aims[n][8], AC_pairwise_aims[n][8], BC_pairwise_aims[n][8])
			pop1_aims.append([snp[1], positions[snp[1]][0], positions[snp[1]][1], float(pop1_frq[i][4]), threeway_Fst[0], threeway_In[0], populations[0]])
			pop2_aims.append([snp[1], positions[snp[1]][0], positions[snp[1]][1], float(pop2_frq[i][4]), threeway_Fst[1], threeway_In[1], populations[1]])
			pop3_aims.append([snp[1], positions[snp[1]][0], positions[snp[1]][1], float(pop3_frq[i][4]), threeway_Fst[2], threeway_In[2], populations[2]])
			n += 1
	print '\n'

	sorted_pop1_aims = sort_pop_aims(pop1_aims, 'In')
	outfile = outstem + '_' + pops[0] + '.aims'
	output_pop_aims(sorted_pop1_aims, pops[0], outfile, n)

	sorted_pop2_aims = sort_pop_aims(pop2_aims, 'In')
	outfile = outstem + '_' + pops[1] + '.aims'
	output_pop_aims(sorted_pop2_aims, pops[1], outfile, n)

	sorted_pop3_aims = sort_pop_aims(pop3_aims, 'In')
	outfile = outstem + '_' + pops[2] + '.aims'
	output_pop_aims(sorted_pop3_aims, pops[2], outfile, n)

	return (sorted_pop1_aims, sorted_pop2_aims, sorted_pop3_aims)



def sort_pop_aims(unsorted_aims, stat = 'In'):
	'''
	Sorts the locus specific branch length aims statistics given a set of aims.  It can sort on the basis
	of In or Fst.  Defaults to In.
	'''
	if stat == 'In':
		print 'sorting output...'
		aims = sorted(unsorted_aims, key = lambda order: order[5])
		aims.reverse()
		return aims
	elif stat == 'Fst':
		print 'sorting output...'
		aims = sorted(unsorted_aims, key = lambda order: order[4])
		aims.reverse()
		return aims
	else:
		print 'Invalid sort command, returning unsorted aims.'
		return aims



def output_pop_aims(aims, population, outfilename, n = -1):
	'''
	Writes the AIMs statistics to a file, given the filename.  n is an optional
	parameter that specifies how many AIMs to output.  It defaults to -1 (all)
	'''
	print 'Writing AIMs statistics for %s population to file %s' %(population, outfilename)
	outfile = open(outfilename, 'w')
	outfile.write('snp\tchr\tposition\tallele_frequency\tLSBL(Fst)\tLSBL(In)\n')
	lines = 0
	for line in aims:
		if lines == n:
			return
		else:
			outfile.write('%s\n' % ('\t'.join([str(val) for val in line[0:6]])))
			lines += 1
	outfile.close()
	return



def calc_paired_aims(positions, freq_pops, pops, outstem, n = -1):
	'''
	Calculates pairwise AIMs statistics, and returns a population - frequency 
	dictionary. The population key is represented by the combination of the 
	conifugred positions of the populations, e.g. '12', '13', '23' for three 
	populations.
	'''
	aims_paired_pops={}
	for i in range(0, len(pops)-1):
		pop1 = pops[i]
		for j in range(i+1, len(pops)):	
			pop2 = pops[j]
			aims = calc_pairwise_aims(positions, freq_pops[pop1], freq_pops[pop2])
			sorted_aims = sort_pairwise_aims(aims)
			aims_paired_pops[str(i+1) + str(j+1)] = sorted_aims
			outfile = outstem + '_' + pop1 + '_' + pop2 + '.aims'
			output_pairwise_aims(sorted_aims, pop1, pop2, outfile, n)

	return aims_paired_pops



def too_close(index_snp, snp_list, distance):
	'''
	Determines whether the index SNP is close to a SNP in snp_list, 
	that has already been selected as an AIM
	'''
	if snp_list == [None]:
		print 'returning none'
		return False
	else:
		for snp in snp_list:
			if snp == None:
				print snp
				print 'snp is empty'
				return False
			else:
				if snp[1] == index_snp[1]:
					if abs(int(snp[2]) - int(index_snp[2])) <= distance:
						return True
	return False



def get_lsbl_aims(positions, lddict, alleles, populations, pop1aims, pop2aims, pop3aims, excluded, distance, n):
	aimslist = []
	print 'Generating informativeness dictionary for %s' %(populations[0])
	pop1stat = dict([[aims[0],aims[5]] for aims in pop1aims])
	print 'Generating informativeness dictionary for %s' %(populations[1])
	pop2stat = dict([[aims[0],aims[5]] for aims in pop2aims])
	print 'Generating informativeness dictionary for %s' %(populations[2])
	pop3stat = dict([[aims[0],aims[5]] for aims in pop3aims])
	print 'Initializing statistics'
	pop1info = 0.
	pop2info = 0.
	pop3info = 0.
	pop1pos = 0
	pop2pos = 0
	pop3pos = 0
	pop1_numaims = 0
	pop2_numaims = 0
	pop3_numaims = 0
	ldex = set()
	pop1ex = set()
	pop2ex = set()
	pop3ex = set()
	n_aims= 0
	n_het_excluded = 0
	while n_aims < n:
		found = False
		if (pop1info < pop2info) and (pop1info < pop3info):
			aim = pop1aims[pop1pos]
			print(aim)
			pop1pos += 1
			if (aim[0] not in (ldex | excluded)) and not too_close(aim, aimslist, distance):
				print 'Found an AIM for %s, %s; %s AIMs found so far' %(populations[0], aim[0], len(aimslist) + 1)
				found = True
				aimslist.append(aim)
				pop1_numaims += 1
		elif pop2info < pop3info:
			aim = pop2aims[pop2pos]
			print(aim)
			pop2pos += 1
			if (aim[0] not in (ldex | excluded)) and not too_close(aim, aimslist, distance):
				print 'Found an AIM for %s, %s; %s AIMs found so far' %(populations[1], aim[0], len(aimslist) + 1)
				found = True
				aimslist.append(aim)
				pop2_numaims += 1
		else:
			aim = pop3aims[pop3pos]
			print(aim)
			pop3pos += 1
			if (aim[0] not in (ldex | excluded)) and not too_close(aim, aimslist, distance):
				print 'Found an AIM for %s, %s; %s AIMs found so far' %(populations[2], aim[0], len(aimslist) + 1)
				found = True
				aimslist.append(aim)
				pop3_numaims += 1
		if found:
			n_aims += 1
			try:
				ldex = ldex | set(lddict[aim[0]])
			except:
				print 'snp %s is not in the ld dictionary.' %(aim[0])
			pop1info += float(pop1stat[aim[0]])
			pop2info += float(pop2stat[aim[0]])			
			pop3info += float(pop3stat[aim[0]])
		else:
			if aim[0] in excluded:
				n_het_excluded += 1
	print 'The total locus specific In for the three populations are:'
	print 'For population %s, found %s aims out of %s evaluated, for a total LSBL of %s' %(populations[0], pop1_numaims, pop1pos, pop1info)
	print 'For population %s, found %s aims out of %s evaluated, for a total LSBL of %s' %(populations[1], pop2_numaims, pop2pos, pop2info)
	print 'For population %s, found %s aims out of %s evaluated, for a total LSBL of %s' %(populations[2], pop3_numaims, pop3pos, pop3info)
	print 'A total of %s AIMs were excluded due to heterogeneity' %(n_het_excluded)
	return aimslist



def print_lsbl_aims(aims, filename, pop1frq, pop2frq, pop3frq, populations):
	'''
	Prints output file for the total LSBL strategy for 3 populations
	'''
	pop1dict = dict([[snp[1],snp[4]] for snp in pop1frq])
	pop2dict = dict([[snp[1],snp[4]] for snp in pop2frq])
	pop3dict = dict([[snp[1],snp[4]] for snp in pop3frq])
	outfile = open(filename, 'w')
	outfile.write('snp\tchr\tposition\t%s_AF\t%s_AF\t%s_AF\tpopulation\tLSBL(Fst)\tLSBL(In)\n' %(populations[0], populations[1], populations[2]))		
	for i, aim in enumerate(aims):
		outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (aim[0], aim[1], aim[2], pop1dict[aim[0]], pop2dict[aim[0]], pop3dict[aim[0]], aim[6], aim[4], aim[5]))
	outfile.close()


def chi(table):
	'''
	Does a chi square test on the given table
	'''
	observed = array(table)
	rowsum = observed.sum(axis = 1)
	colsum = observed.sum(axis = 0)
	expected = (rowsum[:, newaxis] * colsum) / sum(rowsum)
	return stats.chisquare(expected.reshape(-1), observed.reshape(-1))



def calc_het(populations, filenames, pop_freq, mafdict, posdict, threshold = 0.01):
	print 'Calculating heterogeneity for the following populations:',
	for i in populations:
		print i,
	print
	allele_freqs = []
	heterogeneity = []
	for i, filename in enumerate(filenames):
		allele_freqs.append(sort_snps(correct_freq(get_freq(populations[i], len(posdict), filename), mafdict), posdict))
	for i, snp in enumerate(pop_freq):
		table = []
		for j, filename in enumerate(filenames):
			minor_allele_count = float(allele_freqs[j][i][4]) * float(allele_freqs[j][i][5])
			major_allele_count = (1. - float(allele_freqs[j][i][4])) * float(allele_freqs[j][i][5])
			#Counts of zero results in an expected cell count of
			#0 in the contingency tables. This results in the following
			#type of error when the chi() method is called:
			#	Warning: divide by zero encountered in divide
			#	Warning: invalid value encountered in chdtrc
			#To suppress this warning, a count of 0 is set to 1. If
			#an allele has a 0 count in one subpopulation, but has a large count
			#in another subpopulation, setting the count to 0.1 should still
			#result in a small P value, i.e. the SNP will be rejected because of
			#heterogeneity
			if minor_allele_count == 0:
				minor_allele_count = 0.1
			if major_allele_count == 0:
				major_allele_count = 0.1
			table.append([minor_allele_count, major_allele_count])	
		heterogeneity.append(chi(table))
	return heterogeneity



def exclude_het(hetfile, frq, threshold):
	print 'Finding SNPs to exclude on the basis of heterogeneity...'
	exclude = set()
	for i, snp in enumerate(frq):
		if hetfile[i][1] < threshold:
			exclude.add(frq[i][1])
	print 'Excluding %s SNPs on the basis of heterogeneity...' %(len(exclude))
	return exclude



def calc_log(x):
	'''Calculates the natural logarithm. log(0) is 0'''
	if x == 0:
		return 0
	else:
		return log(x)



def calc_multi_In(freq_pops):
	'''
	Calculate Rosenberg et al.'s In statistic for all SNPs,  across all 
	populations. Returns a SNP - In statistic list
	'''
	In_aims = []
	pop_frq = freq_pops[freq_pops.keys()[0]]

	for i,snp in enumerate(pop_frq):
		if int(snp[0]) > 22:
			# ignore non-autosomal SNPs
			continue
		else:
			print_progress('Calculating In statistic', 61, i, len(pop_frq))
		af1_list = []
		af2_list = []
		laf1_list = []
		laf2_list = []
		for pop in freq_pops:
			af1 = float(freq_pops[pop][i][4])
			af2 = 1 - af1
			af1_list.append(af1)
			af2_list.append(af2)
			laf1_list.append(af1 * calc_log(af1))
			laf2_list.append(af2 * calc_log(af2))
		p1 = average(af1_list)
		p2 = average(af2_list)
		lp1 = average(laf1_list)
		lp2 = average(laf2_list)
		In = lp1 + lp2 - p1*calc_log(p1) - p2*calc_log(p2)
		In_aims.append([snp[1], In])

	return In_aims



def get_properties(filename='aims_properties.txt'):
	'''Reads properties from a properties file'''
	properties={}
	try:
		propfile = open(filename, "r")
	except:
		sys.exit('Input error; properties file ' + filename + ' does not exist.')
	for line in propfile:
		data=line.strip().split('=')
		if len(data) == 2:
			properties[data[0].strip()] = data[1].strip()
		else:
			sys.exit('Input error; invalid line found in properties file: ' + line)
	return properties



def get_property(properties, prop, mandatory=False, default='', is_file=False):
	'''Returns the value of a property'''
	val = default
	if prop in properties:
		val = properties[prop]
	if mandatory and (val == ''):
		sys.exit('Input error; mandatory property ' + prop + ' is not specified.')
	if is_file:
		try:
			f = open(val, "r")
		except:
			sys.exit('Input error; file ' + val + ' specified for property ' + \
				prop + ' does not exist.')
	return val



def get_multi_In_aims(positions, lddict, In_aims, excluded, distance,n):
	'''
	Returns a final list of AIMs, according to the best In statistic across
	all the populations. A list consisting of SNP, chr, pos, In is 
	returned.
	'''

	if n == 0:
		return []

	print '\nFinding best multiple In markers'
	
	sorted_In_aims = sorted(In_aims, key = lambda order: order[1])
	sorted_In_aims.reverse()

	aimslist = []
	ldex = set()
	n_aims = 0
	n_het_excluded = 0
	n_ld_excluded = 0
	n_dist_excluded = 0
	for i, snp_in in enumerate(sorted_In_aims):
		snp = snp_in[0]
		aim = [snp, positions[snp][0], positions[snp][1], snp_in[1]]
		if n_aims >= n:
			break
		if (snp not in (ldex | excluded)) and not too_close(aim, aimslist, distance):
			print 'Found an AIM %s; %s AIMs found so far' %(snp, len(aimslist) + 1)
			aimslist.append(aim)
			n_aims += 1
			try:
				ldex = ldex | set(lddict[snp])
			except:
				print 'snp %s is not in the ld dictionary.' %(snp)
		else:
			if snp in excluded:
				n_het_excluded += 1
			elif snp in ldex:
				n_ld_excluded += 1
			else:
				n_dist_excluded += 1
	print 'Totals for best multiple In markers:'
	print 'A total of %s AIMs were excluded due to heterogeneity' %(n_het_excluded)
	print 'A total of %s AIMs were excluded due to LD' %(n_ld_excluded)
	print 'A total of %s AIMs were excluded due to distance' %(n_dist_excluded)

	return aimslist



def get_paired_In_aims(aimslist, positions, lddict, In_aims, aims_paired_pops, excluded, distance, n):
	'''
	Returns a final list of AIMs, by selecting the best markers that distinguish
	between a pair of populations, in such a way that the total In statistic
	is balanced across all the pairs. 
	'''

	if n == 0:
		return aimslist 

	print '\nFinding best pairwise In markers'

	paired_stat = []
	paired_info = []
	paired_pos = []
	paired_numaims = []
	for i, combo in enumerate(aims_paired_pops):
		print 'Generating informativeness dictionary for population pair %s' %(combo)
		paired_stat.append(dict([[aims[0],aims[8]] for aims in aims_paired_pops[combo]]))
		paired_pos.append(0)
		paired_numaims.append(0)
		if len(aimslist) == 0:
			paired_info.append(0)
		else:
			total = 0
			for aim in aimslist:
				total += float(paired_stat[i][aim[0]])
			paired_info.append(total)

	In_aims_dict = dict([[In_aim[0], In_aim[1]] for In_aim in In_aims])
	ldex = set()
	n_aims = 0
	n_het_excluded = 0
	n_ld_excluded = 0
	n_dist_excluded = 0
	while n_aims < n:
		#Find the position of the population pair with the minimum total 
		#pairwise statistic. If there is a tie, select the pair with 
		#the maximum pairwise statistic
		min_val = min(paired_info)
		min_pos_vec = []
		for i in range(0, len(paired_info)):
			if paired_info[i] == min_val:
				min_pos_vec.append(i)
		if (len(min_pos_vec) == 1):
			min_pos = min_pos_vec[0]
		else:
			max_val_pos = []
			for i in range(0, len(min_pos_vec)):
				rec = aims_paired_pops[aims_paired_pops.keys()[min_pos_vec[i]]][paired_pos[min_pos_vec[i]]]
				aim = [rec[0], rec[1], rec[2], In_aims_dict[rec[0]]]
				paired_in =  paired_stat[min_pos_vec[i]][aim[0]]
				max_val_pos.append(paired_in)
			max_min_pos = max_val_pos.index(max(max_val_pos))
			min_pos = min_pos_vec[max_min_pos]
		rec = aims_paired_pops[aims_paired_pops.keys()[min_pos]][paired_pos[min_pos]]
		aim = [rec[0], rec[1], rec[2], In_aims_dict[rec[0]]]
		paired_pos[min_pos] += 1
		if (aim[0] not in (ldex | excluded)) and not too_close(aim, aimslist, distance):
			print 'Found an AIM for pair %s, %s; %s AIMs found so far' \
				%(aims_paired_pops.keys()[min_pos], aim[0], len(aimslist) + 1)
			aimslist.append(aim)
			paired_numaims[min_pos] += 1
			n_aims += 1
			try:
				ldex = ldex | set(lddict[aim[0]])
			except:
				print 'snp %s is not in the ld dictionary.' %(aim[0])
			paired_info[min_pos] += paired_stat[min_pos][aim[0]]
		else:
			if aim[0] in excluded:
				n_het_excluded += 1
			elif aim[0] in ldex:
				n_ld_excluded += 1
			else:
				n_dist_excluded += 1
	print 'Totals for best pairwise In markers:'
	print 'A total of %s AIMs were excluded due to heterogeneity' %(n_het_excluded)
	print 'A total of %s AIMs were excluded due to LD' %(n_ld_excluded)
	print 'A total of %s AIMs were excluded due to distance' %(n_dist_excluded)
	print 'Total pairwise AIMs: %s' %(paired_info)


	return aimslist


def print_aims(aims, populations, freq_pops, aims_paired_pops, filename):
	'''
	Prints output file for the In and pairedIn strategies
	'''
	outfile = open(filename, 'w')
	outfile.write('snp\tchr\tposition')
	pop_stats = []
	for i in range(0, len(populations)):
		outfile.write('\t' + populations[i] + '_AF')
		pop_stats.append(dict([frq[1],frq[4]] for frq in freq_pops[populations[i]]))
	outfile.write('\tIn')
	paired_stats = []
	for combo in aims_paired_pops:
		pop1_index = int(combo[0])-1
		pop2_index = int(combo[1])-1
		pop1_name = populations[pop1_index]
		pop2_name = populations[pop2_index]
		outfile.write('\t' + pop1_name + '_' + pop2_name + '_In')
		paired_stats.append(dict([stat[0],stat[8]] for stat in aims_paired_pops[combo]))
	outfile.write('\n')

	for i, aim in enumerate(aims):
		outfile.write('%s\t%s\t%s' %(aim[0], aim[1], aim[2]))		
		for i in range(0, len(populations)):
			outfile.write('\t' + str(pop_stats[i][aim[0]]))
		outfile.write('\t' + str(aim[3]))
		for i, combo in enumerate(aims_paired_pops):
			outfile.write('\t' + str(paired_stats[i][aim[0]]))
		outfile.write('\n')
	outfile.close()


def print_progress(msg, pos, i, total):
	'''Prints a progress bar'''
	if i % (total/20) == 0:
		nr_eq = i * 20 / total
		nr_spaces = 19 - (i * 20 / total)
		print '%s%s:  [%s%s]' %('\b'*pos, msg, '=' * nr_eq, ' ' * nr_spaces), 
		stdout.flush()



if __name__ == '__main__':
	# bosh 20151103
	print 'started program.' 

	# end bosh
	if len(sys.argv) > 1:
		properties = get_properties(sys.argv[1])
	else:
		properties = get_properties()
	
	ldfile = get_property(properties, 'ldfile', True, '', True)
	try:
		ldthresh = float(get_property(properties, 'ldthresh', False, '0.1'))
		if (ldthresh < 0) or (ldthresh > 1):
			sys.exit('Input error; ldthresh is not between 0 and 1')
	except:
		sys.exit('Input error; ldthresh is not a valid number')
	lddict = calculate_ld(ldfile, ldthresh)
	
	posfile = get_property(properties, 'posfile', True, '', True)
	posdict, alleledict = get_bim(posfile)

	try:
		threshold = float(get_property(properties, 'hetthresh', False, '0.01'))
		if (threshold < 0) or (threshold > 1):
			sys.exit('Input error; hetthresh is not between 0 and 1')
	except:
		sys.exit('Input error; hetthresh is not a valid number')

	distances=[]
	dist_prop = get_property(properties, 'distances', False, '100000').strip().split(',')
	for dist in dist_prop:
		try:
			distances.append(int(dist))
		except:
			sys.exit('Input error; distances contains an invalid number')
	num=[]
	num_prop = get_property(properties, 'nraims', False, '500').strip().split(',')
	for nr in num_prop:
		try:
			num.append(int(nr))
		except:
			sys.exit('Input error; nraims contains an invalid number')

	populations=[]
	freq_pops={}
	exclude=set()
	pops = get_property(properties, 'populations', True).strip().split(',')

	strategy = get_property(properties, 'strategy', False, 'lsbl').strip()
	if strategy not in ('lsbl', 'In'):
		sys.exit('Input error; invalid strategy ' + strategy)
	if (strategy == 'lsbl') and len(pops) != 3:
		sys.exit('Input error; LSBL strategy can only be used for 3 populations')
	if strategy == 'In':
		try:
			proportion_multi_In = float( \
				get_property(properties, 'propmultiIn', False, '0.5').strip())
		except:
			sys.exit('Input error; propmultiIn not a valid number')
		if (proportion_multi_In < 0) or (proportion_multi_In > 1):
			sys.exit('Input error; propmultiIn not in the range [0,1]')
	
	outstem = get_property(properties, 'outstem', False, 'aims_').strip()

	for i, pop in enumerate(pops):
		populations.append(pop)
		freq_file = get_property(properties, pop + '.frq', False, pop + '.frq', True)
		freq = sort_snps(correct_freq(get_freq(populations[i], len(posdict), freq_file), alleledict), posdict)
		freq_pops[pop] = freq
		
		sub_pops = []
		sub_pops_freq_files = []
		if get_property(properties, pop + '.subpopulations') != '':
			sub_pops = get_property(properties, pop + '.subpopulations').strip().split(',')
		for sub_pop in sub_pops:
			sub_pop_freq_file = get_property(properties, sub_pop + '.frq', False, sub_pop + '.frq', True) 
			sub_pops_freq_files.append(sub_pop_freq_file)
		if len(sub_pops) > 0:
			het = calc_het(sub_pops, sub_pops_freq_files, freq, alleledict, posdict)
			exclude = exclude | exclude_het(het, freq, threshold)
	print('A total of %s AIMs have significant heterogeneity' %(len(exclude)))
	print(exclude)


	n_aims = -1
	for i, distance in enumerate(distances):
		for j, n in enumerate(num):
			if (distance < 1000000) and (distance > 1000) and ((distance % 1000)==0):
				dist = str(distance/1000) + 'k'
			elif (distance >= 1000000) and ((distance % 1000000)==0):
				dist = str(distance/1000000) + 'M'
			else:
				dist = str(distance)
			outfile = outstem + dist + '_' + str(n)
			aims_paired_pops = calc_paired_aims(posdict, freq_pops, populations, outfile, n_aims)
			if strategy == 'lsbl':
				pop_aims = calc_lsbl_pop_aims(posdict,
					freq_pops[populations[0]], freq_pops[populations[1]], \
					freq_pops[populations[2]], aims_paired_pops['12'], \
					aims_paired_pops['13'], aims_paired_pops['23'], populations, \
					outstem)
				aims = get_lsbl_aims(posdict, lddict, alleledict, populations, \
					pop_aims[0], pop_aims[1], pop_aims[2], exclude, distance, n)
				print_lsbl_aims(aims, outfile + '.aims', \
					freq_pops[populations[0]], freq_pops[populations[1]], \
					freq_pops[populations[2]], populations)
			elif strategy == 'In':
				nr_multi_aims = round(proportion_multi_In*n)
				nr_paired_aims = n - nr_multi_aims
				In_aims = calc_multi_In(freq_pops)
				aims = get_multi_In_aims(posdict, lddict, In_aims, exclude, \
						distance, nr_multi_aims) 
				aims = get_paired_In_aims(aims, posdict, lddict, In_aims, aims_paired_pops, \
						exclude, distance, nr_paired_aims)
				print_aims(aims, populations, freq_pops, aims_paired_pops, outfile + '.aims')

