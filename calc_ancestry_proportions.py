#!/usr/bin/env python 
##example1: 
##calc_ancestry_proportions.py -p EAS,SAS -l primerPoolList

##example2: 
##calc_ancestry_proportions.py -p AMR,EUR,AFR -l /srv/persistent/bliu2/ancestry/AIMS_selection/multiplexPrimers/AFR.EUR/primerPoolList.txt

######
#setup
######
import os 
import sys 
import argparse 
import subprocess
import timeit
import datetime
######
#const
######
PANEL = "/srv/persistent/bliu2/shared/1000genomes/phase3v5/integrated_call_samples_v3.20130502.ALL.panel"
PLINK = "/srv/persistent/bliu2/tools/plink_1.90_beta3_linux_x86_64/plink"
TGP_VCF = "/srv/persistent/bliu2/ancestry/AIMS_selection/phase3v5_filtered/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.vcf.gz"
TGP_BED_PREFIX = "/srv/persistent/bliu2/ancestry/AIMS_selection/phase3v5_filtered/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05"
# PLINK.PRUNE.IN = "/srv/persistent/bliu2/ancestry/AIMS_selection/prune/plink.prune.in"
ADMIXTURE = "/srv/persistent/bliu2/tools/admixture_linux-1.23/admixture"

######
#funct
######
def gen_pop_file(pops, exe = True):
	'''function to generate pop file'''
	pop_file = ".".join(pops) + ".pop"
	if exe == True:
		with open(PANEL, 'r') as fin, open(pop_file, 'w') as fout:
			for line in fin:
				split_line = line.strip().split()
				if split_line[2] in pops:
					fout.write("%s\t%s\t%s\t%s\n"%(split_line[0], split_line[0], split_line[1], split_line[2]))
		print "generated " + pop_file
	else:
		pass 
	return pop_file

# def make_plink_bed_from_vcf(pops, markers_fname = None, output_prefix = None, log_fname = os.devnull):
# 	'''make plink bed, bim and fam files'''

# 	pop_fname = gen_pop_file(pops, exe = False)

# 	if output_prefix == None: output_prefix = ".".join(pops)

# 	if markers_fname != None:
# 		with open(log_fname, 'w') as out:
# 			subprocess.call([PLINK, "--vcf", TGP_VCF, "--double-id", "--snps-only", "--keep-allele-order", "--keep", pop_fname, "--extract", markers_fname, "--make-bed", "--out", output_prefix], stdout = out)
# 		print "making plink bed file, keeping " + " ".join(pops) + "and markers in " + markers_fname
# 	else:
# 		with open(log_fname, 'w') as out:
# 			subprocess.call([PLINK, "--vcf", TGP_VCF, "--double-id", "--snps-only", "--make-bed", "--keep-allele-order", "--keep", pop_fname, "--out", output_prefix], stdout = out)
# 		print "making plink bed file, keeping " + " ".join(pops) + "and all markers."
	
# 	return  "%s.bed"%output_prefix


def extract_pops_and_markers(pops, markers_fname = None, output_prefix = None, log_fname = os.devnull, bed_prefix = TGP_BED_PREFIX):
	'''make plink bed, bim and fam files, extracting variants specified by markers_fname
	   and populations specified by pops'''

	pop_fname = gen_pop_file(pops, exe = False)
	if output_prefix == None: output_prefix = ".".join(pops)

	if markers_fname != None:
		with open(log_fname, 'w') as out:
			subprocess.call([PLINK, "--bfile", bed_prefix, "--keep-allele-order", "--keep", pop_fname, "--extract", markers_fname, "--make-bed", "--out", output_prefix], stdout = out)
		print "making plink bed file, keeping " + " ".join(pops) + " and markers in " + markers_fname

	else:
		with open(log_fname, 'w') as out:
			subprocess.call([PLINK, "--bfile", bed_prefix, "--keep-allele-order", "--keep", pop_fname, "--make-bed", "--out", output_prefix], stdout = out)
		print "making plink bed file, keeping " + " ".join(pops) + " and all markers."

	return  "%s.bed"%output_prefix




def call_admixture(plink_bed_fname, numpops, log_fname = os.devnull):
	'''call admixture'''
	numpops = str(numpops)
	with open(log_fname, "w") as out: 
		subprocess.call([ADMIXTURE, plink_bed_fname, numpops], stdout = out)
	print "running admixture for %s..."%plink_bed_fname
	Q_fname = plink_bed_fname.replace('bed', '%s.Q'%numpops)
	return Q_fname

def calc_ancestry_proportions(pops, markers_fname = None):
	'''calculate ancestry proportions for the subset of 1000G individuals 
	   specified by pops, using markers in markers_file'''
	##generate pop file: 
	pop_fname = gen_pop_file(pops, exe = False)

	if not os.path.exists("./tmp"): os.mkdir("./tmp")

	if markers_fname != None:
		output_prefix = "%s.%s"%(".".join(pops), os.path.basename(markers_fname))
	else:
		output_prefix = ".".join(pops)


	##make plink bed file for the whole genome: 
	plink_bed_fname = extract_pops_and_markers(pops, markers_fname, output_prefix, './tmp/plink.%s.log'%output_prefix)

	##call admixture:
	numpops = len(pops)
	Q_fname = call_admixture(plink_bed_fname, numpops, "./tmp/admixture.%s.log"%output_prefix)
	return Q_fname

#####
#main
#####
##read command line args: 
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--populations", required = True, type = str, help = "populations, e.g. EAS,EUR,SAS")
parser.add_argument("-l", "--primerPoolList", required = True, type = str, help = "list of primer pools.")
args = parser.parse_args()

##start timer: 
tic = timeit.default_timer()

##make tmp/ to store junk:
if not os.path.exists("./tmp"): os.mkdir("./tmp")

##parse command line args: 
pops = [pop for pop in args.populations.strip().split(",")]
numpops = len(pops)
primerPoolList_fname = args.primerPoolList

##generate pop file: 
pop_fname = gen_pop_file(pops)

##calculate ancestry proportions using whole genome: 
calc_ancestry_proportions(pops, None)

##read list of primer pools:
with open(primerPoolList_fname, 'r') as f:
	primerPoolList = f.readlines()

##calculate ancestry proportions using cumulative primer pools:
poolNum = 0 
cumPrimerPool = []
for primerPool_fname in primerPoolList:
	poolNum += 1
	primerPool_fname = primerPool_fname.strip()

	##read primers:
	with open(primerPool_fname,'r') as f: 
		primerPool = [line.strip().split()[0] for line in f if not line.startswith('snp')]
		print "%i primers in pool %i"%(len(primerPool), poolNum)

	##add current primer pools into cumulative primer pool:
	cumPrimerPool += primerPool 

	##make markers file:
	cumPrimerPool_fname = "./tmp/cumPrimerPool_%i.snpid"%poolNum
	with open(cumPrimerPool_fname,'w') as f:
		f.write("\n".join(cumPrimerPool))

	##calculate ancestry proportions:
	Q_fname = calc_ancestry_proportions(pops, cumPrimerPool_fname)
	print "generated %s!"%Q_fname

##end timer:
toc = timeit.default_timer()
print "Time elapsed: %s"%str(datetime.timedelta(seconds = toc - tic))

