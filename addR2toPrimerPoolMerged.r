#!/usr/bin/env Rscript
##usage: addR2toPrimerPoolMerged.r primerPoolMerged r2
##example: addR2toPrimerPoolMerged.r AIMS_selection/multiplexPrimers/AFR.EUR/primerPoolMerged.txt AIMS_selection/multiplexPrimers/AFR.EUR/AFR.EUR.ancestry_r2_vs_num_AIMs.txt

#######
#setup
#######
library(data.table)
library(dplyr)


#####
#main
#####
##read command line args:
args = commandArgs(trailingOnly = T)
primerPoolMerged_fname = args[1]
r2_fname = args[2]

##read files and merge:: 
primerPoolMerged = fread(primerPoolMerged_fname)

if (ncol(fread(r2_fname)) == 4){
	##if a pair of population, e.g. AFR.EUR:
	r2 = fread(r2_fname) %>% select(poolRank = npool, r2 = pop1) 
	primerPoolR2 = cbind(primerPoolMerged, merge(primerPoolMerged, r2, by = "poolRank") %>% select(r2))
} else if (ncol(fread(r2_fname)) == 5) {
	##for 3 populations like AFR.AMR.EUR
	r2 = fread(r2_fname) %>% select(poolRank = npool, EUR_r2 = pop1, AMR_r2 = pop2, AFR_r2 = pop3)
	primerPoolR2 = cbind(primerPoolMerged, merge(primerPoolMerged, r2, by = "poolRank") %>% select(contains('r2')))
} else {
	message(sprintf('number of columns in %s is neither 4 or 5!',r2_fname))
}

##output:
write.table(primerPoolR2, file = primerPoolMerged_fname, row.names = F, col.names = T, quote = F, sep = '\t')
