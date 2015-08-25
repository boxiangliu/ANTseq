#!/usr/bin/env Rscript
##example Rscript combineAIMsInfoWithPrimers.r ancestry/AIMS_selection/multiplexPrimers/AFR.EUR ancestry/AIMS_selection/AIMs_noSubpop/AFR.EUR
######
#setup
######
library(dplyr)
library(data.table)
library(gtools)

#####
#main
#####
##read command line args:
args = commandArgs(trailingOnly = TRUE)
primer_dir = args[1]
AIMs_fname = args[2]

##read primer pools: 
primerPool_fnames = list.files(path = primer_dir, pattern = 'primerPool_[0-9]*.adaptor.txt', full.names = T) %>% mixedsort()

primerPools = data.table()
n = 0
for (primerPool_fname in primerPool_fnames){ 
	if (nrow(fread(primerPool_fname)) == 0) {next}
	n = n + 1
	primerPool = data.frame(fread(primerPool_fname), pool = n)
	primerPools = rbind(primerPools, primerPool)
}

##read AIMs: 
AIMs = fread(AIMs_fname)

##merge AIMs and primer pools:
merged = merge(primerPools, AIMs, by = 'snp') %>% 
	select(-chr, -position) %>%  
	rename(fwd_adaptor = fwd.adaptor, rev_adaptor = rev.adaptor) %>% 
	arrange(pool) 

##rank pools by cumulative Rosenberg's Informativeness:
if (grepl("AFR_AMR_EUR_446.Galanter.orderedByIn.chrPos.aims", AIMs_fname)) {
	merged = merged %>% group_by(pool) %>% mutate(cum_In = sum(LSBL.In)) %>% arrange(desc(cum_In), desc(LSBL.In)) %>% data.table()
} else {
	merged = merged %>% group_by(pool) %>% mutate(cum_In = sum(In)) %>% arrange(desc(cum_In), desc(In)) %>% data.table()
}


poolRank = 0
previousPool = 0
merged$poolRank = 0
for (i in 1:nrow(merged)){
	if (merged$pool[i] != previousPool) {poolRank = poolRank + 1}
	merged$poolRank[i] = poolRank 
	previousPool = merged$pool[i]
} 

##write merged pools: 
write.table(merged, sprintf('%s/primerPoolMerged.txt',primer_dir), quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
message(sprintf("finished writing '%s/primerPoolMerged.txt",primer_dir))
