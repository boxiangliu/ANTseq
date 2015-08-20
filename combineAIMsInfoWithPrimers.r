#!/usr/bin/env Rscript

######
#setup
######
library(dplyr)
library(data.table)
library(gtools)

#####
#main
#####

##read primer pools: 
primer_dir = "/srv/persistent/bliu2/ancestry/AIMS_selection/multiplexPrimers/AFR.EUR"
primerPool_fnames = list.files(path = primer_dir, pattern = 'primerPool_[0-9]*.adaptor.txt') %>% mixedsort()

primerPools = data.table()
n = 0
for (primerPool_fname in primerPool_fnames){ 
	n = n + 1
	primerPool = data.frame(fread(primerPool_fname), pool = n)
	primerPools = rbind(primerPools, primerPool)
}

##read AIMs: 
AIMs_dir = "/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs_noSubpop/AFR.EUR"
AIMs = fread(list.files(path = AIMs_dir, pattern = ".*_500k_500.aims", full.names = TRUE))

##merge AIMs and primer pools:
merged = merge(primerPools, AIMs, by = 'snp') %>% 
	select(-chr, -position) %>% rename(temp = In) %>% select(-contains('In')) %>% rename(In = temp) %>% 
	rename(fwd_adaptor = fwd.adaptor, rev_adaptor = rev.adaptor) %>% 
	arrange(pool) 

##calculate cumulative Rosenberg's Informativeness:
merged %>% group_by(pool) %>% mutate(cum_In = sum(In)) %>% mutate(rank = order(cum_In, decreasing = T)) %>% select(snp, pool, In, cum_In, rank) %>% data.frame()