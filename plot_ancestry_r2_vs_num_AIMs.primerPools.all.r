#!/usr/bin/env Rscript 

######
#setup
######
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

#####
#main
#####

##read primerPoolMerged: 
primerPoolDir = "/srv/persistent/bliu2/ancestry/AIMS_selection/multiplexPrimers/primerPools"
primerPoolList = list.files(primerPoolDir, pattern = ".*primerPoolMerged.txt", full.names = TRUE)
r2List = list()
for (primerPool_fname in primerPoolList){
	##extract population name:
	populations = str_extract(primerPool_fname, "(?<=primerPools/)(.+)(?=.primerPoolMerged.txt)")
	r2List[[(length(r2List)+1)]] = fread(primerPool_fname) %>% select(contains('r2'))
	names(r2List)[length(r2List)] = populations
}

r2List