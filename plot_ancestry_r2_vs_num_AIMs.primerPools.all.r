#!/usr/bin/env Rscript 

######
#setup
######
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(reshape2))
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

##reshape r2List to long format: 
r2List_long = data.frame()
for (i in length(r2List)) {
	if (ncol(r2List[[i]]) == 1){
		r2List_long = rbind(r2List_long, data.frame(r2 = r2List[[i]], pop = names(r2List)[[1]])
	} else {
		##TODO:melt data.frame
	}
}