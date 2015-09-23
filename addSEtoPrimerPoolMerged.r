#!/usr/bin/env Rscript 
# description:
# add information about SE, bias and RMSE to primerPoolMerged files (as columns).
# usage: 
# Rscript addSEtoPrimerPoolMerged.r [primerPoolMerged] [Q_bias and Q_se directory] [whole-genome Q file] [population file]
# example: 
# Rscript addSEtoPrimerPoolMerged.r \
# 			AIMS_selection/multiplexPrimers/AFR.EUR/primerPoolMerged.txt \ 
# 			AIMS_selection/multiplexPrimers/AFR.EUR/admixture \
# 			AIMS_selection/multiplexPrimers/AFR.EUR/admixture/AFR.EUR.2.Q
#			AIMS_selection/alleleFreq/AFR.EUR.pop
######
#setup
######
library(data.table)
library(dplyr)
library(stringr)
library(lazyeval)

####
#fun
####
guess.number.of.population = function(directory){
	# description:
	# guess the number of populations based on
	# the Q file that exists in [directory]. 
	fileNames = list.files(directory)
	QfileName = fileNames[str_detect(fileNames, pattern = "([A-Z]{3}\\.){2,5}[0-9]\\.Q")]
	numPopulations = as.integer(str_extract(QfileName, pattern = "(?<=([A-Z]{3}\\.){2,5})([0-9])(?=\\.Q)"))
	cat(sprintf("%i populations (based on '%s')\n", numPopulations, QfileName))
	return(numPopulations)
}

order.Q2.by.Q1 <- function(Q1, Q2) {
	# description: 
	# return an ordering Q2 columns to match the ordering of Q1.
	# note:
	# sometimes the ordering is uncertain for some columns. for example,
	# the order can be c(1,1,3). when this happens, an warning is issued. 
	corr = cor(Q1, Q2)
	order = apply(corr, 1, which.max)
	if (sum(duplicated(order)) > 0) warning("ambiguous ordering: ", order)
	return(order)
}

# read.2.populations_ = function(path, pattern){
# 	# description: 
# 	# read files specified by [path] and [pattern] 
# 	# into a data.frame (in long format).
# 	# the file should be Nx2
# 	# use read.multi.populations_ to read files with > 2 columns
# 	fileNames = list.files(pattern = pattern, full.names = T, path = path)
# 	df = data.frame()
# 	for (fileName in fileNames){
# 		file = fread(fileName)
# 		numPool = as.integer(str_extract(fileName, pattern = "(?<=cumPrimerPool_)([0-9]{1,2})(?=\\.)"))
# 		df = rbind(df, data.frame(file[,1, with = F], numPool = numPool))
# 	}
# 	df = df %>% arrange(numPool)
# 	return(df)
# }

read_ = function(path, pattern, numPopulations, wgQfileName){
	# description: 
	# read files specified by [path] and [pattern] 
	# into a data.frame (of long format).
	# 
	# note: 
	# Q files associated with each [Q_bias|Q_se] must be present in the [path] directory too. 
	# whole-genome Q file must also be present. 

	fileNames = list.files(pattern = pattern, full.names = T, path = path)

	# read whole-genome Q file:
	wgQ = fread(wgQfileName)

	# read AIMs (Q_bias|Q_se) and associated Q files:
	df = data.frame()
	for (fileName in fileNames){

		# read (Q_bias|Q_se) file:	
		file = fread(fileName)
		numPool = as.integer(str_extract(fileName, pattern = "(?<=cumPrimerPool_)([0-9]{1,2})(?=\\.)"))
		
		# read Q file:
		QfileName = str_replace(fileName, pattern = str_match(string = pattern, pattern = "Q_se|Q_bias"), replacement = "Q")
		Q = fread(QfileName)
		
		# order columns according to wgQ:
		setcolorder(file, order.Q2.by.Q1(wgQ, Q))

		setnames(file, sprintf("V%i", 1:numPopulations))
		df = rbind(df, data.frame(file, numPool = numPool))
	}
	df = df %>% arrange(numPool)
	return(df)
}

read.Qbias = function(path, numPopulations, wgQfileName){
	# description:
	# read Q_bias files in [path] directory
	df = read_(path, pattern = "Q_bias$", numPopulations, wgQfileName)
	return(df)
}

read.Qse = function(path, numPopulations, wgQfileName){
	# description:
	# read Q_se files in [path] directory
	df = read_(path, pattern = "Q_se$", numPopulations, wgQfileName)
	return(df)
}

merge.Qbias.and.Qse = function(Qbias, Qse){
	# description:
	# merge Qbias and Qse and 
	# create new column rmse = sqrt(se^2 + bias^2) 
	add.key.column = function(x){
		# description: 
		# append row number as a key column
		x$key = seq(nrow(x))
		return(x)
	}

	Qbias = add.key.column(Qbias)
	Qse = add.key.column(Qse)

	Qmerge = merge(Qbias, Qse, by = c("key", "numPool"), suffixes = c(".bias", ".se")) %>% arrange(key)
	return(Qmerge)
}

appendRMSE = function(Qmerge, numPopulations){
	# description:
	# append RMSE (=sqrt(SE^2 + bias^2))
	functionCalls = sapply(X = 1:numPopulations, FUN = function(i) interp("sqrt(x^2 + y^2)", .values = list(x = as.name(sprintf("V%i.se",i)), y = as.name(sprintf("V%i.bias",i)))))
	names(functionCalls) = sprintf("V%i.rmse", 1:numPopulations)
	Qmerge = Qmerge %>% mutate_(.dots = functionCalls)
	return(Qmerge)
}

calculateRMSEmean = function(Qmerge, numPopulations){
	# description:
	# calculate the mean of RMSE of each pool.
	functionCalls = sapply(X = 1:numPopulations, FUN = function(i) interp("mean(x)", x = as.name(sprintf("V%i.rmse", i))))
	names(functionCalls) = sprintf("V%i.rmseMean", 1:numPopulations)
	Qsummary = Qmerge %>% group_by(numPool) %>% summarize_(.dots = functionCalls)
	return(Qsummary)
}

match.columns.and.populations = function(wgQfileName, popFileName){
	# description:
	# match the columns of whole-genome Q file to appropriate populations
	# according to the population file
	# value:
	# character vector of population names with the same order as the columns 
	# of whole-genome Q.  
	wgQ = fread(wgQfileName)
	
	popFile = fread(popFileName, header = F)
	setnames(popFile, c("familyID", "individualID", "population", "superPopulation", "gender"))
	
	populationIndicatorMatrix = model.matrix(~ 0 + superPopulation, popFile)
	colnames(populationIndicatorMatrix) = colnames(populationIndicatorMatrix) %>% str_replace(pattern = "superPopulation", replacement = "")

	match = colnames(populationIndicatorMatrix)[order.Q2.by.Q1(wgQ, populationIndicatorMatrix)]
	names(match) = colnames(wgQ)

	return(match)
}

rename.Qsummary = function(Qsummary, columnNames){
	# description:
	# replace V1, V2, etc with appropriate population names
	for (i in 1:length(columnNames)){
		names(Qsummary) = str_replace(names(Qsummary), names(columnNames)[i], columnNames[i])
	}
	return(Qsummary)
}

merge.primerPoolMerged.and.Qsummary = function(primerPoolMerged, Qsummary){
	# description:
	# merge primerPollMerged and Qsummary by column poolRank.
	Qsummary = Qsummary %>% rename(poolRank = numPool)
	primerPoolMerged = as.data.frame(primerPoolMerged)
	primerPoolMergedAndQsummary = merge(primerPoolMerged, Qsummary, by = "poolRank")
	return(primerPoolMergedAndQsummary)
}

#####
#main
#####
# read [primerPoolMerged] and [Q_bias and Q_se directory]
# args = commandArgs(trailingOnly = T)
# primerPoolMergedFileName = args[1]
# Qdirectory = args[2]
list(primerPoolMergedFileName, Qdirectory, wgQfileName, popFileName) %=% commandArgs(trailingOnly = T)

# read primerPoolMerged:
primerPoolMerged = fread(primerPoolMergedFileName)

# read Qbias and Qse: 
numPopulations = guess.number.of.population(Qdirectory)
Qbias = read.Qbias(Qdirectory, numPopulations, wgQfileName)
Qse = read.Qse(Qdirectory, numPopulations, wgQfileName)

# calculate the mean RMSE of each pool:
Qmerge = merge.Qbias.and.Qse(Qbias, Qse)
Qmerge = appendRMSE(Qmerge, numPopulations)
Qsummary = calculateRMSEmean(Qmerge, numPopulations)

# determine the column names for Qsummary
columnNames = match.columns.and.populations(wgQfileName, popFileName)
Qsummary = rename.Qsummary(Qsummary, columnNames)

# merge Qsummary and primerPoolMerged: 
primerPoolMergedAndQsummary = merge.primerPoolMerged.and.Qsummary(primerPoolMerged, Qsummary)

# write primerPoolMergedAndQsummary to txt
write.table(primerPoolMergedAndQsummary, file = primerPoolMergedFileName, row.names = F, col.names = T, quote = F, sep = '\t')