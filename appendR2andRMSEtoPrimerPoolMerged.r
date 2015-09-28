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
guess.number.of.population = function(QfileName){
	# description:
	# guess the number of populations based on
	# the Q file.
	numPopulations = as.integer(str_extract(QfileName, pattern = "(?<=([A-Z]{3}\\.){2,5}(cumPrimerPool_[0-9]{1,2}\\.){0,1})([0-9])(?=\\.Q)"))
	message(sprintf("%i populations (based on '%s')", numPopulations, QfileName))
	return(numPopulations)
}


guessNumerOfPools = function(fileName){
	# description:
	# guess the number of primer pools based on fileName
	# example: 
	# 27 primer pools based on the filename AFR.EUR.cumPrimerPool_27.2.Q 
	numPools = as.integer(str_extract(fileName, pattern = "(?<=cumPrimerPool_)([0-9]{1,2})(?=\\.)"))
	return(numPools)
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
		numPools = guessNumerOfPools(fileName)

		# read Q file:
		QfileName = str_replace(fileName, pattern = str_match(string = pattern, pattern = "Q_se|Q_bias"), replacement = "Q")
		Q = fread(QfileName)
		
		# order columns according to wgQ:
		setcolorder(file, order.Q2.by.Q1(wgQ, Q))

		setnames(file, sprintf("V%i", 1:numPopulations))
		df = rbind(df, data.frame(file, numPools = numPools))
	}
	df = df %>% arrange(numPools)
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

calculateAncestryCorrelation = function(x, ...){
	# description: 
	# generic method to calculate ancestry correlation
	UseMethod("calculateAncestryCorrelation")
}


calculateAncestryCorrelation.character = function(x, pattern, Q){
	# description: 
	# calculate the Pearson correlation between file Q and 
	# all files specified by [directory, pattern]
	fileNames = list.files(path = x, pattern = pattern, full.names = T)

	r2 = data.frame()

	for (fileName in fileNames){
		# TODO: debug this for loop.

		file = fread(fileName)
		columnOrder = order.Q2.by.Q1(Q, file)
		setcolorder(file, columnOrder)

		numPools = guessNumerOfPools(fileName)
		ancestryCorrelation = calculateAncestryCorrelation(file, Q)

		numPopulations = suppressMessages(guess.number.of.population(fileName))
		colnames(ancestryCorrelation) = sprintf("V%i.r2", 1:numPopulations)

		r2 = rbind(r2, data.frame(numPools = numPools, ancestryCorrelation))
	}

	r2 = r2 %>% arrange(numPools)
	return(r2)
} 

calculateAncestryCorrelation.data.frame = function(x, Q){
	# description: 
	# calculate the Pearson correlation between the columns of x and Q.
	# note: 
	# 1. x and Q should have the same number of columns
	# 2. x and Q should have the same column orders
	# 3. x and Q can be swapped without changing the result (obviously!).
	r2 = t(diag(cor(x,Q)))
	return(r2)
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

	Qmerge = merge(Qbias, Qse, by = c("key", "numPools"), suffixes = c(".bias", ".se")) %>% arrange(key)
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
	Qsummary = Qmerge %>% group_by(numPools) %>% summarize_(.dots = functionCalls)
	return(Qsummary)
}

calculateRMSEmax = function(Qmerge, numPopulations){
	# description:
	# calculate the max of RMSE of each pool.
	functionCalls = sapply(X = 1:numPopulations, FUN = function(i) interp("max(x)", x = as.name(sprintf("V%i.rmse", i))))
	names(functionCalls) = sprintf("V%i.rmseMax", 1:numPopulations)
	Qsummary = Qmerge %>% group_by(numPools) %>% summarize_(.dots = functionCalls)
	return(Qsummary)
}

matchColumnsAndPopulations = function(wgQfileName, popFileName){
	# description:
	# match the columns of whole-genome Q file to appropriate populations
	# according to the population file
	# value:
	# character vector of population names with the same order as the columns 
	# of whole-genome Q.  
	wgQ = fread(wgQfileName)
	
	popFile = fread(popFileName, header = F)
	setnames(popFile, old = c(1,2,3,4), new = c("familyID", "individualID", "population", "superPopulation"))
	
	populationIndicatorMatrix = model.matrix(~ 0 + superPopulation, popFile)
	colnames(populationIndicatorMatrix) = colnames(populationIndicatorMatrix) %>% str_replace(pattern = "superPopulation", replacement = "")

	match = colnames(populationIndicatorMatrix)[order.Q2.by.Q1(wgQ, populationIndicatorMatrix)]
	names(match) = colnames(wgQ)

	return(match)
}

replaceColumnNames = function(Qsummary, columnNames){
	# description:
	# replace V1, V2, etc with appropriate population names
	for (i in 1:length(columnNames)){
		names(Qsummary) = str_replace(names(Qsummary), names(columnNames)[i], columnNames[i])
	}
	return(Qsummary)
}

appendToPrimerPoolMerged = function(primerPoolMerged, addition){
	# description:
	# append addition to primerPollMerged by column poolRank.
	addition = addition %>% rename(poolRank = numPools)
	primerPoolMerged = as.data.frame(primerPoolMerged)
	primerPoolMerged = merge(primerPoolMerged, addition, by = "poolRank")
	return(primerPoolMerged)
}

#####
#main
#####
# read [primerPoolMerged] and [Q_bias and Q_se directory]
list('primerPoolMergedFileName', 'Qdirectory', 'wgQfileName', 'popFileName') %=% commandArgs(trailingOnly = T)

primerPoolMergedFileName="primerPoolMerged.txt"
Qdirectory="admixture/"
wgQfileName="admixture/EAS.SAS.2.Q"
popFileName="plink/EAS.SAS.pop"

# read primerPoolMerged:
primerPoolMerged = fread(primerPoolMergedFileName)

# read whole-genome Q file: 
wgQ = fread(wgQfileName)

# calculate r^2 between AIMs and genomic ancestry: 
r2 = calculateAncestryCorrelation(Qdirectory, pattern = "cumPrimerPool_[0-9]{1,2}\\.[2-5]\\.Q$", Q = wgQ)

# label columns with appropriate names:
columnNames = matchColumnsAndPopulations(wgQfileName, popFileName)
r2 = replaceColumnNames(r2, columnNames) 

# read Qbias and Qse: 
numPopulations = guess.number.of.population(wgQfileName)
Qbias = read.Qbias(Qdirectory, numPopulations, wgQfileName)
Qse = read.Qse(Qdirectory, numPopulations, wgQfileName)

# calculate the root mean squared error (RMSE) of ancestry coefficient by pool:
Qmerge = merge.Qbias.and.Qse(Qbias, Qse)
Qmerge = appendRMSE(Qmerge, numPopulations)

# calculate the mean and max of RMSE grouped by number of AIMs: 
RMSEmean = calculateRMSEmean(Qmerge, numPopulations) %>% replaceColumnNames(columnNames)
RMSEmax = calculateRMSEmax(Qmerge, numPopulations) %>% replaceColumnNames(columnNames)

# append RMSE mean and max to primerPoolMerged: 
primerPoolMerged = appendToPrimerPoolMerged(primerPoolMerged, r2)
primerPoolMerged = appendToPrimerPoolMerged(primerPoolMerged, RMSEmean)
primerPoolMerged = appendToPrimerPoolMerged(primerPoolMerged, RMSEmax)

# write primerPoolMergedAndQsummary to txt
write.table(primerPoolMerged, file = primerPoolMergedFileName, row.names = F, col.names = T, quote = F, sep = '\t')
message(paste('finished', Qdirectory))