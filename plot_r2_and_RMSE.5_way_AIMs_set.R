#--------- setup --------------
setwd('/Volumes//ancestry/')
library(data.table)
library(dplyr)
library(gtools)

#--------- functionns ---------
order.Q2.by.Q1 <- function(Q1, Q2) {
  # description: 
  # return an ordering Q2 columns to match the ordering of Q1.
  # note:
  # sometimes the ordering is uncertain for some columns. for example,
  # the order can be c(1,1,3). when this happens, an warning is issued. 
  corr = cor(Q1, Q2)
  order = apply(corr, 1, which.max)
  if (sum(duplicated(order)) > 0) warning("ambiguous ordering: ", order, immediate. = T)
  return(order)
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

guess.number.of.population = function(QfileName){
  # description:
  # guess the number of populations based on
  # the Q file.
  # numPopulations = as.integer(str_extract(QfileName, pattern = "(?<=([A-Z]{3}\\.){2,5}(cumPrimerPool_[0-9]{1,2}\\.){0,1})([0-9])(?=\\.Q)"))
  numPopulations = as.integer(str_extract(QfileName, pattern = "(?<=\\.)([0-9])(?=\\.Q)"))
  
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
  fileNames = mixedsort(fileNames)
  
  r2 = data.frame()
  for (fileName in fileNames){
    # TODO: debug this for loop.
    
    file = fread(fileName)
    message(fileName) # debug
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


#--------- main --------------
# load WGS ancestry: 
wgs_file = 'AIMS_selection/ADMIXTURE/five_superpopulation/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.5.Q'
wgs_ancestry = read.table(wgs_file)



# calculate RMSE: 
num_cols = dim(wgs_ancestry)[2]
num_pools = 36
rmse = data.frame(matrix(0, nrow = num_pools, ncol = num_cols))
setnames(rmse, col_names)

for (pool_index in 1:num_pools){
  # load AIMs ancestry:
  aims_file = sprintf('AIMS_selection/multiplexPrimers/five_superpopulations/with_heterogeneity_filter/filter_SAS//propmultiIn_0.0//admixture//AFR.AMR.EAS.EUR.SAS.cumPrimerPool_%i.5.Q', pool_index)
  aims_ancestry = read.table(aims_file)
  
  # rearrange aims_ancestry columns to match the order of wgs_ancestry:
  columnOrder = order.Q2.by.Q1(wgs_ancestry, aims_ancestry)
  setcolorder(aims_ancestry, columnOrder)
    
  # name columns of aims_ancestry and wgs_ancestry:
  pop_file = "AIMS_selection/multiplexPrimers/five_superpopulations/with_heterogeneity_filter/filter_SAS/propmultiIn_0.0/plink/AFR.AMR.EAS.EUR.SAS.pop"
  col_names = matchColumnsAndPopulations(wgs_file, pop_file)
  setnames(aims_ancestry, col_names)
  setnames(wgs_ancestry, col_names)
  
  for (col_name in col_names){
    rmse[pool_index,col_name] = sqrt(mean((wgs_ancestry[,col_name] - aims_ancestry[,col_name])^2))
  }  
}



# calculate R2:
Qdirectory = 'AIMS_selection/multiplexPrimers/five_superpopulations/with_heterogeneity_filter/filter_SAS//propmultiIn_0.0//admixture'
r2 = calculateAncestryCorrelation(Qdirectory, pattern = "cumPrimerPool_[0-9]{1,2}\\.[2-5]\\.Q$", Q = wgs_ancestry)
r2 = r2[1:num_pools,] # discard empty pool 37
setnames(r2, c('pool',col_names))

# reshape RMSE and R2 dataframes:
rmse$pool = 1:num_pools
rmse_long = melt(data = rmse, measure.vars = col_names, variable.name = 'population', value.name = 'value')
rmse_long$value_name = 'RMSE'
r2_long = melt(data = r2, measure.vars = col_names, variable.name = 'population', value.name = 'value')
r2_long$value_name = 'R2'


# plot RMSE and R2: 
p = ggplot(rbind(rmse_long, r2_long), aes(x = pool, y = value, color = population)) +
 facet_grid(value_name~., scales = 'free', labeller=label_both) + 
 geom_line() + 
 geom_point(size = 3, alpha = 0.8) + 
 theme_bw() + 
 theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))


# change legend position: 
p = p + theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = c(0.7, 0.7))


# adding x and y axis titles:
p = p + ylab(expression(paste("        RMSE              ","Correlation (", R^{2}, ")"))) + xlab('Number of pools')



# change legend title and label: 
p = p + scale_color_discrete(name="",
                           breaks=c("EUR","AMR","SAS","AFR","EAS"),
                           labels=c("European", "Native American", "South Asian", "African", "East Asian"))


# remove grid and facet strip:
p = p + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank()) + theme(strip.background=element_blank(), strip.text.y = element_blank())

ggsave('AIMS_selection/figures/ancestry_r2_rmse_5way_AIMs_set.SAS_filter.propmultiIn_0.0.real.pdf', plot = p, width = 6, height = 6)