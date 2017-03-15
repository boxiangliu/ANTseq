######
#setup
######
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
rm(list = ls())

##########
#constants
##########
# ANCESTRY_DIR = "/srv/persistent/bliu2/ancestry"
ANCESTRY_DIR = "/Volumes/ancestry"

FIGURES_DIR = paste(ANCESTRY_DIR, "AIMS_selection/figures", sep = "/")

#########
#function
#########
locate = function(character_vector, pattern){
	# description: 
	# find all strings in [character_vector] that contain [pattern] and
	# return their indices.
	return(which(str_detect(character_vector, pattern)))
}

select_r2_columns = function(primer_pool_merged, column_index){
	# description:
	# select r2 columns from [primer_pool_merged] 
	# value: 
	# a data.frame (or data.table) in long format. 
	r2 = melt(primer_pool_merged, id.vars = "poolRank", measure.vars = column_index, variable.name = "population", value.name = "r2")
	r2 = r2 %>% mutate(population = str_replace(population, "\\.r2", ""))
	return(r2)
}

read_r2 = function(population_groups, filename_template){
	# description:
	# read primerPoolMerged.txt specified by [population_group, filename_template], 
	# select r2 columns, and `rbind` them.    
	r2 = data.frame(poolRank = integer(), population = character(), r2 = numeric(), population_group = character())
	for (population_group in population_groups){
		# read primerPoolMerged:
		primer_pool_merged_filename = sprintf(filename_template, population_group)
		primer_pool_merged = fread(primer_pool_merged_filename)

		# plot r2 against number of pools: 
		r2_column_index = locate(colnames(primer_pool_merged), "r2")
		temp = select_r2_columns(primer_pool_merged, r2_column_index)
		temp = data.frame(temp, population_group = population_group)
		r2 = rbind(r2, temp)
	}
	return(r2)
}

plot_r2 = function(r2){
	# description:
	# plot r2 against number of pools.
	p = ggplot(r2, aes(x = poolRank, y = r2, color = population)) + 
		geom_line() + geom_point() + 
		facet_grid(population_group~., scales = 'free') +
		xlab("Number of pools") + ylab(expression(bold(R^"2"))) + scale_color_discrete(name = "Population") + 
		theme_bw() + 
		theme(legend.position = c(1,0), legend.justification = c(1,0), legend.title = element_text(size = 15), legend.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 12)) 
	return(p)
}



select_rmse_columns = function(primer_pool_merged, column_index){
	# description:
	# select rmse columns from [primer_pool_merged] 
	# value: 
	# a data.frame (or data.table) in long format. 
	rmse = melt(primer_pool_merged, id.vars = "poolRank", measure.vars = column_index, variable.name = "population", value.name = "rmse")
	rmse = rmse %>% mutate(population = str_replace(population, "\\.rmse", ""))
	return(rmse)
}



read_rmse = function(population_groups, filename_template){
	# description:
	# read primerPoolMerged.txt specified by [population_group, filename_template], 
	# select rmseMean and rmseMax columns, and `rbind` them.    
	rmse = data.frame(poolRank = integer(), population = character(), rmse = numeric(), population_group = character())
	for (population_group in population_groups){
		# read primerPoolMerged:
		primer_pool_merged_filename = sprintf(filename_template, population_group)
		primer_pool_merged = fread(primer_pool_merged_filename)

		# select r2 columns: 
		rmse_column_index = locate(colnames(primer_pool_merged), "rmse")
		temp = select_rmse_columns(primer_pool_merged, rmse_column_index)
		temp = data.frame(temp, population_group = population_group)
		rmse = rbind(rmse, temp)
	}
	return(rmse)
}

plot_rmse = function(rmse){
	# description: 
	# plot (mean, max) of RMSE against number of pools.
	p = ggplot(rmse, aes(x = poolRank, y = value, color = population, linetype = statistic, group = interaction(population, statistic))) + 
		geom_line() + geom_point() + 
		facet_grid(population_group~.) + 
		xlab("Number of pools") + ylab("RMSE") + scale_color_discrete(name = "Population") + scale_linetype_discrete(name = "Statistic") +
		theme_bw() +
		theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title = element_text(size = 15), legend.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 12)) 
	return(p)
}

#####
#main
#####
population_groups = c("AFR.AMR.EUR", "AFR.EAS", "AFR.EUR", "AFR.SAS", "EAS.EUR", "EAS.SAS", "EUR.SAS")
filename_template = paste(ANCESTRY_DIR, "AIMS_selection/multiplexPrimers/pairwise/%s/primerPoolMerged.txt", sep = "/")

# read ancestry correlation coefficient:
r2 = read_r2(population_groups, filename_template)

# rename population factors:
levels(r2$population_group) = str_replace_all(levels(r2$population_group), '\\.', '+')

# plot correlation coefficient for pairs in {AFR, EAS, EUR, SAS}
p = ggplot(r2 %>% filter(population_group != "AFR+AMR+EUR"), aes(x = poolRank, y = r2)) + 
	geom_line() + geom_point() + 
	facet_grid(population_group~.) +
	xlab("Number of pools") + ylab(expression(R^"2")) + 
	theme_bw() + 
	theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20)) +
	theme(strip.text.y = element_text(size = 20)) + 
	theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank())
ggsave('AIMS_selection/figures/r2.AFR_EAS_EUR_SAS.pdf',p, width = 6, height = 12)
ggsave('AIMS_selection/figures/r2.AFR_EAS_EUR_SAS.png',p, width = 6, height = 12)

# plot correlation coefficient for triplet AFR+EUR+AMR:
p2 = ggplot(r2 %>% filter(population_group == "AFR+AMR+EUR"), aes(x = poolRank, y = r2, color = population)) + 
	geom_line() + geom_point() + 
	xlab("Number of pools") + ylab(expression(R^"2")) + 
	theme_bw() + 
	theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20)) +
	theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank())
ggsave('AIMS_selection/figures/r2.AFR_AMR_EUR.pdf',p2, width = 6, height = 6)
ggsave('AIMS_selection/figures/r2.AFR_AMR_EUR.png',p2, width = 6, height = 6)



# read RMSE:
rmse = read_rmse(population_groups, filename_template)

# rename RMSE factors levels:
levels(rmse$population_group) = str_replace_all(levels(rmse$population_group), '\\.', '+')

# plot RMSE mean and max against number of pools:
# plot correlation coefficient for pairs in {AFR, EAS, EUR, SAS}
p3 = ggplot(rmse %>% filter(population_group != "AFR+AMR+EUR"), aes(x = poolRank, y = rmse)) + 
	geom_line() + geom_point() + 
	facet_grid(population_group~.) +
	xlab("Number of pools") + ylab("RMSE") + 
	theme_bw() + 
	theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20)) +
	theme(strip.text.y = element_text(size = 20)) + 
	theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank())

ggsave('AIMS_selection/figures/rmse.AFR_EAS_EUR_SAS.pdf',p3, width = 6, height = 12)
ggsave('AIMS_selection/figures/rmse.AFR_EAS_EUR_SAS.png',p3, width = 6, height = 12)


# plot RMSE for pairs in {AFR, EAS, EUR, SAS}:
p4 = ggplot(rmse %>% filter(population_group == "AFR+AMR+EUR"), aes(x = poolRank, y = rmse, color = population)) + 
	geom_line() + geom_point() + 
	xlab("Number of pools") + ylab("RMSE") + 
	theme_bw() + 
	theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20)) +
	theme(strip.text.y = element_text(size = 20)) + 
	theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank()) +
	theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = c(0.7, 0.7)) + 
	scale_color_discrete(name="",breaks=c("EUR","AMR","AFR"),labels=c("European", "Native American", "African"))

ggsave('AIMS_selection/figures/rmse.AFR_AMR_EUR.pdf',p4, width = 6, height = 6)
ggsave('AIMS_selection/figures/rmse.AFR_AMR_EUR.png',p4, width = 6, height = 6)

