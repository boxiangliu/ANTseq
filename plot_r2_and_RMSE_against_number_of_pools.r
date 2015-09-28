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

##########
#constants
##########
FIGURES_DIR = "AIMS_selection/figures"

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

read_r2 = function(population_group, filename_template){
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
		facet_grid(.~population_group) +
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
	rmse = melt(primer_pool_merged, id.vars = "poolRank", measure.vars = column_index)
	rmse = rmse %>% mutate(population = str_split_fixed(variable, "\\.", n = 2)[,1], statistic = str_split_fixed(variable, "\\.", n = 2)[,2], variable = NULL)
	rmse = rmse %>% mutate(statistic = str_replace(statistic, "rmse", ""))
	return(rmse)
}

read_rmse = function(population_group, filename_template){
	# description:
	# read primerPoolMerged.txt specified by [population_group, filename_template], 
	# select rmseMean and rmseMax columns, and `rbind` them.    
	rmse = data.frame(poolRank = integer(), value = numeric(), population = character(), statistic = character(), population_group = character())
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
		facet_grid(.~population_group) + 
		xlab("Number of pools") + ylab("RMSE") + scale_color_discrete(name = "Population") + scale_linetype_discrete(name = "Statistic") +
		theme_bw() +
		theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title = element_text(size = 15), legend.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 12)) 
	return(p)
}

#####
#main
#####
population_groups = c("AFR.AMR.EUR", "AFR.EAS", "AFR.EUR", "AFR.SAS", "EAS.EUR", "EAS.SAS", "EUR.SAS")
filename_template = "AIMS_selection/multiplexPrimers/%s/primerPoolMerged.txt"

# read ancestry correlation columns of all population_groups into one data.frame:
r2 = read_r2(population_group, filename_template)
r2_plot = plot_r2(r2)
ggsave(filename = paste(FIGURES_DIR, "r2.pdf", sep = "/"), r2_plot, width = 30, height = 6)

# plot RMSE mean and max against number of pools:
rmse = read_rmse(population_group, filename_template)
rmse_plot = plot_rmse(rmse)
ggsave(filename = paste(FIGURES_DIR, "rmse.pdf", sep = "/"), rmse_plot, width = 30, height = 6)

