#!/usr/bin/env Rscript 

######
#setup
######
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))

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
for (i in 1:length(r2List)) {
	if (ncol(r2List[[i]]) == 1){
		temp = data.frame(pop = names(r2List)[i], r2 = unique(r2List[[i]]))
		temp = temp %>% mutate(npool = seq(1, nrow(temp)))
		r2List_long = rbind(r2List_long, temp)
	} else {
		temp = unique(r2List[[i]])
		temp = temp %>% mutate(npool = seq(1, nrow(temp)))
		temp_long = melt(temp, id.vars = "npool", measure.vars = 1:ncol(r2List[[i]]), value.name = "r2", variable.name = "pop")
		temp_long = temp_long %>% mutate(pop = paste(names(r2List)[i], str_replace(pop, "_r2", ""), sep = ":"))
		r2List_long = rbind(r2List_long, temp_long)
	}
}


##plot r2List_long:
p = ggplot(r2List_long, aes(x = npool, y = r2, color = pop)) + geom_line() + geom_point() + theme_bw() + 
	xlab('Pool number') + ylab(expression(bold(R^"2"))) + 
	scale_color_discrete(name = 'Populations') +
	theme(legend.position = c(1,0), legend.justification = c(1,0), legend.title = element_text(size = 15), legend.text = element_text(size = 12), axis.title = element_text(size = 20, face = 'bold'), axis.text = element_text(size = 12)) 

ggsave(filename = "/srv/persistent/bliu2/ancestry/AIMS_selection/multiplexPrimers/primerPools/ancestry_r2_vs_num_pools.pdf", plot = p)


##tabulate r2 using all pools: 
r2max = r2List_long %>% group_by(pop) %>% summarize(r2 = max(r2), npool = length(r2)) 
write.csv(r2max, "/srv/persistent/bliu2/ancestry/AIMS_selection/multiplexPrimers/primerPools/r2_using_all_pools.csv", row.names = FALSE)