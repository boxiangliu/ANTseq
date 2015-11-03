######
#setup
######
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


#########
#function
#########
detect_changes = function(x, include_first_and_last_element = T){
	# description: 
	# detect changes in value for a vector
	# note: 
	# The vector generally contains stretches of repeated values like
	# [a a a a b b b c c c c]
	x = Q_rearranged$population
	lag_by_1 = lag(x)
	changes = which(lag_by_1 != x)
	if (include_first_and_last_element){
		changes = c(1, changes, length(x))
	}
	return(changes)
}

mean_of_adjacent_elements = function(x){
	# description:
	# calculate the means of adjacent elements of a vector.
	lead_by_1 = lead(x)
	mean_of_adjacent_elements = apply(cbind(x, lead_by_1), 1, mean)
	mean_of_adjacent_elements = mean_of_adjacent_elements[!is.na(mean_of_adjacent_elements)]
}

generate_axis_label = function(x){
	# description: 
	# Make a vector of spaces mixed with labels,
	# useful if there are many repeated labels and you want to 
	# keep one and make others blank. 
	changes = detect_changes(x)
	label_positions = ceiling(mean_of_adjacent_elements(changes))
	labels = rep("", length(x))
	labels[label_positions] = unique(x)
	return(labels)
}

#####
#main
#####
# read 1000 Genome panel file: 
panel_file = "/srv/persistent/bliu2/shared/1000genomes/phase3v5/integrated_call_samples_v3.20130502.ALL.panel"
panel = fread(panel_file, header = F)
setnames(panel, c("ID","subpopulation", "population", "gender"))

# read genomic ancestry (.Q) file:
Q_file = "/srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE/five_superpopulation/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.5.Q"
Q = fread(Q_file)
setnames(Q, c("EUR","AMR","SAS","AFR","EAS"))

# TODO: read AIMs ancestry file: 

# add population and ID to Q 
Q$population = panel$population
Q$ID = panel$ID

# split Q by population: 
Q_split_by_population = split(x = Q, f = Q$population)

# for each population, order individuals by ancestry estimate:
Q_split_by_population_rearranged = list()
for (population in c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')){
	Q_split_by_population_rearranged[[population]] = arrange_(Q_split_by_population[[population]], paste0("-",population))
}

# combine split Q's: 
Q_rearranged = do.call(rbind, Q_split_by_population_rearranged)

# plot ancestry estimates:
x_axis_labels = generate_axis_label(Q_rearranged$population)
pdf_file = "/srv/persistent/bliu2/ancestry/AIMS_selection/figures/genomic_ancestry_estimates.pdf"
pdf(file = pdf_file, width = 10, height = 5)
	barplot(t(as.matrix(Q_rearranged %>% select(1:5))), space = 0, col = rainbow(5), ylab="Genomic ancestry", border=NA, xpd=TRUE, names.arg = x_axis_labels, las = 2)
	# legend will not show up on saved pdf. 
	# legend("bottomright", inset = c(-0.2, 0.1), legend = populations, col = rainbow(5), lty = 1)
dev.off()
