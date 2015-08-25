##plot the correlation of AIMs and true ancestry proportions agains the number of AIMs.  
##usage:
##Rscript plot_ancestry_r2_vs_num_AIMs.primerPools.r pop1,pop2[,pop3,...] numPools
##example: 
##Rscript plot_ancestry_r2_vs_num_AIMs.primerPools.r AFR,EUR 35


#######
#setup
#######
library(ggplot2)
library(reshape2)
#######
#funct
#######
##function to reorder columns of aims df to match the order of wgs df: 
reorder_cols <- function(wgs, aims) {
	corr = cor(wgs, aims)
	aims_reorder = aims
	for (i in 1:nrow(corr)) {
		aims_reorder[,i] = aims[,which.max(corr[i,])]
	}
	return(aims_reorder)
}

#####
#main
#####
##read command line args: 
args <- R.utils::commandArgs(trailingOnly = TRUE)
pops = strsplit(args[1], ",")[[1]]
numPops = length(pops)
numPools = args[2]

message(sprintf("populations: %s\n", args[1]))
message(sprintf("number of primer pools: %s\n", numPools))

##read WGS Q file: 
fname_wgs = sprintf("%s.%i.Q", paste(pops, collapse="."), numPops)
wgs = read.table(fname_wgs)

##calculate the r2 between aims and wgs ancestry proportions: 
r2 = data.frame()
for (n in 1:numPools) {
	aims_Q = read.table(sprintf("%s.cumPrimerPool_%s.%i.Q", paste(pops, collapse = "."), n, numPops), header = F)
	aims_P = read.table(sprintf("%s.cumPrimerPool_%s.%i.P", paste(pops, collapse = "."), n, numPops), header = F)
	numMarkers = nrow(aims_P)

	##reorder columns of aims df to match the order of wgs df: 
	aims_reorder = reorder_cols(wgs, aims_Q)
	##skip if the ordering of aims dataframe is ambiguous:
	if (length(unique(colSums(aims_reorder))) < ncol(aims_reorder)) {next}
	##calculate correlation:
	r2 = rbind(r2, data.frame(nmarkers = numMarkers, t(diag(cor(wgs,aims_reorder)))))
}

##rename columns: 
names(r2)[-1] = paste0("pop", seq(3))
write.table(data.frame(npool = seq(nrow(r2)),r2),sprintf("%s.ancestry_r2_vs_num_AIMs.txt", paste(pops, collapse = ".")), quote = F, col.names = T, row.names = F)
##melt to long format: 
r2_long = melt(r2, id.vars = "nmarkers", variable.name = 'pop', value.name = 'r2')

##plot r2 as a function of number of markers:
p = ggplot(r2_long, aes(x = nmarkers, y = r2, color = pop, label = r2)) + geom_point() + geom_line() + theme_bw() + 
	xlab("Number of AIMs") + ylab("Correlation with true ancestry") + 
	# geom_text() +   
	ggtitle(sprintf("%s Ancestry Proportion Correlation", paste(pops, collapse = "-")))
ggsave(sprintf("%s.ancestry_r2_vs_num_AIMs.pdf", paste(pops, collapse = ".")), p)