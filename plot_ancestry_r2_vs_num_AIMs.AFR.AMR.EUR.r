##plot the correlation of AIMs and true ancestry proportions agains the number of AIMs.  
##example: 
##Rscript plot_ancestry_r2_vs_num_AIMs.r /srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE AFR AMR 500


#######
#setup
#######
library(ggplot2)

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

######
#main
######
#set working directory: 
setwd("/srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE_noSubpop/AFR.AMR.EUR")

##read WGS Q file: 
fname_wgs = "AFR.AMR.EUR.3.Q"
wgs = read.table(fname_wgs)


##calculate the r2 between aims and wgs ancestry proportions: 
r2 = data.frame()
nmarkers = 446
for (n in 1:nmarkers) {
	fname_aims = sprintf("AFR.AMR.EUR.%i.AIMs.3.Q", n)
	aims = read.table(fname_aims)
	
	##reorder columns of aims df to match the order of wgs df: 
	aims_reorder = reorder_cols(wgs, aims)

	##skip if the ordering of aims dataframe is ambiguous:
	if (length(unique(colSums(aims_reorder))) < ncol(aims_reorder)) {next}

	r2 = rbind(r2, data.frame(nmarkers = n, 
							  EUR.r2 = cor(wgs[,1], aims_reorder[,1]), 
							  AMR.r2 = cor(wgs[,2], aims_reorder[,2]), 
							  AFR.r2 = cor(wgs[,3], aims_reorder[,3])
							  ))
}

##remove nmarkers = c(5,6)
r2 = r2[!r2$nmarkers %in% c(5,6),]

##plot r2 as a function of number of markers:
pop = 'EUR'
p1 = ggplot(r2, aes_string(x = 'nmarkers', y = paste0(pop,'.r2'))) + geom_point() + theme_bw() + xlab("Number of AIMs") + ylab("Correlation with true ancestry") + ggtitle(sprintf("%s Ancestry Proportion Correlation", pop))
ggsave(sprintf("AFR_AMR_EUR.%s.ancestry_r2_vs_num_AIMs.pdf", pop), p1)

pop = 'AMR'
p2 = ggplot(r2, aes_string(x = 'nmarkers', y = paste0(pop,'.r2'))) + geom_point() + theme_bw() + xlab("Number of AIMs") + ylab("Correlation with true ancestry") + ggtitle(sprintf("%s Ancestry Proportion Correlation", pop))
ggsave(sprintf("AFR_AMR_EUR.%s.ancestry_r2_vs_num_AIMs.pdf", pop), p2)

pop = 'AFR'
p3 = ggplot(r2, aes_string(x = 'nmarkers', y = paste0(pop,'.r2'))) + geom_point() + theme_bw() + xlab("Number of AIMs") + ylab("Correlation with true ancestry") + ggtitle(sprintf("%s Ancestry Proportion Correlation", pop))
ggsave(sprintf("AFR_AMR_EUR.%s.ancestry_r2_vs_num_AIMs.pdf", pop), p3)
