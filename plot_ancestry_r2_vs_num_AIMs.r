##plot the correlation of AIMs and true ancestry proportions agains the number of AIMs.  
##example: 
##Rscript plot_ancestry_r2_vs_num_AIMs.r /srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE AFR AMR 500


#######
#setup
#######
library(ggplot2)

######
#main
######
##read command line args: 
options(echo=TRUE) # if you want see commands in output file
args <- R.utils::commandArgs(trailingOnly = TRUE)
admixture_dir = args[1]
pop1 = args[2]
pop2 = args[3]
if (!is.na(args[4])) {
	nmarkers = as.integer(args[4])
} else {
	nmarkers = 500 #default.
}
message(sprintf("population 1: %s\n", pop1))
message(sprintf("population 2: %s\n", pop2))
message(sprintf("number of AIMs markers: %i\n", nmarkers))

#set working directory: 
setwd(sprintf("%s/%s.%s", admixture_dir, pop1, pop2))

##read WGS Q file: 
fname_wgs = sprintf("%s.%s.2.Q", pop1, pop2)
wgs = read.table(fname_wgs)

##calculate the r2 between aims and wgs ancestry proportions: 
r2 = data.frame()
for (n in 1:nmarkers) {
	fname_aims = sprintf("%s.%s.%i.AIMs.2.Q", pop1, pop2, n)
	aims = read.table(fname_aims)
	r2 = rbind(r2, data.frame(nmarkers = n, r2 = max(cor(wgs[,1], aims[,1]), cor(wgs[,1], aims[,2]))))
}

##plot r2 as a function of number of markers:
p = ggplot(r2, aes(x = nmarkers, y = r2)) + geom_point() + theme_bw() + 
	xlab("Number of AIMs") + ylab("Correlation with true ancestry") + 
	ggtitle(sprintf("%s-%s Ancestry Proportion Correlation", pop1, pop2))
ggsave(sprintf("%s_%s.ancestry_r2_vs_num_AIMs.pdf", pop1, pop2), p)
