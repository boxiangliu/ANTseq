library(stringr)
library(cowplot)
source("utils.R")

#------ functions -------
#' calculate the mean absolute deviation 
#' @param x a numeric vector
#' @param mean a numeric scalar
#' @return the mean absolute scalar
mean_abs_dev = function(x, mean = NULL){
	if (is.null(mean)){
		mean == mean(x)
	}
	return(mean(abs(x - mean)))
}

#------- main ---------
# read command line args: 
options(echo=TRUE)
args = commandArgs(TRUE)
wd = args[1]
inp = args[2]
out = args[3]
# wd = '../AIMS_selection/allele_freq_unfiltered/'
# inp = 'sample_list.txt'
# out = '../figures/AF_distribution.pdf'

setwd(wd)


# read samples: 
sample_list = scan(inp, what = character())
af = data.table() # allele frequency
for (sample in sample_list){
	temp = fread(sample, header = TRUE)
	af = rbind(af, temp)
}
setnames(af, '#CHROM', 'CHROM')
multi_allelic = which(str_detect(af$ALT, ",")) # remove multi-allelic loci. 
af = af[-multi_allelic, ]
af[,AF := as.numeric(AF)]
af[,EAS_AF := as.numeric(EAS_AF)]
af[,AMR_AF := as.numeric(AMR_AF)]
af[,AFR_AF := as.numeric(AFR_AF)]
af[,EUR_AF := as.numeric(EUR_AF)]
af[,SAS_AF := as.numeric(SAS_AF)]


# calculate minor allele frequencies: 
af[, MAF := ifelse(AF <= 0.5, AF, 1-AF)]


# make histogram of AF distribution: 
pct = sum(af$MAF < 0.05)/length(af$MAF)


# plot_maf = ggplot(af, aes(x = MAF)) + geom_histogram(binwidth = 0.05) + xlab("Minor allele frequency") + geom_vline(xintercept = 0.05, color = 'red', linetype = 3) + annotate('text', x = 0.3, y = 6e7, size = 7, label = sprintf("MAF < 0.05: %0.2f %% ",pct*100)) 
plot_maf = ggplot(af, aes(x = MAF)) + geom_histogram(binwidth = 0.05) + xlab("Minor allele frequency") + geom_vline(xintercept = 0.05, color = 'red', linetype = 3) 
plot_maf1 = add_sub(plot_maf, sprintf("MAF < 0.05: %0.2f %% ",pct*100), y = 15, x = 0.5)


# calculate mean absolute deviation (MAD): 
af[, mad := mean_abs_dev(x = c(EAS_AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF), mean = AF), by = .(CHROM, POS)]


# plot AF vs MAD: 
set.seed(1)
idx = sample(1:nrow(af), 10000, replace = FALSE)
plot_maf_vs_mad = ggplot(af[idx,], aes(x = MAF, y = mad)) + geom_point() + xlab('Minor allele frequency') + ylab('Mean absolute deviation')  + geom_vline(xintercept = 0.05, linetype = 3, color = 'red')



# make paneled plot and save:
plot_1by2 = plot_grid(plot_maf, plot_maf_vs_mad, labels = c('A','B'), align = 'h')
save_plot(out, plot_1by2, base_height = 4, base_aspect_ratio = 2)