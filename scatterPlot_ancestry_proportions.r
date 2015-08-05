##make scatterplot to compare ancestry proportions between
##whole genome SNPs and AIMs.
##Rscript scatterPlot_ancestry_proportion.r WGS.Q AIMs.Q figure.pdf title

######
#setup
######
library(ggplot2)

#########
#function
#########
##return r2 value as an expression:
lm_eqn = function(m) {

  l <- list(a = format(coef(m)[1], digits = 2),
      b = format(abs(coef(m)[2]), digits = 2),
      r2 = format(summary(m)$r.squared, digits = 3));

  if (coef(m)[2] >= 0)  {
    eq <- substitute(~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(~~italic(r)^2~"="~r2,l)    
  }

  as.character(as.expression(eq));
}

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
fname_wgs = args[1]
fname_aims = args[2]
fname_figure = args[3]
title = args[4]
print(paste0("whole-genome SNPs: ", fname_wgs))
print(paste0("AIMs: ", fname_aims))

##read in Q files:
wgs = read.table(fname_wgs)
aims = read.table(fname_aims)

##check if two ancestral clusters match in the two input Q files:
if (cor(wgs[,1],aims[,1]) > cor(wgs[,1], aims[,2])) {
	message('clusters are in the same order in two input files')
	wgs_prop = wgs[,2]
	aims_prop = aims[,2]
} else {
	message('clusters are in different orders in two input files')
	wgs_prop = wgs[,1]
	aims_prop = aims[,2]
}

prop = data.frame(wgs = wgs_prop, aims = aims_prop)
p = ggplot(prop, aes(x = wgs, y = aims)) + geom_point() + theme_bw() + xlab('Genome-wide ancestry') + ylab('Estimated ancestry') + stat_smooth(method = 'lm', col = 'blue') + 
	geom_text(aes(x = 0.3, y = 0.7, label = lm_eqn(lm(aims ~ wgs, prop))), parse = TRUE) + 
	ggtitle(title)
ggsave(filename = fname_figure, plot = p)
