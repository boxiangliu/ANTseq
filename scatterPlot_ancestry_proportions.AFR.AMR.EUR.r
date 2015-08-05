##make scatterplot to compare ancestry proportions between
##whole genome SNPs and AIMs for AFR, AMR and EUR


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

##function to reorder columns of aims df to match the order of wgs df: 
reorder_cols <- function(wgs, aims) {
	corr = cor(wgs, aims)
	aims_reorder = aims
	for (i in 1:nrow(corr)) {
		aims_reorder[,i] = aims[,which.max(corr[i,])]
	}
	return(aims_reorder)
}

##fnames: 
fname_wgs = "/srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE_noSubpop/AFR.AMR.EUR/AFR.AMR.EUR.3.Q"
fname_aims = "/srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE_noSubpop/AFR.AMR.EUR/AFR.AMR.EUR.446.AIMs.3.Q"

##read in Q files:
wgs = read.table(fname_wgs)
aims = read.table(fname_aims)
names(wgs) = paste0('wgs', seq(1,3))
names(aims) = paste0('aims', seq(1,3))

##reorder columns of aims df to match the order of wgs df: 
aims_reorder = reorder_cols(wgs, aims)

##make scatterplot for each pair of populations:
prop = data.frame(wgs, aims_reorder) # prop = ancestry proportion

title = 'EUR ancestry proportion'
p1 = ggplot(prop, aes(x = wgs1, y = aims1)) + geom_point() + theme_bw() + xlab('Genome-wide ancestry') + ylab('Estimated ancestry') + stat_smooth(method = 'lm') + geom_text(aes(x = 0.3, y = 0.7, label = lm_eqn(lm(aims1 ~ wgs1, prop))), parse = TRUE) + ggtitle(title)
fname_figure = "/srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE_noSubpop/AFR.AMR.EUR/AFR_AMR_EUR.446aims.EURancestryProportions.pdf"
ggsave(filename = fname_figure, plot = p1)

title = 'AMR ancestry proportion'
p2 = ggplot(prop, aes(x = wgs2, y = aims2)) + geom_point() + theme_bw() + xlab('Genome-wide ancestry') + ylab('Estimated ancestry') + stat_smooth(method = 'lm') + geom_text(aes(x = 0.3, y = 0.7, label = lm_eqn(lm(aims2 ~ wgs2, prop))), parse = TRUE) + ggtitle(title)
fname_figure = "/srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE_noSubpop/AFR.AMR.EUR/AFR_AMR_EUR.446aims.AMRancestryProportions.pdf"
ggsave(filename = fname_figure, plot = p2)

title = 'AFR ancestry proportion'
p3 = ggplot(prop, aes(x = wgs3, y = aims3)) + geom_point() + theme_bw() + xlab('Genome-wide ancestry') + ylab('Estimated ancestry') + stat_smooth(method = 'lm') + geom_text(aes(x = 0.3, y = 0.7, label = lm_eqn(lm(aims3 ~ wgs3, prop))), parse = TRUE) + ggtitle(title)
fname_figure = "/srv/persistent/bliu2/ancestry/AIMS_selection/ADMIXTURE_noSubpop/AFR.AMR.EUR/AFR_AMR_EUR.446aims.AFRancestryProportions.pdf"
ggsave(filename = fname_figure, plot = p3)

