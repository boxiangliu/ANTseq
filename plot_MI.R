library(cowplot)

#' calculate mutual information between Q and J 
#' see ancestry paper for definition.
#' @param p_j (numeric) overall probability of J = 1
#' @param p_ij (numeric vector) p_ij is the probability that J = 1 in population Q = i. The length of p_ij equals the number of populations.
calcMI = function(p_j, p_ij){
	K = length(p_ij)
	MI = -p_j * log(p_j) + sum(p_ij/K*log(p_ij)) -(1-p_j) * log(1-p_j) + sum((1-p_ij)/K*log(1-p_ij))
	return(MI)
}

options(echo=TRUE)
args = commandArgs(TRUE)
wd = args[1] 
out = args[2]
setwd(wd)

MI = data.frame(p_j = numeric(), d = numeric(), MI = numeric())
for (p_j in seq(0.1, 0.5, 0.1)){
	for (d in seq(-2*p_j, 2*p_j, 0.01)){
		p_0j = p_j + d/2
		p_1j = p_j - d/2
		MI = rbind(MI, data.frame(p_j = p_j, d = d, MI = calcMI(p_j, c(p_0j,p_1j))))
	}
}
p = ggplot(MI, aes(x = d, y = MI, color = as.factor(p_j))) + geom_line() + xlab('Delta') + ylab('Mutual information') + scale_color_discrete(name = expression(paste("  ",p[j])))
save_plot(out,p)