####------ setup -------#####
setwd('/Volumes/ancestry/AIMS_selection/')
library(reshape2)
library(ggplot2)
library(Rmisc)

####------ main --------####
# read aims file:
aims = read.table('AIMs/five_superpopulations/with_heterogeneity_filter/filter_SAS/AFR_AMR_EAS_EUR_SAS_with_heterogeneity_filter500k_500.aims', header = T)

# plot allele frequencies: 
png('figures/allele_frequencies.5_way_AIMs.png', width = 2000, height = 400)
plot(aims$AFR_AF, type = 'l', col = 1, xlab='Number of AIMs', ylab = 'allele frequency or In', main = 'allele frequency at each AIM')
lines(aims$AMR_AF, col = 2)
lines(aims$EAS_AF, col = 3)
lines(aims$EUR_AF, col = 4)
lines(aims$SAS_AF, col = 5)
lines(aims$In, col = 6)
legend('topright', legend = c('AFR','AMR','EAS','EUR','SAS', 'In'), col = 1:6, lty = 1)
dev.off()

# plot pairwise In: 
png('figures/In.5_way_AIMs.png', width = 2000, height = 400)
plot(aims[,9],type = 'l', ylim = c(0,1), main = 'Informativeness Criteria (In)', xlab = 'Number of AIMs', ylab = 'pairwise In')
for (i in 10:19){
  lines(aims[,i], col = i)  
}
legend('topright', legend = colnames(aims)[9:19], col = 9:19, lty = 1)
dev.off()

# plot pairwise In distribution: 
aims$index = 1:nrow(aims)
pairwise_In = melt(aims, id.vars = 'index', measure.vars = colnames(aims)[10:19], value.name = 'In', variable.name = 'pair')

png('figures/In_distribution.5_way_AIMs.png')
plot1 = ggplot(pairwise_In %>% filter(index > 250), aes(x = In, fill = pair)) + geom_histogram(position = 'dodge') + ggtitle('AIMs selected using overall In')
plot2 = ggplot(pairwise_In %>% filter(index <= 250), aes(x = In, fill = pair)) + geom_histogram(position = 'dodge') + ggtitle('AIMs selected using pairwise In') 
multiplot(plot1, plot2)
dev.off()


