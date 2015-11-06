# description: 
# generate mock primer pools using 
# AIMs selected by AIMs_generator.py, 
# effectively bypassing the yamPCR step
# to speed up testing.

####----- setup -----#####
setwd('/srv/persistent/bliu2/ancestry')

####----- main -----####
# read aims: 
aims_file = 'AIMS_selection/AIMs/five_superpopulations/with_heterogeneity_filter/filter_SAS/AFR_AMR_EAS_EUR_SAS_with_heterogeneity_filter500k_500.aims'
aims = read.table(aims_file, header = T)





# write mock primer pools, 
# each with 10 primers. 
num_primers = nrow(aims)
num_primers_per_pool = 10 
num_pools = floor(num_primers/num_primers_per_pool)
out_dir = 'AIMS_selection/AIMs/five_superpopulations/with_heterogeneity_filter/filter_SAS/mock_primer_pools'
for (i in 0:(num_pools-1)){
	start = i*num_primers_per_pool + 1 
	end = (i+1)*num_primers_per_pool
	primer_pool = data.frame(snp = aims$snp[start:end])
	out_file = paste0(out_dir, sprintf('/primerPool_%i.txt', i+1))
	write.table(primer_pool, file = out_file, col.names = T, row.names = F, quote = F)
}
