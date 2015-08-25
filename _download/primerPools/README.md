## The primers pools are names with the following convention: 

{populations}.primerPoolMerged.txt

## Examples: 
AFR.AMR.EUR.primerPoolMerged.txt contains primers for AIMs selected to distinguish African, Native American, and European populations.

## Population Abbreviations: 
AFR = African 
AMR = Native American
EAS = East Asian
EUR = European
SAS = South Asian

## Columns: 
snp = SNP ID by chr_position  
forward = forward primer
fwd_tmp = forward primer Tm
reverse = reverse primer
rev_tmp = reverse primer Tm
range = chr:start-end
fwd_id = forward primer ID
fwd_adaptor = forward primer with adaptor (GCGTTATCGAGGTC)
rev_id = reverse primer ID
rev_adaptor = reverse primer with adaptor (GTGCTCTTCCGATCT)
pool = primer pool index during the design phase
{pop1}_AF = allele frequency for pop1 
{pop2}_AF = allele frequency for pop2 
In = Overall Rosenberg's Informativeness Criteria. 
{pop1}_{pop2}_In = Rosenberg's Informativeness Criteria for pop1 and pop2; the same as In for 2-population comparisons.
cum_In = cumulative In for each primer pool
poolRank = rank of each pool based on cum_In
r2 = correlation between estimated ancestry and true ancestry using all primers up to the current pool 

## Columns for AFR.AMR.EUR.primerPoolMerged.txt: 
LSBL.In = Locus specific branch length based on In used to calculate 3-population comparisons
LSBL.Fst = Locus specific branch length based on Fst used to calculate 3-population comparisons
population = the population for which LSBL.In and LSBL.Fst is calculated. 
EUR_r2, AMR_r2, AFR_r2 = correlation between estimated ancestry and true ancestry using all primers up to the current pool calculated for each population
The rest have the same meanings. 

Interested Readers can refer to: 
*Galanter, J.M., Fernandez-Lopez, J.C., Gignoux, C.R., Barnholtz-Sloan, J., Fernandez-Rozadilla, C., Via, M., Hidalgo-Miranda, A., Contreras, A.V., Figueroa, L.U., Raska, P., et al. (2012). Development of a Panel of Genome-Wide Ancestry Informative Markers to Study Admixture Throughout the Americas. PLoS Genet. 8, e1002554.*

## Usage:
In general, the correlation increases as the number of the genotyped AIMs. We find that 120 AIMs give satisfying estimate because in most cases they reach 0.99 correlation. When ordering primers, you should choose the top-ranking pools.

## FAQ:
1. How many primers does each pool have? 
In general, each primer pools contains 10 primers. However, one pool has 11 and several other has less than 10. We make 10-plex pools for ease of calculation and pipetting so the number of primer pools in each pool should not really matter. 


