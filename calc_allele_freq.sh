#!/bin/bash

##calculate the allele frequency for 5 super populations in 1000 Genomes catalog.
##AFR contains two admixed subpopulations (ASW and ACB0 that needs to be removed. 
##In addition, calculate global LD r2 using all individuals in 1000 Genomes catalog. 

wd="/srv/persistent/bliu2/ancestry/AIMS_selection"
input="$wd/phase3v5_filtered/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.vcf.gz"
panel="/srv/persistent/bliu2/shared/1000genomes/phase3v5/integrated_call_samples_v3.20130502.ALL.panel"
output_dir="$wd/alleleFreq"

##array of superpopulations:
##awk '{print $3}' integrated_call_samples_v3.20130502.ALL.panel | sort | uniq | grep -v 'pop' | awk '{print $1" \\"}'
populations=(AFR \
AMR \
EAS \
EUR \
SAS)

##calculate allele frequency for each population: 
pids=""
for pop in ${populations[*]}; do
	pop_file="$output_dir/$pop.pop"
	cat $panel | grep $pop | awk '{print $1, $1, $2}' > $pop_file
	plink --vcf $input --double-id --snps-only --freq --make-bed --keep-allele-order --keep $pop_file --out $output_dir/$pop &
	echo "calculating allele frequency of population: $pop; pid: $!"
	pids="$pids $!"
done 

##calculate global LD r2:
plink --vcf $input --double-id --snps-only --r2 --ld-window-kb 2000 --ld-window 99999999 --ld-window-r2 0.2 --make-bed --out $output_dir/global_r2_0.2_window_2000k &
pids="$pids $!"
echo "calculating LD; pid $!"
wait $pids

##calculate the allele frequency for AFR without ACB (African Caribbean in Barbados) and ASW (African South West).
pop='AFR-ACB-ASW'
pop_file="$output_dir/$pop.pop"
cat $panel | grep AFR | grep -v ASW | grep -v ACB | awk '{print $1, $1, $2}' > $pop_file
plink --vcf $input --double-id --snps-only --freq --make-bed --keep-allele-order --keep $pop_file --out $output_dir/$pop &
echo "calculating allele frequency of population: $pop; pid: $!"
wait $!
