#!/bin/bash
wd="/srv/persistent/bliu2/ancestry/AIMS_selection"

input="$wd/phase3v5_filtered/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.vcf.gz"
panel="/srv/persistent/bliu2/shared/1000genomes/phase3v5/integrated_call_samples_v3.20130502.ALL.panel"
output_dir="$wd/alleleFreq"

##array of all populations: 
## awk '{print $2}' integrated_call_samples_v3.20130502.ALL.panel | sort | uniq | grep -v 'pop' | awk '{print $1" \\"}'
# populations=(ACB \
# ASW \
# BEB \
# CDX \
# CEU \
# CHB \
# CHS \
# CLM \
# ESN \
# FIN \
# GBR \
# GIH \
# GWD \
# IBS \
# ITU \
# JPT \
# KHV \
# LWK \
# MSL \
# MXL \
# PEL \
# PJL \
# PUR \
# STU \
# TSI \
# YRI)

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

# plink --vcf $input --double-id --snps-only --r2 with-freqs --ld-window-kb 2000000 --ld-window 99999999 --ld-window-r2 0.1 --make-bed --out $output_dir/global &
# echo "calculating LD; pid $!"
# pids="$pids $!"

wait $pids
