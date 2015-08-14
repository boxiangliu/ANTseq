#!/bin/bash 

##Filter out variants with global allele frequency less than 0.05. 
##These variants are unlikely to be informative for ancestry inference
##so we filter them out to speed up downstream computation.

wd="/srv/persistent/bliu2/shared/1000genomes/phase3v5"
cd $wd

input_list=(ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr2.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr3.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr5.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr6.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr7.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr8.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr9.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr10.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr11.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr13.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr15.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr16.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr17.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr18.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr19.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr21.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz)

##filter out SNPs with MAF < 0.05:
pids=""
for input in ${input_list[*]}; do
	output=${input/vcf/maf05.vcf}
	(vcftools --gzvcf $input --recode --recode-INFO-all --maf 0.05 --stdout | gzip -c > $output) &
	pids="$pids $!"
done
wait $pids
echo 'VCF MAF filtering finished!'

##concatenate all filtered VCFs:
vcf-concat ${input_list[*]/vcf/maf05.vcf} | gzip -c > ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.vcf.gz
echo "Finished concatenating VCFs!"

##remove intermediate product: 
# rm ${input_list[*]/vcf/maf05.vcf}
