#!/bin/bash
input="/srv/persistent/bliu2/ancestry/AIMS_selection/phase3v5_filtered/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.vcf.gz"
cd "/srv/persistent/bliu2/ancestry/AIMS_selection/phase3v5_filtered/"
plink --vcf $input --double-id --snps-only --make-bed --keep-allele-order --out ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05 &> ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.log

cp ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.bim ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.bim.bak
awk '{$2=$1"_"$4}1' ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.bim.bak > ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.bim
