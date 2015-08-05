#bin/bash 

##directories: 
ancestry="/srv/persistent/bliu2/ancestry"
alleleFreq_dir="$ancestry/AIMS_selection/alleleFreq"; cd $alleleFreq_dir
AIMs_dir="$ancestry/AIMS_selection/AIMs_noSubpop"
panel="/srv/persistent/bliu2/shared/1000genomes/phase3v5/integrated_call_samples_v3.20130502.ALL.panel"
input="/srv/persistent/bliu2/ancestry/AIMS_selection/phase3v5_filtered/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.vcf.gz"
admixture_dir="$ancestry/AIMS_selection/ADMIXTURE_noSubpop"
log_dir="$ancestry/log"

##array of all populations: 
populations=(AFR \
# AMR \
EAS \
EUR \
SAS)

##run ADMIXTURE on pairs of populations:
npops=${#populations[*]}
nmarkers=500
pids=""
for i in $(seq 1 $((npops-1))); do
	for j in $(seq $((i+1)) $npops); do 

		pop1=${populations[$((i-1))]}
		pop2=${populations[$((j-1))]}

		##make plink bed, keeping individuals in pop1 and pop2:
		pop_file="$pop1.$pop2.pop"
		cat $panel | grep -e $pop1 -e $pop2 | awk '{print $1, $1, $2}' > $pop_file
		echo "making plink bed file, keeping $pop1 and $pop2..."
		plink --vcf $input --double-id --snps-only --make-bed --keep-allele-order --keep $pop_file --out $pop1.$pop2 &> /dev/null

		##change SNP_ID from rs_id to chr_pos within bim file:
		awk '{$2=$1"_"$4}1' $pop1.$pop2.bim > $pop1.$pop2.chrPos.bim 

		##calculate ancestry proportion using genome-wide SNPs: 
		(admixture $alleleFreq_dir/$pop1.$pop2.bed 2 &> $log_dir/admixture.$pop1.$pop2.log) &
		echo "running ADMIXTURE for $pop1 and $pop2 in bg...; pid: $!"
		pids="$pids $!"

		##extract SNP ID from .aims file:
		snpid_dir=$AIMs_dir/$pop1.$pop2/snpid
		if [[ ! -d $snpid_dir ]]; then 
			mkdir $snpid_dir
		fi
		awk 'NR!=1{print $1}' $AIMs_dir/$pop1.$pop2/${pop1}_${pop2}_500k_500.aims > $snpid_dir/${pop1}_${pop2}_500k_500.aims.snpid

		##calculate ancestry proportions using variable number of AIMs: 
		for n in $(seq 1 $nmarkers); do
			if [[ $n -lt $nmarkers ]]; then 
				cat $snpid_dir/${pop1}_${pop2}_500k_500.aims.snpid | head -$n > $snpid_dir/${pop1}_${pop2}_500k_$n.aims.snpid
			fi
			plink --bfile $alleleFreq_dir/$pop1.$pop2 --bim $alleleFreq_dir/$pop1.$pop2.chrPos.bim --extract $snpid_dir/${pop1}_${pop2}_500k_$n.aims.snpid --make-bed --out $pop1.$pop2.${n}.AIMs
			admixture $alleleFreq_dir/$pop1.$pop2.${n}.AIMs.bed 2 &> $log_dir/admixture.$pop1.$pop2.$n.AIMs.log
		done 
	done 
done 
wait $pids

##move to P Q files to ADMIXTURE directories: 
for i in $(seq 1 $((npops-1))); do
	for j in $(seq $((i+1)) $npops); do 

		pop1=${populations[$((i-1))]}
		pop2=${populations[$((j-1))]}

		##move P and Q files to ADMIXTURE directory: 
		admixture_wd=$admixture_dir/$pop1.$pop2
		[[ ! -d  $admixture_wd ]] && mkdir $admixture_wd
		mv $pop1.$pop2*2.? $admixture_wd

	done 
done 

##remove intermediate files: 
trashcan=$ancestry/trashcan
[[ ! -d $trashcan ]] && mkdir $trashcan
for i in $(seq 1 $((npops-1))); do
	for j in $(seq $((i+1)) $npops); do 

		pop1=${populations[$((i-1))]}
		pop2=${populations[$((j-1))]}

		mv $pop1.$pop2*.* $trashcan

	done 
done 


##AMR, EUR, AFR comparison: 
##make plink bed, keeping individuals in pop1 and pop2:
pops="AFR.AMR.EUR"
pop_file="$pops.pop"
cat $panel | grep -e AFR -e AMR -e EUR | awk '{print $1, $1, $2}' > $pop_file
echo "making plink bed file, keeping AFR, AMR and EUR..."
plink --vcf $input --double-id --snps-only --make-bed --keep-allele-order --keep $pop_file --out $pops &> /dev/null

##change SNP_ID from rs_id to chr_pos within bim file:
awk '{$2=$1"_"$4}1' $pops.bim > $pops.chrPos.bim 

##calculate ancestry proportion using genome-wide SNPs: 
(admixture $alleleFreq_dir/$pops.bed 3 &> $log_dir/admixture.$pops.log) &
echo "running ADMIXTURE for $pops in bg...; pid: $!"

##extract SNP ID from .aims file:
snpid_dir=$AIMs_dir/$pops/snpid
if [[ ! -d $snpid_dir ]]; then 
	mkdir $snpid_dir
fi

aims_fname="AFR_AMR_EUR_446.Galanter.ordered.txt"
awk 'NR!=1{print $1}' $AIMs_dir/$pops/$aims_fname > $snpid_dir/$aims_fname.snpid

nmarkers=446
##calculate ancestry proportions using variable number of AIMs: 
for n in $(seq 1 $nmarkers); do
	[[ $n -lt $nmarkers ]] && cat $snpid_dir/$aims_fname.snpid | head -$n > $snpid_dir/${aims_fname/$nmarkers/$n}.snpid

	plink --bfile $alleleFreq_dir/$pops --bim $alleleFreq_dir/$pops.bim --extract $snpid_dir/${aims_fname/$nmarkers/$n}.snpid --make-bed --out $pops.${n}.AIMs
	admixture -s 123 $alleleFreq_dir/$pops.${n}.AIMs.bed 3 &> $log_dir/admixture.$pops.$n.AIMs.log
done 

##move to P Q files to ADMIXTURE directories: 
admixture_wd=$admixture_dir/$pops
[[ ! -d  $admixture_wd ]] && mkdir $admixture_wd
mv $pops*3.? $admixture_wd

##remove intermediate files: 
mv $pops*.* $trashcan



##cross validation to determine the best K:
##EAS.SAS
for K in $(seq 1 10); do
	admixture --cv $alleleFreq_dir/EAS.SAS.500.AIMs.bed $K | tee $alleleFreq_dir/log${K}.out
done
##EUR.SAS
for K in $(seq 1 5); do
	admixture --cv $alleleFreq_dir/EUR.SAS.500.AIMs.bed $K | tee $alleleFreq_dir/EUR.SAS.log${K}.out
done