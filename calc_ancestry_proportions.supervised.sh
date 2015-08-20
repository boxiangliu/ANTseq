######
#const
######
##directories:
ancestry="/srv/persistent/bliu2/ancestry"
AIMs_dir="$ancestry/AIMS_selection/AIMs_noSubpop"
panel="/srv/persistent/bliu2/shared/1000genomes/phase3v5/integrated_call_samples_v3.20130502.ALL.panel"
input="/srv/persistent/bliu2/ancestry/AIMS_selection/phase3v5_filtered/ALL.autosome.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.maf05.vcf.gz"
admixture_dir="$ancestry/AIMS_selection/ADMIXTURE_supervised"
log_dir="$ancestry/log"

##array of all populations: 
populations=( # AFR \
# AMR \
EAS \
EUR \
SAS)

######
#funct
######
##function to create pop file for leave-one-out CV: 
gen_pop_loo () {
	local panel=$1
	local pop1=$2
	local pop2=$3
	local pop_file=$4
	local leaveOut=$5
	cat $panel | grep -e $pop1 -e $pop2 | awk '{print $3}' > $pop_file

	##replace test individual with "-": 
	awk -v i=$leaveOut 'NR==i{$1="-"}1' $pop_file > $pop_file.loo; mv $pop_file.loo $pop_file
}


##generate prune.in file 
# plink --vcf $input --double-id --snps-only --indep-pairwise 50 10 0.1 --out $admixture_dir/plink

##run ADMIXTURE on pairs of populations:
npops=${#populations[*]}
nmarkers=500
nancestor=2
pids=""
for i in $(seq 1 $((npops-1))); do
	for j in $(seq $((i+1)) $npops); do 
		echo "running supervised ADMIXTURE for $pop1 and $pop2..."

		pop1=${populations[$((i-1))]}
		pop2=${populations[$((j-1))]}

		#change working directory:
		cwd=$admixture_dir/$pop1.$pop2
		[[ ! -d $cwd ]] && mkdir $cwd; cd $cwd 

		##make plink bed, keeping individuals in pop1 and pop2:
		pop_file="$pop1.$pop2.all.pop"
		cat $panel | grep -e $pop1 -e $pop2 | awk '{print $1, $1, $2, $3}' > $pop_file
		echo "making plink bed file, keeping $pop1 and $pop2..."
		plink --vcf $input --double-id --snps-only --make-bed --keep-allele-order --keep $pop_file --extract $admixture_dir/plink.prune.in --out $pop1.$pop2 # &> /dev/null

		##change SNP_ID from rs_id to chr_pos within bim file:
		awk '{$2=$1"_"$4}1' $pop1.$pop2.bim > $pop1.$pop2.chrPos.bim 

		##extract SNP ID from .aims file:
		snpid_dir=$AIMs_dir/$pop1.$pop2/snpid
		if [[ ! -d $snpid_dir ]]; then 
			mkdir $snpid_dir
		fi
		n=500 ##number of AIMs:
		awk 'NR!=1{print $1}' $AIMs_dir/$pop1.$pop2/${pop1}_${pop2}_500k_500.aims > $snpid_dir/${pop1}_${pop2}_500k_500.aims.snpid
		plink --bfile $pop1.$pop2 --bim $pop1.$pop2.chrPos.bim --extract $snpid_dir/${pop1}_${pop2}_500k_$n.aims.snpid --make-bed --out $pop1.$pop2.$n.AIMs

		numInd=$(wc -l $pop_file | awk '{print $1}')
		for k in $(seq 1 $numInd); do
			echo "running supervised ADMIXTURE for individual $k..."

			##generate pop file: 
			pop_file_loo="$pop1.$pop2.pop"
			gen_pop_loo $panel $pop1 $pop2 $pop_file_loo $k
			cp $pop_file_loo ${pop_file_loo/pop/$n.AIMs.pop}

			##calculate ancestry proportion using genome-wide SNPs: 
			admixture -j8 --supervised $pop1.$pop2.bed $nancestor &> $log_dir/admixture.$pop1.$pop2.log
			mv $pop1.$pop2.$nancestor.P $pop1.$pop2.loo.$k.$nancestor.P
			mv $pop1.$pop2.$nancestor.Q $pop1.$pop2.loo.$k.$nancestor.Q

			##calculate ancestry proportions using 500 AIMs: 
			admixture --supervised $pop1.$pop2.$n.AIMs.bed $nancestor &> $log_dir/admixture.$pop1.$pop2.$n.AIMs.log
			mv $pop1.$pop2.$n.AIMs.$nancestor.P $pop1.$pop2.$n.AIMs.loo.$k.$nancestor.P
			mv $pop1.$pop2.$n.AIMs.$nancestor.Q $pop1.$pop2.$n.AIMs.loo.$k.$nancestor.Q

			##calculate ancestry proportions using variable number of AIMs: 
			# for n in $(seq 500 $nmarkers); do
			# 	if [[ $n -lt $nmarkers ]]; then 
			# 		cat $snpid_dir/${pop1}_${pop2}_500k_500.aims.snpid | head -$n > $snpid_dir/${pop1}_${pop2}_500k_$n.aims.snpid
			# 	fi
			# 	plink --bfile $pop1.$pop2 --bim $pop1.$pop2.chrPos.bim --extract $snpid_dir/${pop1}_${pop2}_500k_$n.aims.snpid --make-bed --out $pop1.$pop2.$n.AIMs
			# 	cp $pop_file_loo ${pop_file_loo/pop/$n.AIMs.pop}
			# 	admixture --supervised $pop1.$pop2.$n.AIMs.bed $nancestor &> $log_dir/admixture.$pop1.$pop2.$n.AIMs.log
			# 	mv $pop1.$pop2.$n.AIMs.$nancestor.P $pop1.$pop2.$n.AIMs.loo.$k.$nancestor.P
			# 	mv $pop1.$pop2.$n.AIMs.$nancestor.Q $pop1.$pop2.$n.AIMs.loo.$k.$nancestor.Q
			# done 
		done 
	done 
done 