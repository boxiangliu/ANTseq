#!bin/bash 

###########
# functions
###########
add_subpopulations () {
	local pop=$1
	local output_file=$2
	template_line="subpop.frq=/srv/persistent/bliu2/ancestry/AIMS_selection/alleleFreq/subpop.chrPos.NArm.frq"
	subpops=$(eval "echo \${$pop[*]}")
	string1="$pop.subpopulations="
	string2=$(echo ${subpop[*]}| tr " " ",")
	echo "$string1$string2" >> $output_file
	echo "$pop has subpopulations ${subpop[*]}"
	for subpop in ${subpops[*]}; do
		echo $template_line | sed "s/subpop/$subpop/g" >> $output_file
	done
}

AIMs_gen_dir="/srv/persistent/bliu2/ancestry/AIMS_selection/pone.0082224.s018"
output_dir="/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs_noSubpop/"
cd $output_dir

##array of all populations: 
populations=(AFR \
# AMR \
EAS \
EUR \
SAS)

##array of subpopulations for each superpoplation:
AFR=(YRI LWK GWD MSL ESN ASW ACB)
# AMR=(MXL PUR CLM PEL)
EAS=(CHB JPT CHS CDX KHV)
EUR=(CEU TSI FIN GBR IBS)
SAS=(GIH PJL BEB STU ITU)

npops=${#populations[*]}
config_template=$AIMs_gen_dir/aims_properties_cheap_ancestry.template.txt
pids=""
n=0
nproc=2

for i in $(seq 1 $((npops-1))); do
	for j in $(seq $((i+1)) $npops); do 
	n=$(($n+1))

	pop1=${populations[$((i-1))]}
	pop2=${populations[$((j-1))]}
	AIMs_config=$output_dir/aims_properties_cheap_ancestry.${pop1}_${pop2}.txt

	##generate AIMs property configuration file:
	sed -e "s/pop1/$pop1/g" -e "s/pop2/$pop2/g" $config_template > $AIMs_config
	echo "writing AIMs config for $pop1 and $pop2..."

	##add subpopulations: 
	# add_subpopulations $pop1 $AIMs_config
	# add_subpopulations $pop2 $AIMs_config
	
	##run aims_generator: 
	python $AIMs_gen_dir/AIMs_generator.py $AIMs_config &> ${AIMs_config/txt/log} &
	# python $AIMs_gen_dir/AIMs_generator.py $AIMs_config &> ${AIMs_config/txt/log}

	echo "generating AIMs for $pop1 and $pop2; pid: $!"
	pids="$pids $!"

	##if $nproc subprocesses has been spawn, wait for them to finish:
	if [[ $n -ge $nproc ]]; then
		echo "hit $nproc jobs"
		wait $pids 
		n=0
		pids=""
	fi

	done 
done 
wait $pids



##organize files to appropriate directories: 
for i in $(seq 1 $((npops-1))); do
	for j in $(seq $((i+1)) $npops); do 

	pop1=${populations[$((i-1))]}
	pop2=${populations[$((j-1))]}

	[[ ! -d $pop1.$pop2 ]] && mkdir $pop1.$pop2
	mv aims_properties_cheap_ancestry.${pop1}_${pop2}.{log,txt} $pop1.$pop2
	mv ${pop1}_${pop2}_500k_500.aims $pop1.$pop2
	mv ${pop1}_${pop2}_500k_500_${pop1}_${pop2}.aims $pop1.$pop2

	done 
done 