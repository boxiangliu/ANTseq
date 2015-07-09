#!bin/bash 


AIMs_gen_dir="/srv/persistent/bliu2/ancestry/AIMS_selection/pone.0082224.s018"
output_dir="/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs/"
##array of all populations: 
## awk '{print $2}' integrated_call_samples_v3.20130502.ALL.panel | sort | uniq | grep -v 'pop' | awk '{print $1" \\"}'
populations=(ACB \
ASW \
BEB \
CDX \
CEU \
CHB \
CHS \
CLM \
ESN \
FIN \
GBR \
GIH \
GWD \
IBS \
ITU \
JPT \
KHV \
LWK \
MSL \
MXL \
PEL \
PJL \
PUR \
STU \
TSI \
YRI)

npops=${#populations[*]}

pids=""
n=0
nproc=5

for i in $(seq 1 $((npops-1))); do
	for j in $(seq $((i+1)) $npops); do 
	n=$(($n+1))

	pop1=${populations[$((i-1))]}
	pop2=${populations[$((j-1))]}

	# generate AIMs property configuration file:
	sed -e "s/pop1/$pop1/g" -e "s/pop2/$pop2/g" $AIMs_gen_dir/aims_properties_cheap_ancestry.template.txt > $output_dir/aims_properties_cheap_ancestry.${pop1}_${pop2}.txt

	##run aims_generator: 
	python $AIMs_gen_dir/AIMs_generator.py $output_dir/aims_properties_cheap_ancestry.${pop1}_${pop2}.txt &> $output_dir/aims_properties_cheap_ancestry.${pop1}_${pop2}.log &
	
	echo "generating AIMs for $pop1 and $pop2; pid: $!"
	pids="$pids $!"

	##if 100 subprocesses has been spawn, wait for them to finish:
	if [[ $n -ge $nproc ]]; then
		echo "hit $nproc jobs"
		wait $pids 
		n=0
		pids=""
	fi

	done 
done 