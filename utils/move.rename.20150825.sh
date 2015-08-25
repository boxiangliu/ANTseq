#!/bin/bash 

MULTIPLEXPRIMERS_DIR="/srv/persistent/bliu2/ancestry/AIMS_selection/multiplexPrimers"; cd $MULTIPLEXPRIMERS_DIR

pop_contrasts=(AFR.AMR.EUR \
AFR.EAS \
AFR.EUR \
AFR.SAS \
EAS.EUR \
EAS.SAS \
EUR.SAS)

for pop_contrast in ${pop_contrasts[*]}; do

	cd $MULTIPLEXPRIMERS_DIR/$pop_contrast
	[[ ! -d plink ]] && mkdir plink 
	mv $pop_contrast.cumPrimerPool_*.{bed,bim,fam,log,nosex} plink
	mv $pop_contrast.pop plink

	[[ ! -d admixture ]] && mkdir admixture
	mv $pop_contrast.cumPrimerPool_*.?.{P,Q} admixture
	mv $pop_contrast.?.{P,Q} admixture
	mv calc_ancestry_proportions.log admixture


	rm Rplots.pdf

	cp primerPoolMerged.txt ../primerPools/$pop_contrast.primerPoolMerged.txt

done 