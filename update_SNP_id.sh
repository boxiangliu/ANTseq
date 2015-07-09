#!bin/bash
wd=/srv/persistent/bliu2/ancestry/AIMS_selection/alleleFreq
cd $wd
frq_list=(ACB.frq \
ASW.frq \
BEB.frq \
CDX.frq \
CEU.frq \
CHB.frq \
CHS.frq \
CLM.frq \
ESN.frq \
FIN.frq \
GBR.frq \
GIH.frq \
GWD.frq \
IBS.frq \
ITU.frq \
JPT.frq \
KHV.frq \
LWK.frq \
MSL.frq \
MXL.frq \
PEL.frq \
PJL.frq \
PUR.frq \
STU.frq \
TSI.frq \
YRI.frq)


##change the SNP column (2nd column) from rs id to chr_pos:
for frq in ${frq_list[*]}; do

	frq_chrPos=${frq/frq/chrPos.frq}
	bim=${frq/frq/bim}
	bim_chrPos=${bim/bim/chrPos.bim}

	echo "begin processing $(basename $frq_list .frq)..."
	awk 'NR==FNR{a[FNR]=$1"_"$4;next}FNR>1{$2=a[FNR-1]}1' $bim $frq > $frq_chrPos
	echo "generated $frq_chrPos with chr_pos as second column!"

	awk '{$2=$1"_"$4}1' $bim > $bim_chrPos
	echo "generated $bim_chrPos with chr_pos as second column!"

done 

# global_bim=global.bim
# global_bim_chrPos=${global_bim/bim/chrPos.bim}
# global_ld=global.ld 
# global_ld_chrPos=${global_ld/ld/chrPos.ld}

# echo "begin processing $(basename global_bim .bim)..."
# awk 'BEGIN{OFS="\t"}{$2=$1"_"$4}1' $global_bim > $global_bim_chrPos
# echo "generated $global_bim_chrPos with chr_pos as second column!"

# echo "begin processing $(basename global_ld .ld)..."
# awk 'BEGIN{OFS="\t"}NR>1{$3=$1"_"$2;$7=$5"_"$6}1' $global_ld > $global_ld_chrPos
# echo "generated $global_ld_chrPos with chr_pos as second column!"