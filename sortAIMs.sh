#!bin/bash

AIMS_DIR="/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs_noSubpop/"
cd $AIMS_DIR

folders=(AFR.EAS \
AFR.EUR \
AFR.SAS \
EAS.EUR \
EAS.SAS \
EUR.SAS)


for folder in ${folders[*]}; do 
	aims=$(ls $folder/*_500k_500.aims)
	head -1 $aims > ${aims/aims/sorted.aims}
	tail -n +2 $aims | sort -k 7 -r -n  >> ${aims/aims/sorted.aims}
done 


AIMS_DIR="/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs/"
cd $AIMS_DIR

folders=(AFR.AMR \
AFR.EAS \
AFR.EUR \
AFR.SAS \
AMR.EAS \
AMR.EUR \
AMR.SAS \
EAS.EUR \
EAS.SAS \
EUR.SAS)


for folder in ${folders[*]}; do 
	aims=$(ls $folder/*_500k_500.aims)
	head -1 $aims > ${aims/aims/sorted.aims}
	tail -n +2 $aims | sort -k 7 -r -n  >> ${aims/aims/sorted.aims}
done 
