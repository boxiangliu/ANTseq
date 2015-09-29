#!bin/bash 
##########
#constants
##########
AIMS_GENERATOR_DIR="/srv/persistent/bliu2/ancestry/scripts/AIMs_generator"
OUTPUT_DIR="/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs/five_superpopulations" 
AIMS_DIR="/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs/five_superpopulations/"
cd $OUTPUT_DIR

#####
#main
#####
pids=""

# run aims_generator with heterogeneity filter: 
configuration="$AIMS_DIR/with_heterogeneity_filter/filter_AFR_AMR_EAS_EUR_SAS/aims_properties_five_superpopulations_with_heterogeneity_filter.txt"
python $AIMS_GENERATOR_DIR/AIMs_generator.py $configuration &> ${configuration/txt/log} &
pids="$pids $!"
echo "generating AIMs for $configuration; pid: $!"

# run aims_generator with heterogeneity filter for all but AMR: 
configuration="$AIMS_DIR/with_heterogeneity_filter/filter_AFR_EAS_EUR_SAS/aims_properties_five_superpopulations_with_heterogeneity_filter.txt"
python $AIMS_GENERATOR_DIR/AIMs_generator.py $configuration &> ${configuration/txt/log} &
pids="$pids $!"
echo "generating AIMs for $configuration; pid: $!"


# run aims_generator without heterogeneity filter: 
configuration="$AIMS_DIR/no_heterogeneity_filter/aims_properties_five_superpopulations_no_heterogeneity_filter.txt"
python $AIMS_GENERATOR_DIR/AIMs_generator.py $configuration &> ${configuration/txt/log} &
pids="$pids $!"
echo "generating AIMs for $configuration; pid: $!"

wait $pids

echo "Done!"