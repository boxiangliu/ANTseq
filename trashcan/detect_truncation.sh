fastq_dir='/srv/gs1/projects/montgomery/bliu2/ancestry/data/Bosh_ASW_mmPCR/BaseCalls' 
sam_dir="/srv/gs1/projects/montgomery/bliu2/ancestry/data/Bosh_ASW_mmPCR/sam_truncation_test"
index_dir='/srv/gs1/projects/montgomery/bliu2/shared/bowtie2'
log_dir='/srv/gs1/projects/montgomery/bliu2/ancestry/log'

cd $sam_dir
files=(*)
for file in ${files[*]} 
do 
echo $file 
is_truncated=$(samtools view -S $file 2>&1 >/dev/null | grep "truncated")
if  
done 