#### modules #### 
import pysam 
from pdb import set_trace as t
#### TODO ####

#### main #### 

# file paths:
fasta_file_name = "/srv/gs1/projects/montgomery/bliu2/shared/hg19/hg19.fa.hardmasked"
snp_file_name = "/srv/gs1/projects/montgomery/bliu2/ancestry/raw/autosomes/2popaims_wlrld2M_150.aims.snpinfo"
flank_seq_file_name = "/srv/gs1/projects/montgomery/bliu2/ancestry/data/2popaims_wlrld2M_150.aims.snpinfo.fa.up120.down350.01"

FLANK_UP = 120
FLANK_DOWN = 350
fasta_file = pysam.Fastafile(fasta_file_name)
flank_seq_file = open(flank_seq_file_name,'w')
snp_file = open(snp_file_name,'r').readlines()
seq_count = 0
file_count = 1

# retrieve flanking sequences: 
for snp_line in snp_file[1:]:
	# when a flank_seq_file reaches 10 sequences, close it and open a new one: 
	if seq_count == 10:
		file_count += 1 
		flank_seq_file.close()
		flank_seq_file = open('/srv/gs1/projects/montgomery/bliu2/ancestry/data/2popaims_wlrld2M_150.aims.snpinfo.fa.up120.down350.%02i'%file_count,'w')
		
		# reset seq_count: 
		seq_count = 0
	# parse out chr_num, rs_id, and snp_pos:
	snp_items = snp_line.strip().split()
	rs_id = snp_items[0]
	chr_num = "chr%s"%(snp_items[1])
	snp_pos = int(snp_items[2])
	start = snp_pos - FLANK_UP
	end = snp_pos + FLANK_DOWN

	# fetch flanking sequence:
	seq = fasta_file.fetch(chr_num, start, end)

	# fetch a shorter flanking sequence to append to the tag: 
	tag_seq = fasta_file.fetch(chr_num, snp_pos - 5, snp_pos + 5)
	tag = ">%s:%i-%i:%i:%s:%s"%(chr_num,start,end,snp_pos,rs_id, tag_seq)

	# write to output:
	flank_seq_file.write("%s\n%s\n"%(tag,seq))

	# add 1 to seq_count:
	seq_count += 1
flank_seq_file.close()
fasta_file.close()
