import sys 

# input: arg 1
# output: arg 2 

aims_file = sys.argv[1]
gff_file = sys.argv[2]
# aims_file = "/srv/gs1/projects/montgomery/bliu2/ancestry/data/autosomes/2popaims_wlrld2M_150.aims.snpinfo"
with open(aims_file, 'r') as f: 
	aims_lines = f.readlines()

with open(gff_file, 'w') as f: 
	for line in aims_lines:
		if 'snp' in line: continue # skip header line
		split_line = line.strip().split('\t')
		snp = split_line[0]
		chr = split_line[1]
		pos = int(split_line[2])
		gff_line = "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n"%(chr,'aims',snp,pos-1,pos,'.','+','.',' ')
		f.write(gff_line)

