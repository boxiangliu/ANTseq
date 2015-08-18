#!/usr/bin/env python 

##example: cat primerPool_1.txt | python add_adaptor.py > primerPool_1.adaptor.txt

import sys

fwd_adaptor="GCGTTATCGAGGTC";
rev_adaptor="GTGCTCTTCCGATCT"; # bosh 08/07/14 added base T at the end 

for line in sys.stdin:
	split_line = line.strip().split()
	if split_line[0] == "snp": 
		sys.stdout.write("%s\t%s\t%s\t%s\t%s\n"%(line.strip(),'fwd_id', 'fwd+adaptor', 'rev_id','rev+adaptor'))
	else: 
		f_primer = fwd_adaptor + split_line[1]
		r_primer = rev_adaptor + split_line[3]
		f = split_line[0] + "_F"
		r = split_line[0] + "_R"
		sys.stdout.write("%s\t%s\t%s\t%s\t%s\n"%(line.strip(), f, f_primer, r, r_primer))


