#!/usr/bin/env python
from os import chdir 

alleleFreq_dir = "/srv/persistent/bliu2/ancestry/AIMS_selection/alleleFreq"
chdir(alleleFreq_dir)

flist=["ACB.chrPos.frq",
"AFR.chrPos.frq",
"AFR-ACB-ASW.chrPos.frq",
"AMR.chrPos.frq",
"ASW.chrPos.frq",
"BEB.chrPos.frq",
"CDX.chrPos.frq",
"CEU.chrPos.frq",
"CHB.chrPos.frq",
"CHS.chrPos.frq",
"CLM.chrPos.frq",
"EAS.chrPos.frq",
"ESN.chrPos.frq",
"EUR.chrPos.frq",
"FIN.chrPos.frq",
"GBR.chrPos.frq",
"GIH.chrPos.frq",
"GWD.chrPos.frq",
"IBS.chrPos.frq",
"ITU.chrPos.frq",
"JPT.chrPos.frq",
"KHV.chrPos.frq",
"LWK.chrPos.frq",
"MSL.chrPos.frq",
"MXL.chrPos.frq",
"PEL.chrPos.frq",
"PJL.chrPos.frq",
"PUR.chrPos.frq",
"SAS.chrPos.frq",
"STU.chrPos.frq",
"TSI.chrPos.frq",
"YRI.chrPos.frq",
"global_r2_0.2_window_2000k.chrPos.bim"]

##find all SNPs whose freq columns is NA:
NA_snp_list = []
for fname in flist:
	with open(fname,'r') as f:
		for line in f:
			if "NA" in line: 
				NA_snp_list.append(line.strip().split()[1])

##write all non-NA lines to new files: 
for fname in flist:
	out_fname = fname.replace('chrPos','chrPos.NArm')
	with open(fname, 'r') as fin, open(out_fname, 'w') as fout: 
		for line in fin:
			snp = line.strip().split()[1]
			if snp in NA_snp_list:
				continue
			else:
				fout.write(line)
