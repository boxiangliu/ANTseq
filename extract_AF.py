''' extract the allele frequency field (AF) and population specific frequencies (e.g. EAS_AF) from the 1000 genome VCFs'''
import gzip,sys,os
if __name__ == "__main__": 
	wd = sys.argv[1]
	inp = sys.argv[2] 
	out = sys.argv[3]
	os.chdir(wd)

	with gzip.open(inp, 'r') as fin, open(out, 'w') as fout: 
		for line in fin: 
			if line.startswith("##"): continue 
			line = line.strip()
			split_line = line.split()
			if line.startswith("#"):
				chrom_idx = split_line.index("#CHROM")
				pos_idx = split_line.index("POS")
				id_idx = split_line.index("ID")
				ref_idx = split_line.index("REF")
				alt_idx = split_line.index("ALT")
				info_idx = split_line.index("INFO")
				out_line = "#CHROM\tPOS\tID\tREF\tALT\tAF\tEAS_AF\tAMR_AF\tAFR_AF\tEUR_AF\tSAS_AF\n"
			else:
				chrom = split_line[chrom_idx]
				pos = split_line[pos_idx]
				snp_id = split_line[id_idx]
				ref = split_line[ref_idx]
				alt = split_line[alt_idx]
				info = split_line[info_idx]
				split_info = info.split(";")
				info_dict = dict([tuple(field.split("=")) for field in split_info if 'AF' in tuple(field.split("="))[0]])
				out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,snp_id,ref,alt,info_dict['AF'],info_dict['EAS_AF'], info_dict['AMR_AF'], info_dict['AFR_AF'], info_dict['EUR_AF'], info_dict['SAS_AF'])
			fout.write(out_line)