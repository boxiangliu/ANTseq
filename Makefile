# plot reads per marker: 
input1 = plot_reads_per_marker_and_sample.sample_list.txt
figure1 = ../AIMS_selection/figures/reads_per_marker_and_sample.pdf
$(figure1): $(input1) 
	Rscript plot_reads_per_marker_and_sample.R $(input1) $(figure1)

# plot mutual information as function of AF: 
../AIMS_selection/figures/MI_vs_delta_and_pj.pdf: 
	Rscript plot_MI.R ./ ../AIMS_selection/figures/MI_vs_delta_and_pj.pdf

# plot AF distribution: 
../AIMS_selection/figures/AF_distribution.pdf: ../AIMS_selection/allele_freq_unfiltered/sample_list.txt
	Rscript plot_AF.R ../AIMS_selection/allele_freq_unfiltered/ sample_list.txt ../figures/AF_distribution.pdf

# extract AF from 1000 Genomes VCF: 
.PHONY: all
all: ../AIMS_selection/allele_freq_unfiltered/chr1.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr2.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr3.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr4.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr5.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr6.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr7.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr8.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr9.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr10.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr11.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr12.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr13.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr14.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr15.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr16.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr17.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr18.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr19.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr20.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr21.allele_freq.txt ../AIMS_selection/allele_freq_unfiltered/chr22.allele_freq.txt
../AIMS_selection/allele_freq_unfiltered/chr%.allele_freq.txt: ../../shared/1000genomes/phase3v5/ALL.chr%.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
	python extract_AF.py ./ $< $@