library(cowplot)

args = commandArgs(T)

# read list of AIMs:
snpinfo_filename = '../wetlab/autosomes/2popaims_wlrld2M_150.aims.snpinfo'
snpinfo = read.table(snpinfo_filename, header = T, as.is=T)

# read list of primers:
primer_filename = '../wetlab/autosomes/iteration_method_input_10_primers/primer_pools/primer_pool.all'
primers = read.table(primer_filename, as.is = T)

# subset to AIMs with primers:
AIM_with_primers = snpinfo[snpinfo$snp %in% primers$V1, ]
AIM_with_primers$chr = paste0('chr',AIM_with_primers$chr)
AIM_with_primers$ID = paste(AIM_with_primers$chr, AIM_with_primers$pos, sep = '_')


# read mpileup counts: 
inp = args[1]
# inp = 'plot_reads_per_marker_and_sample.sample_list.txt'
filenames = scan(inp, character())

counts = data.frame()
for (filename in filenames){
	ncol = max(count.fields(filename, sep = '\t'))
	col.classes = rep('NULL', ncol) 
	col.classes[c(1,2,5,9,13)] = c('character',rep('integer',4))
	temp = read.table(filename, colClasses = col.classes, fill = T, col.names = 1:ncol)
	colnames(temp) = c('chrom','pos','depth','ref','alt')
	temp$sample = str_split_fixed(basename(filename), "_", 2)[1]
	counts = rbind(counts, temp)
}

# add SNP ID:
counts$ID = paste(counts$chrom, counts$pos, sep = "_")


# subset mpileup counts to AIM with primers: 
counts$primer = ifelse(counts$ID %in% AIM_with_primers$ID, TRUE, FALSE)
counts = counts[counts$primer == TRUE,]

# coerce counts$sample into factor: 
counts$sample = as.factor(counts$sample)

# what percentage of markers have reads larger than 100 (20)? 
mean(counts$depth>100, na.rm=T) # 0.9284703
mean(counts$depth>20, na.rm=T) # 0.973796

# plot reads per sample:
p1 = ggplot(counts, aes(reorder(sample,-depth,FUN=function(x){median(x,na.rm = T)}), depth)) + stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot(outlier.shape = NA) + stat_summary(fun.y=mean,shape=3, col='red',geom='point') + xlab("Individual/sample") + ylab("Reads") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + scale_y_log10(breaks = c(1,10,100,1000))
figure_path = '../AIMS_selection/figures/reads_per_sample.pdf'
save_plot(figure_path, p1)



# find AIMs that failed to amplify: 
temp = with(counts, reorder(ID,depth,FUN=function(x){median(x,na.rm=T)}))
which(is.na(attr(temp,'scores')))
# these AIMs failed to amplify: 
# chr1_1106112  all NAs
# chr14_105953324 all NAs
# chr20_62173925   all NAs
counts[counts$ID %in% levels(temp)[1:3],]
# these AIMs failed to amplify: 
# chr8_145639681   mostly NAs
# chr4_1497319	  mostly NAs
# chr5_100984586   mostly NAs


# plot reads per AIM:
counts_successful = counts[!(counts$ID %in% c('chr8_145639681','chr4_1497319','chr5_100984586','chr1_1106112','chr14_105953324','chr20_62173925')),]
p2 = ggplot(counts_successful, aes(reorder(ID,-depth,FUN=function(x){median(x,na.rm=T)}), depth+1)) + stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot(outlier.shape = NA) + xlab("Ancestry informative marker") + ylab("Reads") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + scale_y_log10(breaks = c(1,10,100,1000))
figure_path2 = '../AIMS_selection/figures/reads_per_marker.pdf'
save_plot(figure_path2,p2,base_width = 16, base_height=8)



# make boxplot with x axis as primers index:  
counts_successful = as.data.table(counts_successful)
counts_successful[,depth := as.numeric(depth)]
counts_successful[,median_count:=median(depth),by='ID']
rank = data.frame(median_count=unique(counts_successful$median_count),rank=rank(-unique(counts_successful$median_count, na.rm=T)))
merged = merge(counts_successful, rank, by='median_count')
p5 = ggplot(merged, aes(x = rank, y = depth, group = as.factor(rank))) + stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot(outlier.shape = NA) + xlab("Ancestry informative marker (ranked)") + ylab("Reads") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + scale_y_log10(breaks = c(1,10,100,1000))
set.seed(1)
counts_subset = counts_successful[counts_successful$ID %in% sample(counts_successful$ID,20),]
p3 = ggplot(counts_subset, aes(reorder(ID,-depth,FUN=function(x){median(x,na.rm=T)}), depth)) + stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot(outlier.shape = NA) + stat_summary(fun.y=mean,shape=3, col='red',geom='point') + xlab("Ancestry informative marker") + ylab("Reads") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
figure_path3 = '../AIMS_selection/figures/reads_per_marker.subset20.pdf'
save_plot(figure_path3,p3)


# make multi-panel plot of reads per sample and AIM: 
p4 = plot_grid(p1, p5, labels = c('A','B'), align = 'h')
figure_path4 = args[2]
# figure_path4 = '../AIMS_selection/figures/reads_per_marker_and_sample.pdf'
save_plot(figure_path4, p4, base_height = 4, base_aspect_ratio = 2)

