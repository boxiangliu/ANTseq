#!usr/bin/env R

##setup:
library(data.table)
library(dplyr)

##read aims:
aims_fname = "/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs_noSubpop/AFR.AMR.EUR/AFR_AMR_EUR_446.Galanter.txt"
aims_df = fread(aims_fname, header = T, sep = '\t')
setnames(aims_df, c('SNP rsID','LSBL(Fst)', 'LSBL(In)'), c('SNP.rsID', 'LSBL.Fst', 'LSBL.In'))

##read global bim file:
bim_fname = "/srv/persistent/bliu2/ancestry/AIMS_selection/alleleFreq/global_r2_0.2_window_2000k.bim"
bim_df = fread(bim_fname, header = T, sep = '\t')
setnames(bim_df, c('chr','SNP.rsID','cM','pos','A1','A2'))

##merge bim and aims:
setkey(aims_df, chr, SNP.rsID)
setkey(bim_df, chr, SNP.rsID)
aims_merge = merge(aims_df, bim_df, all.x = TRUE, all.y = FALSE)

##fill in missing values:
rsID2Pos = list(rs12065716 = c(1,116774045), rs10071261 = c(5,1013694), rs2510719 = c(11,127005791), rs2242865 = c(21,17027031))
for (rsID in names(rsID2Pos)) {
	aims_merge$pos[aims_merge$SNP.rsID == rsID] = rsID2Pos[[rsID]][2]
}

##order by LSBL In: 
aims_merge_ordered = aims_merge[order(aims_merge$LSBL.In, decreasing = T),]

##select columns and write: 
output = aims_merge_ordered %>% select(snp = SNP.rsID, chr, position = pos, NAM_AF, EUR_AF, AFR_AF, population, LSBL.Fst, LSBL.In)
write.table(output, file = "/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs_noSubpop/AFR.AMR.EUR/AFR_AMR_EUR_446.Galanter.orderedByIn.aims", quote = F, sep = '\t', row.names = F, col.names = T)
