#!/usr/bin/env R
aims_fname = "/srv/persistent/bliu2/ancestry/AIMS_selection/AIMs_noSubpop/AFR.AMR.EUR/AFR_AMR_EUR_446.Galanter.txt"
aims_df = read.table(aims_fname, header = T, sep = '\t')
names(aims_df)[10] = 'LSBL.Fst'
names(aims_df)[11] = 'LSBL.In'
aims_df_ordered = aims_df[order(aims_df$LSBL.In, decreasing = T),]
write.table(aims_df_ordered, file = sub("Galanter", "Galanter.ordered", aims_fname), quote = F, sep = '\t', row.names = F, col.names = T)
