#libraries
library(data.table)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)


#laod data
clin.df = fread('data/20221021_DRUP_pembro_LL_final.tsv', data.table = F)

wgs.df = fread('data/pembro_wgs_features.csv', data.table = F)
clin.df$WGS_ID = rep(NA, nrow(clin.df))
clin.df$WGS_ID[clin.df$CPCT_WIDE_CORE %in% wgs.df$patientID] = clin.df$CPCT_WIDE_CORE[clin.df$CPCT_WIDE_CORE %in% wgs.df$patientID]
clin.df$WGS_ID[clin.df$HMFsampleID %in% wgs.df$patientID] = clin.df$HMFsampleID[clin.df$HMFsampleID %in% wgs.df$patientID]
#in wgs data we have same IDs as in clinical now
setequal(wgs.df$patientID, wgs.df$patientID[wgs.df$patientID %in% clin.df$WGS_ID])

rna.df = fread('data/pembro_rna_features.csv', data.table = F)
clin.df$RNA_ID = rep(NA, nrow(clin.df))
clin.df$RNA_ID[clin.df$CPCT_WIDE_CORE %in% rna.df$patientID] = clin.df$CPCT_WIDE_CORE[clin.df$CPCT_WIDE_CORE %in% rna.df$patientID]
clin.df$RNA_ID[clin.df$HMFsampleID %in% rna.df$patientID] = clin.df$HMFsampleID[clin.df$HMFsampleID %in% rna.df$patientID]
#in wgs data we have same IDs as in clinical now
setequal(rna.df$patientID, rna.df$patientID[rna.df$patientID %in% clin.df$RNA_ID])


#merge data
merge.df = full_join(clin.df, wgs.df, by=c('WGS_ID'='patientID'))
final.df = full_join(merge.df, rna.df, by=c('RNA_ID'='patientID'))


#which samples are missing but should be there?
missing_samples.df = final.df %>%  dplyr::select(RNA, DNA, WGS_ID, RNA_ID, CPCT_WIDE_CORE, HMFsampleID)



#write data
write.table(x = final.df, file = "data/20221021_DRUP_pembro_LL_WGS_RNA.tsv", 
                      quote = F, sep = '\t', col.names = T, row.names = F)

