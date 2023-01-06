#libraries
library(data.table)
library(openxlsx)
library(dplyr)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)


#load clin data from DRUP team
clin.df = fread('data/20221021_DRUP_pembro_LL_final.tsv', data.table = F)

#Load metadata from HMF
meta.df = fread("/DATA/share/Voesties/data/HMF/update_10/metadata.tsv") %>% 
  mutate(sampleId = gsub("T$|TI.*","",sampleId))

#Load generated WGS features
wgs.df = fread('data/pembro_wgs_features.csv', data.table = F) %>% 
  rename("sex" = "gender")
clin.df$WGS_ID = rep(NA, nrow(clin.df))
clin.df$WGS_ID[clin.df$CPCT_WIDE_CORE %in% wgs.df$patientID] = clin.df$CPCT_WIDE_CORE[clin.df$CPCT_WIDE_CORE %in% wgs.df$patientID]
clin.df$WGS_ID[clin.df$HMFsampleID %in% wgs.df$patientID] = clin.df$HMFsampleID[clin.df$HMFsampleID %in% wgs.df$patientID]
#in wgs data we have same IDs as in clinical now
setequal(wgs.df$patientID, wgs.df$patientID[wgs.df$patientID %in% clin.df$WGS_ID])


#Load generated RNA features
rna.df = fread('data/pembro_rna_features.csv', data.table = F)
clin.df$RNA_ID = rep(NA, nrow(clin.df))
clin.df$RNA_ID[clin.df$CPCT_WIDE_CORE %in% rna.df$patientID] = clin.df$CPCT_WIDE_CORE[clin.df$CPCT_WIDE_CORE %in% rna.df$patientID]
clin.df$RNA_ID[clin.df$HMFsampleID %in% rna.df$patientID] = clin.df$HMFsampleID[clin.df$HMFsampleID %in% rna.df$patientID]
#in wgs data we have same IDs as in clinical now
setequal(rna.df$patientID, rna.df$patientID[rna.df$patientID %in% clin.df$RNA_ID])

#Load neoantigen preds
neo.df = fread("data/neoantigen/neoantigen_preds.csv") %>% 
  rename("neo_TMB" = "TMB") %>% 
  rename("neo_svTumorMutationalBurden" = "svTumorMutationalBurden") %>% 
  dplyr::select(-c(V27, V28))


#merge data
merge.df = left_join(clin.df, meta.df, by = c("WGS_ID" = "sampleId"))
merge.df = full_join(merge.df, wgs.df, by = c('WGS_ID'='patientID'))
merge.df = left_join(merge.df, neo.df, by = c("hmfSampleId" = "sample_id"))
final.df = full_join(merge.df, rna.df, by=c('RNA_ID'='patientID'))

#prevent mix up of IDs we got from DRUP coordinators and IDs from HMF. 
final.df = final.df %>% 
  rename("bir_CPCT_WIDE_CORE" = "CPCT_WIDE_CORE" ) %>% 
  rename("bir_HMFsampleID" = "HMFsampleID" )

#which samples are missing but should be there?
missing_samples.df = final.df %>%  dplyr::select(RNA, DNA, WGS_ID, RNA_ID, bir_CPCT_WIDE_CORE, bir_HMFsampleID, hmfSampleId)


#write data
write.table(x = final.df, file = "data/20221021_DRUP_pembro_LL_WGS_RNA.tsv", 
                      quote = F, sep = '\t', col.names = T, row.names = F)


write.xlsx(final.df, "data/20221021_DRUP_pembro_LL_WGS_RNA.xlsx")