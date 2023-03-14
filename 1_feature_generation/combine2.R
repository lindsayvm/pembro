#libraries
library(data.table)
library(openxlsx)
library(dplyr)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)

#load clin data from DRUP team
clin.df = fread('data/20221021_DRUP_pembro_LL_final.tsv', data.table = F)
#Load metadata from HMF ~ all IDs from clin have been found in update 10 hmf
meta.df = fread("/DATA/share/Voesties/data/HMF/update_10/metadata.tsv") %>%
  mutate(sampleId = gsub("T$|TI.*","",sampleId))

#Load generated WGS features
wgs.df = fread('data/pembro_wgs_features.csv', data.table = F) 

#Because only WGS is collected based on clin file, all WGS IDs should be present in clin.df
setequal(wgs.df$patientID, wgs.df$patientID[wgs.df$patientID %in% clin.df$patientID])


#Load generated RNA features
# rna.df = fread('data/pembro_rna_features.csv', data.table = F)
# setequal(rna.df$patientID, rna.df$patientID[rna.df$patientID %in% clin.df$patientID])

#Load neoantigen preds
neo.df = fread("data/neoantigen/neoantigen_preds.csv") %>% 
   dplyr::rename("neo_TMB" = "TMB") %>% 
   dplyr::rename("neo_svTumorMutationalBurden" = "svTumorMutationalBurden") %>% 
   dplyr::select(-c(V27, V28))

#merge data
merge.df = left_join(clin.df, meta.df, by = c("patientID" = "sampleId"))
merge.df = full_join(merge.df, wgs.df, by = c('patientID'))
merge.df = left_join(merge.df, neo.df, by = c("hmfSampleId" = "sample_id"))
#final.df = full_join(merge.df, rna.df, by=c('RNA_ID'='patientID'))
final.df = merge.df

#prevent mix up of IDs we got from DRUP coordinators and IDs from HMF. 
final.df = final.df %>% 
  dplyr::rename("bir_CPCT_WIDE_CORE" = "CPCT_WIDE_CORE" ) %>% 
  dplyr::rename("bir_HMFsampleID" = "HMFsampleID" )

#which samples are missing but should be there?
missing_samples.df = final.df %>%  dplyr::select(RNA, DNA, bir_CPCT_WIDE_CORE, bir_HMFsampleID, hmfSampleId)


#write data
write.table(x = final.df, file = "data/20230310_DRUP_pembro_LL_WGS_RNA.tsv", 
            quote = F, sep = '\t', col.names = T, row.names = F)

write.xlsx(final.df, "data/20230310_DRUP_pembro_LL_WGS_RNA.xlsx")
