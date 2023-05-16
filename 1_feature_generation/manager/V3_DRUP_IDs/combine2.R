#libraries
library(data.table)
library(openxlsx)
library(dplyr)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)

#Deprecated
#clin.df = fread('data/20221021_DRUP_pembro_LL_final.tsv', data.table = F) #20221021_DRUP_pembro_LL_final.tsv
#Load metadata from HMF (no DRUP meta...)
#meta.df = fread("/DATA/share/Voesties/data/harmonize/output/metadata.csv") %>%
#  mutate(sampleId = gsub("T$|TI.*","",sampleId))

#load clin data from DRUP team
clin.df = fread("data/20230127_DRUP_pembro_LL_final.tsv", data.table = F)

drup_meta.df =fread("/DATA/share/Voesties/data/DRUP/update_6/metadata.tsv", data.table = F) 
cpct_meta.df = fread("/DATA/share/Voesties/data/HMF/update_10/metadata.tsv", data.table = F)
meta.df = bind_rows(drup_meta.df, cpct_meta.df) %>% 
  mutate(patientID = gsub("T$|TI|TI.*","",sampleId))
#for some reason cannot filter out duplicate IDs, but there are only 2 not overlapping, there are still 14 missing that are in DRUP but not in the database

#Load generated WGS features
wgs.df = fread('data/pembro_wgs_features.csv', data.table = F) 

#Because only WGS is collected based on clin file, all WGS IDs should be present in clin.df
setequal(wgs.df$patientID, wgs.df$patientID[wgs.df$patientID %in% clin.df$patientID])
#run these "CPCT02010936" "CPCT02180009" "CPCT02460001" "DRUP01010249" "CPCT02080234"

#Load generated RNA features
# rna.df = fread('data/pembro_rna_features.csv', data.table = F)
# setequal(rna.df$patientID, rna.df$patientID[rna.df$patientID %in% clin.df$patientID])

#Load neoantigen preds
neo.df = fread("data/neoantigen/neoantigen_preds.csv") %>% 
   dplyr::rename("neo_TMB" = "TMB") %>% 
   dplyr::rename("neo_svTumorMutationalBurden" = "svTumorMutationalBurden") %>% 
   dplyr::select(-c(V27, V28))

#merge data
merge.df = left_join(clin.df, meta.df, by = c("patientID"))
merge.df = left_join(merge.df, wgs.df, by = c('patientID'))
merge.df = left_join(merge.df, neo.df, by = c("hmfSampleId" = "sample_id"))
#final.df = full_join(merge.df, rna.df, by=c('RNA_ID'='patientID'))
final.df = merge.df


#write data
write.table(x = final.df, file = "data/20230310_DRUP_pembro_LL_WGS_RNA.tsv", 
            quote = F, sep = '\t', col.names = T, row.names = F)

write.xlsx(final.df, "data/20230310_DRUP_pembro_LL_WGS_RNA.xlsx")
