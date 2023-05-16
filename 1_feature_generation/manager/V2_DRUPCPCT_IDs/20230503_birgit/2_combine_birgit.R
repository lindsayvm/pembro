#libraries
library(data.table)
library(openxlsx)
library(dplyr)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)

#load clin data from DRUP team
clin.df = fread('data/20230403_DRUP_pembro_LL_final_birgit.tsv', data.table = F)
clin.df = clin.df %>% dplyr::filter(Cohort ==  "Pembro HML>290" )

#Load metadata from HMF ~ all IDs from clin have been found in update 10 hmf
meta.df = fread("/DATA/share/Voesties/data/HMF/update_10/metadata.tsv") %>%
  mutate(sampleId = gsub("T$|TI.*","",sampleId)) 
colnames(meta.df)[-c(1:4)] = paste0("hmf_",colnames(meta.df)[-c(1:4)])

#Load generated WGS features
wgs.df = fread('data/pembro_wgs_features.csv', data.table = F) 

#Because only WGS is collected based on clin file, all WGS IDs should be present in clin.df
wgs.df$patientID[wgs.df$patientID %in% clin.df$patientID]
final.df = dplyr::left_join(clin.df, wgs.df, by = "patientID")

#write data
write.table(x = final.df, file = "data/20230503_DRUP_pembro_LL_WGS_RNA_birgit.tsv", 
            quote = F, sep = '\t', col.names = T, row.names = F)

