#libraries
library(data.table)
library(openxlsx)
library(dplyr)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)

#load clin data from DRUP team
clin.df = fread('data/20230503_DRUP_pembro_LL_final.tsv', data.table = F)

#Load generated WGS features
wgs.df = fread('data/pembro_wgs_features.csv', data.table = F) 

#Load metadata from harmonized
meta.df = fread("/DATA/share/Voesties/data/harmonize/output/metadata.csv") 
colnames(meta.df) = paste0("hmf_",colnames(meta.df))

#Load generated RNA features
# rna.df = fread('data/pembro_rna_features.csv', data.table = F)
# setequal(rna.df$patientID, rna.df$patientID[rna.df$patientID %in% clin.df$patientID])

# #Load neoantigen preds
# neo.df = fread("data/neoantigen/neoantigen_preds.csv") %>% 
#   dplyr::rename("neo_TMB" = "TMB") %>% 
#   dplyr::rename("neo_svTumorMutationalBurden" = "svTumorMutationalBurden") %>% 
#   dplyr::select(-c(V27, V28,cohort))

#aneuploidy scores
aneuploidy.df = fread("/DATA/share/Voesties/data/harmonize/output/aneuploidy-score.csv") %>% 
  mutate(patientID = gsub("T$|TI*$","", sampleId))


maftools.rds = readRDS("/DATA/share/Voesties/data/harmonize/output/hmf-maftools.rda")
maftools.df = maftools.rds@data

#GIE: clonal mutations as those with a subclonal score lower than 0.85.
cTML_pp.df = maftools.df %>% 
  dplyr::filter(SUBCL < 0.85) %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::count() %>% 
  mutate(patientID = gsub("T$|TI*$","",Tumor_Sample_Barcode)) %>% 
  mutate(cTML_maftools = n) %>% 
  dplyr::select(-c(n))

cTML_pp_05.df = maftools.df %>% 
  dplyr::filter(SUBCL < 0.05) %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::count() %>% 
  mutate(patientID = gsub("T$|TI*$","",Tumor_Sample_Barcode)) %>% 
  mutate(cTML_maftools_05 = n) %>% 
  dplyr::select(-c(n))

#merge data
merge.df = left_join(clin.df, meta.df, by = c("patientID" = "hmf_patientId"))
merge.df = left_join(merge.df, wgs.df, by = c('patientID'))
merge.df = left_join(merge.df, aneuploidy.df, by = c('patientID'))
merge.df = left_join(merge.df, cTML_pp.df, by = c('patientID'))
merge.df = left_join(merge.df, cTML_pp_05.df, by = c('patientID'))
#merge.df = left_join(merge.df, neo.df, by = c("hmfSampleId" = "sample_id"))

#final.df = full_join(merge.df, rna.df, by=c('RNA_ID'='patientID'))
final.df = merge.df

#TRUE TML, breast seperated cohort
final.df$Cohort_trueHML = final.df$Cohort
final.df$Cohort_trueHML[final.df$TML>290] = "HML > 290"
final.df$Cohort_trueHML[final.df$TML<=290 & final.df$TML >=140 & final.df$TumorType == "Breast"] = "HML 140-290, breast"
final.df$Cohort_trueHML[final.df$TML<=290 & final.df$TML >=140 & final.df$TumorType != "Breast"] = "HML 140-290, other"


colnames(final.df)



tmp=final.df %>% dplyr::select(subjectkey, Cohort, Cohort_trueHML, TML, TumorType, responders, BOR)
# DRUP01010268
# HML 140-290, breast
# HML > 290
# 355
# 
# DRUP01080039
# HML 140-290, breast
# HML > 290
# 344
# 
# DRUP01390001
# HML 140-290, breast
# HML > 290
# 364

#write data
write.table(x = final.df, file = "data/20230503_DRUP_pembro_LL_WGS_RNA.tsv", 
            quote = F, sep = '\t', col.names = T, row.names = F)

write.xlsx(final.df, "data/20230503_DRUP_pembro_LL_WGS_RNA.xlsx")
