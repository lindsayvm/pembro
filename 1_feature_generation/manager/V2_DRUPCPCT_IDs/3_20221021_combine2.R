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
colnames(meta.df)[-c(1:4)] = paste0("hmf_",colnames(meta.df)[-c(1:4)])

#Load generated WGS features
wgs.df = fread('data/pembro_wgs_features.csv', data.table = F) 

#Because only WGS is collected based on clin file, all WGS IDs should be present in clin.df
setequal(wgs.df$patientID, wgs.df$patientID[wgs.df$patientID %in% clin.df$patientID])


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
merge.df = left_join(clin.df, meta.df, by = c("patientID" = "sampleId"))
merge.df = left_join(merge.df, wgs.df, by = c('patientID'))
merge.df = left_join(merge.df, aneuploidy.df, by = c('patientID'))
merge.df = left_join(merge.df, cTML_pp.df, by = c('patientID'))
merge.df = left_join(merge.df, cTML_pp_05.df, by = c('patientID'))
#merge.df = left_join(merge.df, neo.df, by = c("hmfSampleId" = "sample_id"))

#final.df = full_join(merge.df, rna.df, by=c('RNA_ID'='patientID'))
final.df = merge.df

#prevent mix up of IDs we got from DRUP coordinators and IDs from HMF. 
final.df = final.df %>% 
  dplyr::rename("bir_CPCT_WIDE_CORE" = "CPCT_WIDE_CORE" ) %>% 
  dplyr::rename("bir_HMFsampleID" = "HMFsampleID" )

#which samples are missing but should be there?
missing_samples.df = final.df %>%  dplyr::select(RNA, DNA, bir_CPCT_WIDE_CORE, bir_HMFsampleID, hmfSampleId)

#post processing: SNPeff TML has more patient than HMF TML, only 1 patient that is in HMF but not SNPeff
#x=data.frame(final.df$patientID, final.df$TML_SNPeff, final.df$TML, final.df$Cohort)

#### !!!!!!WHEN SAMPLE HAS SNPEFF THEN USE THAT
final.df$TML_SNPeff[final.df$patientID == "CPCT02180009"] = final.df$TML[final.df$patientID == "CPCT02180009"]
#####


#remove patients where we know that TML is not at least 140.
final.df = final.df %>% dplyr::filter(!patientID %in% c("DRUP01330018", "DRUP01030053"))


final.df$Cohort #as in DRUP admission

#TRUE TML, breast seperated cohort
final.df$Cohort_trueHML = final.df$Cohort
final.df$Cohort_trueHML[final.df$TML>290] = "Pembro HML>290"
final.df$Cohort_trueHML[final.df$TML<=290] = "Pembro HML 140-290" 
final.df$Cohort_trueHML[final.df$Cohort =="Pembro mamma 140-290"] = "Pembro mamma 140-290"

#x=data.frame(final.df$patientID, final.df$TML_SNPeff, final.df$TML, final.df$Cohort, final.df$Cohort_trueHML)


#write data
write.table(x = final.df, file = "data/20221021_DRUP_pembro_LL_WGS_RNA.tsv", 
            quote = F, sep = '\t', col.names = T, row.names = F)

write.xlsx(final.df, "data/20221021_DRUP_pembro_LL_WGS_RNA.xlsx")
