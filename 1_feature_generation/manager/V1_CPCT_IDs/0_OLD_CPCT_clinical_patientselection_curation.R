#Author: Lindsay Leek
#Aim: Pembrolizumab pan-cancer all have high TML, but why are there still so many non-responders?

#libraries
library(data.table)
library(dplyr)

dir_DRUPclinical = "/home/l.leek/pembro/"
setwd(dir_DRUPclinical)

source("src/functions_plots.R")

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

#Load data
#csv has already been curated by Birgit Geurts, Laurien Zeverijn and Lindsay Leek
#(excel has also information about sequence status and failed biopsies)
clin.df = read.csv("data/raw/20221021_DRUP_pembro_LL.csv", header=T, na.strings=c("","NA"), sep = ";")
clin.df = clin.df %>% 
  dplyr::select(CPCT.WIDE.CORD.. , HMFsampleID, Cohort, TumorType, 
                Association, Gender, BOR, MutationalLoad, 
                Start, End, PretreatmentBiopsy, RNA,DNA)  

#' ###########################################################################
#' ###########################################################################
#' Process:  
#' ###########################################################################
#' ###########################################################################
#' 
#Curate
unique(clin.df$Association)

clin.df$Association[clin.df$Association == "BL: ER- PR- HER2-"] = "ER- PR- HER2-"
clin.df$Association[clin.df$Association == "ER- PR+ HER2-  "] = "ER- PR+ HER2-"
clin.df$Association[clin.df$Association == "UV "] = "UV"
clin.df$Association[clin.df$Association == "Hoge tractus digestivus "] = "Hoge tractus digestivus"

#Die AE staat inderdaad voor Adverse Event. Patient heeft behandeling gestaakt ivm ascitis/peritonitis, is gestopt voor <16 weken dus wij boeken die patiënt uit als geen clinical benefit, dus je kan m rekenen onder PD wat mij betreft.
clin.df$BOR = gsub(" |/AE","",clin.df$BOR)

#Ik zag dat er bij 1 ptn die wel een biopt had gehad een 0 stond en bij een andere pt die geen biopt had gehad een 1 stond, dus dat heb ik even voor je aangepast (zie bijlage).
clin.df$PretreatmentBiopsy[clin.df$SubmissionNr == "1690"] = 0
clin.df$PretreatmentBiopsy[clin.df$SubmissionNr == "1111"] = 1

clin.df$Gender = gsub(" ","",clin.df$Gender)

colnames(clin.df)[colnames(clin.df) == "CPCT.WIDE.CORD.."] = "CPCT_WIDE_CORE"
clin.df$CPCT_WIDE_CORE = gsub("-| ","",clin.df$CPCT_WIDE_CORE)
clin.df$CPCT_WIDE_CORE = gsub("T$","",clin.df$CPCT_WIDE_CORE)

clin.df$BOR[clin.df$TumorType == "CUP"]

clin.df$HMFsampleID[clin.df$HMFsampleID == "No blood"] = NA

#"Deze pt is geïncludeerd op basis van TSO500: Totaal TMB is 166 mut/Mb, non synonymous TMB is 121.6 mut/Mb. Op basis van een in silico analyse van HMF leek het ons heel waarschijnlijk dat deze pt een HML van >290 had."
clin.df$MutationalLoad[clin.df$CPCT_WIDE_CORE == "TSO500" & clin.df$MutationalLoad == 166] = NA# "assumed > 290"

clin.df$Cohort[clin.df$CPCT_WIDE_CORE == "WIDE01010534" & clin.df$MutationalLoad == 274] = "Pembro HML>290" # "274 but >290 according to SNPeff and HMF"
clin.df$Cohort[clin.df$HMFsampleID == "DRUP01030053" & clin.df$MutationalLoad == 143] = NA # "143 but low (30) according to SNPeff and HMF"

clin.df$MutationalLoad[clin.df$CPCT_WIDE_CORE == "WIDE01010534" & clin.df$MutationalLoad == 274] = NA# "274 but >290 according to SNPeff and HMF"
clin.df$MutationalLoad[clin.df$HMFsampleID == "DRUP01030053" & clin.df$MutationalLoad == 143] = NA# "143 but low (30) according to SNPeff and HMF"

#When pretreatment biopsy is 0, then the biopsy was still pretreatment but was taken for another study (WIDE/CPCT/COREL). Possibly with another treatment right after biopsy.
table(clin.df$PretreatmentBiopsy)


clin.df = clin.df %>%
  mutate(cohort_manuscript = ifelse(Cohort %in% c("Pembro HML 140-290","Pembro mamma 140-290"), "Pembro 140-290", "Pembro HML>290"))

clin.df$responders = ifelse(clin.df$BOR == "PD", "NR", "R")

#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################
#' 

hmf.df = fread("/DATA/share/Voesties/data/harmonize/output/drivers-data.csv",data.table = F)  %>%
  mutate(sampleId = gsub("T$|TII.*","",sampleId))

#generate a new column where the correct WGS IDs will be
clin.df$patientID = clin.df$CPCT_WIDE_CORE
clin.df$patientID[!clin.df$patientID %in% hmf.df$sampleId] = NA
#if CPCT ID was not found, use the DRUP
clin.df$patientID[is.na(clin.df$patientID)] = clin.df$CPCT_WIDE_CORE[is.na(clin.df$patientID)] 

write.table(clin.df, file='data/20221021_DRUP_pembro_LL_final.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)


