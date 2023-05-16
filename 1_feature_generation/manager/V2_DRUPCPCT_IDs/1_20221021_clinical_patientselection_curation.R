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

#CUPs are probably not melanomas
clin.df$BOR[clin.df$TumorType == "CUP"]

clin.df$HMFsampleID[clin.df$HMFsampleID == "No blood"] = NA

table(clin.df$Cohort)

#When pretreatment biopsy is 0, then the biopsy was still pretreatment but was taken for another study (WIDE/CPCT/COREL). Possibly with another treatment right after biopsy.
table(clin.df$PretreatmentBiopsy)

clin.df$responders = ifelse(clin.df$BOR == "PD", "NR", "R")

clin.df$CPCT_WIDE_CORE[clin.df$CPCT_WIDE_CORE == "CPCT010712"] = "CPCT02010712"

#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################
#' 

write.table(clin.df, file='data/20221021_DRUP_pembro_LL_curated.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)

