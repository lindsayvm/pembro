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

#metadata
# dir_DRUP = "/DATA/share/Voesties/data/DRUP/update_3"
# dir_HMF = "/DATA/share/Voesties/data/HMF/update_10"

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

clin.df$responders = ifelse(clin.df$BOR == "PD", "NR", "R")

clin.df$CPCT_WIDE_CORE[clin.df$CPCT_WIDE_CORE == "CPCT010712"] = "CPCT02010712"
#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################
#' 

#Collect IDs that are in HMF database
###harmonized
hmf.df = fread("/DATA/share/Voesties/data/harmonize/output/metadata-purple.csv",data.table = F)  %>%
  mutate(sampleId = gsub("T$|TII.*","",sampleId))
###DRUP
drup.df = fread("/DATA/share/Voesties/data/DRUP/update_6/metadata.tsv",data.table = F)  %>%
  mutate(sampleId = gsub("T$|TII.*","",sampleId))
###CPCT
cpct.df = fread("/DATA/share/Voesties/data/HMF/update_10/metadata.tsv",data.table = F)  %>%
  mutate(sampleId = gsub("T$|TII.*","",sampleId))

#generate a new column where the correct WGS IDs will be
clin.df$patientID = clin.df$HMFsampleID

#check which HMF database has most DRUP sample IDs available 
clin.df$patientID[!clin.df$patientID %in% hmf.df$sampleId] #less are missing in the harmonized
#clin.df$patientID[!clin.df$patientID %in% c(drup.df$sampleId,cpct.df$sampleId)] 

#If it even cannot be found in hmf.df, than make it NA. Maybe that CPCT_WIDE_CORE is available for those patients
clin.df$patientID[!clin.df$patientID %in% hmf.df$sampleId] = NA

#If DRUP ID was not found, use the CPCT_WIDE_CORE
clin.df$patientID[is.na(clin.df$patientID)] = clin.df$CPCT_WIDE_CORE[is.na(clin.df$patientID)] 

#These IDs can truely not be found 
clin.df$patientID[!clin.df$patientID %in% hmf.df$sampleId] 


#Verder heeft WIDE01010574 helaas geen toestemming gegeven aan HMF om de data op te slaan in de databank, dus dat is de reden dat je deze niet kunt vinden. 
clin.df[clin.df$patientID == "WIDE01010574",] 
#We hebben dus van ongeveer alle patiënten in principe WGS beschikbaar, op ééntje na: namelijk 07-0118, deze pt is geïncludeerd op TSO500 en heeft geen nieuw biopt gehad.
clin.df[clin.df$patientID == "TSO500",] 

write.table(clin.df, file='data/20230403_DRUP_pembro_LL_final_birgit.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)

