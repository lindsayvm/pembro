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

clin.df$BOR = gsub(" |/AE","",clin.df$BOR)

clin.df$Gender = gsub(" ","",clin.df$Gender)

colnames(clin.df)[colnames(clin.df) == "CPCT.WIDE.CORD.."] = "CPCT_WIDE_CORE"
clin.df$CPCT_WIDE_CORE = gsub("-| ","",clin.df$CPCT_WIDE_CORE)
clin.df$CPCT_WIDE_CORE = gsub("T$","",clin.df$CPCT_WIDE_CORE)

clin.df$BOR[clin.df$TumorType == "CUP"]

#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################
#' 

write.table(clin.df, file='data/20221021_DRUP_pembro_LL_final.tsv', 
  quote=TRUE, sep='\t', col.names = T, row.names = F)

