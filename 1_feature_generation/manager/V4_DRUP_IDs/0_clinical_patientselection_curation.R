#Author: Lindsay Leek
#Aim: Pembrolizumab pan-cancer all have high TML, but why are there still so many non-responders?

#libraries
library(data.table)
library(dplyr)
library(readxl)
dir_DRUPclinical = "/home/l.leek/pembro/"
setwd(dir_DRUPclinical)

source("src/functions_plots.R")

#Load data
clin.df = read_excel("data/raw/ExportPembroCohorten_20230503_Lindsay.xlsx")
clin.df$bor[clin.df$subjectkey == "DRUP-01-10-0003"] = "PD"
clin.df = clin.df %>% 
  dplyr::rename(TumorType = location,
                BOR = bor,
                clin_sex = sex,
                Cohort = groep) %>% 
  mutate(responders = ifelse(BOR == "PD", "NR", "R")) %>% 
  mutate(subjectkey = gsub("-","",subjectkey))



#"DRUP_ID" "CPTC_WIDE_CORE" was curated by Birgit and Laurien & Lindsay


write.table(clin.df, file='data/20230503_DRUP_pembro_LL_curated.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)


