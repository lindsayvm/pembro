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
clin.df = read_excel("data/raw/20230127_ExportPembroCohorten.xlsx") %>% 
  dplyr::select(-c("...27","...28")) %>% 
  dplyr::rename(TumorType = location,
                BOR = bor) %>% 
  mutate(responders = ifelse(BOR == "PD", "NR", "R")) %>% 
  mutate(cohort_manuscript = ifelse(groep %in% c("HML 140-290, other","HML 140-290, breast"), "Pembro 140-290", "Pembro HML>290")) %>% 
  mutate(subjectkey = gsub("-","",subjectkey))

##WGS_ID was curated by Birgit and Laurien & Lindsay
clin.df$WGS_ID[clin.df$WGS_ID == "DRUP01080032"] = "DRUP01080022" #TYPO, staat namelijk bij DRUP ID 08-0022, dus moet DRUP01080022 zijn 
clin.df$WGS_ID[clin.df$WGS_ID == "DRUP01080016"] = "DRUP01080023" #TYPO, staat namelijk bij DRUP ID 08-0023, dus moet DRUP01080023 zijn à is deze beschikbaar in de database? Het WGS rapport heeft namelijk een laag TCP.
clin.df$WGS_ID[clin.df$WGS_ID == "DRUP01180002"] #Als het goed is het DRUP rapport gefinished met voldoende TCP, dus zou beschikbaar moeten zijn à dit zou nog kunnen liggen aan het informed consent, want patient heeft geen toestemming gegeven om de informatie op te slaan in de HMF databank. Als deze dus niet beschikbaar is, moeten we denk ik de datamanager van HMF vragen of ze de data vrij willen geven. Je zou dan wel verwachten dat er eerder wel data beschikbaar was van deze patient, want de issues rondom het IC zijn redelijk recent, dus in nieuwe releases blijken inderdaad soms van dit zulk samples te missen.
clin.df$WGS_ID[clin.df$WGS_ID == "NA"] = "TSO500"

write.table(clin.df, file='data/20230127_DRUP_pembro_LL_curated.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)


