clin.df = fread('data/20221021_DRUP_pembro_LL_curated.tsv', data.table = F)

#generate a new column where the correct WGS IDs will be. Start with IDs used in context of DRUP
clin.df$patientID = clin.df$HMFsampleID

#Collect IDs that are in HMF database
###harmonized (DRUP update 3, 5, and 6; CPCT update 10)
hmf.df = fread("/DATA/share/Voesties/data/harmonize/output/metadata-purple.csv",data.table = F)  %>%
  mutate(sampleId = gsub("T$|TII.*","",sampleId))
drup.fn = list.files(path = "/DATA/share/Voesties/data/DRUP/update_6/somatics",
                     pattern = "^DRUP",
                     full.names = F,
                     recursive = F)
drup.fn=gsub("T$|TI$|TI.*$","",drup.fn)
cpct.fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_10/somatics",
                     pattern = "CPCT|WIDE|ACTN",
                     full.names = F,
                     recursive = F)
cpct.fn=gsub("T$|TI$|TI.*$","",cpct.fn)
drupcpct.v = c(drup.fn,cpct.fn)

#No biopsies taken in context of DRPU
clin.df$patientID[!clin.df$patientID %in% hmf.df$sampleId] 
clin.df$patientID[!clin.df$patientID %in% drupcpct.v]

#If it even cannot be found in hmf.df, than make it NA. Maybe that CPCT_WIDE_CORE is available for those patients
clin.df$patientID[!clin.df$patientID %in% drupcpct.v] = NA
clin.df$patientID[is.na(clin.df$patientID)] = clin.df$CPCT_WIDE_CORE[is.na(clin.df$patientID)] 

#These IDs can truely not be found in DRUP or CPCT
clin.df$patientID[!clin.df$patientID %in% hmf.df$sampleId] 
clin.df$WGS_avail = rep("yes", nrow(clin.df))
clin.df$WGS_avail[!clin.df$patientID %in% hmf.df$sampleId] = "no"

#Verder heeft WIDE01010574 helaas geen toestemming gegeven aan HMF om de data op te slaan in de databank, dus dat is de reden dat je deze niet kunt vinden. 
clin.df[clin.df$patientID == "WIDE01010574",] 
#We hebben dus van ongeveer alle patiënten in principe WGS beschikbaar, op ééntje na: namelijk 07-0118, deze pt is geïncludeerd op TSO500 en heeft geen nieuw biopt gehad.
clin.df[clin.df$patientID == "TSO500",] 

#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################
#' 

write.table(clin.df, file='data/20221021_DRUP_pembro_LL_final.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)
