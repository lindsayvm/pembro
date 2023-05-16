dir_DRUPclinical = "/home/l.leek/pembro/"
setwd(dir_DRUPclinical)
clin.df = fread('data/20230503_DRUP_pembro_LL_curated.tsv')

#first start with DRUP IDs
clin.df$patientID = clin.df$subjectkey

#from which samples do we actually have data
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

#31 missing when using DRUP ID
clin.df$patientID[!clin.df$patientID %in% c(drup.fn, cpct.fn)]
#Those missing DRUP ID get CPCT
clin.df$patientID[!clin.df$patientID %in% c(drup.fn, cpct.fn)] = NA
clin.df$patientID[is.na(clin.df$patientID)] = clin.df$CPTC_WIDE_CORE[is.na(clin.df$patientID)]
clin.df$patientID[!clin.df$patientID%in% c(drup.fn, cpct.fn)]
#same as be using WGS ID, conclusion: we need CPCT samples.

clin.df$WGS_available = rep("yes", nrow(clin.df))
clin.df$WGS_available[!clin.df$patientID %in% c(drup.fn, cpct.fn)] = "no"

write.table(clin.df, file='data/20230503_DRUP_pembro_LL_final.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)

