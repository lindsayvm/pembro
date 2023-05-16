#libraries
library(data.table)

dir = "/home/l.leek/pembro/"
setwd(dir)

#clinical file with WGS IDs that are required
clin.df = fread("data/20230503_DRUP_pembro_LL_final.tsv", data.table = F) 
table(clin.df$WGS_avail)


extension = "somatic.vcf.gz"

# Filenames of available PURPLE files (WGS)

drup_raw.fn = list.files(path = "/DATA/share/Voesties/data/DRUP/update_6/somatics",
                         pattern = paste0(extension,"$"),
                         full.names = T,
                         recursive = TRUE)
drup_raw_names = gsub(".*\\/|\\..*","",drup_raw.fn)
drup.id = gsub("T$|TI$|TI.*$","",drup_raw_names)
#Those of DRUP we can find in clin
drup_pembro.fn = drup_raw.fn[drup.id %in% clin.df$patientID]
drup_names = gsub(".*\\/purple/|.purple.somatic.vcf.gz", "", drup_pembro.fn)

write.table(drup_names, file = ("data/pembro_wgs_patientSelection_somatic_DRUP.txt"),  
            quote=F, sep='\n', col.names = F, row.names = F)     
fread("data/pembro_wgs_patientSelection_somatic_DRUP.txt")


hmf_raw.fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_10/somatics",
                        pattern = paste0(extension,"$"),
                        full.names = TRUE,
                        recursive = TRUE)
hmf_raw_names = gsub(".*\\/|\\..*","",hmf_raw.fn)
hmf.id = gsub("T$|TI$|TI.*$","",hmf_raw_names)
#those of HMF we can find in clin
hmf_pembro.fn = hmf_raw.fn[hmf.id %in% clin.df$patientID]
hmf_names = gsub(".*\\/purple/|.purple.somatic.vcf.gz", "", hmf_pembro.fn)



write.table(hmf_names, file = ("data/pembro_wgs_patientSelection_somatic_HMF.txt"),  
            quote=F, sep='\n', col.names = F, row.names = F)     

fread("data/pembro_wgs_patientSelection_somatic_HMF.txt")


head(clin.df)
clin.df$snpeff = NA
clin.df$snpeff[clin.df$patientID %in%  c(gsub("T$|TI$|TI.*$","",hmf_names),
                                         gsub("T$|TI$|TI.*$","",drup_names))] = 1


clin.df$snpeff_HMF = NA
clin.df$snpeff_HMF[clin.df$patientID %in%  gsub("T$|TI$|TI.*$","",hmf_names)] = 1

clin.df$snpeff_DRUP = NA
clin.df$snpeff_DRUP[clin.df$patientID %in%  gsub("T$|TI$|TI.*$","",drup_names)] = 1
View(clin.df)
#No overlap between DRUP and HMF IDs


#check
ann.fn = list.files(path = "data/snpeff_output",
                    pattern = "_ann_filt_oneLine.vcf",
                    full.names = F,  
                    recursive = TRUE)
length(gsub("T.*","",ann.fn))

