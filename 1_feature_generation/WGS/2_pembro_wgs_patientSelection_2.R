#libraries
library(data.table)

dir = "/home/l.leek/pembro/"
setwd(dir)

#clinical file with WGS IDs that are required
clin.df = fread("data/20221021_DRUP_pembro_LL_final.tsv", data.table = F)

#meta data from HMF with all possible WGS IDs in entire database HMF
hmf.df = fread("/DATA/share/Voesties/data/harmonize/output/drivers-data.csv",data.table = F)  %>%
  mutate(sampleId = gsub("T$|TII.*","",sampleId))

#Show IDs that actually can be found back in HMF database: 61 out of 73 have WGS avail
clin.df$patientID[clin.df$patientID %in% hmf.df$sampleId]


select_fn_somatics <- function(extension){
  # Filenames of available PURPLE files (WGS)
  hmf_raw.fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_10/somatics",
                          pattern = paste0(extension,"$"),
                          full.names = TRUE,
                          recursive = TRUE)
  hmf_raw_names = gsub(".*\\/|\\..*","",hmf_raw.fn)
  hmf.id = gsub("T$|TII.*","",hmf_raw_names)
  #those of HMF we can find in clin
  hmf_pembro.fn = hmf_raw.fn[hmf.id %in% clin.df$patientID]
  
  #Load DRUP data to find missing samples
  drup_raw.fn = list.files(path = "/DATA/share/Voesties/data/DRUP/update_3/somatics",
                           pattern = paste0(extension,"$"),
                           full.names = T,
                           recursive = TRUE)
  drup_raw_names = gsub(".*\\/|\\..*","",drup_raw.fn)
  drup.id = gsub("T$|TII.*","",drup_raw_names)
  #Those of DRUP we can find in clin
  drup_pembro.fn = drup_raw.fn[drup.id %in% clin.df$patientID]
  
  
  ######QUALITY CHECK######
  # hmf.id = hmf.id[hmf.id %in% clin.df$patientID]
  # drup.id = drup.id[drup.id %in% clin.df$patientID]
  # total.id = c(hmf.id, drup.id)
  # #check whether there are no duplicates found from DRUP and HMF dirs
  # identical(total.id,unique(total.id))
  # #patients that cannot be matched to either fn in DRUP in HMF based on patient ID column
  # #Because PatientID prioritizes HMFsampleID >> CPCT_WIDE_CORE, check whether CPCT_WIDE_CORE can be found
  # # but that is not the case.
  # tmp = clin.df[!clin.df$patientID %in% total.id, ]
  # ######QUALITY CHECK######
  
  #Merge
  ID = c(hmf_pembro.fn, drup_pembro.fn)
  names = gsub(".*\\/purple/|.purple.purity.tsv|.*\\/linx/|.linx.driver.catalog.tsv|.purple.somatic.vcf.gz|.purple.cnv.somatic.tsv", "", ID)
  if(extension != "somatic.vcf.gz"){
    names = gsub("T$|TI$|TI.*$","",names)
  }
  final.df = data.frame(cbind(names = names, 
                              ID = ID))
  
  return(final.df)
}

select_fn_snpeff <- function(extension, extension2){
  
  clonal = list.files(path = "/home/l.leek/pembro/data/snpeff_output",
                          pattern = "_ann.vcf",
                          full.names = TRUE,
                          recursive = TRUE)
  
  anno = list.files(path = "/home/l.leek/pembro/data/snpeff_output",
                  pattern = "_ann_filt_oneLine.vcf",
                  full.names = TRUE,
                  recursive = TRUE)
  
  names = gsub(".*\\/|_ann.vcf|_ann_filt_oneLine.vcf", "", clonal)
  final.df = data.frame(cbind(names = names, 
                              clonal = clonal,
                              anno = anno))
  
  final.df = final.df[gsub("T$|TI$|TI.*$","",final.df$names) %in% clin.df$patientID, ]
  
  
  return(final.df)
}



purplepurity_fn.df = select_fn_somatics(extension = ".purple.purity.tsv") #x=fread("/DATA/share/Voesties/data/harmonize/output/purple-data.csv")
drivercatalog_fn.df = select_fn_somatics(extension = ".driver.catalog.tsv") #x=fread("/DATA/share/Voesties/data/harmonize/output/drivers-data.csv")
somatic_fn.df = select_fn_somatics(extension = "somatic.vcf.gz")
cnv_fn.df = select_fn_somatics(extension = "purple.cnv.somatic.tsv")
ann_vcf.df = select_fn_snpeff(extension = "_ann.vcf",
                              extension2 = "_ann_filt_oneLine.vcf")


#export
write.table(purplepurity_fn.df, file = ("data/pembro_wgs_patientSelection_purplepurity.tsv"), 
            quote=TRUE, sep='\t', col.names = T, row.names = F)

write.table(drivercatalog_fn.df, file = ("data/pembro_wgs_patientSelection_drivercatalog.tsv"), 
            quote=TRUE, sep='\t', col.names = T, row.names = F)

write.table(cnv_fn.df, file = ("data/pembro_wgs_patientSelection_cnv.txt"),  
            quote=TRUE, sep='\t', col.names = T, row.names = F)     

write.table(somatic_fn.df$names, file = ("data/pembro_wgs_patientSelection_somatic.txt"),  
            quote=F, sep='\n', col.names = F, row.names = F)     

write.table(ann_vcf.df, file = ("data/pembro_wgs_patientSelection_ann_vcf.txt"),
            quote=F, sep='\t', col.names = T, row.names = F)





