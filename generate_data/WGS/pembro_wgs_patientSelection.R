#libraries
library(data.table)

dir = "/home/l.leek/pembro/"
setwd(dir)

clin.df = fread("data/20221021_DRUP_pembro_LL_final.tsv", data.table = F)

select_fn <- function(extension){
  # Filenames of available PURPLE files (WGS)
  hmf_raw.fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_10/somatics",
                          pattern = paste0(extension,"$"),
                          full.names = TRUE,
                          recursive = TRUE)
  hmf_raw_names.fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_10/somatics",
                                pattern = paste0(extension,"$"),
                                full.names = F,
                                recursive = TRUE)
  
  hmf_raw_names = gsub("/.*","",hmf_raw_names.fn)
  hmf.id = gsub("T$|TI$|TI.*$","",hmf_raw_names)
  hmf_drup.id = hmf.id[hmf.id %in% clin.df$CPCT_WIDE_CORE]
  
  #Filter WGS data from HMF database for which clinical data is available 
  hmf_drup.fn = hmf_raw.fn[hmf.id %in% hmf_drup.id]
  
  #remove the later biopsy from the same patient
  hmf_drup.fn = hmf_drup.fn[!str_detect(hmf_drup.fn, "CPCT02010359TII")]
  
  
  #Load DRUP data to find missing samples
  drup_raw.fn = list.files(path = "/DATA/share/Voesties/data/DRUP/update_3/somatics",
                           pattern = paste0(extension,"$"),
                           full.names = T,
                           recursive = TRUE)
  drup_raw_names.fn = list.files(path = "/DATA/share/Voesties/data/DRUP/update_3/somatics",
                                 pattern = paste0(extension,"$"),
                                 full.names = F,
                                 recursive = TRUE)
  drup_raw_names = gsub("/.*","",drup_raw_names.fn)
  drup.id = gsub("T$|TI$|TI.*$","",drup_raw_names)
  
  #Missing IDs (CPCT_WIDE_CORE style)
  missing_drup_CPCT_WIDE_CORE.id = clin.df$CPCT_WIDE_CORE[!clin.df$CPCT_WIDE_CORE %in% hmf_drup.id]
  #Find matching DRUP IDs (DRUP style)
  missing.df = clin.df[clin.df$CPCT_WIDE_CORE %in% missing_drup_CPCT_WIDE_CORE.id, ]
  drup_drup.fn = drup_raw.fn[drup.id %in% missing.df$HMFsampleID]
  
  #Merge
  ID = c(hmf_drup.fn, drup_drup.fn)
  names = gsub(".*\\/purple/|.purple.purity.tsv|.*\\/linx/|.linx.driver.catalog.tsv|.purple.somatic.vcf.gz", "", ID)
  if(extension != "somatic.vcf.gz"){
    names = gsub("T$|TI$|TI.*$","",names)
  }
  final.df = data.frame(cbind(names = names, 
                              ID = ID))
  
  return(final.df)
}

purplepurity_fn.df = select_fn(extension = ".purple.purity.tsv")
drivercatalog_fn.df = select_fn(extension = ".driver.catalog.tsv")
somatic_fn.df = select_fn(extension = "somatic.vcf.gz")

#export
write.table(purplepurity_fn.df, file = ("data/pembro_wgs_patientSelection_purplepurity.tsv"), 
            quote=TRUE, sep='\t', col.names = T, row.names = F)

write.table(drivercatalog_fn.df, file = ("data/pembro_wgs_patientSelection_drivercatalog.tsv"), 
            quote=TRUE, sep='\t', col.names = T, row.names = F)


write.table(somatic_fn.df$names, file = ("data/pembro_wgs_patientSelection_somatic.txt"), 
            quote=F, sep='\n', col.names = F, row.names = F)







