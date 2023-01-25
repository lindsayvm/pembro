#Libraries
library(data.table)
library(ggplot2)
library(fastDummies)
library(vcfR)
library(stringr)

dir = "/home/l.leek/pembro"
setwd(dir)

source("src/1_feature_generation/WGS/0_functions_WGS.R")

#config
path_clinical_data = "data/20221021_DRUP_pembro_LL_final.tsv"
path_cosmicSigs = "data/signatures/cosmicsigs-data.tsv" #UV signatures by Thomas Battaglia
path_purple_patient_fn = "data/pembro_wgs_patientSelection_purplepurity.tsv"
path_driver_patient_fn = "data/pembro_wgs_patientSelection_drivercatalog.tsv"
path_clonalTML_fn = "data/pembro_wgs_patientSelection_ann_vcf.txt"
GOI = c("PIK3CA", "PTEN", "KRAS", "BRAF", "NRAS", "ARID1A", 
             "SMARCB1", "SMARC4A", "PBRM1", "JAK1", "JAK2", "STAT1", 
             "IFNGR1", "IFNGR2","B2M", "TAP1", "TAP2", "TAPB", 
             "CALR", "PDIA3", "CNAX", "HSPA5", "KEAP", "STK11") 
ann.fn = list.files(path = "data/snpeff_output",
                    pattern = "_ann_filt_oneLine.vcf",
                    full.names = TRUE,  
                    recursive = TRUE)

#download data
cosmicSigs.df = fread(path_cosmicSigs, data.table = F)

#read paths to data
pembro_wgs_purple_fn.df = fread(path_purple_patient_fn, data.table = F)
pembro_wgs_driver_fn.df = fread(path_driver_patient_fn, data.table = F)
pembro_wgs_clonalTML_fn.df = fread(path_clonalTML_fn, data.table = F)


#Extract WGS data (TML etc) from PURPLE files
purple.df = my_dataframe(fn = pembro_wgs_purple_fn.df, extension = "purple.purity.tsv")
purple.df = process_purplefiles(df = purple.df)

# #Extract information about TML, clonalTML from SNPeff SNPsift tool
# clonal_and_SNPeffTMLTMB.df = my_clonal_dataframe(fn = pembro_wgs_clonalTML_fn.df)

#Extract Cosmic mutational signatures
#https://www.researchgate.net/figure/Mutational-signatures-of-tumor-samples-Cosmic-mutational-signatures-of-tumor-samples_fig1_336708590
cosmicSigs.df = my_cosmicsignatures(df = cosmicSigs.df) %>% 
                    dplyr::select(matches("sample_id|SBS6|SBS6|SBS15|SBS26|SBS_7|SBS_10")) %>% 
                    dplyr::filter(sample_id %in% purple.df$patientID) # there is info on samples that are not sequenced
#clin.df = fread(path_clinical_data, data.table = F)
####drup_cosmicSigs.df = cosmicSigs.df[cosmicSigs.df$sample_id %in% c(clin.df$CPCT_WIDE_CORE, clin.df$HMFsampleID), ]

#get driver mutations for genes of interest per patient
driver.df = my_dataframe(fn = pembro_wgs_driver_fn.df, extension = ".driver.catalog.tsv") %>% 
            filter(gene %in% GOI)
driver_pp.df = my_GOI_driverMut_pp(df = driver.df, 
                                   driverLikelihood_thres = 0.8,
                                   GOI = GOI,
                                   patientIDs = purple.df$patientID)

#Get snpeff annotated genes
# ann.df = my_ann_snpeff(ann.fn, GOI)
# ann_pp.df = my_ann_pp(ann.df, GOI, ann.fn)



#Merge everything
#merge.df = full_join(x = purple.df, y = clonal_and_SNPeffTMLTMB.df, by = "patientID")
#merge.df = full_join(x = merge.df, y = cosmicSigs.df, by = c("patientID" = "sample_id"))
merge.df = full_join(x = purple.df, y = cosmicSigs.df, by = c("patientID" = "sample_id"))
merge.df = full_join(x = merge.df, y = driver_pp.df, by = "patientID")
#merge.df = full_join(x = merge.df, y = ann_pp.df, by =  "patientID")


write.table(merge.df, file = ("data/pembro_wgs_features.csv"),
            quote=F, sep='\t', col.names = T, row.names = F)

