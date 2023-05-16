#Libraries
library(data.table)
library(ggplot2)
library(fastDummies)
library(vcfR)
library(stringr)
library(dplyr)

dir = "/home/l.leek/pembro"
setwd(dir)

source("src/1_feature_generation/WGS/0_functions_WGS.R")

#config
path_clinical_data = "data/20221021_DRUP_pembro_LL_final.tsv"
path_cosmicSigs = "/DATA/share/Voesties/data/harmonize/output/mutsigs/MutationalPatterns/mutationalpatterns-snv-strict.csv" #UV signatures by Thomas Battaglia
path_purple_patient_fn = "data/pembro_wgs_patientSelection_purplepurity.tsv"
path_driver_patient_fn = "data/pembro_wgs_patientSelection_drivercatalog.tsv"
path_cnv_patient_fn = "data/pembro_wgs_patientSelection_cnv.txt"
path_som_patient_fn = "data/pembro_wgs_patientSelection_somatic.txt"
path_clonalTML_fn = "data/pembro_wgs_patientSelection_ann_vcf.txt"
GOI = c("ARID1A",  "SMARCB1", "SMARC4A", "PBRM1","ATRX","ARID1B","ARID2",
        "JAK1", "JAK2", "STAT1", "IFNGR1", "IFNGR2",
        "B2M", "TAP1", "TAP2", "TAPB", "CALR", "PDIA3", "CNAX", "HSPA5", "KEAP", "STK11",
        "HLA-A","HLA-B", "HLA-C",
        "PIK3CA", "PTEN", "KRAS", "BRAF", "NRAS", "ESR1", "ESR2",
        "APOBEC3B", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H") 


ann.fn = list.files(path = "data/snpeff_output",
                    pattern = "_ann_filt_oneLine.vcf",
                    full.names = TRUE,  
                    recursive = TRUE)

#read paths to data
pembro_wgs_purple_fn.df = fread(path_purple_patient_fn, data.table = F)
pembro_wgs_driver_fn.df = fread(path_driver_patient_fn, data.table = F)
pembro_wgs_cnv.df = fread(path_cnv_patient_fn, data.table = F)
pembro_wgs_clonalTML_fn.df = fread(path_clonalTML_fn, data.table = F)


#Extract WGS data (TML etc) from PURPLE files #is already per patient
purple_pp.df = my_dataframe(fn = pembro_wgs_purple_fn.df, extension = "purple.purity.tsv")
purple_pp.df = process_purplefiles(df = purple_pp.df)

#Extract WGS data from PURPLE CNV files
cnv.df = my_dataframe(fn = pembro_wgs_cnv.df, extension = ".purple.cnv.somatic.tsv")
cnv_pp.df = my_cnv_pp(cnv.df)


#Extract information about TML, clonalTML from SNPeff SNPsift tool
clonal_and_SNPeffTMLTMB.df = my_clonal_dataframe(fn = pembro_wgs_clonalTML_fn.df)

#Extract Cosmic mutational signatures 
cosmicSigs.df = fread(path_cosmicSigs, data.table = F) #(almost 4000 pts in total database)
#https://www.researchgate.net/figure/Mutational-signatures-of-tumor-samples-Cosmic-mutational-signatures-of-tumor-samples_fig1_336708590
colnames(cosmicSigs.df) = gsub("T$|TI$|TI*$","", colnames(cosmicSigs.df))
cosmicSigs.df = cosmicSigs.df[cosmicSigs.df$sig_etiology  %in% c("APOBEC","UV","UV_indirect","HRD"),
                              colnames(cosmicSigs.df) %in% c("SBS", purple_pp.df$patientID)] 
rownames(cosmicSigs.df) = NULL
cosmicSigs.df = cosmicSigs.df %>% tibble::column_to_rownames("SBS")
cosmicSigs.df = as.data.frame(t(cosmicSigs.df))
cosmicSigs.df = cosmicSigs.df %>% tibble::rownames_to_column("patientID")

#get driver mutations for genes of interest per patient
driver.df = my_dataframe(fn = pembro_wgs_driver_fn.df, extension = ".driver.catalog.tsv") %>% 
  filter(gene %in% GOI)
driver_pp.df = my_GOI_driverMut_pp(df = driver.df, 
                                   driverLikelihood_thres = 0.2, #in general 0.8 but then we didnt find any sign genes
                                   GOI = GOI,
                                   patientIDs = purple_pp.df$patientID)

#Get snpeff annotated genes
ann.df = my_ann_snpeff(ann.fn, GOI)
ann_pp.df = my_ann_pp(ann.df, GOI, ann.fn)



#Merge everything
merge.df = full_join(x = purple_pp.df, y = clonal_and_SNPeffTMLTMB.df, by = "patientID")
merge.df = full_join(x = merge.df, y = cosmicSigs.df, by = c("patientID"))
merge.df = full_join(x = merge.df, y = driver_pp.df, by = "patientID")
merge.df = full_join(x = merge.df, y = cnv_pp.df, by = "patientID")
merge.df = full_join(x = merge.df, y = ann_pp.df, by =  "patientID")


write.table(merge.df, file = ("data/pembro_wgs_features.csv"),
            quote=F, sep='\t', col.names = T, row.names = F)

