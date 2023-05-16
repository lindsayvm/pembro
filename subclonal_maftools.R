#Clonal
# 0.85 GIE
# selected everything but silent mutation
# maftools.df = readRDS("/DATA/share/Voesties/data/harmonize/output/hmf-maftools.rda")
# maftools.df@data
# summary(maftools.df@data)
# SUBCL 0 or 1 
# 
# VCF and thus SNPeff has also all non-coding what it needs to be filtered on, and maftools is already filtered out on silent mutations etc.
# 1 Reading all SNPeff annotated files and filtering for nonsyn
# 2 Extracting clonal information (SUBCL) per mutation from vcf
# 3 select clonal (SUBCL <0.05) nonsyn mutations (BIOTYPE == "protein_coding",
#                                               (EFFECT == "missense_variant"))

#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")

clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv")

maftools.rds = readRDS("/DATA/share/Voesties/data/harmonize/output/hmf-maftools.rda")
maftools.df = maftools.rds@data
maftools.df$SUBCL

#GIE: clonal mutations as those with a subclonal score lower than 0.85.
TML_pp.df = maftools.df %>% 
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::count() %>% 
  mutate(TML_maftools = n) %>% 
  dplyr::select(-c(n))

cTML_pp.df = maftools.df %>% 
  dplyr::filter(SUBCL < 0.85) %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::count() %>% 
  mutate(patientID = gsub("T$|TI*$","",Tumor_Sample_Barcode)) %>% 
  mutate(cTML_maftools = n) %>% 
  dplyr::select(-c(n))

#quite some are duplicated
duplicated.v = cTML_pp.df$patientID[duplicated(cTML_TML.df$patientID)]
#but none of the multiple samples per patient are present in our cohort
duplicated.v[duplicated.v %in% clin.df$patientID]

cTML_TML.df = left_join(cTML_pp.df, TML_pp.df, by = "Tumor_Sample_Barcode")


clin.df = left_join(clin.df, cTML_TML.df, by = "patientID")
VAR = "cTML_maftools" #cTML_maftools, TML_maftools, cTML, TML


cTML_all.p = my_wilcoxon(clin.df, VAR)+
  theme(axis.title.x = element_text()) + 
  xlab("all")
cTML_all_nBC.p = my_wilcoxon(clin.df %>% dplyr::filter(TumorType != "Breast cancer"), VAR)+
  theme(axis.title.x = element_text()) + 
  xlab("all (no BC)")
cTML_290.p = my_wilcoxon(clin.df %>% dplyr::filter(TML_SNPeff >= 290), VAR)+
  theme(axis.title.x = element_text()) + 
  xlab(">290 (incl. BC)")
cTML_140_290.p = my_wilcoxon(clin.df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), VAR)+
  theme(axis.title.x = element_text()) + 
  xlab("140-290 (incl. BC)")

cTML_290_nBC.p = my_wilcoxon(clin.df %>% dplyr::filter(TumorType != "Breast cancer" &
                                                               TML_SNPeff >= 290), VAR)+
  theme(axis.title.x = element_text()) + 
  xlab(">290 (no BC)")
cTML_140_290_nBC.p = my_wilcoxon(clin.df %>% dplyr::filter(TumorType != "Breast cancer" &
                                                         TML_SNPeff < 290 & TML_SNPeff >= 140 ), VAR)+
  theme(axis.title.x = element_text()) + 
  xlab("140-290 (no BC)")
cTML_breast.p = my_wilcoxon(clin.df %>% dplyr::filter(TumorType == "Breast cancer"), VAR)+
  theme(axis.title.x = element_text()) + 
  xlab("Breast cancer")

layout="AB#
        CD#
        EFG"
cTML_all.p+cTML_all_nBC.p+
cTML_140_290.p+cTML_290.p+
cTML_140_290_nBC.p+cTML_290_nBC.p+cTML_breast.p+
  plot_layout(design = layout)

