library(data.table)
library(vcfR)
setwd("/home/l.leek/")


# driver = fread("/DATA/share/Voesties/data/HMF/update_10/somatics/CPCT02080234T/purple/CPCT02080234T.driver.catalog.somatic.tsv")
# 17 drivers

vcf.fn = "/DATA/share/Voesties/data/HMF/update_10/somatics/CPCT02080234T/purple/CPCT02080234T.purple.somatic.vcf.gz"
vcf_raw.df = as.data.frame(getFIX(read.vcfR(vcf.fn, verbose = FALSE), getINFO = T))
vcf_raw.df$INFO[1]
#process #SUBCL	1	Subclonal likelihood between 0 and 1
vcf.df = vcf_raw.df %>% 
  mutate(CHROM = as.character(CHROM)) %>% 
  mutate(POS = as.numeric(as.character(POS))) %>% 
  #filter(FILTER == "PASS") %>% 
  mutate(SUBCL = gsub(".*;SUBCL=|;.*","",INFO)) %>% 
  dplyr::select(-INFO)
summary(as.numeric(vcf.df$SUBCL)) #between 0 & 1


ann_vcf.fn = "/home/l.leek/pembro/data/snpeff_output/CPCT02080234T_ann.vcf"
ann_vcf_raw.df = as.data.frame(getFIX(read.vcfR(ann_vcf.fn, verbose = FALSE), getINFO = T))
ann_vcf.df = ann_vcf_raw.df %>% 
  mutate(CHROM = as.character(CHROM)) %>% 
  mutate(POS = as.numeric(as.character(POS))) %>% 
  mutate(SUBCL = gsub(".*;SUBCL=|;.*","",INFO)) %>% 
  dplyr::select(-INFO)
summary(as.numeric(ann_vcf.df$SUBCL)) #between 0 & 1


ann_filt_vcf.fn = "/home/l.leek/pembro/data/snpeff_output/CPCT02080234T_ann_filt.vcf"
ann_filt_vcf_raw.df = as.data.frame(getFIX(read.vcfR(ann_filt_vcf.fn, verbose = FALSE), getINFO = T))
ann_filt_vcf.df = ann_vcf_raw.df %>% 
  mutate(CHROM = as.character(CHROM)) %>% 
  mutate(POS = as.numeric(as.character(POS))) %>% 
  mutate(SUBCL = gsub(".*;SUBCL=|;.*","",INFO)) %>% 
  dplyr::select(-INFO)
summary(as.numeric(ann_filt_vcf.df$SUBCL)) #between 0 & 1
head(ann_filt_vcf.df)


oneline_vcf.df = fread("/home/l.leek/pembro/data/snpeff_output/CPCT02080234T_ann_filt_oneLine.vcf")
#should be between 0 and 1, but max is 1.4089
summary(oneline_vcf.df$PURPLE_AF)


#TML
colnames(oneline_vcf.df) = gsub("EFF\\[\\*\\]\\.","", colnames(oneline_vcf.df))
nonsyn.df = oneline_vcf.df %>% 
  dplyr::filter(IMPACT %in% c("MODERATE","HIGH")) %>% 
  mutate(CHROMPOS = paste0(CHROM, "_", POS)) %>% 
  mutate(CHROMPOSGENE = paste0(CHROMPOS,"_",GENE))


length(unique(nonsyn.df$CHROMPOS)) #2068
length(unique(nonsyn.df$CHROMPOSGENE)) #2083 #what to do when there are multiple genes per position
# for (chrompos in unique(nonsyn.df$CHROMPOS)){
#   tmp = nonsyn.df %>% 
#           dplyr::select(GENE, CHROMPOS) %>% 
#           dplyr::filter(CHROMPOS == chrompos)
#   if(length(unique(tmp$GENE)) > 1){
#     print(chrompos)
#   }
# }
#"5_114956173","5_139905917","6_31910799","9_112918765","12_21032466","12_21176191","12_56693942","16_4463315","16_75485645","16_75485656","19_57929344","20_48713346","X_38546895"
nonsyn.df[nonsyn.df$CHROMPOS == "5_114956173", ]
nonsyn.df[nonsyn.df$CHROMPOS == "5_139905917", ]
nonsyn.df[nonsyn.df$CHROMPOS == "6_31910799", ]
nonsyn.df[nonsyn.df$CHROMPOS == "9_112918765", ]


#
TML = length(unique(nonsyn.df$CHROMPOS)) #2068
subclonal_tml = nonsyn.df %>% 
  filter(PURPLE_AF <= 0.05)
#TMB
prot_nonsyn.df = nonsyn.df %>% 
  dplyr::filter(BIOTYPE %in% c("protein_coding"))   
tmb = length(unique(prot_nonsyn.df$CHROMPOS)) / 2859 

