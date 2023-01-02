#libraries
library(data.table)
library(vcfR)
setwd("/home/l.leek/")

#Ground truth
purity.df = fread("/DATA/share/Voesties/data/HMF/update_10/somatics/CPCT02080234T/purple/CPCT02080234T.purple.purity.tsv")
groundtruth_tml = purity.df$tml
groundtruth_tmb = purity.df$tmbPerMb
print(groundtruth_tml) #1811
print(groundtruth_tmb) #130.43


#Read SNPeff SNPsift output
ann_filt_vcf.fn = "/home/l.leek/pembro/data/snpeff_output/CPCT02080234T_ann_filt.vcf"
ann_filt_vcf_raw.df = as.data.frame(getFIX(read.vcfR(ann_filt_vcf.fn, verbose = FALSE), getINFO = T))

ann_filt_vcf.df = ann_vcf_raw.df %>% 
  mutate(CHROMPOS = paste0(CHROM, "_", POS)) %>%    #counts as 1 mut
  mutate(SUBCL = gsub(".*;SUBCL=|;.*","",INFO)) %>% #extract subclonality info
  dplyr::select(-INFO) %>%                          #filter on quality
  filter(FILTER == "PASS")                          #filter on quality 

#subcl is between 0 and 1
summary(as.numeric(ann_filt_vcf.df$SUBCL)) #between 0 & 1

#Unique mutations that are clonal
clonal_muts = ann_filt_vcf.df$CHROMPOS[as.numeric(ann_filt_vcf.df$SUBCL) < 0.05]
print(length(unique(clonal_muts))) #370463



#Read one line 
oneline_vcf.df = fread("/home/l.leek/pembro/data/snpeff_output/CPCT02080234T_ann_filt_oneLine.vcf")

#should be between 0 and 1, but max is 1.4089
summary(oneline_vcf.df$PURPLE_AF)

#TML
colnames(oneline_vcf.df) = gsub("EFF\\[\\*\\]\\.","", colnames(oneline_vcf.df))
nonsyn.df = oneline_vcf.df %>% 
  mutate(CHROMPOS = paste0(CHROM, "_", POS)) %>%          #counts as 1 mut
  dplyr::filter(FILTER == "PASS") %>%                     #filter on quality 
  dplyr::filter(IMPACT %in% c("MODERATE","HIGH")) %>%     #filter on nonsyn 
  mutate(CHROMPOSGENE = paste0(CHROMPOS,"_",GENE))

#TML
TML = length(unique(nonsyn.df$CHROMPOS)) #2068
print(groundtruth_tml) 
print(TML)
#



#prot_nonsyn
prot_nonsyn.df = nonsyn.df %>% 
  dplyr::filter(BIOTYPE %in% c("protein_coding"))   

#TMB
TMB = length(unique(prot_nonsyn.df$CHROMPOS)) / 2859  #2068
print(groundtruth_tmb) 
print(TMB)
#







