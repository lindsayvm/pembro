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

ann_vcf.fn = "/home/l.leek/pembro/data/snpeff_output/CPCT02080234T_ann.vcf"
ann_vcf_raw.df = as.data.frame(getFIX(read.vcfR(ann_vcf.fn, verbose = FALSE), getINFO = T))

ann_filt_vcf.fn = "/home/l.leek/pembro/data/snpeff_output/CPCT02080234T_ann_filt.vcf"
ann_filt_vcf_raw.df = as.data.frame(getFIX(read.vcfR(ann_filt_vcf.fn, verbose = FALSE), getINFO = T))

oneline_vcf.df = fread("/home/l.leek/pembro/data/snpeff_output/CPCT02080234T_ann_filt_oneLine.vcf")
# oneline_vcf.df$CHROMPOS = paste0(oneline_vcf.df$CHROM, oneline_vcf.df$POS)
# length(unique(oneline_vcf.df$CHROMPOS)) #208 unique locations with mutations
colnames(oneline_vcf.df)
