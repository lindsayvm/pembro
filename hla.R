dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20230310_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F) 
hla_status = fread("/DATA/share/Voesties/data/harmonize/output/mhc_class_i/lilac-hla-wide.tsv", data.table = F)

#
hla.df = fread("/DATA/share/Voesties/data/harmonize/output/mhc_class_i/lilac-hla.csv") %>% 
  dplyr::filter(sampleId %in% clin.df$sampleId)

#TO DO grouped per gene; and then grouped per patient

filtered.df = hla.df[hla.df$TumorCopyNumber <= 0.5 |
                       hla.df$SomaticMissense > 0 |
                       hla.df$SomaticNonsenseOrFrameshift > 0 |
                       hla.df$SomaticSplice > 0 |
                       hla.df$SomaticInframeIndel > 0, ] #dont look at SomaticSynonymous because that wont give a change in peptides
#Can also try 0.5 instead of 0 because most of them go from 0 to 2 alleles; so 0.5 should indicate that at least 1 allele has a mutation

final.df = left_join(clin.df, filtered.df, "sampleId")

#count for each patient whether having a mutation is predictive of resistance and correlates to TMB
# ,, more mutations ..
#.. a specific class of mutations is ... Maybe only B is relevant? 




