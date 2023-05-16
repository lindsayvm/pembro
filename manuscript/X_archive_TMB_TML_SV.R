library(patchwork)
library(data.table)
#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")

####### FILTER BASED ON COHORT NOT TMB
#
#
#
3
#
#
3
#
3
#
3
#
3
#

df = fread("data/20230310_DRUP_pembro_LL_WGS_RNA.tsv")
table(df$cohort_manuscript)
tmp = plyr::count(df$TumorType)
tmp[order(tmp$freq, decreasing = TRUE),]

# TMB high low and response/no response
pval_all = get_fisher_pval(df, "cohort_manuscript")
responders_perCohort.p = my_stacked(df , "cohort_manuscript") +
  ylab("# patients") + 
  xlab("") + 
  annotate("text", x=1.5, y=-1, label= paste0("FisherPval=",pval_all)) 

BOR_perCohort.p = my_stacked_BOR(df , "cohort_manuscript")

pval_bc = get_fisher_pval(df %>% dplyr::filter(TumorType == "Breast cancer"), "TumorType")
responders_breastCohort.p = my_stacked(df %>% dplyr::filter(TumorType == "Breast cancer"), "TumorType") +
  ylab("# patients") +
  xlab("") + 
  theme(legend.position = "right") +
  annotate("text", x=1, y=-1, label= paste0("FisherPval=",pval_bc)) 

png("/home/l.leek/pembro/results/fig_BOR_perCohort_barplot.png",width=2200,height=1600, res=300)
layout = c("AAB")
responders_perCohort.p + responders_breastCohort.p + plot_layout(design = layout)
dev.off()


#The test is called "Wilcoxon-Mann-Whitney Rank Sum test". This means, data has to be ranked. If (in the paired test), the differences between the two samples are not unique, you have ties. There is no algorithm for calculating an exact p-value in that case. 
#####################
TML_clin_290.p = my_wilcoxon2(df %>% dplyr::filter(MutationalLoad >= 290), "MutationalLoad")+
  ylab("HMF_old")
TML_purple_290.p = my_wilcoxon2(df %>% dplyr::filter(TML >= 290), "TML")+
  ylab("TML")
TML_SNPeff_290.p = my_wilcoxon2(df %>% dplyr::filter(TML_SNPeff >= 290), "TML_SNPeff")+
  theme(axis.title.x = element_text()) + 
  xlab(">290")+
  ylab("TML_SNPeff")

TML_clin_140_290.p = my_wilcoxon2(df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "MutationalLoad")+
  ylab("HMF_old")
TML_purple_140_290.p = my_wilcoxon2(df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "TML")+
  ylab("TML")
TML_SNPeff_140_290.p = my_wilcoxon2(df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "TML_SNPeff")+
  theme(axis.title.x = element_text()) + 
  xlab("140-290")+
  ylab("TML_SNPeff")

TML_clinbreast.p = my_wilcoxon2(df %>% dplyr::filter(TumorType == "Breast cancer"), "MutationalLoad")+
  ylab("HMF_old")
TML_purplebreast.p = my_wilcoxon2(df %>% dplyr::filter(TumorType == "Breast cancer"), "TML")+
  ylab("TML")
TML_SNPeffbreast.p = my_wilcoxon2(df %>% dplyr::filter(TumorType == "Breast cancer"), "TML_SNPeff")+
  theme(axis.title.x = element_text()) + 
  xlab("Breast cancer")+
  ylab("TML_SNPeff")

png("/home/l.leek/pembro/results/fig_TML_responders_wilc_boxplot.png",width=3000,height=2000, res=300)
#### Final plot
TML_clin_140_290.p+TML_clin_290.p+TML_clinbreast.p+
TML_purple_140_290.p+TML_purple_290.p+TML_purplebreast.p+
TML_SNPeff_140_290.p +TML_SNPeff_290.p+TML_SNPeffbreast.p+
  plot_layout(ncol = 3)

dev.off()



TML_140_290.p = my_wilcoxon(df %>% dplyr::filter(TML < 290 & TML_SNPeff >= 140 ), "TML")
TMB_140_290.p = my_wilcoxon(df %>% dplyr::filter(TML < 290 & TML_SNPeff >= 140 ), "TMB")
SVL_140_290.p = my_wilcoxon(df %>% dplyr::filter(TML < 290 & TML_SNPeff >= 140 ) %>% 
                              dplyr::rename(SVL = svTumorMutationalBurden), 
                            "SVL") +
  theme(axis.title.x = element_text()) + 
  xlab("140-290")

TML_290.p = my_wilcoxon(df %>% dplyr::filter(TML >= 290), "TML")
TMB_290.p = my_wilcoxon(df %>% dplyr::filter(TML >= 290), "TMB")
SVL_290.p = my_wilcoxon(df %>% dplyr::filter(TML >= 290) %>% 
                          dplyr::rename(SVL = svTumorMutationalBurden), 
                        "SVL") +
  theme(axis.title.x = element_text()) + 
  xlab(">290")


TML_clinbreast.p = my_wilcoxon(df %>% dplyr::filter(TumorType == "Breast cancer"), "TML")+
  ylab("HMF_old")
TMB_purplebreast.p = my_wilcoxon(df %>% dplyr::filter(TumorType == "Breast cancer"), "TMB")
SVL_SNPeffbreast.p = my_wilcoxon(df %>% dplyr::filter(TumorType == "Breast cancer") %>% 
                                    dplyr::rename(SVL = svTumorMutationalBurden), 
                                  "SVL") +
  theme(axis.title.x = element_text()) + 
  xlab("Breast cancer")

png("/home/l.leek/pembro/results/fig_TMB_TML_SV_responders_wilc_barplots.png",width=3000,height=2000, res=300)

TML_140_290.p+TML_290.p+TML_clinbreast.p+
  TMB_140_290.p+TMB_290.p+TMB_purplebreast.p+
  SVL_140_290.p +SVL_290.p+SVL_SNPeffbreast.p+
  plot_layout(ncol = 3)

dev.off()

