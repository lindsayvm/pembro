library(patchwork)
library(data.table)
#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv")


# TMB high low and response/no response
tmb_all.p = my_stacked(df , "TMB_bool") +
  xlab("TMB: Complete cohort")
tmb_cohort140_290.p = my_stacked(df %>% dplyr::filter(Cohort != "Pembro HML>290"), "TMB_bool")  +
  xlab("TMB: Cohort 140-290")
tmb_cohort290.p = my_stacked(df %>% dplyr::filter(Cohort == "Pembro HML>290"), "TMB_bool") +
  xlab("TMB: Cohort >290")

png("/home/l.leek/pembro/results/fig_TMBhighlow_CB_fisher.png",width=2000,height=800, res=300)
tmb_all.p+ylab("# patients")+theme(axis.title.y=element_text(angle = 90))  + tmb_cohort140_290.p + tmb_cohort290.p+theme(legend.position = "right") +
  plot_layout(ncol = 3)
dev.off()

get_fisher_pval(df, "TMB_bool")
get_fisher_pval(df, "TML_bool") 

#####################
tml_all_BOR.p = my_stacked_BOR(df , "TMB_bool") +
  xlab("TMB: Complete cohort")
#####################


# TMB high low and response/no response
tmb_all.p = my_stacked(df , "TMB_bool") +
  xlab("TMB: Complete cohort")
tmb_cohort_BC.p = my_stacked(df %>% dplyr::filter(TumorType == "Breast cancer"), "TMB_bool")  +
  xlab("TMB: Cohort BC")
tmb_cohort_noBC.p = my_stacked(df %>% dplyr::filter(TumorType != "Breast cancer"), "TMB_bool") +
  xlab("TMB: Cohort pan-cancer")

png("/home/l.leek/pembro/results/fig_TML_CB_pertissue.png",width=2000,height=800, res=300)
tmb_all.p +ylab("# patients")+theme(axis.title.y=element_text(angle = 90))  + tmb_cohort_BC.p + tmb_cohort_noBC.p +theme(legend.position = "right") +
  plot_layout(ncol = 3)
dev.off()


#The test is called "Wilcoxon-Mann-Whitney Rank Sum test". This means, data has to be ranked. If (in the paired test), the differences between the two samples are not unique, you have ties. There is no algorithm for calculating an exact p-value in that case. 
#####################
TML_SNPeff_290.p = my_wilcoxon(df %>% dplyr::filter(TML_SNPeff >= 290), "TML_SNPeff")+
  ylim(290,3300)
TML_purple_290.p = my_wilcoxon(df %>% dplyr::filter(TML >= 290), "TML")+
  ylim(290,3300)
TML_clin_290.p = my_wilcoxon(df %>% dplyr::filter(MutationalLoad >= 290), "MutationalLoad")+
  ylim(290,3300)
#####################
TML_SNPeff_140_290.p = my_wilcoxon(df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "TML_SNPeff")+
  ylim(100,290)
TML_purple_140_290.p = my_wilcoxon(df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "TML")+
  ylim(100,290)
TML_clin_140_290.p = my_wilcoxon(df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "MutationalLoad")+
  ylim(100,290)

png("/home/l.leek/pembro/results/fig_TML_CB_wilc.png",width=2000,height=1500, res=300)

#### Final plot
  TML_clin_290.p + TML_purple_290.p + TML_SNPeff_290.p +
  TML_clin_140_290.p +  TML_purple_140_290.p  + TML_SNPeff_140_290.p +
  plot_layout(ncol = 3)

dev.off()

#combined shows its not a power issue
# TML_SNPeff.p = my_wilcoxon(df, "TML_SNPeff")+
#   ylim(0,3300)
# TML_purple.p = my_wilcoxon(df, "TML")+
#   ylim(0,3300)
# TML_clin.p = my_wilcoxon(df, "MutationalLoad")+
#   ylim(0,3300)



png("/home/l.leek/pembro/results/fig_TML_SV_CB_wilc.png",width=2000,height=1500, res=300)

TML_290.p = my_wilcoxon(df %>% dplyr::filter(Cohort == "Pembro HML>290"), "TML")+
  ylim(0,3300)
TMB_290.p = my_wilcoxon(df %>% dplyr::filter(Cohort == "Pembro HML>290"), "TMB")+
  ylim(0,200)
SVL_290.p = my_wilcoxon(df %>% dplyr::filter(Cohort == "Pembro HML>290") %>% 
                          rename(SVL = svTumorMutationalBurden),
                        "SVL") +
  ylim(0,1500)

TML_140_290.p = my_wilcoxon(df %>% dplyr::filter(Cohort != "Pembro HML>290"), "TML")+
  ylim(0,300)
TMB_140_290.p = my_wilcoxon(df %>% dplyr::filter(Cohort != "Pembro HML>290"), "TMB")+
  ylim(0,25)
SVL_140_290.p = my_wilcoxon(df %>% dplyr::filter(Cohort != "Pembro HML>290") %>% 
                                   rename(SVL = svTumorMutationalBurden),
                                    "SVL") +
  ylim(0,1500)

TML_290.p + TMB_290.p + SVL_290.p + 
  TML_140_290.p + TMB_140_290.p + SVL_140_290.p +
  plot_layout(ncol = 3)

dev.off()

