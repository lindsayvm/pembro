library(patchwork)
library(data.table)
library(dplyr)
#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv") %>% 
  dplyr::filter(TumorType == "Breast cancer" )

# TMB high low and response/no response
tmb_all.p = my_stacked(df , "TMB_bool") +
  xlab("TMB: Complete cohort")
tmb_cohort140_290.p = my_stacked(df %>% dplyr::filter(Cohort != "Pembro HML>290"), "TMB_bool")  +
  xlab("TMB: Cohort >290")
tmb_cohort290.p = my_stacked(df %>% dplyr::filter(Cohort == "Pembro HML>290"), "TMB_bool") +
  xlab("TMB: Cohort 140-290")


# TML high low and response/no response
tml_all.p = my_stacked(df , "TML_bool") +
  xlab("TML: Complete cohort")
tml_cohort140_290.p = my_stacked(df %>% dplyr::filter(Cohort != "Pembro HML>290"), "TML_bool")  +
  xlab("TML: Cohort >290")
tml_cohort290.p = my_stacked(df %>% dplyr::filter(Cohort == "Pembro HML>290"), "TML_bool") +
  xlab("TML: Cohort 140-290")

png("/home/l.leek/pembro/results/fig_TML_CB_breast.png",width=2000,height=1500, res=300)

tmb_all.p+ylab("# patients")+theme(axis.title.y=element_text(angle = 90))  + tmb_cohort140_290.p + tmb_cohort290.p+theme(legend.position = "right") + 
tml_all.p + tml_cohort140_290.p + tml_cohort290.p +
  plot_layout(ncol = 3)

dev.off()

get_fisher_pval(df, "TMB_bool")
get_fisher_pval(df, "TML_bool") #0.368

#####################
tml_all_BOR.p = my_stacked_BOR(df , "TML_bool") +
  xlab("TML: Complete cohort")
######???????????
#####################


# TMB high low and response/no response
tmb_all.p = my_stacked(df , "TMB_bool") +
  xlab("TMB: Complete cohort")
tmb_cohort_BC.p = my_stacked(df %>% dplyr::filter(TumorType == "Breast cancer"), "TMB_bool")  +
  xlab("TMB: Cohort BC")
tmb_cohort_noBC.p = my_stacked(df %>% dplyr::filter(TumorType != "Breast cancer"), "TMB_bool") +
  xlab("TMB: Cohort pan-cancer")

tml_all.p = my_stacked(df , "TML_bool") +
  xlab("TML: Complete cohort")
tml_cohort_BC.p = my_stacked(df %>% dplyr::filter(TumorType == "Breast cancer"), "TML_bool")  +
  xlab("TML: Cohort BC")
tml_cohort_noBC.p = my_stacked(df %>% dplyr::filter(TumorType != "Breast cancer"), "TML_bool") +
  xlab("TML: Cohort pan-cancer")

png("/home/l.leek/pembro/results/fig_TML_CB_pertissuu_breast.png",width=2000,height=1500, res=300)

tmb_all.p +ylab("# patients")+theme(axis.title.y=element_text(angle = 90))  + tmb_cohort_BC.p + tmb_cohort_noBC.p +theme(legend.position = "right") + 
tml_all.p + tml_cohort_BC.p + tml_cohort_noBC.p + 
  plot_layout(ncol = 3)

dev.off()

#####################
TML_SNPeff.p = my_wilcoxon(df, "TML_SNPeff")+
  ylim(0,3300)
TML_purple.p = my_wilcoxon(df, "TML")+
  ylim(0,3300)
TML_clin.p = my_wilcoxon(df, "MutationalLoad")+
  ylim(0,3300)
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

png("/home/l.leek/pembro/results/fig_TML_CB_wilc_breast.png",width=2000,height=1500, res=300)

#### Final plot
TML_clin.p + TML_purple.p + TML_SNPeff.p +
  TML_clin_290.p + TML_purple_290.p + TML_SNPeff_290.p +
  TML_clin_140_290.p +  TML_purple_140_290.p  + TML_SNPeff_140_290.p +
  plot_layout(ncol = 3)

dev.off()


png("/home/l.leek/pembro/results/fig_TML_SV_CB_wilc_breast.png",width=2000,height=1500, res=300)

TML.p = my_wilcoxon(df, "TML")+
  ylim(0,3300)
TMB.p = my_wilcoxon(df, "TMB")+
  ylim(0,200)
SVL.p = my_wilcoxon(df , "svTumorMutationalBurden")+
  ylim(0,1500)
TML.p+ TMB.p + SVL.p 

dev.off()

