library(patchwork)
library(data.table)
library(ggplot2)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20230310_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F) 

colSums(!is.na(clin.df))

#clonal_tmb_snv clonal_tmb clonal_tmb_indel clonal_neo neo_TMB indels
vars = c("clonal_tmb_snv", "clonal_tmb", "clonal_tmb_indel", "clonal_neo", "neo_TMB", "indels")
cohort = c("Pembro HML>290", "Pembro 140-290", "Breast cancer")
ls=list()
for(i in vars){
  ls_plots = list()
  print(i)
  for(j in cohort){
    #select cohort
    if(j == "Breast cancer"){
      clin_cohort.df = clin.df %>% dplyr::filter(TumorType %in% "Breast cancer")
    }else if(j %in% c("Pembro 140-290" )){
      clin_cohort.df = clin.df %>% dplyr::filter(TML < 290 & TML >= 140 )
    }else if(j %in% c("Pembro HML>290" )){
      clin_cohort.df = clin.df %>% dplyr::filter(TML >= 290)
    }else if(j == "All"){
      clin_cohort.df = clin.df
    }

      ls_plots[[j]] = my_wilcoxon(clin_cohort.df, i) +   theme(axis.title.x = element_text()) + xlab(j)
  }

  ls[[i]] = ls_plots
}

vars = c("clonal_tmb_snv", "clonal_tmb", "clonal_tmb_indel", "clonal_neo", "neo_TMB", "indels","snvs")

#png("/home/l.leek/pembro/results/fig_clonal_neoantigens_perCohort.png",width=3000,height=2000, res=300)
ls[["clonal_tmb_snv"]][["Pembro 140-290"]]+
  ls[["clonal_tmb_snv"]][["Pembro HML>290"]]+
  ls[["clonal_tmb_snv"]][["Breast cancer"]] + 

  ls[["clonal_tmb"]][["Pembro 140-290"]]+
  ls[["clonal_tmb"]][["Pembro HML>290"]]+
  ls[["clonal_tmb"]][["Breast cancer"]] +
  
  ls[["clonal_neo"]][["Pembro 140-290"]]+
  ls[["clonal_neo"]][["Pembro HML>290"]]+
  ls[["clonal_neo"]][["Breast cancer"]] +
  
  ls[["neo_TMB"]][["Pembro 140-290"]]+
  ls[["neo_TMB"]][["Pembro HML>290"]]+
  ls[["neo_TMB"]][["Breast cancer"]] +
  
  ls[["indels"]][["Pembro 140-290"]]+
  ls[["indels"]][["Pembro HML>290"]]+
  ls[["indels"]][["Breast cancer"]] +
    
  ls[["clonal_tmb_indel"]][["Pembro 140-290"]]+
  ls[["clonal_tmb_indel"]][["Pembro HML>290"]]+
  ls[["clonal_tmb_indel"]][["Breast cancer"]] + plot_layout(ncol=3)


#dev.off()
