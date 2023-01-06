library(patchwork)
library(data.table)
library(plyr)
#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv") %>% 
  filter(TumorType == "Breast cancer")

plyr::count(df$Association)

tnb.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv") %>% 
  filter(Association == "ER- PR- HER2-")
table(tnb.df$BOR)
