library(patchwork)
library(data.table)
library(ggplot2)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F) %>% 
  filter(TumorType == "Breast cancer")

var = "total_neo"
total_neo_aov.p = my_anova_plot(clin.df, var)
total_neo_wilc.p = my_wilcoxon_plot2(clin.df %>%  filter(BOR != "SD") , var)
total_neo_wilc_BOR.p = my_wilcoxon(clin.df, var)

var = "clonal_neo"
clonal_neo_aov.p = my_anova_plot(clin.df, var)
clonal_neo_wilc.p = my_wilcoxon_plot2(clin.df %>%  filter(BOR != "SD") , var)
clonal_neo_wilc_BOR.p = my_wilcoxon(clin.df, var)

var = "subclonal_neo"
subclonal_neo_aov.p = my_anova_plot(clin.df, var)
subclonal_neo_wilc.p = my_wilcoxon_plot2(clin.df %>%  filter(BOR != "SD") , var)
subclonal_neo_wilc_BOR.p = my_wilcoxon(clin.df, var)

var = "fusion_neo"
fusion_neo_aov.p = my_anova_plot(clin.df, var)
fusion_neo_wilc.p = my_wilcoxon_plot2(clin.df %>%  filter(BOR != "SD") , var)
fusion_neo_wilc_BOR.p = my_wilcoxon(clin.df, var)

var = "mut_neo"
mut_neo_aov.p = my_anova_plot(clin.df, var)
mut_neo_wilc.p = my_wilcoxon_plot2(clin.df %>%  filter(BOR != "SD") , var)
mut_neo_wilc_BOR.p = my_wilcoxon(clin.df, var)

png("/home/l.leek/pembro/results/fig_neoantigens_breast.png",width=3000,height=4500, res=300)
total_neo_aov.p + total_neo_wilc.p + total_neo_wilc_BOR.p +
clonal_neo_aov.p + clonal_neo_wilc.p + clonal_neo_wilc_BOR.p + 
subclonal_neo_aov.p + subclonal_neo_wilc.p + subclonal_neo_wilc_BOR.p + 
fusion_neo_aov.p + fusion_neo_wilc.p + fusion_neo_wilc_BOR.p + 
mut_neo_aov.p + mut_neo_wilc.p + mut_neo_wilc_BOR.p +
  plot_layout(ncol = 3)
dev.off()
