library(patchwork)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv") %>% 
  dplyr::filter(TumorType == "Breast cancer" )
perc_clon.df = clin.df[!is.na(clin.df$perc_clon), ] 
perc_clon.df = perc_clon.df %>% 
  arrange(desc(perc_clon)) %>% 
  mutate(id = 1:nrow(perc_clon.df)) %>% 
  dplyr::select(id, perc_clon)


perc_clon_scatter.p = ggplot(perc_clon.df, aes(x=id, y=as.numeric(perc_clon))) +
  geom_point() +
  theme_bw(base_size = 12)+
  theme(legend.position = "None",
        axis.text.x = element_text(size = 12, angle = 33, vjust = 1, hjust=1)) +
  ylab("Clonality percentage") +
  xlab("Total WGS = 62")


aov_TML.p = my_anova_plot(clin.df, "TML_SNPeff", pval_ypos = 3500)
aov_cTML.p = my_anova_plot(clin.df, "cTML", pval_ypos = 3500)
aov_perc_clon.p = my_anova_plot(clin.df, "perc_clon", pval_ypos = 110)

wilc_TML.p = my_wilcoxon_plot2(clin.df %>%  filter(BOR != "SD") , "TML_SNPeff", pval_ypos = 3500)
wilc_cTML.p = my_wilcoxon_plot2(clin.df %>%  filter(BOR != "SD") , "cTML", pval_ypos = 3500)
wilc_perc_clon.p = my_wilcoxon_plot2(clin.df %>%  filter(BOR != "SD") , "perc_clon", pval_ypos = 110)

png("/home/l.leek/pembro/results/fig_clonalTML_breast.png",width=3000,height=4500, res=300)
aov_TML.p + aov_cTML.p + aov_perc_clon.p +
wilc_TML.p + wilc_cTML.p + wilc_perc_clon.p +
   plot_layout(ncol = 3)
dev.off()
