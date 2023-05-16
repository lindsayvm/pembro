library(patchwork)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20230310_DRUP_pembro_LL_WGS_RNA.tsv")
perc_clon.df = clin.df[!is.na(clin.df$perc_clon), ] 
perc_clon.df = perc_clon.df %>% 
  arrange(desc(perc_clon)) %>% 
  mutate(id = 1:nrow(perc_clon.df)) %>% 
  dplyr::select(id, perc_clon, responders)

perc_clon_scatter.p = ggplot(perc_clon.df, aes(x=id, y=perc_clon, color=factor(responders))) + 
  geom_point(size=3) +
  theme_bw(base_size = 12)+
  theme_bw(base_size = 15)+
  theme(legend.position = "None",
        #        axis.text.x = element_text(size = 16, angle = 45, vjust = 1.5, hjust=0.5),
        axis.text.x = element_text(size = 17, angle = 0, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  scale_color_manual(values = c("forestgreen", "brown")) +
  ylab("Clonality percentage") +
  xlab(paste0("IDs (total WGS = ",nrow(perc_clon.df),")"))

perc_clonal_290.p = my_wilcoxon(clin.df %>% dplyr::filter(TML >= 290), "perc_clon",pval_ypos = 110)+
  theme(axis.title.x = element_text()) + 
  xlab(">290")
perc_clonal_140_290.p = my_wilcoxon(clin.df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "perc_clon",pval_ypos = 110)+
  theme(axis.title.x = element_text()) + 
  xlab("140-290")
perc_clonal_breast.p = my_wilcoxon(clin.df %>% dplyr::filter(TumorType == "Breast cancer"), "perc_clon",pval_ypos = 110)+
  theme(axis.title.x = element_text()) + 
  xlab("Breast cancer")

aov_perc_clon_290.p = my_anova_plot(clin.df %>% dplyr::filter(TML >= 290), "perc_clon",pval_ypos = 110)+
  theme(axis.title.x = element_text()) + 
  xlab(">290")
aov_perc_clon_140_290.p = my_anova_plot(clin.df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "perc_clon",pval_ypos = 110)+
  theme(axis.title.x = element_text()) + 
  xlab("140-290")
aov_perc_clonal_breast.p = my_anova_plot(clin.df %>% dplyr::filter(TumorType == "Breast cancer"), "perc_clon",pval_ypos = 110)+
  theme(axis.title.x = element_text()) + 
  xlab("Breast cancer")


cTML_290.p = my_wilcoxon(clin.df %>% dplyr::filter(TML >= 290), "cTML")+
  theme(axis.title.x = element_text()) + 
  xlab(">290")
cTML_140_290.p = my_wilcoxon(clin.df %>% dplyr::filter(TML_SNPeff < 290 & TML_SNPeff >= 140 ), "cTML")+
  theme(axis.title.x = element_text()) + 
  xlab("140-290")
cTMLal_breast.p = my_wilcoxon(clin.df %>% dplyr::filter(TumorType == "Breast cancer"), "cTML")+
  theme(axis.title.x = element_text()) + 
  xlab("Breast cancer")



png("/home/l.leek/pembro/results/fig_perc_clon.png",width=6000,height=3000, res=300)

layout = c("AAABCD
            AAAEFG
            AAAHIJ")
perc_clon_scatter.p+
  perc_clonal_140_290.p +aov_perc_clon_140_290.p +cTML_140_290.p+
  perc_clonal_290.p     +aov_perc_clon_290.p     +cTML_290.p+
    perc_clonal_breast.p +aov_perc_clonal_breast.p+cTMLal_breast.p+
  plot_layout(design = layout)
dev.off()
