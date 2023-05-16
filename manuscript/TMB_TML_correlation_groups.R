library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)


#https://www.nature.com/articles/s41416-020-0762-5/figures/1 
raw.df = fread("data/20230503_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F)
clin.df = raw.df %>% dplyr::select(subjectkey, Cohort, Cohort_trueHML,
                                   TML, TMB) %>% 
  mutate(Cohort = factor(Cohort, levels = c("HML > 290",
                                            "HML 140-290, breast",
                                            "HML 140-290, other",
                                            "Excluded")))

p1 = ggplot(clin.df, aes(x=TMB, y=TML, color=Cohort)) + 
  geom_point(size=1) +
  scale_color_manual(values = c("brown",  "orange", "skyblue" ,"grey")) + 
  theme_bw(base_size = 15) +
  theme(legend.position = "none") +
  ggtitle("Inititial cohorts") +
  geom_hline(yintercept = (290), alpha=0.6, linetype = "dashed") +
  geom_hline(yintercept = (140), alpha=0.6, linetype = "dashed")


p2 = ggplot(clin.df, aes(x=log10(TMB), y=log10(TML), color=Cohort)) + 
  geom_point(size=1) +
  theme_bw(base_size = 15)+
  scale_color_manual(values = c("brown",  "orange", "skyblue", "grey")) + 
  geom_hline(yintercept = log10(290), alpha=0.6, linetype = "dashed") +
  geom_hline(yintercept = log10(140), alpha=0.6, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, size = 0.5)  +
  # geom_smooth(data = clin.df, method = "lm", se = FALSE, color = "grey",  span = 0.3,alpha=0.1)
  stat_cor(label.sep = " = ", 
           aes(color = Cohort), method = "pearson", size = 4)
p1+p2
