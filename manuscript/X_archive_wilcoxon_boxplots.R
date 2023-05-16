library(ggplot2)

raw.df = fread("data/20230503_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F)
clin.df = raw.df %>% dplyr::select(subjectkey, Cohort_trueHML, TumorType, responders,
                                   TML, TMB, svTumorMutationalBurden, 
                                   cTML_maftools_05, cTML_maftools, cTML, cTML_05, perc_clon,
                                   aneuploidyScore, aneuploidy_score,ploidy,avg_cnv,
                                   sex,
                                   gene_PBRM1, gene_snpeff_SMARCB1,gene_snpeff_SMARC4A, gene_ARID1A)
clin.df$swisnf = clin.df$gene_PBRM1
clin.df$swisnf[clin.df$gene_PBRM1 == 1 |
               clin.df$gene_snpeff_SMARCB1 == 1 |
               clin.df$gene_snpeff_SMARC4A == 1 |
               clin.df$gene_ARID1A == 1] = 1

all.df = clin.df %>%
  mutate(group = "All")
other.df = clin.df %>% 
  dplyr::filter(TumorType != "Breast") %>% 
  mutate(group = "Other")
breast.df = clin.df %>% 
  dplyr::filter(TumorType == "Breast") %>% 
  mutate(group = "Breast")
HML140290_other.df = clin.df %>% 
  dplyr::filter(Cohort_trueHML == "HML 140-290, other") %>% 
  mutate(group = "HML 140-290 other")
HML140290_breast.df = clin.df %>% 
  dplyr::filter(Cohort_trueHML == "HML 140-290, breast") %>% 
  mutate(group = "HML 140-290 breast")
HML290.df = clin.df %>% 
  dplyr::filter(Cohort_trueHML == "HML > 290") %>% 
  mutate(group = "HML > 290")

combined.df = rbind(all.df, other.df, breast.df, HML140290_other.df, HML140290_breast.df, HML290.df)
combined.df = combined.df %>% 
  mutate(group = factor(group, 
                        levels = c("All","Other","Breast",
                                   "HML 140-290 other",
                                   "HML 140-290 breast",
                                   "HML > 290")))

my_wilc_grouped <- function(df, VAR, ylab){
  
  count_df = combined.df %>%
    group_by(group, responders) %>%
    summarise(count = n())
  
  ggplot(df, aes(x = group, y = df[[VAR]], fill = responders)) +
    geom_boxplot() +
    labs(x = "",  y=ylab, fill = "") +
    scale_fill_manual(values = c("brown", "forestgreen")) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_bw(base_size = 15)+
    theme(legend.position = "None",
          axis.text.x = element_text(size = 17, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust=1),
          axis.title.x = element_blank(),
          axis.ticks = element_blank()) +
    stat_compare_means(aes(group = responders), method = "wilcox.test",
                       label = "p.format") + #size = 6
    annotate("text",
             #size = 7,
             x = c(0.85, 1.15, 1.85, 2.15, 2.85, 3.15, 3.85, 4.15, 4.85, 5.15, 5.85, 6.15),
             y = 0,
             label = count_df$count,
             vjust = 2 ) +
    scale_y_continuous(expand=expansion(mult = c(0.1,0.1)))
}



tml.p=my_wilc_grouped(combined.df, VAR = "TML", ylab = "TML")
tmb.p=my_wilc_grouped(combined.df, VAR = "TMB", ylab = "TMB")
svl.p=my_wilc_grouped(combined.df, VAR = "svTumorMutationalBurden", ylab = "SVL")
ctml.p=my_wilc_grouped(combined.df, VAR = "cTML_maftools", ylab = "Clonal TML")
percClon.p=my_wilc_grouped(combined.df, VAR = "perc_clon", ylab = "Clonal %")
aneup.p=my_wilc_grouped(combined.df, VAR = "aneuploidyScore", ylab = "Aneuploidy")
ploidy.p=my_wilc_grouped(combined.df, VAR = "ploidy", ylab = "ploidy")

pdf(file = "/home/l.leek/pembro/results/fig_groups_TML_SV_Aneuploidy_etc.pdf",   
    width = 14, # The width of the plot in inches
    height = 30)
tml.p + tmb.p + svl.p + ctml.p + aneup.p + plot_layout(ncol=2)
dev.off()

table(clin.df$sex, clin.df$responders)
table(clin.df$swisnf, clin.df$responders)
table(clin.df$gene_PBRM1, clin.df$responders)
table(clin.df$gene_ARID1A, clin.df$responders)



########Quality check

my_wilcoxon2 <- function(df, VAR = VAR, pval_ypos = log2(max(df[[VAR]], na.rm = T))){
  
  df$responders = factor(df$responders, levels = c("R","NR"))
  df = df[!is.na(df[[VAR]]), ]
  
  log2Var = log2(df[[VAR]])
  wilc.test = pairwise.wilcox.test(log2Var, df$responders,
                                   p.adjust.method="none")
  pvalue = data.frame(signif(wilc.test$p.value,2))
  
  p = df %>%  ggplot( aes(x=responders, 
                          y=log2Var, 
                          fill=responders)) +
    # geom_violin(width=1.4) +
    geom_boxplot() +
    scale_fill_manual(values = c("forestgreen", "brown")) +
    geom_jitter(color="black", size=0.5, alpha=0.3, width = 0.2) +
    ylab(paste0("Log2", "VAR")) +
    theme_bw(base_size = 15)+
    theme(legend.position = "None",
          #        axis.text.x = element_text(size = 16, angle = 45, vjust = 1.5, hjust=0.5),
          axis.text.x = element_text(size = 17, angle = 0, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust=1),
          axis.title.x = element_blank(),
          axis.ticks = element_blank()) +
    geom_text(data = pvalue,
              aes(x = 1.5, y = pval_ypos, label = paste0("pval_wilc=",pvalue)),  
              inherit.aes = FALSE, hjust = "inward", vjust = "inward", size = 3.5) +
    ylim(6.9,12)
  return(p)
}  


my_wilcoxon2(all.df, "TML") #0.13
my_wilcoxon2(other.df, "TML") #0.095
my_wilcoxon2(breast.df, "TML") #0.97
my_wilcoxon2(HML140290_other.df, "TML") #0.26
my_wilcoxon2(HML140290_breast.df, "TML") #1
my_wilcoxon2(HML290.df, "TML") #0.0045
