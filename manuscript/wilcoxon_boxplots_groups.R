#libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
setwd("/home/l.leek/pembro/")

#Params
COHORT = "Cohort_trueHML" #"Cohort_initial"#"Cohort_trueHML"
GROUPS = "140290together" #"140290splitted #140290together
RemoveLowerThan140 = "KeepBelow140" #"RemoveBelow140" "KeepBelow140" #only relevant for cohort_trueHML

#functions
my_wilc_grouped <- function(df, VAR, ylab){
  
  ggplot(df, aes(x = group, y = df[[VAR]], fill = responders)) +
    geom_boxplot() +
    labs(x = "",  y=ylab, fill = "") +
    scale_fill_manual(values = c("brown", "forestgreen")) +
    geom_jitter(color="black", size=0.4, alpha=0.9,width = 0.25) +
    theme_bw(base_size = 15)+
    theme(legend.position = "None",
          axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust=1),
          axis.title.x = element_blank(),
          axis.ticks = element_blank()) +
    stat_compare_means(aes(group = responders), method = "wilcox.test",
                       label = "p.format",
                       label.y = max(combined.df[[VAR]], na.rm = T) + max(combined.df[[VAR]], na.rm = T)*0.1) + #size = 6
    scale_y_continuous(expand=expansion(mult = c(0.1,0.2)))
}

#Load data
raw.df = fread("data/20230503_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F)
clin.df = raw.df %>% dplyr::select(subjectkey, Cohort, Cohort_trueHML, TumorType, responders,
                                   TML, TMB, svTumorMutationalBurden, 
                                   cTML_maftools_05, cTML_maftools, cTML, cTML_05, perc_clon,
                                   aneuploidy_score, ploidy,avg_cnv,
                                   sex, SBS7a, SBS7b,
                                   gene_PBRM1, gene_snpeff_SMARCB1,gene_snpeff_SMARC4A, gene_ARID1A)
tmp = clin.df %>% dplyr::select(SBS7a, responders,TumorType)

#Cohort refers to cohorts used in trial, Cohort_trueHML is defined by true TML and excludes 3 patients with TML<140
if(COHORT == "Cohort_trueHML"){
  table(clin.df$Cohort_trueHML)
  if(RemoveLowerThan140 == "RemoveBelow140"){
    clin.df$Cohort_trueHML[clin.df$TML<100] = "Excluded"
    clin.df = clin.df %>% dplyr::filter(Cohort_trueHML != "Excluded")
  }
  clin.df = clin.df %>% mutate(Cohort = Cohort_trueHML)
  print(COHORT)
}else{
  print(COHORT)
}

#Final boxplot will have groups with overlapping patients ==> df for each cohort and add group to distinguish same patients belonging to >1 cohorts
all.df = clin.df %>%
  mutate(group = "All")
other.df = clin.df %>% 
  dplyr::filter(TumorType != "Breast") %>% 
  mutate(group = "Other")
breast.df = clin.df %>% 
  dplyr::filter(TumorType == "Breast") %>% 
  mutate(group = "Breast")
HML290.df = clin.df %>% 
  dplyr::filter(Cohort == "HML > 290") %>% 
  mutate(group = "HML > 290")
HML140290.df = clin.df %>% 
  dplyr::filter(Cohort %in% c("HML 140-290, other","HML 140-290, breast")) %>% 
  mutate(group = "HML < 290")

HML140290_breast.df = clin.df %>% 
  dplyr::filter(Cohort %in% c("HML 140-290, breast")) %>% 
  mutate(group = "HML 140-290, breast")
HML140290_other.df = clin.df %>% 
  dplyr::filter(Cohort %in% c("HML 140-290, other")) %>% 
  mutate(group = "HML 140-290, other")

if(GROUPS == "140290together"){
  combined.df = rbind(all.df, other.df, breast.df, HML140290.df, HML290.df)
  combined.df = combined.df %>% 
    mutate(group = factor(group, 
                          levels = c("All","Other","Breast",
                                     "HML < 290",
                                     "HML > 290")))
  anno_xpos=c(0.85, 1.15, 1.85, 2.15, 2.85, 3.15, 3.85, 4.15, 4.85, 5.15)
  anno_ypos = 9
}else if(GROUPS == "140290splitted"){
  combined.df = rbind(all.df, other.df, breast.df, HML140290_other.df, HML140290_breast.df, HML290.df)
  combined.df = combined.df %>% 
    mutate(group = factor(group, 
                          levels = c("All","Other","Breast",
                                     "HML 140-290, breast","HML 140-290, other",
                                     "HML > 290")))
  anno_xpos=c(0.85, 1.15, 1.85, 2.15, 2.85, 3.15, 3.85, 4.15, 4.85, 5.15, 5.85, 6.15)
  anno_ypos = 8
}



#Figures
tml.p=my_wilc_grouped(combined.df, VAR = "TML", ylab = "TML") + 
  theme(axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0(COHORT," & ",GROUPS, "(",RemoveLowerThan140, ")")) 
tmb.p=my_wilc_grouped(combined.df, VAR = "TMB", ylab = "TMB") + theme(axis.text.x=element_blank())

#cTML_maftools_05, cTML_maftools, cTML, cTML_05, perc_clon,
ctml.p=my_wilc_grouped(combined.df, VAR = "cTML_maftools_05", ylab = "Clonal TML") + theme(axis.text.x=element_blank())  
svl.p=my_wilc_grouped(combined.df, VAR = "svTumorMutationalBurden", ylab = "SVL") + theme(axis.text.x=element_blank())
uv.p=my_wilc_grouped(combined.df, VAR = "SBS7a", ylab = "UV sig") + theme(axis.text.x=element_blank())

count_df = combined.df %>%
  group_by(group, responders) %>%
  summarise(count = n())
aneup.p=my_wilc_grouped(combined.df, VAR = "aneuploidy_score", ylab = "Aneuploidy") +
  annotate("text",
           size = 3,
           x = anno_xpos,
           y = 70,
           label = count_df$count,
           vjust = anno_ypos)

pdf(file = paste0("/home/l.leek/pembro/results/fig_",COHORT,"_",GROUPS,"_",RemoveLowerThan140,"_TML_TMB_cTML_SV_Aneuploidy.pdf"),   
    width = 6, # The width of the plot in inches
    height = 10)
tml.p + tmb.p + ctml.p + uv.p+ svl.p  + aneup.p  + plot_layout(ncol=1)
dev.off()
tml.p + tmb.p + ctml.p + uv.p + svl.p  + aneup.p  + plot_layout(ncol=1)


exit


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
my_wilcoxon2(HML140290.df, "TML") #0.26
my_wilcoxon2(HML290.df, "TML") #0.0045
