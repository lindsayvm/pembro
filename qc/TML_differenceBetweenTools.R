# Libraries
library(hrbrthemes)
library(GGally)
library(viridis)
library(plotly)

clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv")%>% 
  mutate(Cohort = as.factor(Cohort)) %>% 
  mutate(TML_SNPeff_log2 = log2(TML_SNPeff)) %>% 
  mutate(TML_HMF_log2 = log2(TML)) %>% 
  mutate(TML_DRUPreports_log2 = log2(MutationalLoad))
df = clin.df %>%
  select(c("TML_DRUPreports_log2", "TML_HMF_log2", "TML_SNPeff_log2", "Cohort"))  %>% 
  na.omit()

p= ggparcoord(df,
           columns = c("TML_DRUPreports_log2","TML_HMF_log2","TML_SNPeff_log2"), 
           groupColumn = "Cohort",
           order = "anyClass", 
           scale = "globalminmax",
           showPoints = TRUE, 
           alphaLines = 0.3) + 
  scale_color_viridis(discrete=TRUE) +
  theme_bw(base_size = 15)+
  theme(legend.position = "None",
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Log2 of TML") +
  geom_hline(yintercept = log2(140)) +
  geom_hline(yintercept = log2(290))

ggplotly(p)
idx = 62#38#,62
df[idx,c("TML_SNPeff", "TML", "MutationalLoad", "Cohort")]

clin.df[clin.df$TML == 34, c("bir_CPCT_WIDE_CORE","bir_HMFsampleID","TML_SNPeff", "TML", "MutationalLoad", "TML_DRUPreports_log2","TML_HMF_log2","TML_SNPeff_log2","Cohort") ]
clin.df[clin.df$TML == 364, c("bir_CPCT_WIDE_CORE","bir_HMFsampleID","TML_SNPeff", "TML", "MutationalLoad","TML_DRUPreports_log2","TML_HMF_log2","TML_SNPeff_log2", "Cohort") ]
