# Libraries
library(hrbrthemes)
library(GGally)
library(viridis)
library(plotly)

df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv") %>% 
  mutate(TMB_SNPeff = TMB_SNPeff_pass_nonsyn_protcoding/30) %>%   #30 MB
  mutate(TMB_SNPeff_log = log2(TMB_SNPeff)) %>% 
  mutate(TMB_log = log2(TMB)) %>% 
  mutate(Cohort = as.factor(Cohort)) %>% 
  dplyr::select(c("TMB", "TMB_SNPeff","TMB_SNPeff_log","TMB_log" ,"Cohort"))  %>% 
  na.omit()

p_TMB= ggparcoord(df,
                  columns = c("TMB","TMB_SNPeff"), 
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
  ylab("TMB")
ggplotly(p_TMB)


p_TMB_log= ggparcoord(df,
           columns = c("TMB_log","TMB_SNPeff_log"), 
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
  ylab("TMB")

ggplotly(p_TMB_log)
