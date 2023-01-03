library(patchwork)
library(data.table)
#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv")
perc_clon.df = clin.df %>% drop_na(perc_clon) 
perc_clon.df = perc_clon.df %>% 
  arrange(desc(perc_clon)) %>% 
  mutate(id = 1:nrow(perc_clon.df)) %>% 
  dplyr::select(id, perc_clon)


ggplot(perc_clon.df, aes(x=id, y=as.numeric(perc_clon))) +
  geom_point() +
  theme_bw(base_size = 12)+
  theme(legend.position = "None",
        axis.text.x = element_text(size = 12, angle = 33, vjust = 1, hjust=1)) +
  ylab("Clonality percentage") +
  xlab("Total WGS = 62")

