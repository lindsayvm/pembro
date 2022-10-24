#libraries
library(data.table)
library(randomcoloR)
library(plotly)
library(dplyr)
library(patchwork)

dir = "/home/l.leek/"
setwd(dir)

source("src/pembro/functions_plots.R")

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

clin.df = fread("data/pembro/20221021_DRUP_pembro_LL_final.tsv", data.table = F)


#' ###########################################################################
#' ###########################################################################
#' process:  
#' ###########################################################################
#' ###########################################################################

#balanced cohorts
cohort.p = my_barplot(df = clin.df,
           inputCol = clin.df$Cohort,
           sort = T)

#almost half is breast
tumortype.p = my_barplot(df = clin.df,
                         inputCol = clin.df$TumorType,
                         sort = T)
#ggplotly(tumortype.p)

#subtypes: lots of tumors dont have subtype
subtype.p = my_barplot(df = clin.df,
                       inputCol = clin.df$Association,
                       sort = T)
subtype.p
#All breast are subtyped, all skin is UV typed; and those missing subtypes are derived from 
tmp = clin.df[is.na(clin.df$Association), ]
TumorWithoutSubtypesAnnotated.p = my_barplot(df = tmp ,
           inputCol = tmp$TumorType,
           sort = T)
  


#breast: Heterogeneous (but in context of response to immuno: TNBC is quit prevalent (n=8))
breast.df = clin.df[(clin.df$TumorType == "Breast cancer"), ]
breast.p = my_barplot(df = breast.df,
                      inputCol = breast.df$Association,
                      sort = T)
#ggplotly(breast.p)

#most (55) biopsies are pretreatment; WGS will be stable, but RNA might be influenced in those other 18 if pretreatment was immuno --> CHECK ??? 
table(clin.df$PretreatmentBiopsy)

#more females than males; generally have more active immune system.
table(clin.df$Gender)

#response per mutational load
tml_BOR.p = clin.df %>% 
  select(BOR, MutationalLoad) %>% 
  na.omit() %>% 
  ggplot(aes(x=MutationalLoad, fill=BOR)) +
  geom_density(alpha=0.55)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        legend.position = "bottom") +
    geom_vline(xintercept = 500,
             color = "black", lwd =1, lty = "dashed") +
    geom_vline(xintercept = 150,
             color = "black", lwd =1)

length(response.df$MutationalLoad[response.df$MutationalLoad > 500])

#response
response.df = clin.df %>% 
  mutate(value = rep(1, nrow(clin.df))) %>% 
  mutate(BOR = factor(clin.df$BOR, 
                      levels=c("PR","SD","PD"))) %>% 
  mutate(subtype = factor(clin.df$Association, 
                          levels = names(sort(table(clin.df$Association), decreasing=TRUE)))) %>% 
  mutate(TumorType = factor(clin.df$TumorType, 
                            levels = names(sort(table(clin.df$TumorType), decreasing=TRUE))))


tumortypeResponse.p = response.df %>%  
  select(BOR, value, TumorType) %>% 
  na.omit() %>% 
  ggplot( aes(fill=BOR, y=value, x=TumorType)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, direction = -1) +
  theme_minimal(base_size = 6)+
  theme(axis.text.x = element_text(size = 4, angle = 33, vjust = 0.5, hjust=1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 12)) +
  ylab("BOR")
tumortypeResponse.p
                    

# All 4 patients with UV signature respond relatively well
subtypeResponse.p = response.df %>%  
  select(BOR, value, subtype) %>% 
  na.omit() %>% 
  ggplot( aes(fill=BOR, y=value, x=subtype)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, direction = -1) +
  theme_minimal(base_size = 6)+
  theme(axis.text.x = element_text(size = 4, angle = 33, vjust = 0.5, hjust=1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 12)) +
  ylab("BOR")
subtypeResponse.p


#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################

cohort.p + tumortype.p + plot_layout(widths = c(1,2))

subtype.p + breast.p  + plot_layout(widths = c(2,1))

tml_BOR.p / 
tumortypeResponse.p /
subtypeResponse.p



