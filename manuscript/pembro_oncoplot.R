library(forcats)
library(ggplot2)
library(data.table)
library(dplyr)

raw.df = fread("data/20230503_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F) %>% 
  dplyr::select(gene_snpeff_B2M,gene_BRAF)

raw.df = fread("data/20230503_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F) %>% 
  tidyr::drop_na(gene_snpeff_B2M) 

tml.df = raw.df %>% 
  dplyr::select(subjectkey,Cohort_trueHML, TML) %>% 
  arrange(desc(TML))

cohorts.df = tml.df %>% dplyr::select(subjectkey,Cohort_trueHML)

var.m = raw.df %>% dplyr::select(subjectkey,starts_with("gene_snpeff_")) %>% 
  dplyr::select(subjectkey,!starts_with("gene_snpeff_APOBEC")) 
colnames(var.m) = gsub("gene_snpeff_","",colnames(var.m))


# Long format
var.df = melt(var.m)
var.df = left_join(var.df, cohorts.df, by = "subjectkey")
colnames(var.df) = c("Patients", "y", "value", "cohortGroup")

IFNg = c("JAK1", "JAK2", "STAT1", "IFNGR1", "IFNGR2")
AntigenPresentation = c("B2M", "TAP1", "TAP2", "TAPB", "CALR", "PDIA3", "CNAX", "HSPA5", "KEAP", "STK11")
HLA = c("HLA-A", "HLA-B","HLA-C")
SWISNF = c("PBRM1", "ARID1A", "ARID1B","ARID2", "SMARCB1", "SMARC4A","ATRX")
Other = c( "PIK3CA", "PTEN", "BRAF", "KRAS", "NRAS", "ESR1",  "ESR2")

#' ###########################################################################
#' ###########################################################################
#' Process:  
#' ###########################################################################
#' ###########################################################################

p0 = ggplot(data = tml.df, aes(x = factor(subjectkey, levels = tml.df$subjectkey), y = TML)) + 
  geom_bar(stat = "identity",width=0.8, position = position_dodge(width=0.2)) +
  
  facet_grid(cols = vars(Cohort_trueHML),
             scales = "free",space = "free") +
  theme_classic()+   
  theme(#axis.text.y = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, angle=90),
    strip.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 0.2),
    axis.text.x = element_blank(),
    axis.line.x=element_line(size = 0.2),
    axis.line.y=element_line(size = 0.2),
    strip.text.x = element_text(size=10, angle=0))+
  labs(y = "TML")

#Plot tile
var.df = var.df %>% 
  mutate(groups = forcats::fct_collapse(y, 
                                        IFNg=IFNg,
                                        AntigenPresentation=AntigenPresentation,                      
                                        HLA=HLA,
                                        SWISNF=SWISNF,
                                        Other=Other)) %>% 
  #order
  mutate(groups = factor(groups, levels = c("IFNg", "AntigenPresentation", "HLA", "SWISNF","Other")))

p1= ggplot(var.df, 
           aes(x = factor(Patients, levels = tml.df$subjectkey), 
               y = factor(y, levels = c(rev(IFNg),
                                        rev(AntigenPresentation),                      
                                        rev(HLA),
                                        rev(SWISNF),
                                        rev(Other))), fill = value))+
  
  geom_tile(color = "white",
            lwd = 0.1,
            linetype = 1) +
  
  scale_fill_gradient(
    low = "skyblue", 
    high = "seagreen",
    na.value = "lightgrey") +
  
  facet_grid(rows = vars(groups),
             cols = vars(cohortGroup),
             scales = "free",space = "free", 
             switch = "y") + 

    theme(legend.position = "none",
        axis.title = element_blank(),       # Change y axis title only
        axis.text.x=element_blank(),
        axis.text=element_text(size=12),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.placement = "outside",
        axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))


var.df2 = raw.df %>%   tidyr::drop_na(gene_snpeff_B2M) %>% dplyr::select(subjectkey,responders) %>% 
  mutate(y = "responders") %>% 
  mutate(value = ifelse(responders=="R",1,0)) %>% 
  dplyr::select(-responders)
var.df2 = left_join(var.df2, cohorts.df, by = "subjectkey")
colnames(var.df2) = c("Patients", "y", "value", "cohortGroup")
length(unique(var.df$Patients))
my_y_title = expression(paste(italic(n[total]), " = ", 55))
#plot
p2 =  ggplot(var.df2, 
             aes(x = factor(Patients, levels = tml.df$subjectkey), 
                 y = y, 
                 fill = factor(value))) +
  
  geom_tile(color = "white",
            lwd = 0.1,
            linetype = 1) +
  
  scale_fill_manual(values = c("salmon", "seagreen"),na.value="#EEEEEE", name = "CB") +
  
  facet_grid(cols = vars(cohortGroup),
             scales = "free",space = "free") +
  
  theme_void()+   
  theme(#axis.text.y = element_blank(),
    legend.position = "none",
    axis.title.x=element_text(size=14,face="bold"),
    strip.text.x = element_blank()) + 
  labs(x = my_y_title)

p0+p1+p2 + plot_layout(ncol=1, heights = c(3,30,1))
