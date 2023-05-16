#libraries
library(data.table)
library(dplyr)
source("src/functions_plots.R")

dir = "/home/l.leek/pembro/"
setwd(dir)

#Load data
clin.df = fread("data/20230503_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F) 
hla_raw.df = fread("/DATA/share/Voesties/data/harmonize/output/mhc_class_i/lilac-hla.csv", data.table = F)
#hla_status = fread("/DATA/share/Voesties/data/harmonize/output/mhc_class_i/lilac-hla-wide.tsv", data.table = F)

#Somatic Missense and SomaticNonsenseOrFrameshift are most prevalent
colSums(hla_raw.df %>% dplyr::select(SomaticMissense,SomaticNonsenseOrFrameshift,SomaticSplice,SomaticInframeIndel))

#Filter on pembro patients
hla.df = hla_raw.df %>% 
  mutate(patientID = gsub("T$|TI$|TI.*$","",sampleId)) %>% 
  dplyr::filter(patientID %in% clin.df$patientID)

#How many patients from pembro? 59 (not 62 for which we have WGS)
length(unique(hla.df$patientID))
colSums(hla.df %>% dplyr::select(TumorCopyNumber,SomaticMissense,SomaticNonsenseOrFrameshift,SomaticSplice,SomaticInframeIndel))

#Loss of mutation/splice
hla_filtered.df = hla.df[hla.df$TumorCopyNumber <= 0.5 |
                       hla.df$SomaticMissense > 0 |
                       hla.df$SomaticNonsenseOrFrameshift > 0 |
                       hla.df$SomaticSplice > 0 |
                       hla.df$SomaticInframeIndel > 0, ]  
#dont look at SomaticSynonymous because that wont give a change in peptides
#Can also try 0.5 instead of 0 because most of them go from 0 to 2 alleles; so 0.5 should indicate that at least 1 allele has a mutation

#Only a third of the patients has a mutation
length(unique(hla_filtered.df$patientID)) #only 22


#Is there HLA mut?
clin.df = clin.df %>% mutate(HLA_LILAC = factor(ifelse(patientID %in% hla_filtered.df$patientID, 1, 0)))
#How many HLA mutations are there?
clin.df = clin.df %>% mutate(HLA_LILAC_counts = table(factor(hla_filtered.df$patientID, levels = clin.df$patientID))[as.character(clin.df$patientID)])
#Maybe only HLA-B relevant?
hla_B_filtered.df = hla_filtered.df %>%  
  mutate(Locus = substr(Allele, 1, 1)) %>% 
  dplyr::filter(Locus == "B")
clin.df = clin.df %>% mutate(HLA_B_LILAC = factor(ifelse(patientID %in% hla_B_filtered.df$patientID, 1, 0)))
clin.df = clin.df %>% mutate(HLA_B_LILAC_counts = table(factor(hla_B_filtered.df$patientID, levels = clin.df$patientID))[as.character(clin.df$patientID)])

#####################
######RESULTS########
#####################

#barplots for different cohorts with overlapping patients
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
HML140290_breast.df = clin.df %>% 
  dplyr::filter(Cohort %in% c("HML 140-290, breast")) %>% 
  mutate(group = "HML 140-290, breast")
HML140290_other.df = clin.df %>% 
  dplyr::filter(Cohort %in% c("HML 140-290, other")) %>% 
  mutate(group = "HML 140-290, other")

combined.df = rbind(all.df, other.df, breast.df, HML140290_other.df, HML140290_breast.df, HML290.df)
combined.df = combined.df %>% 
  mutate(group = factor(group, 
                        levels = c("All","Other","Breast",
                                   "HML 140-290, breast","HML 140-290, other",
                                   "HML > 290")))

#Are mutations in HLA class I genes associated to high TML? Little bit but not sign
p1 = ggplot(combined.df, aes(x = group, y = TML, fill = HLA_LILAC)) +
  geom_boxplot() +
  scale_fill_manual(values = c("brown", "forestgreen")) +
  geom_jitter(color="black", size=0.4, alpha=0.9,width = 0.25) +
  theme_bw(base_size = 15)+
  theme(legend.position = "None",
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust=1),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  stat_compare_means(aes(group = HLA_LILAC), method = "wilcox.test",
                     label = "p.format",
                     label.y = max(combined.df$TML, na.rm = T) + max(combined.df$TML, na.rm = T)*0.1) + #size = 6
  scale_y_continuous(expand=expansion(mult = c(0.1,0.2)))
p1



# #Maybe only allele B relevant?
# p1B = ggplot(combined.df, aes(x = group, y = TML, fill = HLA_B_LILAC)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("brown", "forestgreen")) +
#   geom_jitter(color="black", size=0.4, alpha=0.9,width = 0.25) +
#   theme_bw(base_size = 15)+
#   theme(legend.position = "None",
#         axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
#         axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust=1),
#         axis.title.x = element_blank(),
#         axis.ticks = element_blank()) +
#   stat_compare_means(aes(group = HLA_B_LILAC), method = "wilcox.test",
#                      label = "p.format",
#                      label.y = max(combined.df$TML, na.rm = T) + max(combined.df$TML, na.rm = T)*0.1) + #size = 6
#   scale_y_continuous(expand=expansion(mult = c(0.1,0.2)))
# p1B




# #Not grouped
# wilc.test = pairwise.wilcox.test(clin.df$TML, clin.df$HLA_B_LILAC,
#                                  p.adjust.method="none")
# pvalue = data.frame(signif(wilc.test$p.value,4))
# p1Xtra = clin.df %>% ggplot(aes(x=HLA_LILAC,
#                             y=log2(TML),
#                             fill=HLA_LILAC)) +
#   geom_violin(width = 0.6) +
#   geom_jitter(color="black", size=0.5, alpha=0.3, width = 0.2) +
#   ylab("TML") +
#   theme_bw(base_size = 15)+
#   theme(legend.position = "None",
#         axis.text.x = element_text(size = 17, angle = 0, vjust = 1, hjust=1),
#         axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust=1)) +
#   geom_text(data = pvalue,
#             aes(x = 1.5, y = max(log2(clin.df$TML), na.rm = T), label = paste0("pval_wilc=",pvalue)),
#             inherit.aes = FALSE, hjust = "inward", vjust = "inward", size = 3.5)  +
#   scale_fill_manual(values = c("forestgreen", "brown"))
# p1Xtra




#Are MORE mutations in HLA class I genes associated to higher TML?
p2 = ggplot(clin.df, aes(x=HLA_LILAC_counts, y=log2(TML), color=Cohort_trueHML)) + 
  geom_point(size=1) +
  scale_color_manual(values = c("brown",  "orange", "skyblue" ,"grey")) + 
  theme_bw(base_size = 15) +
  theme(legend.position = "none") 
p2



# p2B = ggplot(clin.df, aes(x=HLA_B_LILAC_counts, y=log2(TML), color=Cohort_trueHML)) + 
#   geom_point(size=1) +
#   scale_color_manual(values = c("brown",  "orange", "skyblue" ,"grey")) + 
#   theme_bw(base_size = 15) +
#   theme(legend.position = "none") 
# p2B




#Do NR have more HLA muts in NR - yes but not sign
p3 = clin.df %>% ggplot(aes(x=responders, 
                           y=HLA_LILAC_counts, 
                           fill=factor(responders, levels = c("R","NR")))) +
  geom_violin(width = 0.6) +
  geom_jitter(color="black", size=0.5, alpha=0.3, width = 0.2) +
  theme_bw(base_size = 15)+
  theme(legend.position = "None",
        axis.text.x = element_text(size = 17, angle = 0, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust=1),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(data = pvalue,
            aes(x = 1.5, y = max(clin.df$HLA_LILAC_counts, na.rm = T), label = paste0("pval_wilc=",pvalue)),
            inherit.aes = FALSE, hjust = "inward", vjust = "inward", size = 3.5)  +
  scale_fill_manual(values = c("forestgreen", "brown")) 
p3



# p3B = clin.df %>% ggplot(aes(x=responders, 
#                             y=HLA_B_LILAC_counts, 
#                             fill=factor(responders, levels = c("R","NR")))) +
#   geom_violin(width = 0.6) +
#   geom_jitter(color="black", size=0.5, alpha=0.3, width = 0.2) +
#   theme_bw(base_size = 15)+
#   theme(legend.position = "None",
#         axis.text.x = element_text(size = 17, angle = 0, vjust = 1, hjust=1),
#         axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust=1),
#         axis.title.x = element_blank(),
#         axis.ticks = element_blank()) +
#   geom_text(data = pvalue,
#             aes(x = 1.5, y = max(clin.df$HLA_B_LILAC_counts, na.rm = T), label = paste0("pval_wilc=",pvalue)),
#             inherit.aes = FALSE, hjust = "inward", vjust = "inward", size = 3.5)  +
#   scale_fill_manual(values = c("forestgreen", "brown")) 
# p3B




