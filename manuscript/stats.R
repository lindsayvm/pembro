# STATS TABLE

library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(ggpubr)
library(dplyr)
library(broom)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")

#load data
clin.df = fread("data/20230503_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F)
##clin.df$Cohort #initial DRUP admission (140-290 and >290 and seperated breast cohort); were not gonna use this
clin.df$Cohort_trueHML #TRUE TML: 140-290 and >290 (incl breast in >290)


missing_count = clin.df %>% 
  summarize_all(~ sum(is.na(.)))

missing_count_tidy <- missing_count %>% 
  tidyr::gather(variable, missing_count) %>%
  arrange(missing_count)


# Create the barplot using ggplot2
ggplot(missing_count_tidy  %>% dplyr::filter(!str_detect(missing_count_tidy$variable,"hmf_|DRUP_ID|Tumor_Sample_Barcode"),), aes(x = reorder(variable, missing_count), y = missing_count)) +
  geom_bar(stat = "identity",width=0.6, position = position_dodge(width=0.2))  +
  xlab("Variable") +
  ylab("Number of missing values") +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = missing_count), vjust = -0.5, size = 2) + 
  theme_classic()+   
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))



cohorts = c("All",
            "Other",  
            "Breast", 
            "HML 140-290, other",
            "HML 140-290, breast",
            "HML < 290",
            "HML > 290")

#select vars of interest
names1 = colnames(clin.df)[str_detect(colnames(clin.df),pattern = "rna_|gene|SBS")]
names2 = c("purity"                    ,     "normFactor"          ,                      "score"    ,                                
           "diploidProportion"         ,               "ploidy"              ,                      "sex"      ,                                
           "polyclonalProportion" ,                     "minPurity" ,                               
           "maxPurity"                    ,             "minPloidy"                ,                 "maxPloidy"   ,                             
           "minDiploidProportion"          ,            "maxDiploidProportion"   ,                                           
           "somaticPenalty"                 ,           "wholeGenomeDuplication"  ,                  "msIndelsPerMb" , 
           "MSI_bool"                       ,           "TML"                      ,                 "TML_bool"      ,                           
           "TMB"                             ,          "TMB_bool"                ,                  "svTumorMutationalBurden"     ,             
           "TML_SNPeff"    , 
           "neo_TMB"                            ,       "frameshift"            ,                    "missense"                    ,             
           "neo_svTumorMutationalBurden"          ,     "clonal_tmb_dbs"          ,                  "clonal_tmb_indel"              ,           
           "clonal_tmb_snv"                       ,     "subclonal_tmb_dbs"       ,                  "subclonal_tmb_indel"           ,           
           "subclonal_tmb_snv"           ,              "total_fusions"              ,               "total_neo"                       ,         
           "clonal_neo"                  ,              "subclonal_neo"            ,                 "fusion_neo"                      ,         
           "mut_neo"                      ,             "clonal_tmb"                 ,               "subclonal_tmb"                    ,        
           "snvs"                          ,            "indels"                  ,                  "dbs"            ,
           "aneuploidy_score", "avg_cnv", "aneuploidyScore",
           "cTML",           "perc_clon"     , "cTML_05", "perc_clon_05", 
           "cTML_maftools", "cTML_maftools_05" )     
variables = c(names2, names1)

#Framework for output
stats.ls = list()
cohort.ls = list()


for(cohort in cohorts){
  print(cohort)

  if(cohort == "All"){
    clin_cohort.df = clin.df
  }else if(cohort == "Other"){
    clin_cohort.df = clin.df %>% dplyr::filter(TumorType != "Breast")
  }else if(cohort == "Breast"){
    clin_cohort.df = clin.df %>% dplyr::filter(TumorType == "Breast")
  }else if(cohort == "HML 140-290, other"){
    clin_cohort.df = clin.df %>% dplyr::filter(Cohort_trueHML == "HML 140-290, other")
  }else if(cohort == "HML 140-290, breast"){
    clin_cohort.df = clin.df %>% dplyr::filter(Cohort_trueHML == "HML 140-290, breast")
  }else if(cohort == "HML > 290"){
    clin_cohort.df = clin.df %>% dplyr::filter(Cohort_trueHML == "HML > 290")
  }else if(cohort == "HML < 290"){
    clin_cohort.df = clin.df %>% dplyr::filter(Cohort_trueHML %in% c("HML 140-290, other","HML 140-290, breast"))
  }
    
  data = data.frame(varname=NA,
                    n_tot=NA, 
                    wilcoxon_responders=NA, 
                    mean_responder=NA,
                    mean_nonresponders=NA,
                    wilcoxon_PDvsPR=NA, 
                    wilcoxon_PDvsSD=NA,
                    anova_BOR=NA, 
                    glm_pval=NA,glm_pval_correctedByTMB=NA,
                    fisher_responders=NA, fish_OR=NA, fish_conf.low=NA,   fish_conf.high=NA  ) 
  
  for (i in 1:length(variables)){   
    
    varname = variables[i]
    var = clin_cohort.df[[varname]]
    
    #print(varname)
    
    var_complete = var[!is.na(var)]
    n_tot = length(var_complete)
    
    # skip if there are no groups possible because there is just 1 output
    if (length(unique(var[!is.na(var)])) <= 1){
      print(paste("     skipped",varname))
      
      new_row = c(varname, n_tot, NA, NA, NA, NA, NA, NA, NA, NA,NA,NA,NA,NA ) 
      data[i, ] = new_row            
      rownames(data)[i] = paste0(varname) 
      
      next
    }
    
    #FISHER
    if(is.logical(clin_cohort.df[[varname]]) || str_detect(varname, "_bool$|gene_|wholeGenomeDuplication|sex") ){
      # categorical data
      dat = table(var, clin_cohort.df$responders) #useNA = "ifany"
      
      fish.test = fisher.test(dat)
      fisher_responders = signif(fish.test$p.value,2) 
      fish_OR = fish.test$estimate
      fish_conf.low = fish.test$conf.int[1]
      fish_conf.high = fish.test$conf.int[2]
      
    } else {
      #continuous data
      fisher_responders = NA
      fish_OR=NA
      fish_conf.low=NA
      fish_conf.high=NA
    }
    
    #Wilcoxon
    if(is.logical(clin_cohort.df[[varname]]) || str_detect(varname, "_bool$|gene_|wholeGenomeDuplication|sex") ){
      wilcoxon_responders=NA
      wilcoxon_PDvsPR=NA
      wilcoxon_PDvsSD=NA
      anova_BOR=NA
      mean_responders=NA
      mean_nonresponders=NA
      glm_pval=NA
      glm_pval_correctedByTMB=NA
    } else {
      wilc.test = pairwise.wilcox.test(clin_cohort.df[[varname]], clin_cohort.df$responders,
                                       p.adjust.method="none")
      wilcoxon_responders = signif(wilc.test$p.value,2)
      
      responders = clin_cohort.df %>% dplyr::filter(responders == "R")
      nr = clin_cohort.df %>% dplyr::filter(responders == "NR")
      mean_responders = mean(responders[[varname]], na.rm = T)
      mean_nonresponders = mean(nr[[varname]], na.rm = T)
      
      
      PDvsPR.df = clin_cohort.df %>% dplyr::filter(BOR != "SD") %>% tidyr::drop_na(varname)
      if(length(unique(PDvsPR.df$responders)) == 1){
        wilc.test_PDvsPR = NA
      }else{
        wilc.test_PDvsPR = pairwise.wilcox.test(PDvsPR.df[[varname]], PDvsPR.df$responders,
                                                p.adjust.method="none")
        wilcoxon_PDvsPR = signif(wilc.test_PDvsPR$p.value,2)
      }
      
      
      PDvsSD.df = clin_cohort.df %>% dplyr::filter(BOR != "PR") %>% tidyr::drop_na(varname) 
      if(length(unique(PDvsSD.df$responders)) == 1){
        wilcoxon_PDvsSD = NA
      }else{
        wilc.test_PDvsSD = pairwise.wilcox.test(PDvsSD.df[[varname]], PDvsSD.df$responders,
                                                p.adjust.method="none")
        wilcoxon_PDvsSD = signif(wilc.test_PDvsSD$p.value,2)
      }
      
      aov.test = aov(clin_cohort.df[[varname]] ~ BOR, data = clin_cohort.df)
      anova_BOR = data.frame(signif(unlist(summary(aov.test))["Pr(>F)1"], 4))[[1]]
      
      
      vars = c(varname)
      clin_cohort.df = clin_cohort.df %>% mutate(responders = as.factor(responders))
      fit = glm(reformulate(termlabels = paste(vars, "- 1"), response = "responders"), 
                data = clin_cohort.df,
                family = binomial(link = "logit"), 
                na.action = na.omit) 
      glm.df = tidy(fit, exponentiate = T, conf.int = T) %>% 
        mutate(OR = estimate) %>% 
        mutate(Biomarkers = reorder(term, estimate))
      glm_pval = glm.df[glm.df$term == varname, "p.value"][[1]]
      
      vars = c("TML_SNPeff",varname)
      fit = glm(reformulate(termlabels = paste(vars, "- 1"), response = "responders"), 
                data = clin_cohort.df,
                family = binomial(link = "logit"), 
                na.action = na.omit) 
      glm.df = tidy(fit, exponentiate = T, conf.int = T) %>% 
        mutate(OR = estimate) %>% 
        mutate(Biomarkers = reorder(term, estimate))
      glm_pval_correctedByTMB = glm.df[glm.df$term == varname, "p.value"][[1]]
    }
    
    new_row = c(varname,n_tot, wilcoxon_responders, mean_responders, mean_nonresponders, 
                wilcoxon_PDvsPR, wilcoxon_PDvsSD, anova_BOR, 
                glm_pval,glm_pval_correctedByTMB,
                fisher_responders, fish_OR, fish_conf.low, fish_conf.high ) 
    data[i, ] = new_row            
    rownames(data)[i] = paste0(varname) 

  }
  
  cohort.ls[[cohort]] = clin_cohort.df
  stats.ls[[cohort]] = data
}


cohorts
cohort.ls[["All"]]
all.stats = stats.ls[["All"]]

cohort.ls[["Other"]]
other.stats = stats.ls[["Other"]]

cohort.ls[["Breast"]]
breast.stats = stats.ls[["Breast"]]

cohort.ls[["HML 140-290, other"]]
other140290.stats = stats.ls[["HML 140-290, other"]]

cohort.ls[["HML 140-290, breast"]]
breast140290.stats = stats.ls[["HML 140-290, breast"]]

cohort.ls[["HML > 290"]]
HML290.stats = stats.ls[["HML > 290"]]

cohort.ls[["HML < 290"]]
HMLlower290.stats = stats.ls[["HML < 290"]]

tmp = cohort.ls[["All"]]
table(tmp %>% dplyr::select('gene_snpeff_HLA-A', responders))
#HLA-A, PBRM1, ARID1A, ARID1B, BRAF, STK11, SMARCB1, IFNGR2, APOBEC3G


####
VAR = "perc_clon"
cohort = "HML > 290"
pval_ypos = 110
perc_clonal_290.p = my_wilcoxon(clin.df %>% dplyr::filter(Cohort_trueHML == cohort), VAR ,pval_ypos = pval_ypos)+
  theme(axis.title.x = element_text()) + 
  xlab("")
perc_clonal_290.p
#
##
##
##
##
#

#many non-responders are breast cancer, and these are mostly women
table(clin.df$sex, clin.df$responders)

#rna_sig_signature_batf3, rna_sig_krt,rna_sig_mesench SBS_7_UV_sum, rna_sig_ESTROGEN_RESPONSE_EARLY, rna_sig_MYC_TARGETS_V2
var = "rna_sig_krt"


wilc.p = my_wilcoxon(clin.df, var)
wilc_PD_PR.p = my_wilcoxon_plot2(clin.df %>% filter(BOR != "SD"), var)
anova.p = my_anova_plot(clin.df, var)

wilc.p+wilc_PD_PR.p+ theme(axis.title.y = element_blank()) + 
  anova.p + theme(axis.title.y = element_blank())

write.table(x = data, file = "/home/l.leek/pembro/results/20221021_pembro_stats.csv", 
            quote = F, sep = "\t", col.names = T, row.names = T)
x=fread("/home/l.leek/pembro/results/20221021_pembro_stats.csv", data.table = F)

