# STATS TABLE

library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(ggpubr)

#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv", data.table = F)



names1 = colnames(clin.df)[str_detect(colnames(clin.df),pattern = "rna_|gene|SBS")]
names2 = c("purity"                    ,                "normFactor"          ,                      "score"    ,                                
 "diploidProportion"         ,                "ploidy"              ,                      "sex"      ,                                
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
 "cTML"                           ,           "perc_clon"                   )     

variables = c(names2, names1)

#Framework for output
data = data.frame(n_tot=NA, 
                  wilcoxon_responders=NA, 
                  wilcoxon_PDvsPR=NA, 
                  wilcoxon_PDvsSD=NA,
                  anova_BOR=NA, 
                  fisher_responders=NA, fish_OR=NA, fish_conf.low=NA,   fish_conf.high=NA  ) 



for (i in 1:length(variables)){   

  varname = variables[i]
  var = clin.df[[varname]]
  
  #print(table(var))
  var_complete = var[!is.na(var)]
  n_tot = length(var_complete)
  
  # skip if there are no groups possible because there is just 1 output
  if (length(unique(var[!is.na(var)])) <= 1){
      print(paste("skipped",varname))
      
      new_row = c(n_tot, NA, NA, NA, NA, NA, NA, NA, NA ) 
      data[i, ] = new_row            
      rownames(data)[i] = paste0(varname) 
      
      next
  }
  
  #FISHER
  if(is.logical(clin.df[[varname]]) || str_detect(varname, "_bool$|gene_|wholeGenomeDuplication|sex") ){
    # categorical data
    dat = table(var, clin.df$responders) #useNA = "ifany"
    
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
  if(is.logical(clin.df[[varname]]) || str_detect(varname, "_bool$|gene_|wholeGenomeDuplication|sex") ){
    wilcoxon_responders=NA
    wilcoxon_PDvsPR=NA
    wilcoxon_PDvsSD=NA
    anova_BOR=NA
  } else {
    wilc.test = pairwise.wilcox.test(clin.df[[varname]], clin.df$responders,
                                     p.adjust.method="none")
    wilcoxon_responders = signif(wilc.test$p.value,2)
    
    PDvsPR.df = clin.df %>% filter(BOR != "SD")
    wilc.test_PDvsPR = pairwise.wilcox.test(PDvsPR.df[[varname]], PDvsPR.df$responders,
                                     p.adjust.method="none")
    wilcoxon_PDvsPR = signif(wilc.test_PDvsPR$p.value,2)
    
    PDvsSD.df = clin.df %>% filter(BOR != "PR")
    wilc.test_PDvsSD = pairwise.wilcox.test(PDvsSD.df[[varname]], PDvsSD.df$responders,
                                            p.adjust.method="none")
    wilcoxon_PDvsSD = signif(wilc.test_PDvsSD$p.value,2)
    
    aov.test = aov(clin.df[[varname]] ~ BOR, data = clin.df)
    anova_BOR = data.frame(signif(unlist(summary(aov.test))["Pr(>F)1"], 4))[[1]]
    
  }
  
  
  new_row = c(n_tot, wilcoxon_responders, wilcoxon_PDvsPR, wilcoxon_PDvsSD, anova_BOR, fisher_responders, fish_OR, fish_conf.low, fish_conf.high ) 
  data[i, ] = new_row            
  rownames(data)[i] = paste0(varname) 
  
}

#many non-responders are breast cancer, and these are mostly women
table(clin.df$sex, clin.df$responders)

#rna_sig_signature_batf3, rna_sig_krt,rna_sig_mesench SBS_7_UV_sum, rna_sig_ESTROGEN_RESPONSE_EARLY, rna_sig_MYC_TARGETS_V2
var = "rna_sig_krt"

clin.df = clin.df %>%  filter(TumorType != "Breast cancer")

wilc.p = my_wilcoxon(clin.df, var)
wilc_PD_PR.p = my_wilcoxon_plot2(clin.df %>% filter(BOR != "SD"), var)
anova.p = my_anova_plot(clin.df, var)

wilc.p+wilc_PD_PR.p+ theme(axis.title.y = element_blank()) + 
  anova.p + theme(axis.title.y = element_blank())



