#Libraries
library(data.table)
library(ggplot2)
library(fastDummies)
library(vcfR)
library(stringr)
library(patchwork)

dir = "/home/l.leek/pembro"
setwd(dir)

source("src/1_feature_generation/WGS/0_functions_WGS.R")
source("src/functions_plots.R")

clin.df= fread("~/TONIC/TONIC_decodingfile_LOCKEDFORSHARING.csv")

#config
GOI = c("PIK3CA", "PTEN", "KRAS", "BRAF", "NRAS", "ARID1A",
        "SMARCB1", "SMARC4A", "PBRM1", "JAK1", "JAK2", "STAT1",
        "IFNGR1", "IFNGR2","B2M", "TAP1", "TAP2", "TAPB",
        "CALR", "PDIA3", "CNAX", "HSPA5", "KEAP", "STK11","TP53")
ann.fn = list.files(path = "~/TONIC/vcf_snpeffsift/",
                    pattern = "_ann_filt_oneLine.vcf",
                    full.names = TRUE,  
                    recursive = TRUE)

#Get snpeff annotated genes
ann.df = my_ann_snpeff(ann.fn, GOI)
ann_pp.df = my_ann_pp(ann.df, GOI, ann.fn)
ann_pp.df$'Identifier DNA tumor baseline (WES)' = paste0("CF",gsub("m2.*CF|_.*","",gsub("_vs_.*","",ann_pp.df$patientID)))
#There are duplicates with a diff identifier, these are removed. 
ann_pp.df$`Identifier DNA tumor baseline (WES)`[duplicated(ann_pp.df$`Identifier DNA tumor baseline (WES)`)]
ann_pp.df$`Identifier DNA tumor baseline (WES)`[str_detect(ann_pp.df$patientID, "5077")] = NA


clin.df = dplyr::left_join(clin.df, ann_pp.df, by = "Identifier DNA tumor baseline (WES)")
clin.df$`Best overall response iRECIST`[clin.df$`Best overall response iRECIST` == "Non CR/non PD (or SD > 24 weeks)"] = NA
clin.df = clin.df %>% mutate(responders = ifelse(`Best overall response iRECIST` %in% c("PR","CR"), "R","NR"))


genes.v = colnames(clin.df)[str_detect(colnames(clin.df),"gene")] #_snpeff
final.df = as.data.frame(matrix(rep(NA, 3*0), ncol=3))

for (gene in genes.v){
  complete.df = clin.df %>%  dplyr::select(gene, responders) %>%  na.omit()
  pval_fish = get_fisher_pval(df = complete.df, var = gene)
  n_mut = length(complete.df[[gene]][complete.df[[gene]] == 1])
  new_row = c(gene,round(pval_fish,4), n_mut)
  final.df = rbind(final.df, new_row)
}
colnames(final.df) = c("gene", "fisher", "n_patients_mut")

#png("/home/l.leek/pembro/results/fig_mutations.png",width=2000,height=2000, res=300)
ggplot(final.df, aes(x=as.integer(n_patients_mut), y=fisher)) +
  geom_point() + 
  geom_text(
    label=final.df$gene,#gsub("gene_snpeff_","",final.df$gene), 
    size = 3,
    aes(colour = factor(gene)),
    nudge_x = 0.25, nudge_y = 0.25, angle = 33,
    check_overlap = F) +
  theme_bw(base_size = 12)+
  theme(legend.position = "None",
        axis.text.x = element_text(size = 12, angle = 33, vjust = 1, hjust=1)) +
  xlab("number of patients with mutation (total WGS = 62)")
#dev.off()

final.df = clin.df %>% dplyr::select(c(bir_CPCT_WIDE_CORE,gene_snpeff_STK11, gene_snpeff_ARID1A, gene_snpeff_BRAF, gene_snpeff_IFNGR2, gene_snpeff_PBRM1, BOR, responders, TumorType)) %>% 
  na.omit()

table(final.df$gene_snpeff_PBRM1, final.df$responders)
table(final.df$gene_snpeff_ARID1A, final.df$responders)
#table(final.df$gene_snpeff_BRAF, final.df$responders)

swisli.df = clin.df %>% 
  select(c(gene_snpeff_ARID1A, gene_snpeff_SMARCB1, gene_snpeff_SMARC4A, gene_snpeff_PBRM1, BOR, responders, TumorType)) %>% 
  na.omit()

swisli_breast.df = clin.df %>% 
  filter(TumorType == "Breast cancer") %>% 
  select(c(bir_CPCT_WIDE_CORE,gene_ARID1A,gene_PBRM1, gene_snpeff_ARID1A, gene_snpeff_SMARCB1, gene_snpeff_SMARC4A, gene_snpeff_PBRM1, BOR, responders, TumorType)) %>% 
  na.omit()


table(swisli_breast.df$gene_snpeff_PBRM1, swisli_breast.df$responders)
table(swisli_breast.df$gene_snpeff_ARID1A, swisli_breast.df$responders)
#table(swisli_breast.df$gene_snpeff_BRAF, swisli_breast.df$responders)

