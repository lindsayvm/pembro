library(patchwork)
library(data.table)
#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv") #%>% 
  #filter(TumorType == "Breast cancer")
p1= clin.df[clin.df$patientID=="CPCT02260029", c("patientID", "responders","TumorType","gene_PBRM1")]
p2= clin.df[clin.df$patientID=="DRUP01220002", c("patientID", "responders","TumorType","gene_PBRM1")] #CPCT02220026
p2_cpct=fread("/DATA/share/Voesties/data/HMF/update_10/somatics/CPCT02220026T/linx/CPCT02220026T.linx.driver.catalog.tsv")
p2_drup=fread("/DATA/share/Voesties/data/DRUP/update_3/somatics/DRUP01220002T/linx/DRUP01220002T.linx.driver.catalog.tsv")


genes.v = colnames(clin.df)[str_detect(colnames(clin.df),"gene_snpeff")]
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
    label=gsub("gene_snpeff_","",final.df$gene), 
    size = 3,
    aes(colour = factor(gene)),
    nudge_x = 0.25, nudge_y = 0.25, angle = 33,
    check_overlap = F) +
  theme_bw(base_size = 12)+
  theme(legend.position = "None",
        axis.text.x = element_text(size = 12, angle = 33, vjust = 1, hjust=1)) +
  xlab("number of patients with mutation (total WGS = 62)")
#dev.off()

final.df = clin.df %>%  select(c(bir_CPCT_WIDE_CORE,gene_snpeff_STK11, gene_snpeff_ARID1A, gene_snpeff_BRAF, gene_snpeff_IFNGR2, gene_snpeff_PBRM1, BOR, responders, TumorType)) %>% 
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
