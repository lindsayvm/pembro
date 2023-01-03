library(patchwork)
library(data.table)
#env
dir = "/home/l.leek/pembro/"
setwd(dir)
source("src/functions_plots.R")


clin.df = fread("data/20221021_DRUP_pembro_LL_WGS_RNA.tsv")


genes.v = colnames(clin.df)[str_detect(colnames(clin.df),"gene_")]
final.df = as.data.frame(matrix(rep(NA, 3*0), ncol=3))

for (gene in genes.v){
  complete.df = clin.df %>%  dplyr::select(gene, responders) %>%  na.omit()
  pval_fish = get_fisher_pval(df = complete.df, var = gene)
  n_mut = length(complete.df[[gene]][complete.df[[gene]] == 1])
  new_row = c(gene,round(pval_fish,4), n_mut)
  final.df = rbind(final.df, new_row)
}
colnames(final.df) = c("gene", "fisher", "n_patients_mut")


ggplot(final.df, aes(x=as.integer(n_patients_mut), y=fisher)) +
  geom_point() + 
  geom_text(
    label=gsub("gene_","",final.df$gene), 
    size = 3,
    aes(colour = factor(gene)),
    nudge_x = 0.25, nudge_y = 0.25, angle = 33,
    check_overlap = F) +
  theme_bw(base_size = 12)+
  theme(legend.position = "None",
        axis.text.x = element_text(size = 12, angle = 33, vjust = 1, hjust=1)) +
  xlab("number of patients with mutation (total WGS = 62)")

