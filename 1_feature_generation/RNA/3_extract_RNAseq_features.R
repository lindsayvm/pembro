#' Author: Lindsay Leek
#' Project started: begin 2021
#' Most recent update: 19/12/2021
#' Project: Neoadjuvant immunotherapy AI
#' RQ: Can we integrate different data modalities (CT, H&E, RNAseq, WES) to 
#'    predict which patients do not benefit from neoadjuvant immunotherapy?
#' Aim: add RNAseq data; high level mechanisms


#Libraries
library(stringr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(plotly)
library(dplyr)
library(plyr)
require(gridExtra)
library(survminer)
library(survival)
library(randomcoloR)
library(org.Hs.eg.db) 
library(RColorBrewer)

dir = "/home/l.leek/"
setwd(dir)

source("src/pembro/functions_plots.R")
source("src/pembro/functions.R")

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

rna.df = fread("/home/l.leek/data/pembro/pembro_RNAseq_normalized_cts.tsv", data.table = F)
colnames(rna.df)

clin.df = fread("data/pembro/20221021_DRUP_pembro_LL_final_1_WGS.tsv", data.table = F)
clin.df$patientID = paste0("ID_",gsub("CPCT0|DRUP0|WIDE0|CORELR0|COREDB0","",clin.df$CPCT_WIDE_CORE))
clin.df = clin.df[clin.df$patientID %in% colnames(rna.df), ]

#check whether all RNA can be found in clin.df
colnames(rna.df)[colnames(rna.df) %in% clin.df$patientID ]

#' ###########################################################################
#' ###########################################################################
#' Process: 
#' ###########################################################################
#' ###########################################################################


genes.v = c("CD274","TP53","VEGFA", "ALB", "ACTA1",
            "KRT10", "KRT14", "KRT17", "KRT4" , "KRT5" , "KRT6A" ,"KRT6B", "KRT80","SBSN","DSG3","SFN",
            "MYH6", "ACTC1","MB","CLDN11","TNNC1","BEST3","XIRP2","MYOG","MYH7",
            "CYP21A2", "CYP2C8", "CYP4A11", "CYP4F11", "CYP4F2", "CYP8B1", 
            "SLC17A2", "SLC22A7", "SLC2A2")



for (i in genes.v){
  
  #counts
  gene.df = rna.df[rna.df$gene == i, ]
  
  # transpose all but the first column (name)
  gene.df = as.data.frame(t(gene.df[,-1]))
  
  #assign colname gene with rna prefix
  colnames(gene.df) = paste0("rna_", i)
  gene.df$patientID = rownames(gene.df)
  
  #merge to main clinical file
  clin.df = join(clin.df, gene.df, by = "patientID")
}



# Immune checkpoints
ic.v = c("CD27", "CD28", "CD40", "CD122", "CD137", "CD134", "GITR", "ICOS", "A2AR", "CD272", "CD152", "IDO", "LAG3", "NOX2")

#retrieve ENSEMBL symbols for those that can be retrieved through os Hs eg db package
symbols = mapIds(org.Hs.eg.db, keys = ic.v, keytype = "SYMBOL", column="ENSEMBL")
ensembl.df = data.frame(ic.v, symbols)
ensembl.v = ensembl.df$symbols[!is.na(ensembl.df$symbols)]

#add missing
ensemble2.v = c("ENSG00000276977","ENSG00000120217", "ENSG00000197646", "ENSG00000163599", "ENSG00000135077")

#merge
ensembl.v = c(ensembl.v, ensemble2.v)


#give correct gene names from ENSEMBL
symbols = mapIds(org.Hs.eg.db, keys = ensembl.v, keytype = "ENSEMBL", column="SYMBOL")
ensembl.df = data.frame(ensembl.v, symbols)

#select IC genes from rna df
ic.df = rna.df[rna.df$gene %in% ensembl.df$symbols, ]

#make heatmap
ic.m = as.matrix(ic.df %>% dplyr::select(-gene))
my_heatmap(clin.df, ic.df, ic.m, "BOR" )


##### SIGNATURES
sig.fn = list.files(path = "/home/l.leek/data/tango/signatures",
                    pattern = "^HALLMARK", full.names = TRUE, recursive = TRUE)
names = gsub(".*\\/|.txt$", "", sig.fn)

# read as list
sig.ls = lapply(sig.fn, function(x) {
  df = read.delim(x, col.names = "gene", skip = 1)
  return(df)})

#Retrieve IDs
names(sig.ls) = names


#Additional signatures from literature
sig2.fn = list.files(path = "/home/l.leek/data/tango/signatures",
                     pattern = "^signature", full.names = TRUE, recursive = TRUE)
names2 = gsub(".*\\/|.txt$", "", sig2.fn)

# read as list
sig.ls2 = lapply(sig2.fn, function(x) {
  df = read.delim(x,col.names = "gene")
  return(df)})
#Retrieve IDs
names(sig.ls2) = names2


#append lists
sig.ls = append(sig.ls, sig.ls2)


#starting df: rna counts for genes in first signature
finalSum.df = rna.df[rna.df$gene %in% sig.ls[[1]]$gene, ]
#make matrix (numeric)
finalSum.m = finalSum.df %>% dplyr::select(-gene)
finalSum.m = sapply(finalSum.m, as.numeric)
#calc sum of each gene across all patients
finalSum.df = as.data.frame(colSums(finalSum.m))
#name of signature
colnames(finalSum.df) = c(names(sig.ls)[1])

#add all other signatures
for (i in 2:length(sig.ls)){
  print(i)
  
  tmp.df = rna.df[rna.df$gene %in% sig.ls[[i]]$gene, ]
  finalSum.m = tmp.df%>% dplyr::select(-gene)
  finalSum.m = sapply(finalSum.m, as.numeric)
  sum.df = as.data.frame(colSums(finalSum.m))
  
  #merge
  finalSum.df = cbind(finalSum.df, sum.df)
  
  #add name of signature
  colnames(finalSum.df)[length(finalSum.df)] = names(sig.ls)[i]
}

colnames(finalSum.df) = paste0("rna_sig_", gsub("HALLMARK_","",names(sig.ls)) )

# transpose
finalSum.m = as.matrix((finalSum.df))
finalSum.m = t(finalSum.m)

#make heatmap
my_heatmap2(clin.df, finalSum.m, "BOR")

#as data frame
finalSum.df$patientID = rownames(finalSum.df)
cols = colnames(finalSum.df)[-ncol(finalSum.df)]

#Merge
clin.df = merge(clin.df, finalSum.df, by = "patientID")

clin.df$BOR = factor(clin.df$BOR)


######## IMMUNE CELL TYPES

#60 markers for immune cell types proposed by Danaher et al.
#tumor infiltrating immune cells
tii.df = fread("/home/l.leek/data/tango/signatures/danahercelltypes.csv")
colnames(tii.df) = c("gene", colnames(tii.df)[-1])
tii_rna.df = rna.df[rna.df$gene %in% tii.df$gene, ]
tii.df = merge(tii.df, tii_rna.df, by = "gene")


tii_prognosis = function(tii.df, PROGNOSIS) {
  #' leave out DC, Treg, Mast cells 
  #' because they didn’t correlate to CD45 (leukocyte common antigen (LCA)) in Danaher
  if (PROGNOSIS == "goodprognosis"){
    rm.v = c("DC", "Neutrophils", "Macrophages", "Mast cells", "CD45")
    tii.df = tii.df[!tii.df$`Cell Type` %in% rm.v, ]
    
    #CD45 is a leukocyte common antigen (LCA), more like a summary statistic
  }else if (PROGNOSIS == "poorprognosis"){
    keep.v = c("DC", "Neutrophils", "Macrophages")
    tii.df = tii.df[tii.df$`Cell Type` %in% keep.v, ]
  }
  
  tii_counts.df = tii.df[ ,-c("gene", "Cell Type")]  
  
  tii.v = unique(tii.df$`Cell Type`)
  
  #make df
  ncol = ncol(tii_counts.df)
  nrow = 1
  final.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
  colnames(final.df) = colnames(tii_counts.df)
  rownames(final.df) = 'Celltype'
  
  for (i in tii.v){
    #select rows from matrix based on rows of df with clinical data that match cell type
    tmp.df = tii_counts.df[tii.df$`Cell Type` == i, ]
    tmpSum.df = as.data.frame(colSums(tmp.df))
    colnames(tmpSum.df) = i
    tmpSum.df = t(tmpSum.df)
    final.df = rbind(tmpSum.df, final.df)
  }
  colnames(final.df)
  final.df = final.df[-nrow(final.df), ]
  final.df$celltype = rownames(final.df) 
  m = as.matrix(final.df[ ,-ncol(final.df)])
  
  return(m)
}

# tumor infiltrating lymphocytes for good prognosis
til.m = tii_prognosis(tii.df, "goodprognosis")
#make heatmap
my_heatmap2(clin.df, til.m, "BOR")

#tumor infiltrating immune cells for bad prognosis
tii.m = tii_prognosis(tii.df, "poorprognosis")
#make heatmap
my_heatmap2(clin.df, tii.m, "BOR")


#add to clinical file
toMerge.df = as.data.frame((data.frame(rna_IC = colSums(ic.m),
                                       rna_TIL = colSums(til.m),
                                       rna_TII_poorprognosis = colSums(tii.m))))
toMerge.df$patientID = rownames(toMerge.df)
rownames(toMerge.df) = NULL

t_tii.m = t(tii.m)
t_til.m = t(til.m)
toMerge.df2 = merge(t_tii.m,
                    t_til.m, by = 0)
colnames(toMerge.df2)[1] = "patientID"
final.df = merge(toMerge.df2, toMerge.df, by="patientID")


#merge
clin.df = merge(clin.df, final.df, by = "patientID")


add_sig_from_vector <- function(df, input_vector, name){
  
  #select genes from rna 
  selection.df = df[df$gene %in% input_vector, ]
  
  #order cell types specific order as specified in the keep vector (later for visualization)
  selection.df$gene = factor(selection.df$gene, levels=input_vector)
  selection.df = selection.df[order(selection.df$gene), ]
  
  #numeric
  m = selection.df[ ,-1]
  m = sapply(m, as.numeric)
  rownames(m) = selection.df$gene
  
  
  mean.df = as.data.frame(colSums(m))
  colnames(mean.df) = name
  mean.df$patientID = rownames(mean.df)
  
  return(mean.df)
}  

#' Rozeman et al 2021: IFN-γ score of each patient is the average z-score
#' of all the genes within the IFN-γ signature
clin.df = clin.df %>%  dplyr::select(-"rna_sig_INTERFERON_GAMMA_RESPONSE")
# that a high IFN-γ score (average z-score of all genes within the IFN-γ signature described by Ayers et al
# Ayers, M. et al. IFN-γ-related mRNA profile predicts clinical response
# to PD-1 blockade. J. Clin. Invest. 127, 2930–2940 (2017).
ifng.v = c("STAT1", "CXCL9","CXCL10","HLA-DRA", "GZMA","PRF1","IDO1", "CXCL11","CCR5", "IFNG")
ifng_mean.df = add_sig_from_vector(rna.df, ifng.v, "rna_sig_INTERFERON_GAMMA_RESPONSE")
clin.df = merge(clin.df, ifng_mean.df, by = "patientID")

krt.v = c("KRT4","KRT5","KRT6A","KRT6B","KRT10","KRT14","KRT17","KRT80")
krt_mean.df = add_sig_from_vector(rna.df, krt.v, "rna_sig_krt")
clin.df = merge(clin.df, krt_mean.df, by = "patientID")

epith.v = c("KRT4","KRT5","KRT6A","KRT6B","KRT10","KRT14","KRT17","KRT80","SBSN","DSG3","SFN")
epith_mean.df = add_sig_from_vector(rna.df, epith.v, "rna_sig_epith")
clin.df = merge(clin.df, epith_mean.df, by = "patientID")

mesench.v = c("MYH6", "ACTC1","MB","CLDN11","TNNC1","BEST3","XIRP2","MYOG","MYH7")
mesench_mean.df = add_sig_from_vector(rna.df, mesench.v, "rna_sig_mesench")
clin.df = merge(clin.df, mesench_mean.df, by = "patientID")




clin.df$BOR_bool = ifelse(clin.df$BOR == "PD",1,0)
#Make BOOLEANS
#all rna columns
rna.v = colnames(clin.df)[str_detect(string = colnames(clin.df), pattern = "rna_")]

for (i in rna.v){
  #make densities of input column for clinical and non clinical benefit
  LOW_LIM = min(clin.df[[i]])
  UPPER_LIM = max(clin.df[[i]])
  cb.density = density(subset(clin.df, BOR_bool == 1)[[i]], 
                       from = LOW_LIM, to = UPPER_LIM)
  no_cb.density = density(subset(clin.df, BOR_bool == 0)[[i]], 
                          from = LOW_LIM, to = UPPER_LIM)
  #substract values of input col of non cb from cb
  diff.v = cb.density$y - no_cb.density$y
  #find intersect
  INTERSECT = cb.density$x[which(diff(diff.v > 0) != 0) + 1]
  
  #when multiple intersects
  if(length(INTERSECT) > 1){
    print(paste(length(INTERSECT),"intersects for",i))
    
    #take intersect that is closest to median of input column
    INTERSECT = INTERSECT[which.min(abs(INTERSECT - median(clin.df[[i]])))]
  } 
  
  #add new boolean columns
  clin.df$new = clin.df[[i]] > INTERSECT
  clin.df$new = ifelse(clin.df$new == TRUE, 1, 0)
  colnames(clin.df)[ncol(clin.df)] = paste0(i, "_bool")
}


# ggplot(clin.df, aes(clin.df[[i]], fill = BOR_bool)) + 
#   geom_density(alpha = 0.2) + 
#   geom_vline(xintercept = INTERSECT, color = "red") 


clin.df$rna_combiTmlTil = rep(NA, nrow(clin.df))
clin.df$rna_combiTmlTil[clin.df$wgs_tml_bool == 1 & clin.df$rna_sig_INTERFERON_GAMMA_RESPONSE_bool == 1] = "TML-IFNGhigh"
clin.df$rna_combiTmlTil[clin.df$wgs_tml_bool == 1 & clin.df$rna_sig_INTERFERON_GAMMA_RESPONSE_bool == 0] = "TMLhigh_IFNGlow"
clin.df$rna_combiTmlTil[clin.df$wgs_tml_bool == 0 & clin.df$rna_sig_INTERFERON_GAMMA_RESPONSE_bool == 1] = "TMLlow_IFNGhigh"
clin.df$rna_combiTmlTil[clin.df$wgs_tml_bool == 0 & clin.df$rna_sig_INTERFERON_GAMMA_RESPONSE_bool == 0] = "TML-IFNGlow"


#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################


new_cols = colnames(clin.df)[str_detect(colnames(clin.df), pattern = "rna")]
final.df = clin.df %>% dplyr::select(CPCT_WIDE_CORE, new_cols)
raw.df = fread('data/pembro/20221021_DRUP_pembro_LL_final_1_WGS.tsv')
df = join(raw.df, final.df, by = "CPCT_WIDE_CORE")

#Save output
write.table(df, file='data/pembro/20221021_DRUP_pembro_LL_final_2_WGS_RNA.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)



