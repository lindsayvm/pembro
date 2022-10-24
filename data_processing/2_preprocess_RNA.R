library(data.table)
library(edgeR)
library(dplyr)
# install.packages("devtools")
# devtools::install_github("stephenturner/annotables")
library(annotables)

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

#Load data
df = fread("data/pembro/20221021_DRUP_pembro_LL_final_1_WGS.tsv", data.table = F)

#Load RNAseq all of HMF
extension = "-featurecounts.txt.gz"

fn = list.files("/DATA/share/Voesties/data/resources/rnaseq/results/quantification", 
                pattern = extension, recursive = T, full.names = T)
names = list.files("/DATA/share/Voesties/data/resources/rnaseq/results/quantification", 
                        pattern = extension, recursive = T, full.names = F)

#match RNAseq with data
rna.names = gsub("/.*","",names)
rna.names = gsub("T$|TI$|TI.*$","",rna.names)

#??! there should be 41 but there are 34
fn = fn[rna.names %in% df$HMFsampleID]
df = df[df$HMFsampleID %in% rna.names, ]

#read gene data
ls = lapply(fn, function(x) {
  
  print(x)
  
  #read rnaseq file of each patient
  df = data.table::fread(x, data.table = F)
  
  #for every 1,000,000 RNA molecules in the RNA-seq sample, x came from this gene
  df = df[,c(1, ncol(df))]
  
  names(df)
  return(df)
})
#Retrieve IDs
names(ls) = gsub('\\D+|CPCT0|DRUP0', "", gsub(pattern = paste(".*\\/|",extension), "", fn))

#check whether rownames are same
df1=ls[[1]]
df2=ls[[20]]
setequal(df1$GeneId, df2$GeneId)

#concatenate all files
rna.df = ls[[1]]
colnames(rna.df) = c("GeneId",names(ls)[1])
for(i in 2:length(fn)){
  print(i)
  tmp.df = ls[[i]]
  tmp.df = tmp.df[ ,-1]
  rna.df = cbind(rna.df, tmp.df)
  colnames(rna.df)[ncol(rna.df)] = names(ls[i])
}


#genes are denoted as ENSG00000159023.21 --> gene names
rna.df = rna.df %>% 
  mutate(GeneId = gsub("\\.[^.]*$", "", .$GeneId)) %>%   #rm everything after dot to make ENSEMBL names
  dplyr::left_join(grch38, by = c("GeneId" = "ensgene")) %>%  #get info about genes --> gene name
  distinct(symbol, .keep_all = T) %>% # keep 1 for duplicated symbol genes (only 40 genes not unique)
  mutate(GeneId = symbol) %>% 
  dplyr::select(-c(symbol, entrez, chr, start, end, strand, biotype, description)) %>% 
  filter(GeneId != "") #not all missing values are reported as NA therefore empty ""

write.table(x = as.data.frame(rna.df),
            file = "/home/l.leek/data/pembro/pembro_RNAseq_raw_cts.tsv",
            sep = "\t", quote = F,na = "NA")


###########
rna.df = rna.df[!is.na(rna.df$GeneId),] 
rna.m = rna.df %>% 
  remove_rownames %>% 
  column_to_rownames(var = "GeneId") 

#counts per million (CPM), per library/per sample calc fraction if total was 1e6
cpm.df = cpm(rna.m)
colSums(cpm.df)

#Select (not remove) genes (rows) with < 10 counts per million per library/patient
# Dont remove because you need all count for proper normalization
cutoff = 10 #30
drop = which(apply(cpm.df, 1, max) < cutoff)

# you do this on the cpm data because that is not filtered
#Joris: log2(sample * 10e6 / library size + 1)
cpm_normalized.df = log10(cpm.df +1)

#Filter genes out that needed to be removed
cpm_normalized_filtered.df = cpm_normalized.df[-drop, ] 

final.df = as.data.frame(dplyr::as_tibble(cpm_normalized_filtered.df, rownames = "gene"))
colnames(final.df)[-1] = paste0("ID_", colnames(final.df)[-1])

write.table(x = as.data.frame(final.df),
            file = "/home/l.leek/data/pembro/pembro_RNAseq_normalized_cts.tsv",
            sep = "\t", quote = F,na = "NA", row.names = FALSE)
