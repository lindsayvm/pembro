library(data.table)
library(edgeR)
library(dplyr)
# install.packages("devtools")
# devtools::install_github("stephenturner/annotables")
library(annotables)

dir = "/home/l.leek/pembro"
setwd(dir)

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

#Load data
df = fread("data/20221021_DRUP_pembro_LL_final.tsv", data.table = F)

#Load RNAseq all of HMF
extension = "-featurecounts.txt.gz"

fn = list.files("/DATA/share/Voesties/data/resources/rnaseq/results/quantification", 
                pattern = extension, recursive = T, full.names = T)
names = list.files("/DATA/share/Voesties/data/resources/rnaseq/results/quantification", 
                        pattern = extension, recursive = T, full.names = F)
#match RNAseq with data
rna.names = gsub("/.*","",names)
rna.names = gsub("T$|TI$|TI.*$","",rna.names)

#because we dont know if names are used from CPCT_WIDE_CORE or HMFsampleID we use both. Some patients are submitted under both names I guess which is why we get 66 in total. (32 form CPCT... and 34 from HMF...)
#its not thatttt many patients, so for now just run on all and select later for which one we want the features. 
PATIENTS = unique(c(df$CPCT_WIDE_CORE, df$HMFsampleID))
PATIENTS = rna.names[rna.names %in% PATIENTS]
fn = fn[rna.names %in% PATIENTS]

# df1 = df[df$CPCT_WIDE_CORE %in% PATIENTS, ]
# df2 = df[df$HMFsampleID %in% PATIENTS, ]
# df = rbind(df1, df2)
# #

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
names(ls) = gsub(pattern = ".*\\/|T-featurecounts.txt.gz", "", fn)

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
            file = "data/pembro_RNAseq_raw_cts.tsv",
            sep = "\t", quote = F,na = "NA", row.names = FALSE)


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

write.table(x = as.data.frame(final.df),
            file = "data/pembro_RNAseq_normalized_cts.tsv",
            sep = "\t", quote = F,na = "NA", row.names = FALSE)
