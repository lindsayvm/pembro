#libraries
library(data.table)
library(randomcoloR)
library(plotly)
library(dplyr)
library(patchwork)
library(ggpubr)

dir = "/home/l.leek/"
setwd(dir)

source("src/pembro/functions_plots.R")
source("src/pembro/functions.R")

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

clin.df = fread("data/pembro/20221021_DRUP_pembro_LL_final.tsv", data.table = F)

# Filenames of available PURPLE files (WGS)
hmf_raw.fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_10/somatics",
                         pattern = ".purple.purity.tsv$",
                         full.names = TRUE,
                         recursive = TRUE)
hmf_raw_names.fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_10/somatics",
                              pattern = ".purple.purity.tsv$",
                              full.names = F,
                              recursive = TRUE)

hmf_raw_names = gsub("/.*","",hmf_raw_names.fn)
hmf.id = gsub("T$|TI$|TI.*$","",hmf_raw_names)
hmf_drup.id = hmf.id[hmf.id %in% clin.df$CPCT_WIDE_CORE]

#Filter WGS data from HMF database for which clinical data is available 
hmf_drup.fn = hmf_raw.fn[hmf.id %in% hmf_drup.id]

#remove the later biopsy from the same patient
hmf_drup.fn = hmf_drup.fn[hmf_drup.fn != "/DATA/share/Voesties/data/HMF/update_10/somatics/CPCT02010359TII/purple/CPCT02010359TII.purple.purity.tsv"]


#Load DRUP data to find missing samples
drup_raw.fn = list.files(path = "/DATA/share/Voesties/data/DRUP/update_3/somatics",
                               pattern = ".purple.purity.tsv$",
                               full.names = T,
                               recursive = TRUE)
drup_raw_names.fn = list.files(path = "/DATA/share/Voesties/data/DRUP/update_3/somatics",
                               pattern = ".purple.purity.tsv$",
                               full.names = F,
                               recursive = TRUE)
drup_raw_names = gsub("/.*","",drup_raw_names.fn)
drup.id = gsub("T$|TI$|TI.*$","",drup_raw_names)

#Missing IDs (CPCT_WIDE_CORE style)
missing_drup_CPCT_WIDE_CORE.id = clin.df$CPCT_WIDE_CORE[!clin.df$CPCT_WIDE_CORE %in% hmf_drup.id]
#Find matching DRUP IDs (DRUP style)
missing.df = clin.df[clin.df$CPCT_WIDE_CORE %in% missing_drup_CPCT_WIDE_CORE.id, ]
drup_drup.fn = drup_raw.fn[drup.id %in% missing.df$HMFsampleID]


#Merge
final.fn = c(hmf_drup.fn, drup_drup.fn)

#' ###########################################################################
#' ###########################################################################
#' Process:  
#' ###########################################################################
#' ###########################################################################

#Extract WGS data from PURPLE files
extension = "purple.purity.tsv"
purple.df = my_dataframe(final.fn, extension)

purple.df = purple.df %>% 
  dplyr::select(patientID, purity, ploidy, 
                wholeGenomeDuplication, 
                msIndelsPerMb, msStatus, 
                tml, tmlStatus, tmbPerMb, tmbStatus, svTumorMutationalBurden)

#Add wgs prefix to all colnames
colnames(purple.df) = c("patientID", paste0("wgs_", colnames(purple.df)[-1]))

#Rename & make clear which columns are booleans
names(purple.df)[names(purple.df) == "wgs_gender"] = "wgs_gender_bool"
names(purple.df)[names(purple.df) == "wgs_wholeGenomeDuplication"] = "wgs_wholeGenomeDuplication_bool"
names(purple.df)[names(purple.df) == "wgs_msStatus"] = "wgs_msStatus_bool"
names(purple.df)[names(purple.df) == "wgs_tmlStatus"] = "wgs_tml_bool"
names(purple.df)[names(purple.df) == "wgs_tmbStatus"] = "wgs_tmb_bool"
names(purple.df)[names(purple.df) == "wgs_tmbPerMb"] = "wgs_tmb"

#Rename continous wgs_ploidy into categorical high/low
purple.df$wgs_ploidy_bool = round(purple.df$wgs_ploidy)
purple.df$wgs_ploidy_bool = ifelse(purple.df$wgs_ploidy_bool > 2, 1, 0)

#also log of tml/tmb
purple.df$wgs_tml_log = log(purple.df$wgs_tml)
purple.df$wgs_tmb_log = log(purple.df$wgs_tmb)

#to change values make factors --> characters
purple.df = as.data.frame(lapply(purple.df,as.character))

#to numbers
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "LOW"] = 0
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "HIGH"] = 1
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "MALE"] = 0
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "FEMALE"] = 1
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "false"] = 0
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "true"] = 1
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "NORMAL"] = 0
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "MSS"] = 0
purple.df[ ,1:ncol(purple.df)][purple.df[ ,1:ncol(purple.df)] == "MSI"] = 1

#as.numeric
purple.df = as.data.frame(lapply(purple.df,as.numeric))


# Add cols of WGS info to clinical file according to patientID
clin.df$patientID = gsub('\\D+','', clin.df$CPCT_WIDE_CORE)
clin.df$patientID = as.integer(gsub("^0","",clin.df$patientID))

#MERGE
final.df = left_join(clin.df, purple.df, by = "patientID")
head(final.df)


#UV signatures by Thomas Battaglia
mutSig.df = fread("/home/l.leek/data/tango/signatures/cosmicsigs-data.tsv")
colnames(mutSig.df)[1]  = "patientID"
mutSig.df$patientID = gsub("T$","",mutSig.df$patientID)
colnames(mutSig.df)[colnames(mutSig.df) == "unknown"] = "SBS_unknown"
colnames(mutSig.df)[-1] = paste0("wgs_", colnames(mutSig.df)[-1])

#calc the cumulative sum for UV signatures
mutSig.df$wgs_uv_sum = rowSums(mutSig.df[ ,c("wgs_SBS7a",  
                                             "wgs_SBS7b",   
                                             "wgs_SBS7c", 
                                             "wgs_SBS7d")])
#merge
final.df = inner_join(final.df, mutSig.df, by = c("CPCT_WIDE_CORE"="patientID"), keep =F)


#Threshold UV sig
dens = ggplot(final.df, aes(x = wgs_uv_sum)) + 
  geom_density()
ggplotly(dens)

threshold = 0.5
final.df$wgs_uv_bool = ifelse(final.df$wgs_uv_sum > threshold, 1, 0)

#Sometimes patients with low UV have high UV signature
uv = ifelse(final.df$wgs_uv_bool == 1, "high", "low")
tml = ifelse(final.df$wgs_tml_bool == 1, "high", "low")
table(uv, tml) #useNA = "ifany"

final.df$wgs_tml_scaled = log(final.df$wgs_tml)
final.df$wgs_uv_sum_scaled = scale(final.df$wgs_uv_sum)
ggplot(final.df, aes(x = wgs_tml_scaled, y = wgs_uv_sum_scaled)) +
  geom_point()+
  geom_smooth(method = "loess") +
  stat_cor(method = "pearson") 



#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################

final.df= final.df %>% dplyr::select(-c("Cohort","TumorType",                      
                       "Association","Gender","BOR",                            
                       "MutationalLoad","Start","End",                           
                       "PretreatmentBiopsy","RNA","DNA","patientID"    ))
raw.df = fread('data/pembro/20221021_DRUP_pembro_LL_final.tsv')
df = join(raw.df, final.df, by = "CPCT_WIDE_CORE")

#Save output
write.table(df, file='data/pembro/20221021_DRUP_pembro_LL_final_1_WGS.tsv', 
            quote=TRUE, sep='\t', col.names = T, row.names = F)







