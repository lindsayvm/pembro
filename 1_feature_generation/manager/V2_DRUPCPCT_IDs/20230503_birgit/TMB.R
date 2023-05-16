source("src/functions_plots.R")

clin.df = fread('data/20230503_DRUP_pembro_LL_WGS_RNA_birgit.tsv')
clin.df = clin.df %>% dplyr::select(c("CPCT_WIDE_CORE","HMFsampleID","patientID","Cohort",
                                      "MutationalLoad","TML","TML_SNPeff","responders"))
summary(clin.df$TML_SNPeff)
# Min: 355
# Median: 549.5
# Mean: 976.2
# Max: 3192
# Missing IDs: 4
# Wilcoxon pval = 0.081
my_wilcoxon2(clin.df, "TML_SNPeff") + ylab("TML (SNPeff)")
clin.df$patientID[is.na(clin.df$TML_SNPeff)]


summary(clin.df$TML)
# Min: 302
# Median: 350.0
# Mean: 798.9
# Max: 2880
# Missing: 3
# Wilcoxon pval = 0.056
clin.df$patientID[is.na(clin.df$TML)]
clin.df$CPCT_WIDE_CORE[is.na(clin.df$TML)]
clin.df$HMFsampleID[is.na(clin.df$TML)]
my_wilcoxon2(clin.df, "TML") + ylab("TML (HMF)")

clin.df %>%
  group_by(responders) %>%
  summarise(median_value = median(TML, na.rm=T))
