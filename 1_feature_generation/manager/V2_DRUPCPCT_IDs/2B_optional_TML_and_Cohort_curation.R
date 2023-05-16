###!!!!! TML curation stuff 

#"Deze pt is geïncludeerd op basis van TSO500: Totaal TMB is 166 mut/Mb, non synonymous TMB is 121.6 mut/Mb. Op basis van een in silico analyse van HMF leek het ons heel waarschijnlijk dat deze pt een HML van >290 had."
#this one is still included because we dont have WGS anyway
clin.df$MutationalLoad[clin.df$CPCT_WIDE_CORE == "TSO500" & clin.df$MutationalLoad == 166] = NA# "assumed > 290"

clin.df$Cohort[clin.df$CPCT_WIDE_CORE == "WIDE01010534" & clin.df$MutationalLoad == 274] = "Pembro HML>290" # "274 but >290 according to SNPeff and HMF"
clin.df$Cohort[clin.df$HMFsampleID == "DRUP01030053" & clin.df$MutationalLoad == 143] = NA # "143 but low (30) according to SNPeff and HMF"

clin.df$MutationalLoad[clin.df$CPCT_WIDE_CORE == "WIDE01010534" & clin.df$MutationalLoad == 274] = NA# "274 but >290 according to SNPeff and HMF"
clin.df$MutationalLoad[clin.df$HMFsampleID == "DRUP01030053" & clin.df$MutationalLoad == 143] = NA# "143 but low (30) according to SNPeff and HMF"
