# - Enrichment mutaties in patiënten met CB (SD/PR) versus patiënten zonder CB (PD), zoals bijv
# -        PI3K pathway (PIK3CA, PTEN)
# -        RAS/MAPK pathway (KRAS, BRAF, NRAS)
# -        SWI/SNF pathway  (ARID1A, SMARCB1, SMARC4A, PBRM1)
# -        IFN pathway (JAK1/JAK2STAT1, IFNGR1, IFNGR2)
# -        MHC antigen presentation (B2M, TAP1, TAP2, TAPB, CALR, PDIA3, CNAX, HSPA5,KEAP/STK11)
# 
# - Clonal vs. subclonal mutaties in patiënten met CB (SD/PR) versus patiënten zonder CB (PD)
# 
# - CUPPA analyses van de cancer of unknown primary’s à komt hun tumor toevallig overeen met primaire tumoren waarvan we weten dat ze goed reageren op anti-PD1?
#   
# Heeft geen zin:   
# - Vergelijking TMB/TML etc in patiënten met CB (SD/PR) versus patiënten zonder CB (PD) (al zijn de verschillende afkapwaarden binnen de cohorten misschien wel wat lastig om dit eerlijk te vergelijken)

# Is niet ready:
# - Immunogenecity score à dit wordt door HMF ontworpen, maar is volgens mij nog steeds niet af helaas. Tom gebruikt wel een neoantigen prediction pipeline voor het analyses van 3rd stage nivolumab MSI cohort, wellicht kunnen we zoiets ook gebruiken voor deze data?


# If more from RNAseq is needed; we need quantify these samples too. 
#missing
# [1] "CPCT02350008" "WIDE01011073" "COREDB010053" "WIDE01010451" "CPCT02080234" "CPCT02120203" "WIDE01010859"
# [8] "CPCT02130063" "WIDE01010760" "WIDE01010416" "CPCT02010936" "WIDE01010489" "WIDE01010093" "CPCT02460001"
# [15] "CPCT02080290" "CPCT02120137" "CPCT02210059" "CPCT02260029"

# can be found in list.files("/DATA/share/Voesties/data/HMF/update_10/somatics/", full.names = F)
# [1] "CPCT02350008" "WIDE01011073" "WIDE01010451" "CPCT02080234" "CPCT02120203" "CPCT02130063" "WIDE01010760"
# [8] "WIDE01010416" "CPCT02010936" "WIDE01010489" "WIDE01010093" "CPCT02460001" "CPCT02080290" "CPCT02120137"
# [15] "CPCT02210059" "CPCT02260029"
