clin1.df = fread('data/20221021_DRUP_pembro_LL_final.tsv', data.table = F)
clin1.df$patientID

clin.df = read_excel("data/20230127_ExportPembroCohorten.xlsx") %>% 
  dplyr::select(-c("...27","...28"))
clin.df$WGS_ID

clin1.df$patientID[clin1.df$patientID %in% clin.df$WGS_ID]
clin1.df$patientID[!clin1.df$patientID %in% clin.df$WGS_ID]

clin.df$WGS_ID[clin.df$WGS_ID %in% clin1.df$patientID]
clin.df$WGS_ID[!clin.df$WGS_ID %in% clin1.df$patientID]

#1 wrong selection DRUP patients by looking at HMF and not DRUP ID
#2 some still had CPCT but we needed DRUP ID, most patients are the same but 7 are different ..., but we should at least find 65 in meta.df

drup_meta.df =fread("/DATA/share/Voesties/data/DRUP/update_6/metadata.tsv", data.table = F) 
cpct_meta.df = fread("/DATA/share/Voesties/data/HMF/update_10/metadata.tsv", data.table = F) %>% 
  dplyr::filter(!hmfSampleId %in% drup_meta.df$hmfSampleId)
x=fread("/DATA/share/Voesties/data/harmonize/output/metadata.csv", data.table = F)

length(cpct_meta.df$hmfSampleId) + length(drup_meta.df$hmfSampleId) #5066
length(unique(c(cpct_meta.df$hmfSampleId, drup_meta.df$hmfSampleId)))#5064 #filter doesnt wrok
length(unique(c(cpct_meta.df$hmfSampleId, drup_meta.df$hmfSampleId, x$hmfPatientId)))#5064 #filter doesnt wrok

#meta.df = bind_rows(drup_meta.df, cpct_meta.df) %>% dplyr::distinct("hmfSampleId", .keep_all = TRUE) #just keeps 1
meta.df = bind_rows(drup_meta.df, cpct_meta.df,x) %>% 
  mutate(WGS_ID = gsub("T$|TI|TI.*","",sampleId))

merge.df = left_join(clin.df, meta.df, by = "WGS_ID")
colSums(is.na(merge.df)) #14 NAs in sampleId. ==> There are 14 patients in WGS_ID that do have WGS but are not in the meta data
