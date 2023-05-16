final.df = data.frame(id = character(), sampleId = character(), update = numeric())
for(i in 6:1){ #start with the newest update
  print(i)
  
  drup.df = fread(paste0("/DATA/share/Voesties/data/DRUP/update_",i,"/metadata.tsv"),data.table = F)  
  
  id.v = gsub("T$|TII.*","",drup.df$sampleId)
  sampleId.v = drup.df$sampleId
  update.v = i
  new.df = data.frame(id = id.v, sampleId = sampleId.v, update = update.v)
  new.df = new.df %>% dplyr::filter(!sampleId %in% final.df$sampleId)
  
  final.df = rbind(final.df, new.df)
}
table(final.df$update)
length(unique(final.df$id))
length(unique(final.df$sampleId))

#but metadata is not always complete but the files can still be in there

fn = list.files(path = "/DATA/share/Voesties/data/DRUP/update_6/somatics",
                pattern = "^DRUP",
                full.names = F,
                recursive = F)
length(fn)

#only when looking at the actual files (and not metadata) we find the same.

#HOWEVER
#CPCT02060143, CPCT02120203, CPCT02260029 were found in HMF harmonized but not in CPCT samples
"/DATA/share/Voesties/data/HMF/update_10/somatics/CPCT02060143T/purple/"

