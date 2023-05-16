
read_files <- function(fn){
  lapply(fn, function(x) {
    df = read.table(file = x, sep = '\t', header = TRUE, stringsAsFactors = TRUE)
    return(df)})
}

concatenate_files <- function(ls){
  nrow = 0
  ncol = ncol(ls[[1]]) +1 # +1 because want to add patientID
  final.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
  colnames(final.df) = colnames(ls[[1]])
  for(i in 1:length(ls)){
    print(i)
    tmp.df = ls[[i]]
    if(nrow(tmp.df) != 0){
      tmp.df = cbind(patientID = names(ls[i]), tmp.df)
      final.df = rbind(final.df, tmp.df)
    }
  }
  return(final.df)
}



my_dataframe <- function(fn, extension) {
  #' make df from HMF data ordered per patient
  
  #read gene data
  ls = read_files(fn = fn$ID)
  names(ls) = fn$names
  
  #concatenate all files
  purple.df = concatenate_files(ls = ls)
  
  return(purple.df)
}

process_purplefiles <- function(df){
  df = df %>% 
    dplyr::rename(MSI_bool = msStatus,
                  TML_bool = tmlStatus,
                  TMB_bool = tmbStatus,
                  sex = gender,
                  TMB = tmbPerMb,
                  TML = tml) 
  
  # #to change values make factors --> characters
  # df = as.data.frame(lapply(df,as.character))
  # #...as.numeric
  # df = as.data.frame(lapply(df,as.numeric))
  
  return(df)
  
}





my_cosmicsignatures <-function(df) {
  df %>% 
    mutate(sample_id = gsub("T$","", sample_id)) %>% 
    dplyr::rename(SBS_unknown = unknown) %>% 
    mutate(SBS_7_UV_sum = rowSums(df %>% dplyr::select(starts_with("SBS7")))) %>% 
    mutate(SBS_10_Poly_e_sum = rowSums(df %>% dplyr::select(starts_with("SBS10")))) %>% 
    mutate(SBS_7_UV_bool = ifelse(SBS_7_UV_sum > 0.5, 1, 0)) %>% 
    mutate(SBS_10_Poly_e_bool = ifelse(SBS_10_Poly_e_sum > 0.5, 1, 0))
}



my_driverMut_pp <- function(df, driverLikelihood_thres){
  nrow = 0
  ncol = ncol(df) 
  pp.df = as.data.frame(matrix(rep(NA, ncol * nrow), 
                               ncol = ncol, 
                               nrow = nrow))
  colnames(pp.df) = colnames(df)
  #add info
  for (i in unique(df$patientID)){
    tmp.df = df[df$patientID == i, ]
    
    tmp.df = tmp.df %>% 
      filter(driverLikelihood >= driverLikelihood_thres)  %>% 
      distinct(gene, .keep_all = TRUE)
    
    pp.df = rbind(pp.df, tmp.df)
  }
  return(pp.df)
}

process_genes_to_dummy <- function(driverMut_pp.df, GOI, patientIDs){
  GOI_driver.df = driverMut_pp.df %>% 
    dplyr::filter(gene %in% GOI) %>% 
    mutate(gene = droplevels(gene)) %>% 
    dplyr::select(patientID, gene)
  GOI_dummy_driver.df = dummy_cols(GOI_driver.df, select_columns = "gene")
  GOI_dummy_driver.df = GOI_dummy_driver.df %>% 
    dplyr::select(-gene) %>% 
    group_by(patientID) %>% 
    dplyr::summarise(across(everything(), sum))
  
  #add zero for the samples that are sequenced but dont have driver muts
  noDriverMutPatients.v = patientIDs[!patientIDs %in% GOI_dummy_driver.df$patientID]
  nrow = length(noDriverMutPatients.v)
  ncol = ncol(GOI_dummy_driver.df)-1
  noDriverMutPatients.df = as.data.frame(matrix(rep(0, ncol * nrow), 
                                                ncol = ncol, 
                                                nrow = nrow))
  noDriverMutPatients.df = cbind(noDriverMutPatients.v, noDriverMutPatients.df)
  colnames(noDriverMutPatients.df) = colnames(GOI_dummy_driver.df)
  GOI_dummy_driver.df = rbind(GOI_dummy_driver.df,noDriverMutPatients.df)
  
  return(GOI_dummy_driver.df)
}

my_GOI_driverMut_pp <- function(df, driverLikelihood_thres, GOI, patientIDs){
  driverMut_pp.df = my_driverMut_pp(df, driverLikelihood_thres)
  GOI_driverMut_pp.df = process_genes_to_dummy(driverMut_pp.df, GOI, patientIDs)
  return(GOI_driverMut_pp.df)
}






my_cnv_pp <- function(df){
  nrow = 0
  ncol = 3
  pp.df = as.data.frame(matrix(rep(NA, ncol * nrow), 
                               ncol = ncol, 
                               nrow = nrow))
  colnames(pp.df) = c("patientID", "avg_cnv", "aneuploidy_score")
  #add info
  for (i in unique(df$patientID)){
    tmp.df = df[df$patientID == i, ]
    
    tmp.df = tmp.df %>% 
      mutate(cnv_len = (end-start)) %>% 
      mutate(cnv_weights = cnv_len * copyNumber)
    
    #average cnv per entire genome
    avg_cnv = sum(tmp.df$cnv_weights)/sum(tmp.df$cnv_len)
    
    #calc per chrom
    agg_tbl_perChrom = tmp.df %>% group_by(chromosome) %>% 
      dplyr::summarise(sum_cnv_weights = sum(cnv_weights),
                       sum_cnv_len = sum(cnv_len),
                       .groups = 'drop') %>% 
      mutate(avg_cnv_perChrom = sum_cnv_weights/sum_cnv_len)
    #sum of average cnv per chrom
    aneuploidy_score = sum(agg_tbl_perChrom$avg_cnv_perChrom)
    
    tmp.df = data.frame(cbind(patientID = i,
                  avg_cnv = as.numeric(round(avg_cnv,5)),
                  aneuploidy_score = as.numeric(round(aneuploidy_score,5))))
    tmp.df[-1] = as.data.frame(lapply(tmp.df[-1],as.numeric))
    pp.df = rbind(pp.df, tmp.df)
  }
  return(pp.df)
}




read_anno_vcf_files <- function(fn) {
  lapply(fn, function(x) {
    df = read.table(file = x, sep = '\t', header = TRUE, stringsAsFactors = TRUE) %>% 
      mutate(CHROMPOS = paste0(CHROM, "_", POS)) %>%          #counts as 1 mut
      dplyr::filter(FILTER == "PASS") %>%                     #filter on quality 
      dplyr::filter(EFF....IMPACT %in% c("MODERATE","HIGH")) %>%   #filter on nonsyn 
      rename_with(~ gsub("EFF....", "", .x, fixed = TRUE)) %>% 
      dplyr::select(CHROMPOS, GENE, FILTER, IMPACT, EFFECT, BIOTYPE)
    return(df)})
}


read_clonal_vcf_files <- function(fn) {
  lapply(fn, function(x) {
    df = as.data.frame(getFIX(read.vcfR(x, verbose = FALSE), getINFO = T)) %>% 
      mutate(CHROMPOS = paste0(CHROM, "_", POS)) %>%    #counts as 1 mut
      mutate(SUBCL = gsub(".*;SUBCL=|;.*","",INFO)) %>% #extract subclonality info
      dplyr::select(CHROMPOS, SUBCL)                    #speed up by rm info
    return(df)})
  
}


get_TML_cTML <- function(nonsyn.ls, clonal.ls) {
  
  df = as.data.frame(matrix(rep(NA, 0), nrow = 0, ncol = 8))
  colnames(df) = c("patientID","TML_SNPeff", "TMB_SNPeff_pass_nonsyn_protcoding", 
                   "TMB_SNPeff_pass_nonsyn_protcoding_snps", "cTML","perc_clon", "cTML_05", "perc_clon_05")
  for(i in names(nonsyn.ls)){
    print(i)
    
    TML = length(unique(nonsyn.ls[[i]]$CHROMPOS))
    
    TMB.df = nonsyn.ls[[i]] %>% 
      dplyr::filter(BIOTYPE == "protein_coding")
    TMB = length(unique(TMB.df$CHROMPOS))

    TMB_snps.df = nonsyn.ls[[i]] %>% 
      dplyr::filter(BIOTYPE == "protein_coding") %>%  
      dplyr::filter(EFFECT == "missense_variant")
    TMB_snps = length(unique(TMB_snps.df$CHROMPOS))

    clonal_muts = clonal.ls[[i]]$CHROMPOS[as.numeric(clonal.ls[[i]]$SUBCL) < 0.85]
    clonal_nonsyn.df = nonsyn.ls[[i]] %>% 
      dplyr::filter(CHROMPOS %in% clonal_muts)
    cTML = length(unique(clonal_nonsyn.df$CHROMPOS))
    
    perc_clon = round((cTML / TML ) * 100,1)

    #different cutoff
    clonal_muts_05 = clonal.ls[[i]]$CHROMPOS[as.numeric(clonal.ls[[i]]$SUBCL) < 0.05]
    clonal_nonsyn_05.df = nonsyn.ls[[i]] %>% 
      dplyr::filter(CHROMPOS %in% clonal_muts_05)
    cTML_05 = length(unique(clonal_nonsyn_05.df$CHROMPOS))
    
    perc_clon_05 = round((cTML / TML ) * 100,1)
    
        
    patientID = gsub("T$|TI$|TI.*$","",i)
    
    new_row = c(patientID, as.integer(TML), as.integer(TMB), as.integer(TMB_snps),  as.integer(cTML), perc_clon, cTML_05, perc_clon_05)
    
    df[nrow(df) + 1, ] = new_row
  }
  return(df)
}

my_clonal_dataframe <- function(fn) {
  
  print("Reading all SNPeff annotated files and filtering for nonsyn")
  nonsyn.ls = read_anno_vcf_files(fn$anno)
  names(nonsyn.ls) = fn$names

  print("Extracting clonal information per mutation from vcf. This might take a while.")  
  clonal.ls = read_clonal_vcf_files(fn$clonal)
  names(clonal.ls) = fn$names
  
  print("Select clonal nonsyn mutations")
  final.df = get_TML_cTML(nonsyn.ls, clonal.ls)
  final.df = final.df %>% 
    mutate_at(vars("TML_SNPeff", "TMB_SNPeff_pass_nonsyn_protcoding", "TMB_SNPeff_pass_nonsyn_protcoding_snps", "cTML", "perc_clon","cTML_05", "perc_clon_05"), as.numeric)
  
  return(final.df)
}




my_ann_snpeff <- function(ann.fn, GOI){
  #read annotated gene data
  ls = lapply(ann.fn, function(x) {
    df = fread(x, data.table = F) %>% 
      dplyr::rename_with(~ gsub("[*]", "", .x, fixed = TRUE)) %>%
      dplyr::rename_with(~ gsub("EFF.", "", .x, fixed = TRUE)) %>% 
      mutate(chrompos = paste0(CHROM,"_", POS)) %>% 
      filter(GENE %in% GOI)
    
    return(df)})
  #Retrieve IDs
  names(ls) = gsub("T$|T.$|T..$","",gsub('.*\\/|_ann_filt_oneLine.vcf', "", ann.fn))
  
  #concatenate all files
  #make starting df
  nrow = 0
  ncol = ncol(ls[[1]]) +1 # +1 because want to add patientID
  ann.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
  colnames(ann.df) = c("patientID",colnames(ls[[1]]))
  #add info
  for(i in 1:length(ann.fn)){
    tmp.df = ls[[i]]
    if(nrow(tmp.df) != 0){
      tmp.df = cbind(patientID = names(ls[i]), tmp.df)
      ann.df = rbind(ann.df, tmp.df)
    }
  }
  return(ann.df)
}

my_ann_pp <- function(ann.df, GOI, ann.fn){
  
  total_pts = gsub("T$|T.$|T..$","",gsub(".*/|_ann_filt_oneLine.vcf","",ann.fn))
  clin.df = data.frame(patientID = total_pts)
  #add new cols of each mutation as WT/MUT
  for (i in GOI){
    
    #make new col
    clin.df$new = rep(0)
    
    #select only moderate and high impact genes
    impact.df = ann.df[ann.df$GENE == i ,]
#                         ann.df$IMPACT %in% c("MODERATE", "HIGH"), ]
    mut.id = unique(impact.df$patientID)
    clin.df$new[clin.df$patientID %in% mut.id] = 1
    colnames(clin.df)[colnames(clin.df) == "new"] = paste0("gene_snpeff_", i)
  }
  return(clin.df)
}















