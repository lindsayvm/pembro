
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
  df = df %>% dplyr::select(patientID, 
                            purity, 
                            msIndelsPerMb, msStatus,
                            tml, tmlStatus, tmbPerMb, tmbStatus, svTumorMutationalBurden) %>% 
    
    dplyr::rename(MSI_bool = msStatus,
                  TML_bool = tmlStatus,
                  TMB_bool = tmbStatus,
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
    filter(gene %in% GOI) %>% 
    mutate(gene = droplevels(gene)) %>% 
    dplyr::select(patientID, gene)
  GOI_dummy_driver.df = dummy_cols(GOI_driver.df, select_columns = "gene")
  GOI_dummy_driver.df = GOI_dummy_driver.df %>% 
    dplyr::select(-gene) %>% 
    group_by(patientID) %>% 
    summarise(across(everything(), sum))
  
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


