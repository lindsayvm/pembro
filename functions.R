
my_dataframe <- function(fn, extension) {
  #' make df from HMF data ordered per patient
  
  #read gene data
  ls = lapply(fn, function(x) {
    df = read.table(file = x, sep = '\t', header = TRUE, stringsAsFactors = TRUE)
    return(df)})
  #Retrieve IDs
  names(ls) = gsub('\\D+', "", gsub(pattern = paste(".*\\/|",extension), "", fn))
  
  #concatenate all files
  nrow = 0
  ncol = ncol(ls[[1]]) +1 # +1 because want to add patientID
  purple.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
  colnames(purple.df) = colnames(ls[[1]])
  for(i in 1:length(fn)){
    print(i)
    tmp.df = ls[[i]]
    if(nrow(tmp.df) != 0){
      tmp.df = cbind(patientID = names(ls[i]), tmp.df)
      purple.df = rbind(purple.df, tmp.df)
    }
  }
  return(purple.df)
}
