#####################################################
# Define a function to evaluate correction accuracy #
#####################################################
compareGeno <- function(geno, true, out_list){
  if(is.null(geno)){
    geno <- out_list$best_geno[, out_list$param_list$samples_index]
    geno <- t(abs(geno - 2))
  }
  ref2het <- true == 0 & geno == 1
  alt2het <- true == 2 & geno == 1
  het2ref <- true == 1 & geno == 0
  het2alt <- true == 1 & geno == 2
  ref2alt <- true == 0 & geno == 2
  alt2ref <- true == 2 & geno == 0
  misscall <- ref2het | alt2het | het2ref | het2alt | ref2alt | alt2ref
  missing <- is.na(geno)
  correct <- !misscall & !missing
  n_mar <- ncol(true)
  df_ind <- data.frame(
    correct = apply(correct, 1, sum, na.rm = T) / n_mar,
    miscall = apply(misscall, 1, sum, na.rm = T) / n_mar,
    missing = apply(missing, 1, sum, na.rm = T) / n_mar
  )
  n_ind <- nrow(true)
  df_mar <- data.frame(
    correct = apply(correct, 2, sum, na.rm = T) / n_ind,
    missing = apply(missing, 2, sum, na.rm = T) / n_ind,
    miscall = apply(misscall, 2, sum, na.rm = T) / n_ind
  )
  return(list(ind = df_ind, mar = df_mar))
 }
