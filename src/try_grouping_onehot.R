# Author: T.Y.
# Date created: 01.03.2021
# Script contents:
# Usage

## one-hot encode cats
onehot_enc_categoricals <- function(df1){
  cols_to_encode <- cns(df1)
  for(col_nm in cols_to_encode){
    if(col_nm %in% cns(df1))
      df1 = set_onehot_cols_in_df (df1, col_nm, onehot_enc_data(df1[,col_nm], col_nm = col_nm))
  }
  df1
}

set_onehot_cols_in_df <- function(df1, col_nm, new_columns, n_cats = len(uniq(df1[,col_nm]))){
  var_indx <- which(colnames(df1) == col_nm)
  for(j in 1:n_cats){
    if(n_cats == 2 && j == 2)
      next  #  for 2 cats only we need 1 iteration (i.e. binary var)
    if(j != 1)
      df1 <- cbind(  df1[,1:(var_indx-1)], new_columns[,j],  df1[,(var_indx:dim(df1)[2]) ] )
    colnames(df1)[var_indx] = sprintf("%s_cat_%d_", col_nm, j)
    if(n_cats > 2)
      df1[,var_indx] = new_columns[,j]
    else
      df1[,var_indx] = new_columns

    var_indx = var_indx + 1
  }
  df1
}

onehot_enc_data <- function(data_in, n_cats = -1, col_nm = ""){
  n_rows <- len(data_in)
  if(n_cats == 2)
    n_cats = 1 # binary , dont need to variables..
  if(n_cats == -1) {
    n_cats <- len(table(data_in))
    try_log_warn(0 ,"@onehot_enc_data- called onehot_enc_data without defining n_cats")
  }
  data_out <- data.frame(matrix(nrow = n_rows, ncol = n_cats))
  cat_vals <- as.numeric(names(table(data_in)))

  try_log_debug("@onehot_enc_data- **** COLUMN %s; Going to one-hot encode values ****", col_nm)
  try_log_debug("[ %s ];", paste( cat_vals, collapse = ","))
  try_log_debug(0,"**** into _cat_ ****")
  try_log_debug("[ %s ];", paste( 1: n_cats, collapse = ","))

  for(i in 1:n_cats){
    data_out[,i] = as.numeric(data_in == cat_vals[i])
    if(len(is.na(data_out[,i])) == 1 && is.na(data_out[,i]))
      data_out[,i] = 0
  }

  if(n_cats == 2) # if its binary, no need to split into more cats
    data_in
  else
    data_out[,1:n_cats]
}
