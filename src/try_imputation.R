# Author: T. Yordanov
# Date: 18.10.2021

remove_non_imputable_columns <-  function(df1) {
  col_to_remove <- c()
  for(cn in cns(df1)){
    if(len(uniq(df1[,cn])) > MAX_N_CATEGORIES_IN_DATASET || cn == "outcome_of_interest")
      next
    if(any(table(df1[,cn], df1$outcome_of_interest) <= 8) || nrow(table(df1[,cn])) < 2)
      col_to_remove <- c(col_to_remove, cn)
  }
  try_log_info("@remove_non_imputable_columns- removing following columns (reason: non-NA vals seen only when outcome_of_interest = 0 or 1): \n[ %s ];",
               paste(col_to_remove, collapse = ", "))
  if(len(col_to_remove) > 0)
    df1 <- df1[,setdiff(cns(df1),col_to_remove)]
  df1
}

get_imputed_dfs_list <- function(df1, n_copies) {
  imputed_dfs_list <- impute_n_copies(df1, n_copies = n_copies)
  imputed_dfs_list <- map(1:n_copies , function(x) {
    categorical_cols = cns(imputed_dfs_list[[x]]) [which(   apply(imputed_dfs_list[[x]], 2, function(x) len(uniq(x))) <= MAX_N_CATEGORIES_IN_DATASET)]
    cur_df <-  try_factorize_df(imputed_dfs_list[[x]], categorical_cols)
    return(cur_df)
  })
  return(imputed_dfs_list)
}

impute_n_copies <- function(df1, n_copies = N_IMPUTED_DATASETS, is_test_set = F){
  imputed_dfs = NULL
  imputed_dfs <- tryCatch( {
    mice::mice(df1, m = n_copies, meth="pmm", seed = random.seed,  maxit = 5 , printFlag = F)
  },
  error = function(e){
    try_log_error("@impute_n_copies mice::mice - failed, going to skip current dataset...")
  },
  finally = {})
  if(len(imputed_dfs) == 0){
    return(NULL)
  }
  imputed_dfs_list <- map(1:n_copies, function(x) {mice::complete(imputed_dfs, x)})

  nas_present_after_mice <- are_NAs_present_df_list(imputed_dfs_list)
  # if(nas_present_after_mice)
  #   try_log_warn(0,"@impute_n_copies-s NAS PRESENT AFTER IMPUTING (with mice)!  **** ")

  nas_present <- nas_present_after_mice
  try_imp_cntr <- 0
  unimputable_columns = NULL
  while(nas_present && try_imp_cntr < 3) {
    unimputable_columns = print_NAs_and_counts_df(imputed_dfs_list[[1]])
    try_imp_cntr <- try_imp_cntr + 1
    # try_log_debug("@impute_n_copies- trying to impute hard columns.. **** -> try_imp_cntr = %d", try_imp_cntr)
    # if(try_imp_cntr > 1)
    #   try_log_warn("@impute_n_copies- trying to impute hard columns. try_imp_cntr greater than 1")
    imputed_dfs_list <- try_impute_hard_columns_list(imputed_dfs_list)
    nas_present <- are_NAs_present_df_list(imputed_dfs_list)
  }
  if(nas_present) {
    # try_log_error("@impute_n_copies- *NAS PRESENT AFTER MICE AND COULD NOT BE FIXED!! Removing columns: \n\t\t-> : [ %s ]; <- ****",
    #               paste(unimputable_columns, collapse = ","))
    imputed_dfs_list = map(1:len(imputed_dfs_list), function(x) {
      imputed_dfs_list[[x]] = imputed_dfs_list[[x]][,-which(cns(imputed_dfs_list[[1]]) %in% unimputable_columns)]})

    # TODO make it set to mean rather than removing! what about those with just 1 repeating value?
    # imputed_dfs_list = map(1:len(imputed_dfs_list), function(x) {
    #   c_df = imputed_dfs_list[[x]]
    #   for(c_col in unimputable_columns){ # c_col  = unimputable_columns[1]
    #     c_na_rows = which(is.na(c_df[,c_col]))
    #
    #   }
    #   # imputed_dfs_list[[x]][,which(cns(imputed_dfs_list[[1]]) %in% unimputable_columns)] = imputed_dfs_list[[x]]
    #   })
  }
  imputed_dfs_list
}

# remove columns not present in first two dfs interstect
clean_imputed_dfs <- function(test_df, imputed_dfs_train_list){
  ## intersect columns between test & train dfs
  #prints
  cols_missing_in_test <- c()
  cols_missing_in_train <- c()
  for(i in 1:len(imputed_dfs_train_list)){
    cols_missing_in_test <- uniq(c(cols_missing_in_test, setdiff(cns(imputed_dfs_train_list[[i]]), cns(test_df)) ))
    cols_missing_in_train <- uniq(c(cols_missing_in_train, setdiff(cns(test_df) ,cns(imputed_dfs_train_list[[i]]))))
  }

  if(len(cols_missing_in_test) > 0)
    try_log_debug("@clean_imputed_dfs- Removing following columns [reason: could not be imputed in test dataset of current fold]: %s", paste(cols_missing_in_test))
  if(len(cols_missing_in_train) > 0)
    try_log_debug("@clean_imputed_dfs - Removing following columns [reason: could not be imputed in train datasets of current fold]: %s", paste(cols_missing_in_train))
  #intersect
  cols_to_remove <- union(cols_missing_in_test, cols_missing_in_train)
  if(len(which(cns(test_df) %in% cols_to_remove)) > 0 ){
    test_df <- test_df[,-which(cns(test_df) %in% cols_to_remove)]
  }
  cols_to_remove2 <-  c()
  for(cln in cns(test_df)){
    if(len(which(table(test_df[,cln]) == 0)))
      cols_to_remove2 <- c(cols_to_remove2, cln)
    else
      for(i in 1:len(imputed_dfs_train_list)){
        if(len(which(table(imputed_dfs_train_list[[i]][,cln]) == 0))){
          cols_to_remove2 <- c(cols_to_remove2, cln)
          break
        }
      }
  }
  if(len(cols_to_remove2) > 0){
    try_log_debug(0, "GOING TO REMOVE VAR %s; REASON: contains only 1 value observed in fold", cols_to_remove2)
    test_df[, cols_to_remove2] <- NULL
    for(i in 1:len(imputed_dfs_train_list)){
      imputed_dfs_train_list[[i]][, cols_to_remove2] <- NULL
    }
  }

  for(imp_indx in 1:len(imputed_dfs_train_list)){
    if(len(which(cns(imputed_dfs_train_list[[imp_indx]]) %in% cols_to_remove)) > 0)
      imputed_dfs_train_list[[imp_indx]] <- imputed_dfs_train_list[[imp_indx]][,-which(cns(imputed_dfs_train_list[[imp_indx]]) %in% cols_to_remove)]
  }
  return(list(train=imputed_dfs_train_list, test=test_df))
}

# when imputation with mice fails, we try imputing using the confusion matrix
try_impute_hard_columns <- function(df1){
  hard_cols <- c()
  impossible_cols <- c()
  for(x in cns(df1)){
    is_categorical <- len(uniq(df1[,x])) <= MAX_N_CATEGORIES_IN_DATASET
    if(is_categorical){
      if( any(is.na(df1[,x])) && x != "outcome_of_interest" && !all(is.na(df1[,x])))
        hard_cols <- c(hard_cols, x)
      if(all(is.na(df1[,x])))
        impossible_cols <- c(impossible_cols, x)
    }
  }
  if(len(impossible_cols) > 0){
    df1 <- df1[,-which(cns(df1) %in% impossible_cols)]
    try_log_info("@try_impute_hard_columns- removing columns which are impossible to impute:\n\t\t -> [ %s ]; <-",
                 paste(impossible_cols, collapse = ", "))
  }
  if(len(hard_cols) > 0)
    hard_cols <- intersect(hard_cols, cns(df1))
  for(h_col in hard_cols){
    conf_ma3x <- table(df1[,h_col], df1$outcome_of_interest, useNA = "always")
    conf_ma3x <- conf_ma3x[,-which(is.na(cns(conf_ma3x)))]
    ma3x_ncols = len(cns(conf_ma3x))
    na_rows <- which(is.na(df1[,h_col]))

    for(i in 1:len(na_rows)){
      c_outcome <- df1[na_rows[i], "outcome_of_interest"]
      is_vect = is.null(nrow(conf_ma3x))
      c_vect_from_cm <- if(!is_vect) conf_ma3x[,as.character(c_outcome)] else conf_ma3x[as.character(c_outcome)]
      c_vect_from_cm_opp <- if(!is_vect) conf_ma3x[,which(cns(conf_ma3x) != as.character(c_outcome))] else
        conf_ma3x[which(cns(conf_ma3x) != as.character(c_outcome))]

      c_prob <- sum(c_vect_from_cm) / (sum(c_vect_from_cm_opp) + sum(c_vect_from_cm))
      c_outcome_prob = ifelse(c_outcome == 1, c_prob, 1-c_prob)
      new_val <- NA
      if(!is_vect && nrow(conf_ma3x) >= 3)
        new_val = sample(c(0,1), size = 1, replace = T, prob = c(1-c_outcome_prob, c_outcome_prob))
      df1[na_rows[i],h_col] <- new_val
    }

  }
  df1
}

try_impute_hard_columns_list <- function(df_list){
  for(i in 1:len(df_list)){
    df_list[[i]] <- try_impute_hard_columns(df_list[[i]])
  }
  df_list
}

