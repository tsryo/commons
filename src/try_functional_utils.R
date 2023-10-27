# Author: T. Yordanov
# Date: 15.10.2021

map <- function(in_v, in_f) {
  res <- list()
  if(len(in_v) > 0)
    for(i in 1:len(in_v)){
      res[[i]] <- in_f(in_v[i])
    }
  return(res)
}

umap <- function(in_v, in_f) {
  return(unlist(map(in_v, in_f)))
}

strlen <- function(x) { stringr::str_length(x) }

howmany <- function(x) { len(which(x))}

nuniq <- function(x) {len(uniq(x))}

lnames <- function(x) {len(names(x))}

try_round <- function(x, digits = 0) {
  is_df = "data.frame" %in% class(x)
  if(is_df){
    for(c_nm in cns(x)){
      if(is.numeric(x[,c_nm]))
        x[,c_nm] = round(x[,c_nm], digits = digits)
    }
    return(x)
  } else return(round(x, digits = digits))

}

apply_on_all_cells_df <- function(df1, fun1) {
  apply(df1, MARGIN = c(1,2), FUN = fun1)
}

round_up <- function(x) { if(round(x) >= x) round(x) else round(x) + 1 }
round_down <- function(x) { if(round(x) <= x) round(x) else round(x) - 1 }


### Mathematical
get_prevalence <- function(x) { howmany(x == 1) / len(x)}

# z = 1.96
# c_mean + z*sd = ci_upper
# sd = (ci_upper -  c_mean)/z
get_sd_from_ci95 <- function(ci_lower, ci_upper, volume = NA) {
  if(len(volume) == 0 || any(is.na(volume))){
    try_log_error("Calling get_sd_from_ci95 without volume! break")
    quit(-1)
  }
  c_mean = (ci_lower + ci_upper)/2
  c_sd = ((ci_upper - c_mean)/(1.96))*sqrt(volume)
  return(c_sd)
}

get_ci_95_from_mean_and_sd <- function(c_mean, c_sd, volume) {
  return(c(c_mean - (c_sd/sqrt(volume))*1.96,
           c_mean + (c_sd/sqrt(volume))*1.96))
}

###
# save yourself some typing
try_table <- function(x) {
  return(table(x,useNA = 'always'))
}

try_hist <- function(x) {
  return(hist(x, breaks = 1000))
}

logit2prob <- function(v1){ gtools::inv.logit(v1) }
prob2logit <- function(v1) {
  res = gtools::logit(v1)
  if(any(is.infinite(res)) ){
    try_log_warn("@prob2logit -  encountered Infinity after logit transform.... setting value to 999999")
     res = unlist(map(res, function(x) {  if(is.infinite(x)) 999999 else x  }))
  }
  else res
  }

convert_results_list_to_df <- function(metrics_results_list, metrics_nms = NULL) {
  first_non_null_res_list <- NULL
  if(len(metrics_results_list) == 0)
    return(NULL)
  for(i in 1:len(metrics_results_list)) {
    if(len(metrics_results_list[[i]]) != 0){
      first_non_null_res_list <- metrics_results_list[[i]]
      break
    }
  }

  results_df <- as.data.frame(matrix(ncol = len(names(first_non_null_res_list)), nrow =0))
  if(nrow(results_df) == 0)
    return(results_df)

  if(is.null(metrics_nms))
    metrics_nms = names(first_non_null_res_list)

  for(m_nm in metrics_nms)
    if(!m_nm %in% cns(results_df))
      results_df[,m_nm] = -999999999

  colnames(results_df) <- metrics_nms
  for(i in 1:len(metrics_results_list)){
    if(len(metrics_results_list[[i]]) != 0) {
      temp_ma3x <- matrix(ncol = len(metrics_nms), nrow =1)
      colnames(temp_ma3x) <- metrics_nms
      for(metric_nm in metrics_nms){
        #print(metric_nm)
        metric_val <- metrics_results_list[[i]][[metric_nm]]
        temp_ma3x[1, metric_nm] <- ifelse(is.null(metric_val), -999999, metric_val)
      }
      results_df <- rbind(results_df, temp_ma3x)
    }

  }
  colnames(results_df) <- metrics_nms
  results_df = as.data.frame(apply(results_df, c(1,2), function(x) if(len(x) == 0 || x == -999999) NA else x))
  return(results_df)
}

convert_df_list_to_df <- function(df_l) {
  f_df = NULL
  f_cns = cns(df_l[[1]])
  for(c_df in df_l) {
    c_cns = cns(c_df)
    if(!all(intersect(c_cns, f_cns) == union(c_cns, f_cns))){
      fcns_missing = setdiff(c_cns, f_cns) # fcns missing something
      ccns_missing = setdiff(f_cns, c_cns) # curr cns missing something
      if(len(fcns_missing) > 0 ) {
        for(c_missing in fcns_missing) f_df[,c_missing] = NA
      }
      if(len(ccns_missing) > 0) {
        for(c_missing in ccns_missing) c_df[,c_missing] = NA
      }
    }
    f_df = rbind(f_df, c_df)
  }
  colnames(f_df) = f_cns
  return(f_df)
}

grep_on_colnames <- function(df1, str1) { unlist(map(cns(df1), function(x) { grep(str1, x, value = T) })) }

get_lps_from_model_and_data <- function(m1, df1) {
  is_lr = any(class(m1) ==  "lm" )
  predict_fn =  if(is_lr) predict else NULL
  df_in = df1
  if(is_lr) {
    df_in = try_factorize_df(df1)
  }
  m_preds = predict_fn(m1, newdata = df_in)
  m_logits = m_preds #predict method on LR gives bak logits already!
  m_logits
}

numerize_df <- function(df1){
  for(cur_col in cns(df1)) { if(class(df1[,cur_col]) != "numeric") df1[,cur_col] = as.numeric(df1[,cur_col])-1 }
  df1
}

try_factorize_df <- function(df1, categoricals) {
  per_col_type <- map(df1[1,], class)
  types_of_columns <- uniq( per_col_type )
  df_factors <- df1
  for(cur_col in colnames(df1)){
    cur_col_vect <- df1[,cur_col]
    is_col_categorical <- (cur_col %in% categoricals) || is_col_binary(df1[,cur_col])
    if(class(cur_col_vect[1]) %in% c("integer", "double", "numeric") && ( is_col_binary(cur_col_vect) || is_col_categorical ))
      df_factors[,cur_col] <- factor(df1[,cur_col])
  }
  df_factors
}

is_col_binary <- function(col_vect) {
  if(len(uniq(col_vect)) > 3)
    return(F)
  val1 = sort(uniq(col_vect))[1]
  val2 = sort(uniq(col_vect))[2]
  count_0s <- len(which(col_vect == val1))
  count_1s <- len(which(col_vect == val2))
  count_nas <- len(which(is.na(col_vect)))
  count_nulls <- len(which(is.null(col_vect)))
  count_0s + count_1s + count_nas + count_nulls == len(col_vect)
}

is_categorical_column_from_name <- function(col1) {
  len(grep(pattern = "^.*_cat_\\d_\\d?$", x = col1)) == 1
}

is_categorical_column <- function(c_v) {
  return(len(uniq(c_v)) <= MAX_N_CATEGORIES_IN_DATASET)
}

remove_trailing_one <- function(str1) {
  if(len(str1) == 0)
    return(NULL)
  if(len(str1) > 1){
    res = umap(str1, function(x) { remove_trailing_one(x) })
    return(res)
  }

  if(len(grep("^.*1$", str1)) == 1)
    return(remove_last_char(str1))
  return(str1)
}

remove_last_char <- function(str1, n = 1) {
  stringr::str_sub(str1, 1, stringr::str_length(str1) - n)
}

remove_suffix_from_categorical_names <- function(c_var) {
  if(len(c_var) == 0)
    return(NULL)
  if(len(c_var) > 1){
    res = umap(c_var, function(x) { remove_suffix_from_categorical_names(x) })
    return(res)
  }
  if(!is_categorical_column(c_var)) {
    if(len(grep(pattern = "^.*\\d$", x = c_var)) == 1)
      c_var = stringr::str_sub(c_var, 1,stringr::str_length(c_var) -1 )
    return(c_var)
  }
  cat_strt_idx = regexpr(pattern = "_cat_",  c_var )[[1]]
  return(stringr::str_sub(c_var, 1, cat_strt_idx-1))
}

get_all_categoricals_names_like_current <- function(cat_names, df1) {
  if(len(cat_names) == 0)
    return(NULL)
  if(len(cat_names) > 1) {
    res = umap(cat_names, function(x) { get_all_categoricals_names_like_current(x, df1) })
    return(uniq(res))
  }
  col1 = cat_names
  if(is_categorical_column(col1)) {
    cat_name = stringr::str_sub(col1, end  = strlen(col1)-3)
    return(cns(df1)[grep(cns(df1), pattern = cat_name)])
  }
  else
    return(col1)
}

are_NAs_present_df_list <- function(dfl1) { any(unlist(map(1:len(dfl1), function(x){any(is.na(dfl1[[x]]))})) == T) }

sample_class_balanced <- function(df1, balance_on = "outcome_of_interest", balance_ratio = c(0.036), sample_size = NULL, verbose = F) {
  balance_on_remaining = NULL
  balance_ratio_ramaining = NULL
  if(is.null(sample_size)) sample_size = nrow(df1)
  if(len(balance_on) > 1) {
    balance_on_remaining = balance_on[-1]
    balance_on = balance_on[1]

    balance_ratio_ramaining = balance_ratio[-1]
    balance_ratio = balance_ratio[1]
  }
  if(verbose)
    try_log_debug("balancing on %s with ratio %0.2f (current ratio %0.2f)", balance_on, balance_ratio, get_prevalence(df1[,balance_on]))
  cases = which(df1[,balance_on] == 1)
  controls = which(df1[,balance_on] == 0)
  n_cases = len(cases)
  n_controls = len(controls)
  c_cases = sample(cases, round(sample_size*balance_ratio), replace = T)
  c_controls = sample(controls, round(sample_size*(1-balance_ratio)), replace = T)
  df_cases = df1[c_cases,]
  df_controls = df1[c_controls,]
  res = rbind(df_cases, df_controls)
  if(len(balance_on_remaining) > 0) {
    res =rbind(sample_class_balanced(df_cases, balance_on_remaining, balance_ratio_ramaining) ,
               sample_class_balanced(df_controls, balance_on_remaining, balance_ratio_ramaining))
  }

  return(res)
}

try_get_class_for_each_column_df <- function(df1) {
  clss <- c()
  for(c_nm in cns(df1)){
    c_v = df1[,c_nm]
    c_cls = class(c_v)
    clss <- c(clss, c_cls)
  }
  names(clss) = cns(df1)
  return(clss)
}

try_factorize_col <- function(c_v, cateogries_grouping = list()){
  nex_uniq_val = max(c_v, na.rm  = T) + 1

  if(len(cateogries_grouping) == 0)
    return(as.factor(c_v))

  excluded_grp_indxs = c()
  for(i in 1:len(cateogries_grouping)) {
    c_grp = cateogries_grouping[[i]]
    c_grp_nm = names(cateogries_grouping)[i]
    c_indxs = which(c_v %in% c_grp)
    if(len(c_indxs) == 0){
      try_log_warn("No records found for group %s", c_grp_nm)
      excluded_grp_indxs = c(excluded_grp_indxs, i)
      next
    }
    c_v[c_indxs] = nex_uniq_val
    nex_uniq_val = nex_uniq_val + 1
  }
  grp_names = names(cateogries_grouping)
  res = as.factor(c_v)
  if(len(grp_names) > 0)
    if(len(excluded_grp_indxs) > 0)
      levels(res) = grp_names[-excluded_grp_indxs]
    else
      levels(res) = grp_names
  return(res)
}


########## Functions used for k-fold CV  ##########

last <- function(x) {return(x[len(x)])}

# col_pos = only used when cols_to_remove = NULL
try_place_column_in_order_df <- function(df1, cols_to_remove, col_to_add, col_nm, col_pos = -1) {
  orig_cns = cns(df1)
  var_pos = col_pos
  if(!is.null(cols_to_remove))
    var_pos = which(cns(df1) %in% cols_to_remove)[1]

  if(!is.null(cols_to_remove))
    df1[,cols_to_remove] = NULL

  df1[,col_nm] = col_to_add
  new_cns = c (orig_cns[1:var_pos-1], col_nm,  orig_cns[ (var_pos +len(cols_to_remove)) :  (len(orig_cns) ) ] )
  df1 = df1[,new_cns]
  return(df1)
}

### e.g. - treat nyha class 1 + 2 as a single category
try_merge_boolean_variables_into_one <- function(df1, vars_to_merge, merged_var_name) {
  merged_var = umap(1:nrow(df1), function(x) {
    as.numeric(any(df1[x,vars_to_merge]))
  })
  df1 = try_place_column_in_order_df(df1, vars_to_merge, merged_var, merged_var_name)
  return(df1)
}

### e.g. - treat nyha class 1 + 2 as a single category AND make NYHA into continuous var
try_convert_ohe_categorical_to_cont <- function(df1, ohe_vars, cont_var_nm) {
  cont_var = umap(1:nrow(df1),  function(x) {
    if(any(is.na(df1[x, ohe_vars])) ) NA else ( which(df1[x, ohe_vars] == 1) - 1)
  })
  df1 = try_place_column_in_order_df(df1, ohe_vars, cont_var, cont_var_nm)
  return(df1)
}

try_add_removed_category_from_ohe <- function(df1, vars_ohe, removed_var_nm) {
  removed_var = umap(1:nrow(df1),  function(x) {
    as.numeric(all(df1[x, vars_ohe] == F))
  })
  df1 = try_place_column_in_order_df(df1, NULL, removed_var, removed_var_nm, which(cns(df1) == last(vars_ohe) ) + 1 )

  return(df1)
}



get_imputed_dfs_train_and_test <- function(cur_df_train, cur_df_test, logfile_prefix="", cur_fold=-1){
  test_df <- get_imputed_dfs_list(cur_df_test, n_copies = 1)[[1]]
  imputed_dfs_train_list = get_imputed_dfs_list(cur_df_test, n_copies = 5)
  ## intersect columns between test & train dfs
  clean_res <- clean_imputed_dfs(test_df, imputed_dfs_train_list)
  test_df <- clean_res$test
  imputed_dfs_train_list <- clean_res$train
  rm(clean_res)
  return(list(train=imputed_dfs_train_list, test=test_df))
}

# stratified by outcome only
create_kfolds <- function(df1, n_folds = CV_N_FOLDS, override.excluded.predictors = NULL, verbose = T, add_final_fold = T) {
  folds <- splitTools::create_folds(df1$outcome_of_interest, k = n_folds, type = "stratified")
  train_dfs <- list()
  test_dfs <- list()
  cntr = 0
  df1 = remove_excluded_predictors(df1, override.excluded.predictors = override.excluded.predictors, verbose = verbose)
  for(cur_fold in folds){
    cntr = cntr + 1
    cur_tt <- get_test_train_df_for_fold(cur_fold, df1)
    train_dfs[[cntr]] = remove_excluded_predictors(cur_tt$train, override.excluded.predictors = override.excluded.predictors, verbose = verbose)
    test_dfs[[cntr]] = remove_excluded_predictors(cur_tt$test, override.excluded.predictors = override.excluded.predictors, verbose = verbose)

  }
  # add one last fold with all data in training - to generate final model
  if(add_final_fold) {
    train_dfs[[len(train_dfs)+1]] = df1
    test_dfs[[len(test_dfs)+1]] = df1
  }
  return(list(folds = folds, train_dfs =  train_dfs, test_dfs = test_dfs))
}
# stratified by outcome per center
create_kfold_strat_center <- function(df1, n_folds = CV_N_FOLDS, override.excluded.predictors = NULL) {
  train_dfs <- list()
  test_dfs <- list()
  df1 = remove_excluded_predictors(df1, override.excluded.predictors = override.excluded.predictors, verbose = F)
  #note: assume you dont use the folds - you may have to fix it in try_LCO.R l43 ande try_kfold_cv l19
  par_df = partition_df_from_pipeline(df1)
  for(i in names(par_df)) {
    c_train_dfs = list()
    c_test_dfs = list()
    c_df = par_df[[i]]
    tmp = create_kfolds(c_df, n_folds = n_folds, override.excluded.predictors, verbose = F, add_final_fold = F)
    counter = 1
    for(cur_fold in tmp$folds) {
      cur_tt <- get_test_train_df_for_fold(cur_fold, c_df)
      cur_tt$train = remove_excluded_predictors(cur_tt$train, override.excluded.predictors = override.excluded.predictors, verbose = F)
      cur_tt$test = remove_excluded_predictors(cur_tt$test, override.excluded.predictors = override.excluded.predictors, verbose = F )
      c_train_dfs[[counter]] = cur_tt$train
      c_test_dfs[[counter]] = cur_tt$test
      counter = counter + 1
    }
    # first time adding
    if(as.numeric(i) == 1){
      train_dfs = c_train_dfs
      test_dfs = c_test_dfs
    } #not first time
    else {
      for(j in 1:len(c_test_dfs)){
        train_dfs[[j]] = rbind(train_dfs[[j]], c_train_dfs[[j]])
        test_dfs[[j]] = rbind(test_dfs[[j]], c_test_dfs[[j]])
      }
    }

  }

  # add one last fold with all data in training - to generate final model
  train_dfs[[len(train_dfs)+1]] = df1
  test_dfs[[len(test_dfs)+1]] = df1
  return(list(folds = NULL, train_dfs = train_dfs, test_dfs = test_dfs))
}


# train/test pairs
#' @param UOF_VN - unit of federation varname
create_folds_lcoa <- function(df1, n_folds = -1, override.excluded.predictors = NULL, UOF_VN = "centrum") {
  train_dfs <- list()
  test_dfs <- list()
  center_numeric = as.numeric(levels(df1[, UOF_VN]))[df1[, UOF_VN]]
  fake_center = max(center_numeric)+1 # for final model training on all data
  for(i in c(sort(uniq(center_numeric)), fake_center )) {
    c_train_df = df1
    c_test_df = df1
    if(i != fake_center) {
      c_train_df = df1[df1[, UOF_VN] != i,]
      c_test_df = df1[df1[, UOF_VN] == i,]
    }
    train_dfs[[i]] = remove_excluded_predictors(c_train_df, override.excluded.predictors = override.excluded.predictors, verbose = F)
    test_dfs[[i]] = remove_excluded_predictors(c_test_df, override.excluded.predictors = override.excluded.predictors, verbose = F)
  }
  # this is doing it [center][fold] , you want it to become [fold]
  # solution - when it does it for a given center i, you add to each fold those rows

  return(list(folds = NULL, train_dfs = train_dfs, test_dfs = test_dfs))
}

get_test_train_df_for_fold <- function(cur_fold, df1){
  neg_cur_fold <- -1*cur_fold
  df_train <- df1[cur_fold, ]
  df_test <- df1[neg_cur_fold, ]
  return(list(train = df_train, test = df_test))
}



compute_prediction_errors <- function(predictions, true_vals) { abs(true_vals - predictions) }

rename_column_df <- function(df1, old_nm, new_nm) {
  if(!old_nm %in% cns(df1)) {
    try_log_error("rename_column_df - no column with name %s", old_nm)
    exit(-1)
  }
  old_cns = cns(df1)
  old_cm_idx = which(old_cns == old_nm)

  cols_infront = if(old_cm_idx > 1) old_cns[1:(old_cm_idx-1)] else NULL
  cols_inback = if(old_cm_idx < len(old_cns)) old_cns[(old_cm_idx+1):len(old_cns)] else NULL

  colnames(df1) = c( cols_infront, new_nm,  cols_inback)
  return(df1)
}

remove_strings_matching_regex <- function(v_strings, rgx) {
  to_remove = umap(v_strings, function(x) { grep(rgx, x, value = T) })
  res = v_strings
  if(len(to_remove) > 0 )
    res = setdiff(v_strings, to_remove)

  if(len(res) == 0)
    res = c()

  return(res)
}

remove_excluded_predictors <- function(df1, override.excluded.predictors = NULL, verbose = T){
  if(len(override.excluded.predictors) > 0)
    excluded.predictors = override.excluded.predictors
  if(len(excluded.predictors) > 0 && any(cns(df1) %in% excluded.predictors) ){
    vars_to_exclude = cns(df1)[which(cns(df1) %in% excluded.predictors)]
    vars_to_keep = cns(df1)[-which(cns(df1) %in% excluded.predictors)]
    if(verbose) try_log_info("removing following vars from use in model derivation: \n\t\t -> [ %s ]; <- \n",
                 paste(vars_to_exclude, collapse = ","))
    df1 = df1[, vars_to_keep]
  }
  df1
}
