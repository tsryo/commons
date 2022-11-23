# Author: T.Y.
# Date created: 19.10.2021
# Script contents:
### 7 -> perform backwards step-wise elimination based on AIC on each of the

### STEP 7

create_formula_allvars <- function(df1) {
  create_formula_somevars(df1)
}


create_formula_somevars <- function(df1, vars = c()) {
  if(len(vars) == 0){
    vars <- colnames(df1)
    vars <- vars[-which(vars == "outcome_of_interest")]
  }
  eq_rhs <- paste(vars, collapse = '+')
  eq_whole_str <- paste( c("outcome_of_interest", eq_rhs), collapse = '~')
  eq_whole_str
}

# @param agree_stren (between-set agreement ratio)  -> tuning param
#             : when set to 1, returns vars selected in all sets from imputed_dfs_train_list
#             : when set to 0.5, returns vars selected in at least 50% of all sets from imputed_dfs_train_list
#             : when set to 0, returns vars selected in at least 1 of the sets from imputed_dfs_train_list
# rename: select_variables_AIC_imp_sets
select_variables_AIC_imp_sets <- function(imputed_dfs_train_list, fold_cntr=0, agree_stren = 1.0, verbose = T){
  if(agree_stren < 0 || agree_stren > 1.0){
    try_log_error(0,"@select_variables_AIC_imp_sets- invalid parameter value for agree_stren")
    return(NULL)
  }
  cur_vars_selected <- c()
  aic_tm_start <- as.numeric(Sys.time())
  aic_tm_end <- as.numeric(Sys.time())
  vars_occur_cntr <- list()
  for(cn in cns(imputed_dfs_train_list[[1]])) {
    if(is.factor(imputed_dfs_train_list[[1]][,cn])) {
      for(cn_lv in levels(imputed_dfs_train_list[[1]][,cn])) {
        cn_val = sprintf("%s%s", cn, cn_lv)
        vars_occur_cntr[[cn_val]] = 0
      }
    } else
      vars_occur_cntr[[cn]] = 0
    }

  for(cur_imp_ds in 1:len(imputed_dfs_train_list)){
    aic_tm_end <- as.numeric(Sys.time())
    if(cur_imp_ds != 1 && verbose)
      try_log_debug("@select_variables_AIC_imp_sets- AIC for imputed DS #%d Elapsed time: %0.1fs", cur_imp_ds-1, aic_tm_end - aic_tm_start)
    aic_tm_start <- as.numeric(Sys.time())

    try_log_info("@select_variables_AIC_imp_sets- [FOLD #%d] Running AIC with imputed ds #%d", fold_cntr, cur_imp_ds)

    xx <- imputed_dfs_train_list[[cur_imp_ds]]

    if(cur_imp_ds == 1)
      cur_vars_selected <- colnames(xx)
    if(len(cur_vars_selected) > 3)
      cur_vars_selected <- perform_backwards_AIC(xx)
    for(cv in cur_vars_selected) vars_occur_cntr[[cv]] = unlist(vars_occur_cntr[[cv]]) + 1

    if(verbose)
      try_log_debug("@select_variables_AIC_imp_sets- cur vars selected: \n ['%s']\n", paste(cur_vars_selected, collapse =   "', '"))
  }
  n_occurs_thresh = round_up(len(imputed_dfs_train_list) * agree_stren)
  n_occurs_thresh = ifelse(n_occurs_thresh == 0, 1, n_occurs_thresh)

  vars_selected <- c()
  for(i in 1:len(vars_occur_cntr)) {
    if(len(vars_occur_cntr) == 0)
      break
    if(vars_occur_cntr[[i]] >= n_occurs_thresh )
      vars_selected <- c(vars_selected, names(vars_occur_cntr)[i])
  }

  return(vars_selected)
}

fix_naming_issue_in_step_output <- function(sout, valid_col_nms) {
  coef_labs <- names(sout$coefficients)
  # remove "(Intercept)"
  coef_labs <- coef_labs[which(coef_labs != "(Intercept)")]
  is_col_nm_original <- function(col_nm) { col_nm %in% valid_col_nms }
  is_col_nm_char_suffixed <- function(col_nm, char1) { substr(col_nm, nchar(col_nm), nchar(col_nm)) == char1 }
  is_col_nm_num_suffixed <- function(col_nm) {
    res = F
    for(i in 0:9) res = res || all(is_col_nm_char_suffixed(col_nm, as.character(i)))
    return(res)
    }
  is_col_nm_1_valid <- function(col_nm) {substr(col_nm, 1, nchar(col_nm)-1) %in% valid_col_nms }
  for(col_nm in coef_labs){
    if(is_col_nm_original(col_nm))
      nullvoid <- sprintf("Orginal found, %s", col_nm)
    else{
      if( is_col_nm_num_suffixed(col_nm)   && is_col_nm_1_valid(col_nm)){
        #          print(sprintf("Found column to fix, %s", col_nm))
        coef_labs[which(coef_labs == col_nm)] <- substr(col_nm, 1, nchar(col_nm)-1)
      }
      else
        try_log_error("@fix_naming_issue_in_step_output - found column that can not fix! (please investigate!), %s", col_nm)
    }
  }
  uniq(coef_labs)
}


## step 7 - backwards AIC
perform_backwards_AIC <- function(df1, out_var = "outcome_of_interest", return_step_out = F, n_steps = NULL) {

  vars_to_remove <- unlist(map(cns(df1), function(cn){
    if(class(df1[,cn]) == "factor"){
      if(len(levels(df1[,cn])) < 2 ||
         (len(table(df1[,cn])) < 7 && any(table(df1[,cn]) == 0) )
         )
        return(cn)
    }
    return(NULL)
  }))
  vars_to_remove = remove_strings_matching_regex(vars_to_remove, "center_column_cat_\\d+_")
  if(len(vars_to_remove) > 0){
    try_log_info("@perform_backwards_AIC- REMOVING VARS where len(levels(df1[1,var])) < 2 [ %s ];",
                 paste(vars_to_remove, collapse = ","))
    df1 = df1[, -which(cns(df1) %in% vars_to_remove)]
  }


  eq_whole_str <- create_formula_allvars(df1)
  allvar_model <- create_logreg_model(df1, eq_whole_str, hardcode_control = T)

  step_out = NULL
  if(!is.null(n_steps))
    step_out <- step(allvar_model, direction = "backward" , trace=0, steps = n_steps)
  else
    step_out <- step(allvar_model, direction = "backward" , trace=0)

  if(return_step_out)
    return(step_out)

  variables_selected <- names(step_out$coefficients)[-1]

  variables_selected
}

# Impute only if NAs present in @df1 # used for debug
PRIVATE_maybe_impute_and_run_backwards_AIC <- function(df1, vselc_bs_agree_stren = 1.0) {
  df_list = list(df1)
  if(any(is.na(df1)))
    df_list = get_imputed_dfs_list(df1, n_copies = 1)
  select_variables_AIC_imp_sets(df_list, fold_cntr = -1, agree_stren = vselc_bs_agree_stren)
}
