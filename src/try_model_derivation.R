# Author: T.Y.
# Date created: 17.11.2021

# Wrapper functions

derive_model_multi_imputation <- function(cur_df_train, cur_df_test, vselc_bs_agree_stren = 1.0,
                                       logfile_prefix = "", cur_fold = -1) {
  PRIVATE__derive_model_multi_imputation(cur_df_train, cur_df_test, vselc_bs_agree_stren = vselc_bs_agree_stren,
                                      logfile_prefix = logfile_prefix, cur_fold = cur_fold)
}

# Logic functions

#' @param vselc_bs_agree_stren =  AIC between-set agreement strength (1.0 = vars appearing in 5/5 sets (intersect), 0.4 = vars appearing in 2/5 sets, 0.0 = appearing in at least 1/5 sets (union))
PRIVATE__derive_model_multi_imputation <- function(cur_df_train, cur_df_test, vselc_bs_agree_stren = 1.0,
                                                logfile_prefix = "", cur_fold = -1) {
  imputed_dfs = get_imputed_dfs_train_and_test(cur_df_train, cur_df_test, logfile_prefix=logfile_prefix, cur_fold=cur_fold)
  test_df = imputed_dfs$test
  imputed_dfs_train_list = imputed_dfs$train

  vars_selected <- select_variables_AIC_imp_sets(imputed_dfs_train_list, fold_cntr = cur_fold, agree_stren = vselc_bs_agree_stren)

  if(len(vars_selected) == 0){
    try_log_crit(0, "@PRIVATE__derive_model_multi_imputation- No variables selected from backwards AIC")
    stop(-2)
  }

  try_log_debug("|%s| FOLD #%d AIC selected vars : \n\t\t -> [ %s ]; <- ", logfile_prefix, cur_fold, paste(vars_selected, collapse = ","))
  lr_res <- NULL

  lr_res <- train_pool_test_pp_logreg(imputed_dfs_train_list = imputed_dfs_train_list,
                                      test_df = test_df, vars_selected =  vars_selected, fold_cntr = cur_fold, logfile_prefix = logfile_prefix)


  return(list(lr_res = lr_res, aic_res = vars_selected ))
}
