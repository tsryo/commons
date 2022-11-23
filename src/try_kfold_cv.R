# Author: T.Y.
# Date created: 19.10.2021

# @vselc_bs_agree_stren =  AIC between-set agreement strength (1.0 = vars appearing in 5/5 sets (intersect), 0.4 = vars appearing in 2/5 sets, 0.0 = appearing in at least 1/5 sets (union))
try_kfold_cv <- function(df1, log_to_file = T, logfile_prefix = "", vselc_bs_agree_stren = 1.0) {
  logfile_prefix = paste0(logfile_prefix, "-10CV")

  if(log_to_file)
    start_logging_to_file(logfile_prefix)

  folds_created = create_kfolds(df1)

  vars_selected_per_fold <- list()
  lr_CV_restuls_list <- list()

  fold_tm_start <- as.numeric(Sys.time())
  fold_tm_end <- as.numeric(Sys.time())

  for(cur_fold in 1:len(folds_created$folds)){
    cur_df_train <- folds_created$train_dfs[[cur_fold]]
    cur_df_test <- folds_created$test_dfs[[cur_fold]]

    # timing & printing
    fold_tm_end <- as.numeric(Sys.time())
    if(cur_fold != 1)
      try_log_info("|%s| FOLD #%d Elapsed time: %0.1fs", logfile_prefix,  cur_fold-1, fold_tm_end - fold_tm_start)
    fold_tm_start <- as.numeric(Sys.time())
    try_log_debug("|%s| CV status: Entering fold #%d", logfile_prefix, cur_fold)
    strat_s_res = NULL


    strat_s_res <- derive_model_multi_imputation(cur_df_train, cur_df_test, vselc_bs_agree_stren, logfile_prefix, cur_fold)


    vars_selected_per_fold[[cur_fold]] <- strat_s_res$aic_res
    lr_CV_restuls_list[[cur_fold]] <- strat_s_res$lr_res
  }

  ## derive final model
  try_log_debug(0,"***** START derive_model_multi_imputation final *****")
  df1 = remove_excluded_predictors(df1)
  strat_s_res <- derive_model_multi_imputation(df1, df1, vselc_bs_agree_stren, logfile_prefix, -1)
  vars_selected_per_fold[["final"]] <- strat_s_res$aic_res
  lr_CV_restuls_list[["final"]] <- strat_s_res$lr_res


  # convert metric lists into dfs to compare results
  logreg_CV_results_df <- NULL
  metric_nms <- c("aucroc", "aucpr", "bss", "bs", "cali_slope", "cali_intercept",
                  "cali_slope_ci_lo", "cali_slope_ci_hi", "cali_intercept_ci_lo",
                  "cali_intercept_ci_hi", "cali_p_val")

  logreg_CV_results_df <- convert_results_list_to_df(lr_CV_restuls_list, metric_nms)
  if(nrow(logreg_CV_results_df) == 0)
    return(logreg_CV_results_df)
  logreg_CV_results_df$model_type = "lr"
  try_log_info("@try_kfold_cv- mean AUC-ROC logreg = %0.5f", mean(logreg_CV_results_df$aucroc))


  try_log_debug(0, "@try_kfold_cv- COMPLETE")
  print(Sys.time())
  if(log_to_file)
    sink()
  try_log_debug(0, "@try_kfold_cv- COMPLETE")

  return(list(lr=logreg_CV_results_df))
}



######   Run code #######

REPEAT_ANALYSIS_N_TIMES = 1
RUN_CODE = F

if(RUN_CODE) {
  df3 = datasets::infert
  for(c_iter in 1:REPEAT_ANALYSIS_N_TIMES) {
    c_logfile_prefix =  sprintf("my-cv-experiment-%d", c_iter)
    kfold_cv_res = try_kfold_cv(df3, log_to_file = T, logfile_prefix = c_logfile_prefix, vselc_bs_agree_stren = 1.0)
    sink()
  }
}
