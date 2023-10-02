# Author: T.Y.
# Date created: 19.10.2021
# Script contents:
###
##############  code for building the logreg models
#@hardcode_control - for some reason calling step (AIC) with a control list, takes the value as a string literal (i.e. object verbose/maxit not found error...)
create_logreg_model <- function(df1, eq_str, starting_coefficients = NULL, maxit = 500, verbose = F, hardcode_control = F){
  verbose = if(verbose == F) 0 else 1
  control_l = list(maxit = maxit, trace = verbose)
  logreg_model = NULL
  if(hardcode_control)
    logreg_model <- glm(data = df1, formula = as.formula(eq_str), family = binomial, control = list(maxit = 500, trace = 0))
  else
    logreg_model <- glm(data = df1, formula = as.formula(eq_str), family = binomial, control = control_l, start = starting_coefficients)
  logreg_model
}

train_pool_test_pp_logreg <- function(imputed_dfs_train_list, test_df, vars_selected, fold_cntr, logfile_prefix ="") {

  pooled_lrm = train_and_pool_logreg(imputed_dfs_train_list, vars_selected, fold_cntr, logfile_prefix)

  lrm_test_res = test_logreg(pooled_lrm, test_df)

  lrm_test_pred_probs = lrm_test_res$preds
  lrm_test_outcomes = lrm_test_res$true_vals

  perf_metrics <- NULL
  tryCatch( {
    perf_metrics <- compute_perf_metrics(predicted_probs = lrm_test_pred_probs,
                                         actual_outcomes = lrm_test_outcomes)

  },
  error = function(cond){
    try_log_warn("[LR] Could not compute perf metrics on %s ", logfile_prefix)
  },
  finally = {}
  )
  ##
  ## compute metrics
  ret_list <- perf_metrics
  ret_list[["coeffs"]] <- pooled_lrm$coefficients
  ret_list[["coef_nms"]] <- names(pooled_lrm$coefficients)
  ret_list[["test_predictions"]] <- lrm_test_pred_probs
  ret_list[["test_outcomes"]] <- lrm_test_outcomes

  try_log_debug(0, "*****  @train_pool_test_pp_logreg complete  *****")
  return(ret_list)
}

train_and_pool_logreg <- function(imputed_dfs_train_list, vars_selected, fold_cntr=-1, logfile_prefix, with_pooling = T,
                                  starting_coefficients = NULL, maxit = 500, verbose = F) {
  # derive 5 logregs using vars selected on each imputed ds
  list_lrms <- list()
  cur_lrm <- NULL
  for(cur_imp_ds_idx in 1:len(imputed_dfs_train_list)){
    cur_imp_df <-  imputed_dfs_train_list[[cur_imp_ds_idx]]

    #go over factor columns , remove levels not selected
    # compute simple name for leveled vars selected
    vars_selected_simple_names = c()
    vars_selected_simple_names = vars_selected

    eq_str <- create_formula_somevars(cur_imp_df, vars_selected_simple_names)
    nvars = len(vars_selected_simple_names)
    if(nvars == 0)
      try_log_error("@train_and_pool_logreg - 0 vars selected", nvars)

    if(ncol(cur_imp_df) > 2 && nvars > 0){
      cur_lrm  <- create_logreg_model(cur_imp_df, eq_str, starting_coefficients = starting_coefficients , maxit = maxit, verbose = verbose)
      list_lrms[[cur_imp_ds_idx]] <- cur_lrm
    }
    else
      list_lrms[[cur_imp_ds_idx]] <- NULL
  }
  # pooled lrm models (Rubin's Rules)
  if(!with_pooling)
    return(list_lrms)

  rubins_pooled_lrms <- mice::pool(list_lrms)

  try_log_debug("@train_logreg- FOLD #%d RUBIN's pooled coefficients of logreg", fold_cntr)
  try_log_debug("@train_logreg- VAR = %s; COEFF = %0.4f;", rubins_pooled_lrms$pooled$term ,rubins_pooled_lrms$pooled$estimate)

  # use body of lrm, fill in coefficients from pooled
  pooled_lrm <- derive_lrm_from_pooled_obj(rubins_pooled_lrms, cur_lrm)

  return(pooled_lrm)
}

### MODEL DERIVATION
derive_lrm_from_pooled_obj <- function(pooled_obj, lrm_obj) {
  n_terms <- len(pooled_obj$pooled$term)
  if(len(lrm_obj$coefficients) != n_terms)
    try_log_error(0, "@derive_lrm_from_pooled_obj- len(cur_lrm$coefficients) != n_terms; please investigate ")
  else
    for(i in 1:n_terms){
      lrm_obj$coefficients[pooled_obj$pooled$term[i]] <- pooled_obj$pooled$estimate[i]
    }
  lrm_obj
}

test_logreg <- function(lrm1, test_df) {
  lrm_test_preds <- predict(lrm1, newdata = test_df)

  lrm_test_pred_probs <- logit2prob(lrm_test_preds)
  lrm_test_outcomes <- as.numeric(test_df[,"outcome_of_interest"])-1
  return(list(true_vals = lrm_test_outcomes, preds=lrm_test_pred_probs))
}


