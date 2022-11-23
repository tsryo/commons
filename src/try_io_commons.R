# Author: T. Yordanov
# Date: 17.11.2021


load_with_assign <- function(fp){
  return(eval(parse(text=load(fp))))
}

save_incremental_name <- function(objToSave, saveFilepath) {
  pathElems = unlist(stringr::str_split(saveFilepath,  "/"))
  fn1 = pathElems[len(pathElems)]
  pathElems = pathElems[-len(pathElems)]
  f_name_and_ext = unlist(stringr::str_split(fn1, stringr::fixed( ".") ))
  cur_index = 1
  fn_new = sprintf("%s-%s.%s", f_name_and_ext[1], cur_index, f_name_and_ext[2])
  while(file.exists(saveFilepath)) {
    saveFilepath = paste(c(pathElems,fn_new), collapse  = "/")
    cur_index =  cur_index + 1
    fn_new = sprintf("%s-%s.%s", f_name_and_ext[1], cur_index, f_name_and_ext[2])
  }
  try_log_debug("@save_incremental_name - Saving file %s", saveFilepath)
  save(objToSave, file =  saveFilepath)
}

# @results_type - oneof(LCO, 10CV)
load_model_predictions_and_outcomes <- function(sink_prefix = "proba11", model_type = "lr", results_type = "LCO",
                                                excluding_final_model_10CV = T, load_train_results_also = F,
                                                root_dir = getwd()) {
  root_dir = "G:/divjk/kik/NHR-Onderzoek/Abu Hanna/R scripts(TY)/PhD/SRP-rewrite"
  fn_template_map = list("LCO" = SAVED_FILENAME_TEMPLATES$results_list$lco,
                         "10CV" = SAVED_FILENAME_TEMPLATES$results_list$kfoldcv)
  m_res = load_with_assign(paste0(root_dir, "/out/",sprintf(fn_template_map[[results_type]], sink_prefix, model_type)))
  test_preds_n_outcomes = NULL
  for(i in 1:len(m_res)) {
    if(results_type == "10CV" && i == 11)
      break
    cur_preds <- m_res[[i]]$test_predictions
    #cur_preds_train <- m_res[[i]]$train_predictions todo: gotta do it the hard way..
    test_preds_n_outcomes <- rbind(test_preds_n_outcomes,as.data.frame(list( predictions = cur_preds, true_val = m_res[[i]]$test_outcomes, fold = rep(i, len(cur_preds)))))
    if(load_train_results_also)
      try_log_error("@load_model_predictions_and_outcomes: called with load_train_results_also = T, but this is not implemented yet!")
  }
  return(test_preds_n_outcomes)
}


load_train_and_test_dfs <- function(sink_prefix = "proba11", folds = 1:CV_N_FOLDS, root_dir = getwd()) {
  train_dfs <- list()
  test_dfs <- list()
  for(c_fold in folds){
    c_train_dfs <- load_with_assign(paste0(root_dir, "/out/",sprintf(SAVED_FILENAME_TEMPLATES$train_df, sink_prefix, c_fold)))
    c_test_df <- load_with_assign(paste0(root_dir, "/out/",sprintf(SAVED_FILENAME_TEMPLATES$test_df, sink_prefix, c_fold)))

    train_dfs[[c_fold]] = c_train_dfs
    test_dfs[[c_fold]] = c_test_df
  }
  # train = list of lists of dfs (multiple-imputations), test = list of dfs (1-imputation)
  return(list(train = train_dfs, test= test_dfs))
}