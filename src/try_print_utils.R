# Author: T. Yordanov
# Date: 15.10.2021


get_current_date_printable <- function(){
  format(Sys.time(),"%b %d %X %Y")
}

###### LOG functions
try_log <- function(severity, format_str, p1 = "", p2 = "", p3 = "", p4 = "", p5 = "", p6= "", p7 = "") {
  if(which(LOGGING_SEVERITIES == severity) < which(LOGGING_SEVERITIES == logging.level))
    return()
  format_str = if(class(format_str) == "numeric" && format_str <= 0) "%s" else format_str

  cur_filepath = rstudioapi::getActiveDocumentContext()$path
  cur_filename = tail(strsplit(cur_filepath, "/")[[1]], n = 1)
  cur_filename = if(len(cur_filename) == 0) "CONSOLE" else cur_filename
  cur_date = get_current_date_printable()

  format_str = sprintf("[%s] |%s| [%s] %s", severity, cur_date, cur_filename, format_str)

  n_params = stringr::str_count(format_str, "%") - 2*stringr::str_count(format_str, "%%")

  eval_str = ""
  if(n_params > 0)
  {
    eval_str = ", p1"
    if(n_params > 1)
      for(i in 2:n_params)
        eval_str = sprintf("%s, p%d", eval_str, i)
  }
  eval_str = sprintf("writeLines(sprintf(format_str %s))", eval_str)

  eval(parse(text = eval_str))
}

try_log_trace <- function(format_str, p1 = "", p2 = "", p3 = "", p4 = "", p5 = "", p6= "", p7 = ""){
  try_log(severity = "TRACE", format_str, p1, p2, p3, p4, p5, p6, p7)
}

try_log_debug <- function(format_str, p1 = "", p2 = "", p3 = "", p4 = "", p5 = "", p6= "", p7 = ""){
  try_log(severity = "DEBUG", format_str, p1, p2, p3, p4, p5, p6, p7)
}

try_log_info <- function(format_str, p1 = "", p2 = "", p3 = "", p4 = "", p5 = "", p6= "", p7 = ""){
  try_log(severity = "INFO", format_str, p1, p2, p3, p4, p5, p6, p7)
}

try_log_warn <- function(format_str, p1 = "", p2 = "", p3 = "", p4 = "", p5 = "", p6= "", p7 = ""){
  try_log(severity = "WARN", format_str, p1, p2, p3, p4, p5, p6, p7)
}

try_log_error <- function(format_str, p1 = "", p2 = "", p3 = "", p4 = "", p5 = "", p6= "", p7 = ""){
  try_log(severity = "ERROR", format_str, p1, p2, p3, p4, p5, p6, p7)
}

try_log_crit <- function(format_str, p1 = "", p2 = "", p3 = "", p4 = "", p5 = "", p6= "", p7 = ""){
  try_log(severity = "CRITICAL", format_str, p1, p2, p3, p4, p5, p6, p7)
}

start_logging_to_file <- function(logfile_prefix){
  # create logs dir if not present
  dir.create("./logs", showWarnings = F)
  sink_fn = sprintf("logs/%srunlog%s.log", logfile_prefix ,stringr::str_replace_all(Sys.time(), ":", "-"))
  try_log_info("Logging output to %s", sink_fn)
  sink(sink_fn)
  print(Sys.time())
}

######### df functions

print_df_with_rounded_numbers <- function(df1, round_to = 2) {
  for(c1 in cns(df1)){
    if(class(df1[,c1]) == "numeric") {
      df1[,c1] = round(df1[,c1], round_to)
    }
  }
  df1
}

print_NAs_and_counts_df <- function(df1) {
  offending_columns = c()
  for(i in cns(df1)) {
    na_count = tail(table(df1[,i], useNA = "always"), n = 1)
    if(na_count > 0) {
      try_log_warn(0, i)
      try_log_warn(0, na_count)
      offending_columns = c(offending_columns, i)
    }
  }
  offending_columns
}


######### results print functions
print_simple_LCO_pp_pooled <- function(PP_reses, model_type = "lr") {
  res_accessor = sprintf("PP_%s", model_type)
  c_res = PP_reses[[res_accessor]]
  for(i in 1:len(c_res)) cat(sprintf("\t-> %s = %s;\n",toupper(names(c_res)[i]),  round(c_res[[names(c_res)[i]]], 2)   ))
}

print_simple_LCO_pp_per_center <- function(PP_reses, model_type = "lr") {
  res_accessor = sprintf("PP_per_center_%s", model_type)
  c_res = PP_reses[[res_accessor]]
  for(i in names(c_res)) {
    c_c_res = c_res[[i]]
    print(i)
    for(j in 1:len(c_c_res))
      cat(sprintf("\t-> %s = %s;\n",toupper(names(c_c_res)[j]),  round(c_c_res[[names(c_c_res)[j]]], 2)   ))
  }
}

print_simple_LCO_pp_per_center_compare <- function(PP_reses, PP_reses2, model_type = "lr", decimal_precision = 2) {
  res_accessor = sprintf("PP_per_center_%s", model_type)
  c_res = PP_reses[[res_accessor]]
  c_res2 = PP_reses2[[res_accessor]]
  for(i in names(c_res)) {
    c_c_res = c_res[[i]]
    c_c_res2 = c_res2[[i]]
    print(i)
    for(j in names(c_c_res)) {
      diff_f =  round(c_c_res[[j]] - c_c_res2[[j]], decimal_precision)
      precision_str = sprintf("+%%0.%df", decimal_precision)
      diff_str = umap(diff_f, function(x) {if(is.na(x) || is.infinite(x) || is.null(x) || x >= 0) sprintf(precision_str, x) else as.character( x)})
      cat(sprintf("\t-> %s = %s (%s);\n",toupper(j),  round(c_c_res[[j]], 2),diff_str
      ))
    }
  }
}

######### misc
print_vect_byline <- function(c_v) {
  cat(sprintf("%s\n", c_v))
}
map_column_shortname_to_longname <- function(c_nm) {
  COLUMN_NAME_MAPPINGS[COLUMN_NAME_MAPPINGS$shortname %in% c_nm,]$longname
}

try_table_per_column <- function(df1, warn_thresh = 17, quiet =F) {
  if(!quiet) try_log_debug("@try_table_per_column for df with dim = [ %s ]", paste0(dim(df1),collapse = " : "))
  thresh_vars = c()
  for(c_col in cns(df1)) {
    if(!quiet) cat("\n##########\n\n")
    c_tbl = table(df1[,c_col], useNA = "always")
    if(!quiet) cat(c_col, "\n")
    if(len(c_tbl) >= 15) {
      if(!quiet) print(try_round(quantile(df1[,c_col], probs = seq(0,1,0.1), na.rm = T),3) )
    }
    else {
      if(!quiet) cat(paste0(names(c_tbl),collapse = "\t"))
      if(!quiet) cat("\n")
      if(!quiet) cat(paste0(c_tbl,collapse = "\t"))

      out_tbl = table(df1[,c(c_col, "outcome_of_interest")])
      if(any(out_tbl  <= warn_thresh)) {
        if(!quiet) print(out_tbl)
        thresh_vars = c(thresh_vars, c_col)
        if(!quiet) cat("\n^^^warn_thresh^^^\n")
      }
    }
    if(!quiet) cat("\n***********\n")
  }
  return(thresh_vars)
}


try_print_class_for_each_column_df <- function(df1, quiet = F){
  should_be_categoricals = c()
  should_be_binaries = c()
  for(c_nm in cns(df1)){
    c_v = df1[,c_nm]
    c_cls = class(c_v)
    n_uniq_vls = nuniq(c_v)
    should_be_binary = n_uniq_vls <= 3
    should_be_categorical = n_uniq_vls <= 17
    c_txt = if(should_be_binary || should_be_categorical) "factor" else "numeric"
    if(!quiet) cat(paste0(c_nm, " = " , c_cls, " [  suggested = ", c_txt, "  ]\n"))
    if(should_be_binary)
      should_be_binaries = c(should_be_binaries, c_nm)
    else if(should_be_categorical)
      should_be_categoricals = c(should_be_categoricals, c_nm)
  }
  return(list(should_be_categoricals = should_be_categoricals, should_be_binaries = should_be_binaries))
}

