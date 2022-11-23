# Author: T.Y.
# Date created: 26.09.2022

is_col_ohed_categorical_ohed <- function(x, col) {
  f_cat = paste0(col, "_cat_1_")
  if(!f_cat %in% cns(x))
    return(NULL)
  a_cats = c(f_cat, get_all_categoricals_names_like_current(f_cat, x))
  if(all(a_cats %in% cns(x)))
    return(a_cats)
  else return(NULL)
}


# - easy subgroup data analysis
# ::> inputs - dataframe, model predictions, subgroups definition (e.g. "geslacht" = splits into 2 groups male/female;
#                                                                  "leeftijd" + 30 + 5 = splits into X groups from ages < 30, 30-35, 35-40, etc. )
# "NYHA" + 2 = splits into 2 groups with grp1 = NYHA1&2 , grp2 = NYHA3&4 # default would split it per category, i.e. 4 in NYHA
# NOTE: handling missing ...
# make 5 imputed copies of the dataset
# compute analytics for each copy
# report mean metric results
#
# ::> outputs - per group - calibration [graph, Cox intercept & slope], AUC-ROC, AUC-PR, prevalence, BS, BSS, MSE
try_easy_subgroup_analysis <- function(x, y, preds, sg_def) {
  if(is.factor(y)) y = as.numeric(y) - max(as.numeric(levels(y)))
  sg_col = sg_def[1]
  sg_par1 = sg_def[2]
  sg_par2 = sg_def[3]

  sg_cols = is_col_ohed_categorical_ohed(x, sg_col)
  if(! is.null(sg_cols) ) {
    sg_cols =     c(sg_cols, paste0(sg_col,"_cat_", len(sg_cols), "_" ) )
    x[,last(sg_cols)] = 9999999999
  }
  is_categorical = len(sg_cols) > 0 || is_categorical_column( x[,sg_col] )


  is_binary = is.null(sg_cols) && is_col_binary(  x[,sg_col] )
  is_numeric = !is_categorical
  Xs = list()
  Ys = list()
  preds_l = list()

  if(is_binary) {
    for(unq_v in uniq( x[,sg_col] )) {
      c_indxs = which(x[,sg_col] == unq_v)
      Xs[[unq_v]] = x[c_indxs,]
      Ys[[unq_v]] = y[c_indxs]
      preds_l[[unq_v]] = preds[c_indxs]
    }
  }
  else if(is_categorical) {
    a_indxs = c()
    for(s_col in sg_cols) {
      for(unq_v in uniq( x[,s_col] )) {
        if(as.numeric(unq_v) == 0)
          next
        if(s_col == last(sg_cols)) {
          c_indxs = -a_indxs
        }
        else {
          c_indxs = which(x[,s_col] == unq_v)
          a_indxs = c(a_indxs, c_indxs)
        }
        if(unq_v >= 999999) # magic numbers are bad, todo!
          unq_v = 1

        c_lab = paste0(s_col,"=",unq_v)
        Xs[[c_lab]] = x[c_indxs,]
        Ys[[c_lab]] = y[c_indxs]
        preds_l[[c_lab]] = preds[c_indxs]
      }
    }
  }
  else if(is_numeric) {
    max_v = max(x[,sg_col])
    start_interval = as.numeric(sg_par1)
    interval_step = as.numeric(sg_par2)
    for(i in seq(start_interval, max_v, interval_step)){
      int_start = i
      int_end =  i+interval_step
      if(i + interval_step > max_v)
        int_end = max_v
      c_lab = paste0(sg_col, int_start, "-",int_end)
      c_indxs = intersect(which(x[,sg_col] >= int_start), which(x[,sg_col] < int_end))
      Xs[[c_lab]] = x[c_indxs,]
      Ys[[c_lab]] = y[c_indxs]
      preds_l[[c_lab]] = preds[c_indxs]
    }
  }


  metrics_df = NULL
  plots_l = list()
  for(c_nm in names(Xs)) {
    c_sgn = sprintf("%s-%s", sg_col, c_nm)
    c_x = Xs[[c_nm]]
    c_y = Ys[[c_nm]]
    c_preds = preds_l[[c_nm]]
# TODO: wrap this in try catch
    c_cali = NULL
    tryCatch(
      {
        c_cali = try_get_calibration(c_preds, c_y, verbose = F)
      },error=function(cond){
        try_log_error("@try_easy_subgroup_analysis - failed to compute calibration metrics for %s", c_nm)
      }
      , finally = { }
    )

    c_aucroc = NULL
    tryCatch(
      {
        c_aucroc = try_get_aucroc_ci(c_preds, c_y)
      },error=function(cond){
        try_log_error("@try_easy_subgroup_analysis - failed to compute auc-roc for %s", c_nm)
      }
      , finally = { }
    )

    c_aucpr = NULL
    tryCatch(
      {
        c_aucpr = try_get_aucpr(c_preds, c_y)
      },error=function(cond){
        try_log_error("@try_easy_subgroup_analysis - failed to compute auc-pr for %s", c_nm)
      }
      , finally = { }
    )

    c_bs = NULL
    c_bss = NULL
    tryCatch(
      {
        c_bs = try_get_bs(c_preds, c_y)
        c_bss = try_get_bss(c_preds, c_y)
      },error=function(cond){
        try_log_error("@try_easy_subgroup_analysis - failed to compute BS / BSS for %s", c_nm)
      }
      , finally = { }
    )


    plot.new()
    par(mfrow= c(1,2))
    plot_calibration_graph(c_preds, c_y, smoothed = F, title = sprintf("%s all-preds", c_sgn))

    # dev.off()
    plot_calibration_graph(c_preds, c_y, smoothed = F, excluded_percentiles = c(0,95), title = sprintf("%s bottom-95%%-preds", c_sgn))
    c_plot = recordPlot()
    dev.off()
    plots_l[[c_sgn]] = c_plot

    c_row = as.data.frame(list(c_slope = if(len(c_cali) == 0 || len(c_cali$cal.slope) == 0) NA else  c_cali$cal.slope,
                               c_int = if(len(c_cali) == 0 ||  len(c_cali$cal.intercept) == 0) NA else  c_cali$cal.intercept,
                               aucroc = if(len(c_aucroc[2]) == 0) NA else  c_aucroc[2],
                               aucpr = if(len(c_aucpr) == 0) NA else  c_aucpr,
                               bs = if(len(c_bs) == 0) NA else  c_bs,
                               bss = if(len(c_bss) == 0) NA else  c_bss,
                               subgroup = c_sgn,
                               nrow = len(c_y),
                               prevalence = howmany(c_y == 1)/len(c_y)))
    metrics_df = rbind(metrics_df, c_row)
  }
  return(list(metrics_df = metrics_df, plots = plots_l ))
}

all(df1$outcome_of_interest == df_preds$true_val)

# :: test for binary group
xxx = try_easy_subgroup_analysis(x = df1, y = df1$outcome_of_interest, pred = df_preds$predictions, sg_def = "geslacht" )

# 1 = female, 0 = male
try_round(xxx$metrics_df, 3)

xxx$plots$`geslacht-0`
xxx$plots$`geslacht-1`






# # :: test for categorical group
#
#
# xxx = try_easy_subgroup_analysis(x = df1, y = df1$outcome_of_interest, pred = df_preds$predictions, sg_def = "NYHA" )
#
# # 1 = female, 0 = male
# try_round(xxx$metrics_df, 3)
# xxx$plots$`NYHA-NYHA_cat_1_=1`



# # :: test for integer group
#
#
xxx = try_easy_subgroup_analysis(x = df1, y = df1$outcome_of_interest, pred = df_preds$predictions, sg_def = c("leeftijd", 30, 5) )

# 1 = female, 0 = male
try_round(xxx$metrics_df, 3)
xxx$plots$`leeftijd-leeftijd70-101`




# age_group_increment = 5
# age_group_start = 40
# plt_df = NULL
# df3 = df1
# df3$preds = df_preds$predictions
# for(c_age in seq(age_group_start, 95, age_group_increment)) {
#   c_rows = intersect( which(df3$leeftijd >= c_age),  which(df3$leeftijd < c_age + age_group_increment) )
#   try_log_debug("age range %0.1f - %0.1f , n records = %d", c_age, c_age + age_group_increment, len(c_rows))
#   m_mort = mean(as.numeric(df3[c_rows,]$outcome_of_interest)-1)
#   m_preds = mean(df3[c_rows,]$preds)
#
#   plt_df = rbind(plt_df, as.data.frame(list(age_start = c_age, m_mort = mean(m_mort)*100, m_pred = mean(m_preds)*100 )))
# }
# plot(plt_df$age_start, plt_df$m_mort, pch=19)
# points(plt_df$age_start, plt_df$m_pred, col ="red", pch=19)
# for(i in 1:nrow(plt_df)){
#   c_row = plt_df[i,]
#   l_col = if(c_row$m_mort - c_row$m_pred > 0 ) "yellow" else "red"
#   lines(c(c_row$age_start,c_row$age_start) , c(c_row$m_mort, c_row$m_pred) , col = l_col)
# }
#
#
#
# hist(df3$leeftijd, breaks = 1000)
# df3_cpy = df3
# pats_55_65 = df3[intersect( which(df3$leeftijd >= 55),  which(df3$leeftijd < 65) ),]
# pats_65_95 = df3[intersect( which(df3$leeftijd >= 65),  which(df3$leeftijd < 95) ), ]
# nrow(pats_55_65)
# nrow(pats_65_95)

# federated_lasso
df2 = df1[federated_lasso$row_id,]

xxx = try_easy_subgroup_analysis(x = df2, y = federated_lasso$outcome, pred = federated_lasso$prediction, sg_def = "geslacht" )
# 1 = female, 0 = male
try_round(xxx$metrics_df, 3)

xxx$plots$`geslacht-0`
xxx$plots$`geslacht-1`


