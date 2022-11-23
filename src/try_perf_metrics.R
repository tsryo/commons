# Author: T.Y.
# Date created: 19.10.2021

compute_perf_metrics <- function(predicted_probs, actual_outcomes){
  ## compute metrics
  c_aucroc <- try_get_aucroc(predictions = predicted_probs, true_vals = actual_outcomes)

  c_aucroc_95ci = try_get_aucroc_ci(predictions = predicted_probs, true_vals = actual_outcomes)

  c_aucpr <- try_get_aucpr(predictions = predicted_probs,
                           true_vals = actual_outcomes)

  # BSS + BS
  c_bs <- verification::brier(obs = actual_outcomes, pred = predicted_probs)$bs
  c_bss <- verification::brier(obs = actual_outcomes, pred = predicted_probs)$ss

  # CALI
  cxr <- try_get_calibration(predicted_probs, actual_outcomes, verbose = F)

  try_log_debug("AUC-ROC = %0.5f", c_aucroc)
  try_log_debug("AUC-PR = %0.5f", c_aucpr)
  try_log_debug("BS = %0.5f", c_bs)
  try_log_debug("BSS = %0.5f", c_bss)
  try_log_debug("CALI intercept = %0.5f", cxr$cal.intercept)
  try_log_debug("CALI slope = %0.5f", cxr$cal.slope)

  return(list(aucroc=c_aucroc, aucroc_95ci = c_aucroc_95ci,
              aucpr=c_aucpr, bs=c_bs, bss=c_bss, cali_slope=cxr$cal.slope, cali_intercept=cxr$cal.intercept,
              cali_slope_ci_lo = cxr$ci.cal.slope[1], cali_slope_ci_hi = cxr$ci.cal.slope[2],
              cali_intercept_ci_lo = cxr$ci.cal.intercept[1], cali_intercept_ci_hi = cxr$ci.cal.intercept[2],
              cali_eo_ratio = cxr$eo_ratio,
              cali_p_val = cxr$p.value.groups
  ))
}

try_get_bs <- function(predictions, true_vals) {
  return(verification::brier(obs = true_vals, pred = predictions)$bs)
}

try_get_bss <- function(predictions, true_vals) {
  return(verification::brier(obs = true_vals, pred = predictions)$bss)
}

try_get_aucroc <- function(predictions, true_vals, return_auc_only = T) {
  if(return_auc_only)
    return(PRROC::roc.curve(
                            scores.class0 = predictions,
                            weights.class0 =  true_vals,
                            )$auc
           )
  else
    return(pROC::roc(true_vals, predictions))
}

try_get_aucroc_ci <- function(predictions, true_vals) {
  roc1 = suppressMessages(pROC::roc( true_vals, predictions ))
  pROC::ci(roc1)
}

## AUC-PR
#https://www.biostat.wisc.edu/~page/rocpr.pdf TODO: cite if using this method!
try_get_aucpr <- function(predictions, true_vals) {
  PRROC::pr.curve(
    scores.class0 = predictions,
    weights.class0 =  true_vals,
  )$auc.davis.goadrich
}

# p is predicted probability and y is the observed event (0/1). This would be usually applied on a test set.
## understanding of outputs:
##  original.intercept   -> this is what they call  the "intercept" in rms::val.prob plot
##  calmodel.intercept   -> this is what we get from fitting a glm on the offset
cox.recalibration <- function(p, y, use_logits_in_glms = T, verbose = T) {
  if(!verbose) try_log_debug = function(...) {}
  # Preprocess
  p <- ifelse(p <= 0, 0, ifelse(p >= 1, 1, p))
  eo_ratio = sum(y) / sum(p)
  logits = p
  if(use_logits_in_glms)
    logits <- log(p/(1-p))
  if(any(logits == -Inf)) logits[which(logits == -Inf)] = -99
  transformer_fn = function(x){
    if(use_logits_in_glms)
      x
    else
      inv.logit(x)
  }
  # make decile groups

  if(len(table(quantile(p,(1:9)/10))) < 9){
    try_log_warn(0, "@cox.recalibration- Can not calculate calibration - non-unique breaks in quantiles!")
    return( list(inter= NULL,
                 ci.inter= NULL,
                 cal.intercept= NULL,
                 ci.cal.intercept =  NULL,
                 cal.slope =  NULL,
                 ci.cal.slope =  NULL,
                 cal.test.slope =  NULL,
                 ci.cal.test.slope =  NULL,
                 sig =  NULL,
                 sig.test =  NULL,
                 cal.groups =  NULL,
                 cal.nogroups =  NULL,
                 p.value.groups= NULL))
  }
  group <- cut(p, c(-Inf, quantile(p,(1:9)/10), Inf))

  # Coarse view
  try_log_debug(0, cat("Total events: ", sum(y), "total predictions: ", round(sum(p)), "\n"))

  calmodel <- glm(y ~ logits, family=binomial)
  original.intercept <- transformer_fn(calmodel$coefficients[1])
  suppressMessages(ci.original.intercept <- transformer_fn(confint(calmodel, verbose=F)[1,]))
  try_log_debug(0,cat("Original intercept is: ", round(original.intercept, 3), "With 95%CI: ",
      round(ci.original.intercept[1], 3), " to: ", round(ci.original.intercept[2], 3), "\n"))

  # what is intercept when we have slope = 1?
  # https://stats.stackexchange.com/questions/175144/what-does-the-1-in-the-r-formula-y-1-mean
  # In most R regression packages, y ~ 1 means "fit an intercept only".
  # so you want to get the
  calmodel.intercept <- glm(y ~ 1 + offset(logits), family=binomial, control = glm.control(maxit=25))
  cal.intercept <- transformer_fn(calmodel.intercept$coefficients)

  xxx =summary(calmodel.intercept)
  intercept_ci = c(xxx$coefficients[1,1]-1.96*xxx$coefficients[1,2],
                   xxx$coefficients[1,1]+1.96*xxx$coefficients[1,2])

  #suppressMessages(ci.cal.intercept <- transformer_fn(as.numeric(confint(calmodel.intercept))))
  suppressMessages(ci.cal.intercept <- transformer_fn(as.numeric(intercept_ci)))
  # does zero appear in conf interval?
  significance <- ifelse(ci.cal.intercept[1] < 0 & ci.cal.intercept[2] > 0, "", "*")
  try_log_debug(0, cat("Calibration intercept is: ", round(cal.intercept, 3), "With 95%CI: ",
      round(ci.cal.intercept[1], 3), " to: ", round(ci.cal.intercept[2], 3), significance, "\n"))

  # this should be original.slope?
  cal.slope <- transformer_fn(calmodel$coefficients[2])
  suppressMessages(ci.cal.slope <- transformer_fn(confint(calmodel)[2,]))
  # To test for significance of calibration slope deviating from 1 we can fit another glm and
  # see if coefficient of the free variable logits is different than 0 after introducing offset
  calmodel.testslope <- glm(y ~ logits + offset(logits), family=binomial)
  cal.test.slope <- transformer_fn(calmodel.testslope$coefficients[2])
  suppressMessages(ci.cal.test.slope <- confint(calmodel.testslope)[2,])
  significance.test.slope <- ifelse(ci.cal.test.slope[1] < 0 & ci.cal.test.slope[2] > 0, "", "*")
  try_log_debug(0, cat("Calibration slope is: ", round(cal.slope, 3), "With 95%CI: ", round(ci.cal.slope[1], 3),
      " to: ", round(ci.cal.slope[2], 3), significance.test.slope, "\n"))

  # Working with decile groups, due to [2]
  # A test for group effect (i.e., ??1 = ??2 = . = ??k =0) in model 3 is asymptotically equivalent
  # to the Hosmer-Lemeshow goodness-of-fit test, though not precisely identical for finite samples.
  #  A score test for group effect in the following model is the goodness-of-fit test that was first
  #   proposed by Tsiatis.

  cal.groups <- glm(y ~ -1 + group +offset(logits), family= binomial)
  cal.nogroups <- glm(y ~ -1 + offset(logits), family= binomial)
  test.groups <- anova(cal.groups, cal.nogroups, test="Chisq")
  p.value.groups <-  test.groups$"Pr(>Chi)"[2]
  try_log_debug(0, cat("P-valueof decile groups is: ", p.value.groups, " [Low p-value is associated with miscalibration]\n"))

  return(list(inter=original.intercept, ci.inter=ci.original.intercept, cal.intercept=cal.intercept, ci.cal.intercept = ci.cal.intercept,
              cal.slope = cal.slope, ci.cal.slope = ci.cal.slope, cal.test.slope = cal.test.slope, ci.cal.test.slope = ci.cal.test.slope,
              sig = significance, sig.test = significance.test.slope, cal.groups = cal.groups, cal.nogroups = cal.nogroups,
              eo_ratio = eo_ratio,
              p.value.groups=p.value.groups))
}
# aliasing
try_get_calibration <- cox.recalibration

## plot metrics
# @replace - only applicable when smoothed = F
plot_calibration_graph <- function(preds, outcomes, excluded_percentiles = c(0,100), title = "", smoothed = T, df = 4, replace = T,
                                   shade_col = "lightyellow", shade_density = NULL, line_par = list(col = "black"),
                                   rug_par = list(side = 1), xlim = "auto", y_lim  = NULL, add_quantiles_and_mean_text = T) {
  if(smoothed) # deprecared
    return(PRIVATE_plot_smooth_calibration_graph_loess(preds, outcomes, excluded_percentiles, title))
  return(PRIVATE_plot_calibration_graph_gbm(preds, outcomes, excluded_percentiles, title, df = df, replace =  replace,
                                            shade_col = shade_col, shade_density = shade_density, line_par = line_par, rug_par = rug_par,
                                            xlim = xlim, y_lim  = y_lim, add_quantiles_and_mean_text = add_quantiles_and_mean_text))

}

PRIVATE_plot_calibration_graph_gbm <- function(preds, outcomes, excluded_percentiles = c(0,100), title = "",
                                               df = 6, replace = T, shade_col = "lightyellow", shade_density = NULL,
                                               line_par = list(col = "black"), rug_par = list(side = 1),
                                               xlim = "auto", y_lim = NULL, add_quantiles_and_mean_text = T) {
  c_plot = NULL
  c_df <- as.data.frame(list(preds = preds, outs = outcomes))
  c_df = c_df[order(c_df$preds),]
  excl_percentile_lower = excluded_percentiles[1]
  excl_percentile_upper = excluded_percentiles[2]
  n_obs = nrow(c_df)
  c_df =  c_df[(n_obs*excl_percentile_lower/100 + 1):(n_obs*excl_percentile_upper/100), ]
  tryCatch(
    {
      xlim = if(xlim == "auto") c(0,max(c_df$preds)*2) else xlim
      y_lim = if(is.null(y_lim)) xlim else y_lim

      calibrate.plot(y = c_df$outs, p = c_df$preds,
                          xlab = "Predicted value", main = title, xlim = xlim, ylim = y_lim,
                          knots = NULL, df =df, replace = replace, shade.col = shade_col, shade.density = shade_density,
                          line.par = line_par, rug.par = rug_par, add_quantiles_and_mean_text = add_quantiles_and_mean_text)
    },error=function(cond){
      try_log_error(0,"@plot_calibration_graph - failed to plot calibration graph") # this usually happens when your preds are really not well calibrated
    }
    , finally = { }
  )
}

#' Deprecated
PRIVATE_plot_smooth_calibration_graph_loess <- function(preds , outcomes, excluded_percentiles = c(0,100),
                                                  label = "") {
  excl_percentile_lower = excluded_percentiles[1]
  excl_percentile_upper = excluded_percentiles[2]

  c_df = as.data.frame(list(preds = preds, true_vals = outcomes))
  c_df = c_df[order(c_df$preds),]
  n_obs = nrow(c_df)
  c_df =  c_df[(n_obs*excl_percentile_lower/100 + 1):(n_obs*excl_percentile_upper/100), ]


  plx = predict(loess(c_df$true_vals ~ c_df$preds, span = 1, degree = 1), se=T)
  c_df$tracer = 1:nrow(c_df)

  CI_Xs = c_df$preds
  CI_Xs = unlist(map(CI_Xs, function(x){c(x,x)}))
  CI_Ys = unlist(map(1:(len(plx$fit)),
                     function(x){
                       c(plx$fit[c_df$tracer[x]]  -qt(0.900,plx$df)*plx$se[c_df$tracer[x]],
                         plx$fit[c_df$tracer[x]] +qt(0.900,plx$df)*plx$se[c_df$tracer[x]])
                     }))

  were_lines_added = T
  iter_counter = 0
  while(were_lines_added) {
    ci_xs = c()
    ci_ys = c()
    if(iter_counter > 0)
      try_log_debug("@plot_LCO_cali_grid -  ### Adding lines... %d", iter_counter)
    iter_counter = iter_counter + 1
    were_lines_added = F
    for(i in 1:(len(CI_Xs)-1) ) {
      cx = CI_Xs[i]
      cxn = CI_Xs[i+1]
      cy = CI_Ys[i]
      cyn = CI_Ys[i+1]
      should_add_more_lines = abs(cx - cxn) > 0.02

      if(should_add_more_lines)
        were_lines_added =  T

      if(should_add_more_lines) {
        closest_x_indx = which( abs(c_df$preds - cx) == min(abs(c_df$preds - cx)))
        closest_y_val = plx$fit[c_df$tracer][closest_x_indx] # you want the CI to be centered aounr this y val

        ci_xs = c(ci_xs, cx, (cx+(cxn-cx)/3), (cx+2*(cxn-cx)/3), cxn)

        cy_range = abs(cy - cyn)
        cy_is_upper = cy > cyn

        if(cy_is_upper)
          ci_ys = c(ci_ys, closest_y_val + cy_range, closest_y_val - cy_range, closest_y_val + cy_range, closest_y_val - cy_range)
        else
          ci_ys = c(ci_ys, closest_y_val - cy_range, closest_y_val + cy_range, closest_y_val - cy_range, closest_y_val + cy_range)
      }
      else {
        ci_xs = c(ci_xs, cx, cxn)
        ci_ys = c(ci_ys, cy, cyn)
      }
    }
    ci_xs = c(ci_xs, CI_Xs[len(CI_Xs)])
    ci_ys = c(ci_ys, CI_Ys[len(CI_Ys)])

    CI_Xs = ci_xs
    CI_Ys = ci_ys
  }


  plot(c(0,max(c_df$preds, plx$fit)),c(0,max(c_df$preds, plx$fit)), col = "white", main = label,
       xlab= "Predicted value", ylab = "Observed average")
  lines(CI_Xs, CI_Ys,
        col = rgb(red = 255/255, green = 255/255, blue =102/255, alpha = 0.3), lty= 1, lwd = 5)
  lines(c_df$preds, plx$fit[c_df$tracer], col = "black", lty= 2, lwd = 3)
  lines(x = c(0,100), y = c(0,100), col ="red", lwd = 2)
}


# @results_df = preds, true_vals, center
plot_LCO_cali_grid <- function(logfile_prefix = "", model_type = "lr", excluded_percentiles = c(0,100), center_indexes = CENTER_FOLDS_NO68,
                               smooth_plot = T, df = 6) {
  results_df = load_model_predictions_and_outcomes(logfile_prefix, model_type = model_type, results_type = "LCO")
  colnames(results_df) = c("preds", "true_vals", "center")

  plot_cali_grid <- function() {
    do_plot_single_center = function(center) {
      centr_label = hidden_center_mapping[[center]]
      centr_indx = which(sort(center_indexes) == center)
      c_df = results_df[results_df$center == centr_indx,]
      c_df = c_df[order(c_df$preds),]
      n_obs = nrow(c_df)
      title = sprintf("Center %s", centr_label)
      plot_calibration_graph(preds = c_df$preds, outcomes = c_df$true_vals, title = title, smoothed = smooth_plot,
                             df = df, excluded_percentiles = excluded_percentiles, add_quantiles_and_mean_text = F)
    }
    par(mfrow=c(4,4))
    for(i in center_indexes){
      do_plot_single_center(i)
    }
  }
  plot_cali_grid()
}

#################################################

########################################################
plot_perf_metrics <- function(predicted_probs, actual_outcomes, fold_cntr = 0, output_path = "out/tmp",
                              with_aucroc = F, with_aucpr = F, with_cali = T,
                              filename_prefix = "tmp") {

  # AUC-ROC
  if(with_aucroc){
    png(sprintf("%s/fold-%d-%s-auc-roc.png", output_path, fold_cntr, filename_prefix), width = 640, height = 480)
    plot(try_plot_aucroc(predictions = predicted_probs,
                         true_vals = actual_outcomes,
                         label = sprintf("logreg fold #%d", fold_cntr)))
    dev.off()
  }

  if(with_aucpr){
    # AUC-PR
    png(sprintf("%s/fold-%d-%s-auc-pr.png", output_path, fold_cntr, filename_prefix), width = 640, height = 480)
    plot(try_plot_aucpr(predictions = predicted_probs,
                        true_vals = actual_outcomes,
                        label = sprintf("logreg fold #%d", fold_cntr)))
    dev.off()
  }

  if(with_cali){
    # CALI
    png(sprintf("%s/fold-%d-%s-cali.png", output_path, fold_cntr, filename_prefix), width = 640, height = 480)
    tryCatch(
      {
        gbm::calibrate.plot(y = actual_outcomes, p = predicted_probs,
                            xlab = sprintf("Predicted value\n logreg fold #%d", fold_cntr))
      },error=function(cond){
        print("error in plot cali!") # this usually happens when your preds are really not well calibrated
      }
      , finally = {}

    )
    dev.off()
  }

}
