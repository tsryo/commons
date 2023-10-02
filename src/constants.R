# Author: T. Yordanov
# Date: 15.10.2021



##### ALIASES
len = length
cns = colnames
uniq = unique
e = exp(1)
beep = function() { beepr::beep(11) }

######## Operational constants
DIR_ROOT = getwd()
LOGGING_SEVERITIES <- c("TRACE", "DEBUG", "INFO", "WARN", "ERROR", "CRITICAL")

SAVED_FILENAME_TEMPLATES = list( vars_selected = "AIC-vars-%s-fold-%d.RData",
                                 model = "model-%s-fold_%d-%s.RData",
                                 results_list = list( kfoldcv = "results-list-%s-10CV-%s.RData" , lco = "results-list-%s-%s.RData"),
                                 test_df = "%s-test_df_fold%d.RData",
                                 train_df = "%s-train_df_list_fold%d.RData")





