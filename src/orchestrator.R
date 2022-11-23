# Author: T. Yordanov
# Date: 15.10.2021


setwd(sprintf("%s/..",dirname(rstudioapi::getActiveDocumentContext()$path)))
print(getwd())
if(length(unlist(strsplit(getwd(), split = "/"))) < 6)
  stop(-1)

# read & set flags from properties file
FLAGS_FILE = "flags.properties"
flags_l = properties::read.properties(FLAGS_FILE)
mvbutils::extract.named(flags_l)
MAX_N_CATEGORIES_IN_DATASET = as.numeric(MAX_N_CATEGORIES_IN_DATASET)
MIN_N_POS_OUTCOMES_IN_DATASET = as.numeric(MIN_N_POS_OUTCOMES_IN_DATASET)
MIN_N_VARS_SELECTED_PER_MODEL_ENSEMBLE = as.numeric(MIN_N_VARS_SELECTED_PER_MODEL_ENSEMBLE)
excluded.predictors = unlist(strsplit(excluded.predictors, ","))
random.seed = as.numeric(random.seed)

set.seed(random.seed)

source("../commons/src/constants.R")
source("../commons/src/gbm_calibrate_plot_fork.R")
source("../commons/src/try_functional_utils.R")
source("../commons/src/try_grouping_onehot.R")
source("../commons/src/try_imputation.R")
source("../commons/src/try_io_commons.R")
source("../commons/src/try_kfold_cv.R")
source("../commons/src/try_logreg.R")
source("../commons/src/try_model_derivation.R")
source("../commons/src/try_print_utils.R")
source("../commons/src/try_perf_metrics.R")
source("../commons/src/try_variable_selection.R")


# cast flags
try_log_info(0, "******************************* *** *******************************")
try_log_info(0, "*** BEGIN PRINT FLAGS VALUES ***")
try_log_info("\t\t - %s = %s\n", names(flags_l), flags_l)
try_log_info(0, "*** END PRINT FLAGS VALUES ***")
try_log_info(0, "******************************* *** *******************************")


