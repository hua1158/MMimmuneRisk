##############################################################################
# prepare_model.R
#
# Run this AFTER training your RSF + plsRcox model.
#
# Prerequisites:
#   1. Install the package: devtools::install("path/to/MMimmuneRisk")
#   2. Have these objects in your R session:
#      - fit        : fitted plsRcox model
#      - cox        : fitted rms::cph model (Surv ~ rs + Age)
#      - exp_train  : training expression matrix (genes x samples)
#      - pre_var    : character vector of genes used in plsRcox training
##############################################################################

library(MMimmuneRisk)

# Example — adapt object names to match your session:
setup_model(
  plsrcox_model   = fit,         # plsRcox model object
  cox_model        = cox,         # rms::cph model object
  train_expr       = exp_train,   # gene expression matrix (genes x samples)
  model_genes      = pre_var,     # all genes used by plsRcox
  risk_cutoff      = 0.8527786,   # surv_cutpoint optimal cutoff
  highlight_genes  = c("FABP5", "BIRC5", "RASGRP3", "FGF13", "RAC3",
                        "BLNK", "MBL2", "CSF2", "SOD1")
)

# After running this, the model is embedded in the installed package.
# Test it:
#   run_dynnom()
