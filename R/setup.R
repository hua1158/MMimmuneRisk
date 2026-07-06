#' Save trained model objects into the package
#'
#' After training your RSF + plsRcox model and Cox nomogram model,
#' call this function to embed them into the installed package so that
#' \code{predict_mmrisk()} and \code{run_dynnom()} work out of the box.
#'
#' @param plsrcox_model The fitted plsRcox model object.
#' @param cox_model The fitted \code{rms::cph} model (Surv ~ rs + Age).
#' @param train_expr Training expression matrix (genes in rows, samples in
#'   columns). Used for computing scaling parameters for single-sample
#'   prediction. Should contain at least the model genes.
#' @param model_genes Character vector of gene names used in plsRcox training
#'   (i.e., the \code{pre_var} vector). If NULL, extracted from column names
#'   of the plsRcox training data.
#' @param risk_cutoff Optimal risk score cutoff (default: 0.8527786).
#' @param highlight_genes Character vector of key signature gene names
#'   for display (default: the 9 PCDI genes).
#' @return Invisible NULL.
#' @export
setup_model <- function(plsrcox_model,
                        cox_model,
                        train_expr,
                        model_genes = NULL,
                        risk_cutoff = 0.8527786,
                        highlight_genes = c("FABP5", "BIRC5", "RASGRP3",
                                            "FGF13", "RAC3", "BLNK",
                                            "MBL2", "CSF2", "SOD1")) {

  if (is.null(model_genes)) {
    model_genes <- tryCatch({
      colnames(plsrcox_model$dataX)
    }, error = function(e) {
      stop("Cannot extract gene names from model. ",
           "Please provide model_genes explicitly.", call. = FALSE)
    })
  }

  available <- intersect(model_genes, rownames(train_expr))
  if (length(available) < length(model_genes)) {
    warning(length(model_genes) - length(available),
            " model genes not found in train_expr. ",
            "Scale params computed for available genes only.")
  }

  train_sub <- as.data.frame(t(train_expr[available, , drop = FALSE]))
  scale_params <- list(
    mean = colMeans(train_sub),
    sd   = apply(train_sub, 2, sd)
  )

  pkg_dir <- system.file(package = "MMimmuneRisk")
  if (pkg_dir == "") pkg_dir <- find.package("MMimmuneRisk")

  model_dir <- file.path(pkg_dir, "model")
  if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)

  save(plsrcox_model, cox_model, scale_params,
       risk_cutoff, model_genes, highlight_genes,
       file = file.path(model_dir, "model_objects.rda"))

  .mmrisk_env[["plsrcox_model"]]    <- plsrcox_model
  .mmrisk_env[["cox_model"]]         <- cox_model
  .mmrisk_env[["scale_params"]]      <- scale_params
  .mmrisk_env[["risk_cutoff"]]       <- risk_cutoff
  .mmrisk_env[["model_genes"]]       <- model_genes
  .mmrisk_env[["highlight_genes"]]   <- highlight_genes

  message("Model saved successfully.")
  message("  Model genes: ", length(model_genes))
  message("  Risk cutoff: ", risk_cutoff)
  message("  Location: ", file.path(model_dir, "model_objects.rda"))
  invisible(NULL)
}
