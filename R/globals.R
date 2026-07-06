#' @title MMimmuneRisk package constants
#' @description Internal constants for the MMimmuneRisk classifier.
#' @keywords internal

.mmrisk_env <- new.env(parent = emptyenv())

HIGHLIGHT_GENES <- c("FABP5", "BIRC5", "RASGRP3", "FGF13", "RAC3",
                      "BLNK", "MBL2", "CSF2", "SOD1")

DEFAULT_CUTOFF <- 0.8527786

.onLoad <- function(libname, pkgname) {
  model_path <- system.file("model", "model_objects.rda", package = pkgname)
  if (file.exists(model_path)) {
    load(model_path, envir = .mmrisk_env)
  }
}

get_model <- function(name) {
  obj <- .mmrisk_env[[name]]
  if (is.null(obj)) {
    stop(
      "Model not found. Run setup_model() first to embed your trained model.\n",
      "See ?setup_model for details.",
      call. = FALSE
    )
  }
  obj
}

has_model <- function() {
  !is.null(.mmrisk_env[["plsrcox_model"]])
}
