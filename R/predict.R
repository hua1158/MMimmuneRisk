#' Predict immune risk score for multiple myeloma patients
#'
#' Computes the PCDI risk score from bulk RNA-seq expression data using
#' the pre-trained RSF + plsRcox model.
#'
#' @param expr A numeric matrix or data.frame of gene expression.
#'   Accepts either genes-in-rows (gene x sample) or genes-in-columns
#'   (sample x gene) format. Must contain all genes used in training
#'   (stored by \code{setup_model}).
#' @param genes_in_rows Logical. TRUE if rows are genes and columns are
#'   samples (default). Set FALSE if rows are samples.
#' @return A named numeric vector of risk scores, one per sample.
#' @export
predict_mmrisk <- function(expr, genes_in_rows = TRUE) {

  model     <- get_model("plsrcox_model")
  model_genes <- get_model("model_genes")

  expr <- as.data.frame(expr)

  if (genes_in_rows) {
    missing <- setdiff(model_genes, rownames(expr))
    if (length(missing) > 0)
      stop("Missing ", length(missing), " model genes. First 10: ",
           paste(head(missing, 10), collapse = ", "))
    expr <- as.data.frame(t(expr[model_genes, , drop = FALSE]))
  } else {
    missing <- setdiff(model_genes, colnames(expr))
    if (length(missing) > 0)
      stop("Missing ", length(missing), " model genes. First 10: ",
           paste(head(missing, 10), collapse = ", "))
    expr <- expr[, model_genes, drop = FALSE]
  }

  expr_scaled <- as.data.frame(scale(expr))

  if (nrow(expr_scaled) == 1) {
    train_params <- get_model("scale_params")
    for (g in model_genes) {
      expr_scaled[[g]] <- (expr[[g]] - train_params$mean[g]) / train_params$sd[g]
    }
  }

  rs <- as.numeric(predict(model, type = "lp", newdata = expr_scaled))
  names(rs) <- rownames(expr_scaled)
  rs
}


#' Classify patients into High/Low risk groups
#'
#' @param expr Gene expression matrix (see \code{\link{predict_mmrisk}}).
#' @param genes_in_rows Logical. TRUE if rows are genes (default).
#' @param age Optional numeric vector of patient ages (years). If provided,
#'   survival probabilities at 1/2/3 years are estimated via the Cox nomogram.
#' @param cutoff Risk score cutoff. Default uses the training-set optimal value.
#' @return A data.frame with columns: sample, risk_score, risk_group, and
#'   optionally surv_1y, surv_2y, surv_3y if age is provided.
#' @export
classify_mmrisk <- function(expr, genes_in_rows = TRUE, age = NULL,
                            cutoff = NULL) {

  rs <- predict_mmrisk(expr, genes_in_rows = genes_in_rows)

  if (is.null(cutoff)) {
    cutoff <- tryCatch(get_model("risk_cutoff"), error = function(e) DEFAULT_CUTOFF)
  }

  result <- data.frame(
    sample     = names(rs),
    risk_score = rs,
    risk_group = ifelse(rs > cutoff, "High", "Low"),
    stringsAsFactors = FALSE
  )

  if (!is.null(age)) {
    if (length(age) != nrow(result))
      stop("Length of 'age' must match number of samples.")

    cox_mod <- tryCatch(get_model("cox_model"), error = function(e) NULL)
    if (!is.null(cox_mod)) {
      nd <- data.frame(rs = rs, Age = age)
      sf <- rms::survest(cox_mod, newdata = nd, times = c(12, 24, 36))
      if (nrow(result) == 1) {
        result$surv_1y <- sf$surv[1]
        result$surv_2y <- sf$surv[2]
        result$surv_3y <- sf$surv[3]
      } else {
        result$surv_1y <- sf$surv[, 1]
        result$surv_2y <- sf$surv[, 2]
        result$surv_3y <- sf$surv[, 3]
      }
    } else {
      message("Cox model not found; survival probabilities not computed.")
    }
  }

  rownames(result) <- NULL
  result
}


#' List the genes required by the model
#'
#' @return Character vector of gene names needed for prediction.
#' @export
model_genes <- function() {
  get_model("model_genes")
}

#' List the 9 key signature genes
#'
#' @return Character vector of the 9 PCDI highlight genes.
#' @export
signature_genes <- function() {
  HIGHLIGHT_GENES
}
