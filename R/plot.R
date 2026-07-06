#' Plot Kaplan-Meier survival curves by risk group
#'
#' @param surv_time Numeric vector of survival times (months).
#' @param surv_status Numeric vector of event status (1 = event, 0 = censored).
#' @param risk_group Character vector of risk groups ("High" / "Low").
#' @param title Plot title (default: "MMimmuneRisk").
#' @param ylab Y-axis label (default: "Overall survival").
#' @param palette Colors for High and Low risk curves.
#' @param break_x X-axis break interval in months.
#' @return A ggsurvplot object (invisibly).
#' @export
plot_km_risk <- function(surv_time, surv_status, risk_group,
                         title = "MMimmuneRisk",
                         ylab = "Overall survival",
                         palette = c("#74a3b2", "#c873a8"),
                         break_x = 12) {

  df <- data.frame(
    os_time   = surv_time,
    os_status = as.numeric(surv_status),
    group     = factor(risk_group, levels = c("High", "Low"))
  )

  fit  <- survival::survfit(survival::Surv(os_time, os_status) ~ group, data = df)
  diff <- survival::survdiff(survival::Surv(os_time, os_status) ~ group, data = df)

  pval <- 1 - stats::pchisq(diff$chisq, df = 1)
  pval_txt <- if (pval < 0.001) "P<0.001" else paste0("P=", sprintf("%.03f", pval))

  p <- survminer::ggsurvplot(
    fit, data = df,
    pval          = pval_txt,
    risk.table    = TRUE,
    surv.median.line = "hv",
    palette       = palette,
    legend.labs   = c("High risk", "Low risk"),
    legend.title  = title,
    ylab          = ylab,
    xlab          = "Time (months)",
    censor.shape  = 124, censor.size = 2,
    conf.int      = TRUE,
    break.x.by    = break_x
  )

  p$plot <- p$plot +
    ggplot2::annotate("text", x = Inf, y = 0.15, label = "Log-rank",
                      size = 5, hjust = 1.1, vjust = -1)

  p$plot$theme$axis.title.x$size <- 16
  p$plot$theme$axis.title.y$size <- 16
  p$plot$theme$axis.text.x$size  <- 14
  p$plot$theme$axis.text.y$size  <- 14
  p$plot$theme$legend.text$size  <- 14

  p$table <- p$table +
    ggplot2::xlab(NULL) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line    = ggplot2::element_blank()
    )

  pt <- p$plot / p$table + patchwork::plot_layout(nrow = 2, heights = c(2, 0.8))
  print(pt)
  invisible(p)
}
