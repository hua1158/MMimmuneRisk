#' Launch the dynamic nomogram Shiny app
#'
#' Opens an interactive web application where users can:
#' \itemize{
#'   \item Input gene expression + age for single-patient prediction
#'   \item Upload a bulk expression matrix for batch classification
#'   \item View dynamic nomogram, survival curves, and risk distributions
#' }
#'
#' @param host Host address (default "127.0.0.1").
#' @param port Port number (default: random available port).
#' @param launch.browser Logical, open in browser (default TRUE).
#' @return Does not return; runs the Shiny app.
#' @export
run_dynnom <- function(host = "127.0.0.1", port = NULL,
                       launch.browser = TRUE) {
  if (!has_model()) {
    stop("No model found. Run setup_model() first.", call. = FALSE)
  }

  app <- build_app()
  shiny::runApp(app, host = host, port = port,
                launch.browser = launch.browser)
}


build_app <- function() {

  mg        <- get_model("model_genes")
  cutoff    <- get_model("risk_cutoff")
  hi_genes  <- tryCatch(get_model("highlight_genes"), error = function(e) HIGHLIGHT_GENES)

  ui <- shiny::navbarPage(
    title = "MMimmuneRisk - Dynamic Nomogram",
    theme = NULL,

    # ---- Tab 1: Single Patient ----
    shiny::tabPanel(
      "Single Patient",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          width = 4,
          shiny::h4("Patient Information"),
          shiny::numericInput("age", "Age (years)", value = 65,
                              min = 20, max = 100, step = 1),
          shiny::hr(),
          shiny::h4("Gene Expression"),
          shiny::helpText(
            "Upload a single-column CSV/TSV with gene names as row names ",
            "and expression values (log2 FPKM). Or paste values below."
          ),
          shiny::fileInput("single_file", "Upload expression file",
                           accept = c(".csv", ".tsv", ".txt")),
          shiny::hr(),
          shiny::h4("Or enter key gene values manually"),
          shiny::helpText("9 signature genes (remaining genes set to 0):"),
          lapply(hi_genes, function(g) {
            shiny::numericInput(
              inputId = paste0("gene_", g),
              label   = g,
              value   = 0, step = 0.1
            )
          }),
          shiny::hr(),
          shiny::radioButtons("input_mode", "Input mode:",
                              choices = c("File upload" = "file",
                                          "Manual entry" = "manual"),
                              selected = "file"),
          shiny::actionButton("predict_btn", "Predict Risk",
                              class = "btn-primary",
                              style = "width:100%; font-size:16px;")
        ),
        shiny::mainPanel(
          width = 8,
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::wellPanel(
                shiny::h4("Risk Score"),
                shiny::verbatimTextOutput("risk_score_txt"),
                shiny::h4("Risk Group"),
                shiny::uiOutput("risk_group_ui")
              )
            ),
            shiny::column(
              6,
              shiny::wellPanel(
                shiny::h4("Survival Probability"),
                shiny::tableOutput("surv_table")
              )
            )
          ),
          shiny::hr(),
          shiny::h4("Nomogram"),
          shiny::plotOutput("nomogram_plot", height = "350px"),
          shiny::hr(),
          shiny::h4("Individualized Survival Curve"),
          shiny::plotOutput("surv_curve_plot", height = "400px")
        )
      )
    ),

    # ---- Tab 2: Batch Analysis ----
    shiny::tabPanel(
      "Batch Analysis",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          width = 3,
          shiny::h4("Upload Expression Data"),
          shiny::fileInput("expr_file", "CSV / TSV file",
                           accept = c(".csv", ".tsv", ".txt")),
          shiny::helpText(
            "Format: genes in rows, samples in columns. ",
            "First column = gene names. Values = log2(FPKM)."
          ),
          shiny::hr(),
          shiny::h4("Optional: Age Data"),
          shiny::fileInput("age_file", "CSV with sample & age columns",
                           accept = c(".csv", ".tsv", ".txt")),
          shiny::hr(),
          shiny::actionButton("batch_btn", "Run Classification",
                              class = "btn-primary",
                              style = "width:100%; font-size:16px;")
        ),
        shiny::mainPanel(
          width = 9,
          shiny::h4("Classification Results"),
          DT::DTOutput("batch_table"),
          shiny::hr(),
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::h4("Risk Score Distribution"),
              shiny::plotOutput("score_dist_plot", height = "350px")
            ),
            shiny::column(
              6,
              shiny::h4("Group Summary"),
              shiny::verbatimTextOutput("group_summary")
            )
          ),
          shiny::hr(),
          shiny::h4("Download Results"),
          shiny::downloadButton("download_results", "Download CSV")
        )
      )
    ),

    # ---- Tab 3: About ----
    shiny::tabPanel(
      "About",
      shiny::fluidRow(
        shiny::column(
          8, offset = 2,
          shiny::wellPanel(
            shiny::h3("MMimmuneRisk"),
            shiny::p("A machine-learning-based risk stratification tool for ",
                     "Multiple Myeloma using immune-related gene expression."),
            shiny::h4("Key Signature Genes (PCDI)"),
            shiny::tags$ul(
              lapply(hi_genes, function(g) shiny::tags$li(g))
            ),
            shiny::h4("Method"),
            shiny::p("The classifier uses a two-step pipeline:"),
            shiny::tags$ol(
              shiny::tags$li("Random Survival Forest (RSF) for feature selection"),
              shiny::tags$li("Partial Least Squares Cox Regression (plsRcox) for risk scoring")
            ),
            shiny::p(paste0("Total model genes: ", length(mg))),
            shiny::p(paste0("Optimal risk score cutoff: ", round(cutoff, 4))),
            shiny::h4("Training & Validation"),
            shiny::tags$ul(
              shiny::tags$li("Training: MMRF-CoMMpass"),
              shiny::tags$li("Internal validation: MMRF-test"),
              shiny::tags$li("External validation: GSE2658, JSPH cohort")
            )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {

    # ---- Helper: read single-patient expression ----
    get_single_expr <- shiny::reactive({
      if (input$input_mode == "file" && !is.null(input$single_file)) {
        fpath <- input$single_file$datapath
        fname <- input$single_file$name
        if (grepl("\\.tsv$|\\.txt$", fname)) {
          expr <- read.delim(fpath, row.names = 1, check.names = FALSE)
        } else {
          expr <- read.csv(fpath, row.names = 1, check.names = FALSE)
        }
        return(as.matrix(expr))
      }

      vals <- sapply(hi_genes, function(g) input[[paste0("gene_", g)]])
      full_vals <- rep(0, length(mg))
      names(full_vals) <- mg
      for (g in hi_genes) {
        if (g %in% mg) full_vals[g] <- vals[g]
      }
      matrix(full_vals, ncol = 1, dimnames = list(mg, "Patient"))
    })

    # ---- Single Patient ----
    single_result <- shiny::eventReactive(input$predict_btn, {
      expr_mat <- get_single_expr()
      rs <- predict_mmrisk(expr_mat, genes_in_rows = TRUE)
      cl <- classify_mmrisk(expr_mat, genes_in_rows = TRUE,
                            age = input$age, cutoff = cutoff)
      list(rs = rs, cl = cl)
    })

    output$risk_score_txt <- shiny::renderText({
      req(single_result())
      sprintf("%.4f", single_result()$rs)
    })

    output$risk_group_ui <- shiny::renderUI({
      req(single_result())
      grp <- single_result()$cl$risk_group
      col <- if (grp == "High") "#c873a8" else "#74a3b2"
      shiny::tags$span(
        style = paste0("font-size:24px; font-weight:bold; color:", col),
        paste0(grp, " Risk")
      )
    })

    output$surv_table <- shiny::renderTable({
      req(single_result())
      cl <- single_result()$cl
      if ("surv_1y" %in% names(cl)) {
        data.frame(
          Timepoint = c("1-Year", "2-Year", "3-Year"),
          `Survival Probability` = sprintf("%.1f%%",
            c(cl$surv_1y, cl$surv_2y, cl$surv_3y) * 100),
          check.names = FALSE
        )
      } else {
        data.frame(Note = "Cox model with Age needed for survival estimation.")
      }
    })

    output$nomogram_plot <- shiny::renderPlot({
      req(single_result())
      cox_mod <- tryCatch(get_model("cox_model"), error = function(e) NULL)
      if (is.null(cox_mod)) {
        plot.new(); text(0.5, 0.5, "Cox model not available", cex = 1.5)
        return()
      }
      tryCatch({
        nom <- rms::nomogram(cox_mod,
                             fun = list(
                               function(x) 1 - rms::Survival(cox_mod)(12, x),
                               function(x) 1 - rms::Survival(cox_mod)(24, x),
                               function(x) 1 - rms::Survival(cox_mod)(36, x)
                             ),
                             funlabel = c("1-Year Mortality",
                                          "2-Year Mortality",
                                          "3-Year Mortality"))
        plot(nom)
      }, error = function(e) {
        plot.new(); text(0.5, 0.5, paste("Nomogram error:", e$message), cex = 1)
      })
    })

    output$surv_curve_plot <- shiny::renderPlot({
      req(single_result())
      cox_mod <- tryCatch(get_model("cox_model"), error = function(e) NULL)
      if (is.null(cox_mod)) {
        plot.new(); text(0.5, 0.5, "Cox model not available", cex = 1.5)
        return()
      }
      nd <- data.frame(rs = single_result()$rs, Age = input$age)
      tryCatch({
        sf <- rms::survest(cox_mod, newdata = nd,
                           times = seq(0, 60, by = 1))
        surv_vals <- if (is.matrix(sf$surv)) sf$surv[1, ] else sf$surv
        low_vals  <- if (is.matrix(sf$lower)) sf$lower[1, ] else sf$lower
        up_vals   <- if (is.matrix(sf$upper)) sf$upper[1, ] else sf$upper
        plot_df <- data.frame(time = sf$time, surv = surv_vals,
                              lower = low_vals, upper = up_vals)
        grp_col <- if (single_result()$cl$risk_group == "High") "#c873a8" else "#74a3b2"

        ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = surv)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                               fill = grp_col, alpha = 0.2) +
          ggplot2::geom_line(color = grp_col, linewidth = 1.2) +
          ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                              color = "grey50") +
          ggplot2::labs(x = "Time (months)", y = "Survival Probability",
                        title = paste0("Predicted Survival - ",
                                       single_result()$cl$risk_group, " Risk")) +
          ggplot2::scale_y_continuous(limits = c(0, 1),
                                     labels = scales::label_percent()) +
          ggplot2::theme_minimal(base_size = 14)
      }, error = function(e) {
        plot.new(); text(0.5, 0.5, paste("Error:", e$message), cex = 1)
      })
    })

    # ---- Batch Analysis ----
    batch_result <- shiny::eventReactive(input$batch_btn, {
      req(input$expr_file)
      fpath <- input$expr_file$datapath
      fname <- input$expr_file$name
      if (grepl("\\.tsv$|\\.txt$", fname)) {
        expr <- read.delim(fpath, row.names = 1, check.names = FALSE)
      } else {
        expr <- read.csv(fpath, row.names = 1, check.names = FALSE)
      }

      age_vec <- NULL
      if (!is.null(input$age_file)) {
        age_df <- read.csv(input$age_file$datapath, check.names = FALSE)
        sample_col <- intersect(c("sample", "Sample", "ID", "id"),
                                colnames(age_df))[1]
        age_col <- intersect(c("age", "Age", "AGE"), colnames(age_df))[1]
        if (!is.na(sample_col) && !is.na(age_col)) {
          age_matched <- age_df[[age_col]][match(colnames(expr),
                                                 age_df[[sample_col]])]
          if (!all(is.na(age_matched))) age_vec <- age_matched
        }
      }

      classify_mmrisk(expr, genes_in_rows = TRUE, age = age_vec, cutoff = cutoff)
    })

    output$batch_table <- DT::renderDT({
      req(batch_result())
      DT::datatable(batch_result(), options = list(pageLength = 20),
                    rownames = FALSE) |>
        DT::formatRound(columns = "risk_score", digits = 4)
    })

    output$score_dist_plot <- shiny::renderPlot({
      req(batch_result())
      df <- batch_result()
      ggplot2::ggplot(df, ggplot2::aes(x = risk_score, fill = risk_group)) +
        ggplot2::geom_histogram(bins = 30, alpha = 0.7, color = "white") +
        ggplot2::scale_fill_manual(values = c("High" = "#c873a8",
                                              "Low" = "#74a3b2")) +
        ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed",
                            color = "red", linewidth = 0.8) +
        ggplot2::annotate("text", x = cutoff, y = Inf,
                          label = paste0("Cutoff = ", round(cutoff, 3)),
                          vjust = 2, hjust = -0.1, color = "red") +
        ggplot2::labs(x = "Risk Score", y = "Count", fill = "Risk Group") +
        ggplot2::theme_minimal(base_size = 14)
    })

    output$group_summary <- shiny::renderPrint({
      req(batch_result())
      df <- batch_result()
      cat("Total patients:", nrow(df), "\n")
      cat("High risk:", sum(df$risk_group == "High"),
          sprintf("(%.1f%%)", mean(df$risk_group == "High") * 100), "\n")
      cat("Low risk:", sum(df$risk_group == "Low"),
          sprintf("(%.1f%%)", mean(df$risk_group == "Low") * 100), "\n")
      cat("\nRisk Score Summary:\n")
      print(summary(df$risk_score))
    })

    output$download_results <- shiny::downloadHandler(
      filename = function() {
        paste0("MMimmuneRisk_results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(batch_result(), file, row.names = FALSE)
      }
    )
  }

  shiny::shinyApp(ui = ui, server = server)
}
