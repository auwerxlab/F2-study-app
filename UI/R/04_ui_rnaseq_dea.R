# RNA-seq differential expression module UI.
# Contrast choices come from the precomputed RDS (rnaseq_dea_path in app.R),
# read once and cached. Returns a tagList for the RNA-seq tabset.
.rnaseq_dea_ui_cache <- new.env(parent = emptyenv())

rnaseq_dea_ui <- function(id) {
  ns <- NS(id)

  # Prefer the tiny meta sidecars (rnaseq_dea_meta.rds / rnaseq_ora_meta.rds) so
  # the UI need not decompress the full DEA result tables / ORA bundle just to
  # build dropdowns; fall back to the full bundles if a sidecar is absent.
  if (is.null(.rnaseq_dea_ui_cache$loaded)) {
    dsrc <- if (exists("rnaseq_dea_meta_path") && file.exists(rnaseq_dea_meta_path)) rnaseq_dea_meta_path
            else if (exists("rnaseq_dea_path") && file.exists(rnaseq_dea_path)) rnaseq_dea_path
            else NULL
    osrc <- if (exists("rnaseq_ora_meta_path") && file.exists(rnaseq_ora_meta_path)) rnaseq_ora_meta_path
            else if (exists("rnaseq_ora_path") && file.exists(rnaseq_ora_path)) rnaseq_ora_path
            else NULL
    .rnaseq_dea_ui_cache$data <- if (is.null(dsrc)) NULL else tryCatch(readRDS(dsrc), error = function(e) NULL)
    .rnaseq_dea_ui_cache$ora  <- if (is.null(osrc)) NULL else tryCatch(readRDS(osrc), error = function(e) NULL)
    .rnaseq_dea_ui_cache$loaded <- TRUE
  }
  d   <- .rnaseq_dea_ui_cache$data
  ora <- .rnaseq_dea_ui_cache$ora

  contrast_choices <- if (!is.null(d)) d$contrast_choices else c("Sex: male vs female" = "sex")
  ora_coef_choices <- if (!is.null(ora)) ora$coef_choices
                      else c("Fibrosis (%) [adj.]" = "fibrosis", "Plasma TIMP-1 [adj.]" = "timp1")

  tagList(
    p("Differential expression (limma-voom). Pick a contrast to see its volcano plot and full results table. Significant genes pass the FDR and |logFC| thresholds below."),

    method_panel(
      h6("Differential expression (limma-voom)"),
      p("Analyses used F2 mice with complete measurements for the traits of interest (fibrosis %, plasma TIMP-1, relative liver weight, steatosis %, and body fat % at 22 weeks). Raw counts were TMM-normalised and transformed with voom to log2-CPM with precision weights, gene-wise linear models were fitted with lmFit(), and standard errors moderated by empirical Bayes (eBayes). Genes with a Benjamini-Hochberg adjusted P < 0.05 (and, where set, an absolute log2 fold-change above the threshold) were called significant."),
      h6("Trait-association models"),
      p("For each phenotype or covariate, continuous variables were centred and scaled and sex modelled as categorical; the same edgeR / limma-voom pipeline yielded the phenotype-associated gene sets used across the app."),
      h6("GO enrichment"),
      tags$ul(
        tags$li(tags$b("GSEA:"), " all genes are ranked by the direction and significance of their change (signed -log10 adjusted P) and tested for enriched biological-process terms."),
        tags$li(tags$b("Set-specific ORA:"), " significant genes from the multivariate liver-disease model are split into those specific to one coefficient vs shared between two, then each set is tested for GO over-representation.")
      )
    ),

    bslib::accordion(
      id = ns("dea_params"),
      bslib::accordion_panel(
        title = "Volcano plot",
        value = "volcano",
        open = TRUE,

        # Two aligned rows: four col-3 cells in row 1, and row 2 lines its controls up
        # under them (Height under Contrast, Show non-significant under FDR, downloads
        # grouped on the right).
        fluidRow(
          column(3, selectInput(ns("contrast"), label_tip("Contrast", "Which comparison / model term to display"),
                                choices = contrast_choices,
                                selected = if (length(contrast_choices)) unname(contrast_choices[1]) else NULL)),
          column(3, sliderInput(ns("sig_fdr"), label_tip("FDR threshold", "Adjusted p-value (BH) cutoff for calling a gene significant"),
                                min = 0.01, max = 0.10, value = 0.01, step = 0.01)),
          column(3, sliderInput(ns("lfc"), label_tip("Min |logFC|", "Minimum absolute log2 fold-change for significance"),
                                min = 0, max = 2, value = 1, step = 0.1)),
          column(3, numericInput(ns("top_n"), label_tip("Label top genes", "Number of top significant genes to label on the volcano"),
                                 value = 15, min = 0, max = 50, step = 1))
        ),
        fluidRow(
          column(3, numericInput(ns("plot_height"), label_tip("Height", "Vertical size of the volcano plot in pixels"),
                                 value = 520, min = 300, max = 1200, step = 50)),
          column(3, div(class = "d-flex align-items-end h-100 pb-2",
                        checkboxInput(ns("show_ns"), label_tip("Show non-significant", "Show points below the FDR threshold (grey). Off (default) shows only genes significant by FDR."), value = FALSE))),
          column(6, div(class = "d-flex align-items-end justify-content-end h-100 pb-1 gap-2",
                        downloadButton(ns("download_plot"), "Volcano (PDF)", class = "btn btn-primary")))
        ),

        div(class = "mt-2 text-center", actionButton(ns("refresh"), "Refresh plots", icon = icon("arrows-rotate"),
                                         class = "btn btn-outline-secondary btn-sm px-4")),
        div(class = "mt-2", tags$small(class = "text-muted", textOutput(ns("dea_caption"), inline = TRUE))),
        div(class = "mt-2", uiOutput(ns("plot_container")))
      )
    ),

    div(class = "mt-3",
        bslib::accordion(
          id = ns("dea_overview"),
          open = FALSE,
          bslib::accordion_panel(
            title = "DEG counts across all contrasts",
            value = "overview",
            plotOutput(ns("deg_bar"), height = "300px")
          )
        )
    ),

    div(class = "mt-3",
        bslib::accordion(
          id = ns("gsea_panel"),
          open = FALSE,
          bslib::accordion_panel(
            title = "GO enrichment (GSEA, biological process)",
            value = "gsea",
            p(class = "text-muted",
              "Gene-set enrichment for the contrast above (genes ranked by sign(logFC) × -log10(adj. P)). Bars/points show the top enriched biological-process terms."),
            fluidRow(
              column(3, numericInput(ns("gsea_top_n"), label_tip("Top terms / direction", "How many of the most significant enriched terms to show for each direction (up / down)"),
                                     value = 10, min = 3, max = 30, step = 1)),
              column(3, sliderInput(ns("gsea_padj"), label_tip("Max p.adjust", "Only show GO terms at or below this adjusted p-value"),
                                    min = 0.01, max = 0.25, value = 0.05, step = 0.01)),
              column(6, div(class = "d-flex align-items-end h-100 pb-1 gap-2",
                            actionButton(ns("refresh_gsea"), "Refresh", icon = icon("arrows-rotate"), class = "btn btn-outline-secondary"),
                            downloadButton(ns("download_gsea_plot"), "Dot plot (PDF)", class = "btn btn-outline-secondary"),
                            downloadButton(ns("download_gsea_table"), "Table (CSV)", class = "btn btn-outline-secondary")))
            ),
            plotOutput(ns("gsea_dot"), height = "520px"),
            div(class = "mt-3", DT::dataTableOutput(ns("gsea_table")))
          )
        )
    ),

    div(class = "mt-3",
        bslib::accordion(
          id = ns("ora_panel"),
          open = FALSE,
          bslib::accordion_panel(
            title = "Set-specific GO enrichment (ORA across model coefficients)",
            value = "ora",
            p(class = "text-muted",
              "From the multivariate liver-disease model: split significant genes into those specific to one coefficient vs shared between two, then GO over-representation. Up terms point right, down terms left."),
            fluidRow(
              column(3, selectInput(ns("ora_coef_a"), label_tip("Coefficient A", "First model coefficient to compare"),
                                    choices = ora_coef_choices,
                                    selected = if (length(ora_coef_choices)) unname(ora_coef_choices[1]) else NULL)),
              column(3, selectInput(ns("ora_coef_b"), label_tip("Coefficient B", "Second model coefficient to compare"),
                                    choices = ora_coef_choices,
                                    selected = if (length(ora_coef_choices) > 1) unname(ora_coef_choices[2]) else NULL)),
              column(3, selectInput(ns("ora_segment"), label_tip("Gene set", "Which gene set to enrich: specific to A, specific to B, or shared by both"),
                                    choices = c("Coefficient A - specific" = "a",
                                                "Coefficient B - specific" = "b",
                                                "Shared (both)" = "common"),
                                    selected = "a")),
              column(3, numericInput(ns("ora_top_n"), label_tip("Top terms / direction", "How many top GO terms to show for each direction"),
                                     value = 10, min = 3, max = 30, step = 1))
            ),
            fluidRow(
              column(3, sliderInput(ns("ora_padj"), label_tip("Max p.adjust", "Only show GO terms at or below this adjusted p-value"),
                                    min = 0.01, max = 0.25, value = 0.05, step = 0.01)),
              column(9, div(class = "d-flex align-items-end h-100 pb-1 gap-2",
                            actionButton(ns("refresh_ora"), "Refresh", icon = icon("arrows-rotate"), class = "btn btn-outline-secondary"),
                            downloadButton(ns("download_ora_plot"), "Dot plot (PDF)", class = "btn btn-outline-secondary"),
                            downloadButton(ns("download_ora_table"), "Table (CSV)", class = "btn btn-outline-secondary")))
            ),
            div(class = "mt-1", tags$small(class = "text-muted", textOutput(ns("ora_caption"), inline = TRUE))),
            uiOutput(ns("ora_plot_container")),
            div(class = "mt-3", DT::dataTableOutput(ns("ora_table")))
          )
        )
    ),

    # Full results table for the selected contrast, in its own accordion.
    div(class = "mt-3",
        bslib::accordion(
          id = ns("table_panel"),
          open = FALSE,
          bslib::accordion_panel(
            title = "Results table",
            value = "table",
            div(class = "d-flex justify-content-end mb-2",
                downloadButton(ns("download_table"), "Download CSV", class = "btn btn-outline-secondary btn-sm")),
            DT::dataTableOutput(ns("dea_table"))
          )
        )
    )
  )
}
