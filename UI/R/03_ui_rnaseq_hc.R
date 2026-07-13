# RNA-seq hierarchical clustering module UI.
# Choices come from the precomputed RDS (rnaseq_hc_path defined in app.R),
# read once and cached. Returns a tagList for use inside the RNA-seq tabset.
.rnaseq_hc_ui_cache <- new.env(parent = emptyenv())

rnaseq_hc_ui <- function(id) {
  ns <- NS(id)

  # Prefer the tiny meta sidecar (rnaseq_hc_meta.rds) so the UI need not
  # decompress the full clustering bundle just to build dropdowns; fall back to
  # the full bundle if the sidecar is absent (older data).
  if (is.null(.rnaseq_hc_ui_cache$loaded)) {
    src <- if (exists("rnaseq_hc_meta_path") && file.exists(rnaseq_hc_meta_path)) rnaseq_hc_meta_path
           else if (exists("rnaseq_hc_path") && file.exists(rnaseq_hc_path)) rnaseq_hc_path
           else NULL
    .rnaseq_hc_ui_cache$data <- if (is.null(src)) NULL else tryCatch(readRDS(src), error = function(e) NULL)
    .rnaseq_hc_ui_cache$loaded <- TRUE
  }
  d <- .rnaseq_hc_ui_cache$data

  sample_sets  <- c("All samples" = "all")
  pheno_vars   <- c("Plasma TIMP-1" = "TIMP1_24")
  default_bars <- "TIMP1_24"
  reorder_default <- ""
  strain_levels <- c("C57BL/6J", "129S1/SvImJ", "CAST/EiJ", "PWK/PhJ",
                     "B6CASTF1", "129SPWKF1", "B6CAST-129SPWK-F2")

  if (!is.null(d)) {
    have <- d$set_names %||% names(d$sets)
    sample_sets <- c(
      if ("all"    %in% have) c("All samples"  = "all"),
      if ("male"   %in% have) c("Males only"   = "male"),
      if ("female" %in% have) c("Females only" = "female")
    )
    pheno_vars   <- d$pheno_vars
    default_bars <- d$default_bars
    if (!is.null(d$strain_levels)) strain_levels <- d$strain_levels
    reorder_default <- if ("Fibrosis_perc_24" %in% d$pheno_vars) "Fibrosis_perc_24"
                       else if (length(d$pheno_vars)) unname(d$pheno_vars[1]) else ""
  }

  tagList(
    p("Hierarchical clustering of liver transcriptomes (Euclidean distance, complete linkage). Branches are coloured by strain; the bars below show phenotypes (z-scored), sex and strain."),

    method_panel(
      tags$ul(
        tags$li(tags$b("Input:"), " the same filtered liver RNA-seq matrix as the PCA tab (STAR to GRCm38, RSEM counts, genes with CPM > 1 in at least 50% of samples)."),
        tags$li(tags$b("Clustering:"), " unsupervised hierarchical clustering of whole-liver transcriptomes (here Euclidean distance, complete linkage). In the study, established MASH signature genes clustered together with clinical measures of liver-disease severity, indicating that the most affected mice have a molecular phenotype resembling human MASH."),
        tags$li(tags$b("Phenotype bars:"), " the bars below the dendrogram show z-scored phenotypes (skewed measures such as TIMP-1, ALT, steatosis and fibrosis are log-transformed first), alongside sex and strain.")
      )
    ),

    bslib::accordion(
      id = ns("hc_params"),
      bslib::accordion_panel(
        title = "Plot Parameters",
        value = "params",
        open = TRUE,

        # Two aligned rows: four col-3 cells each, so every row-2 control sits directly
        # under the row-1 field above it. Height shares its cell with the PDF button.
        fluidRow(
          column(3, selectInput(ns("sample_set"), label_tip("Samples", "Cluster all samples, or within one sex only"),
                                choices = sample_sets, selected = "all")),
          column(3, selectInput(ns("reorder_by"), label_tip("Order leaves by", "Rotate branches so samples are ordered by this phenotype (clustering itself is unchanged)"),
                                choices = c("None (cluster order)" = "", pheno_vars), selected = reorder_default)),
          column(3, selectizeInput(ns("strain_filter"), label_tip("Strains", "Show only these strains; the subset is re-clustered. Leave empty for all."),
                                   choices = strain_levels, selected = NULL, multiple = TRUE, width = "100%",
                                   options = list(placeholder = "All strains", plugins = list("remove_button")))),
          column(3, selectizeInput(ns("bar_vars"), label_tip("Phenotype rows", "Phenotypes shown as z-scored heatmap rows under the dendrogram"),
                                   choices = pheno_vars, selected = default_bars, multiple = TRUE, width = "100%",
                                   options = list(plugins = list("remove_button"))))
        ),
        fluidRow(
          column(3, div(class = "d-flex align-items-end h-100",
                        checkboxInput(ns("show_sex"), label_tip("Sex bar", "Show a sex colour bar"), value = TRUE))),
          column(3, div(class = "d-flex align-items-end h-100",
                        checkboxInput(ns("show_strain"), label_tip("Strain bar", "Show a strain colour bar"), value = TRUE))),
          column(3, div(class = "d-flex align-items-end h-100",
                        checkboxInput(ns("log_bars"), label_tip("Log-transform bars", "Apply log1p before z-scoring (the paper logs skewed measures like TIMP-1, ALT, steatosis, fibrosis)"), value = TRUE))),
          column(3, fluidRow(
                      column(6, numericInput(ns("plot_height"), label_tip("Height", "Vertical size of the plot in pixels"),
                                             value = 600, min = 300, max = 1400, step = 50)),
                      column(6, numericInput(ns("tree_height"), label_tip("Tree ht", "Height of the dendrogram in points; increase to enlarge the tree relative to the heatmap body"),
                                             value = 150, min = 40, max = 500, step = 20))))
        )
      )
    ),

    div(class = "mt-2 d-flex justify-content-center align-items-center gap-2",
        actionButton(ns("refresh"), "Refresh plot", icon = icon("arrows-rotate"),
                     class = "btn btn-outline-secondary btn-sm px-4"),
        downloadButton(ns("download_plot"), "PDF", class = "btn btn-primary btn-sm")),
    div(class = "mt-2", tags$small(class = "text-muted", textOutput(ns("hc_caption"), inline = TRUE))),
    div(class = "mt-2", uiOutput(ns("plot_container")))
  )
}
