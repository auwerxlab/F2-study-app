# Protein vs RNA-seq comparison module UI (the "Data comparison" tab).
# Gene list is read once from the precomputed bundle (prot_rna_path in app.R).
.compare_ui_cache <- new.env(parent = emptyenv())

compare_ui <- function(id) {
  ns <- NS(id)

  # The UI only needs the cluster levels (to decide whether to offer the
  # "Liver cluster" colour option). Read the few-KB meta sidecar written by
  # Scripts/preprocess_prot_rna.R rather than decompressing the 82 MB bundle;
  # fall back to the full bundle if the sidecar is absent (older data).
  if (is.null(.compare_ui_cache$loaded)) {
    src <- if (exists("prot_rna_meta_path") && file.exists(prot_rna_meta_path)) prot_rna_meta_path
           else if (exists("prot_rna_path") && file.exists(prot_rna_path)) prot_rna_path
           else NULL
    obj <- if (is.null(src)) NULL else tryCatch(readRDS(src), error = function(e) NULL)
    .compare_ui_cache$cluster_levels <- if (is.null(obj)) character(0) else obj$cluster_levels
    .compare_ui_cache$loaded <- TRUE
  }
  has_cluster <- length(.compare_ui_cache$cluster_levels) > 0

  color_choices <- c("Strain" = "Strain")
  if (has_cluster) color_choices <- c(color_choices, "Liver cluster" = "liver_cluster")
  method_choices <- c("Pearson" = "pearson", "Spearman" = "spearman")

  # Plain-language "Pearson vs Spearman" help, shown in the hover tooltip on the
  # Correlation "?" icon (fresh instance per call). The tooltip is widened in app.R
  # (`.tooltip-inner { max-width: ... }`) so this reads across the screen, not in a
  # narrow column.
  corr_tip <- function() tags$div(
    style = "text-align: left;",
    tags$p(class = "mb-1",
           "Both run from -1 to +1 (+1 = rise together, 0 = no link, -1 = opposite)."),
    tags$p(class = "mb-1",
           tags$b("Pearson"), ": how close the points are to a straight line (sensitive to outliers)."),
    tags$p(class = "mb-0",
           tags$b("Spearman"), ": whether they agree in rank order - robust to outliers and curved trends.")
  )

  tagList(
    method_panel(
      tags$ul(
        tags$li(tags$b("Inputs:"), " matched liver transcriptomes (log2-CPM) and proteomes (log2 abundance) from the same mice."),
        tags$li(tags$b("Correlation:"), " for each gene, a Pearson correlation between transcript and protein level was computed within each liver cluster. The tabs also offer Spearman and per-sample correlations for exploration."),
        tags$li(tags$b("Differential correlation:"), " to compare clusters, Pearson r values were Fisher r-to-z transformed and compared with a Z-test between each test cluster (Clusters 2-3) and the reference (Cluster 1); Cluster 4 was excluded for low sample size. P-values were Benjamini-Hochberg adjusted, and genes with FDR < 0.05 were called differentially correlated.")
      )
    ),

  bslib::navset_card_tab(
    id = ns("tabs"),
    title = tagList(icon("exchange-alt"), " Protein vs RNA"),

    bslib::nav_panel(
      "Per-sample",
      p("RNA-protein correlation within each mouse (across genes). How well does protein abundance track transcript level per sample?"),
      fluidRow(
        column(3, selectInput(ns("smethod"), label_tip("Correlation", corr_tip()),
                              choices = method_choices, selected = "pearson")),
        column(3, selectInput(ns("scolor"), label_tip("Color by", "Color the per-sample bars by strain or liver cluster"),
                              choices = color_choices, selected = "Strain")),
        column(6, div(class = "d-flex align-items-end h-100 pb-1 gap-2",
                      actionButton(ns("refresh_s"), "Refresh", icon = icon("arrows-rotate"), class = "btn btn-outline-secondary"),
                      downloadButton(ns("download_sample"), "Download barplot (PDF)", class = "btn btn-outline-secondary")))
      ),
      div(class = "mt-1", tags$small(class = "text-muted", textOutput(ns("sample_caption"), inline = TRUE))),
      fluidRow(
        column(5, div(class = "mt-2", h6("Distribution"), plotOutput(ns("sample_hist"), height = "320px"))),
        column(7, div(class = "mt-2", h6("Per sample (ranked)"), plotOutput(ns("sample_bar"), height = "320px")))
      )
    ),

    bslib::nav_panel(
      "Per-gene",
      p("RNA-protein correlation for each gene (across mice). Genes with positive correlation have concordant transcript and protein levels."),
      fluidRow(
        column(3, selectInput(ns("gmethod"), label_tip("Correlation", corr_tip()),
                              choices = method_choices, selected = "pearson")),
        column(6, div(class = "d-flex align-items-end h-100 pb-1 gap-2",
                      actionButton(ns("refresh_g"), "Refresh", icon = icon("arrows-rotate"), class = "btn btn-outline-secondary"),
                      downloadButton(ns("download_gene"), "Download histogram (PDF)", class = "btn btn-outline-secondary")))
      ),
      fluidRow(
        column(6, div(class = "mt-2", h6("Distribution (all genes)"), plotOutput(ns("gene_hist"), height = "300px"))),
        column(6, div(class = "mt-2", h6("By liver cluster"), plotOutput(ns("gene_density"), height = "300px")))
      ),
      div(class = "mt-3",
          div(class = "d-flex justify-content-between align-items-center",
              h6("Per-gene correlation table", class = "mb-0"),
              downloadButton(ns("download_gene_table"), "Download CSV", class = "btn btn-outline-secondary btn-sm")),
          DT::dataTableOutput(ns("gene_table")))
    ),

    bslib::nav_panel(
      "Gene scatter",
      p("Protein vs RNA across mice for one gene. Each point is a mouse; the grey line is the ",
        tags$b("linear fit"), " - the straight line that best follows the points. An ",
        tags$b("upward"), " line means mice with more RNA tend to have more protein (they move together); a ",
        tags$b("downward"), " line means more RNA goes with less protein; a roughly ", tags$b("flat"),
        " line means RNA tells you little about that gene's protein level. The steeper the line, the stronger the trend."),
      fluidRow(
        column(4, selectizeInput(ns("gene"), label_tip("Gene", "Choose a gene (by symbol) to plot protein vs RNA"),
                                 choices = NULL, multiple = FALSE,
                                 options = list(placeholder = "Type a gene symbol..."))),
        column(3, selectInput(ns("pcolor"), label_tip("Color by", "Color points by liver cluster or strain"),
                              choices = rev(color_choices),
                              selected = if (has_cluster) "liver_cluster" else "Strain")),
        column(5, div(class = "d-flex align-items-end h-100 pb-1 gap-2",
                      actionButton(ns("refresh_sc"), "Refresh", icon = icon("arrows-rotate"), class = "btn btn-outline-secondary"),
                      downloadButton(ns("download_scatter"), "Download scatter (PDF)", class = "btn btn-outline-secondary")))
      ),
      div(class = "mt-1", tags$small(class = "text-muted", textOutput(ns("scatter_caption"), inline = TRUE))),
      div(class = "mt-2", plotOutput(ns("gene_scatter"), height = "460px"))
    )
  )
  )
}
