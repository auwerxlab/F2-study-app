# RNA-seq PCA module UI.
# Control choices (sample sets, PC axes, colour-by options) are read once from
# the precomputed RDS (rnaseq_pca_path defined in app.R) and reused across
# navigations / sessions. Falls back to a minimal UI if the data is absent.
.rnaseq_ui_cache <- new.env(parent = emptyenv())

rnaseq_pca_ui <- function(id) {
  ns <- NS(id)

  # Prefer the tiny meta sidecar (rnaseq_pca_meta.rds) so the UI need not
  # decompress the full PCA bundle just to build dropdowns; fall back to the full
  # bundle if the sidecar is absent (older data). Both expose dims/discrete_vars/
  # pheno_vars; only the sample-set names differ (set_names vs names(sets)).
  if (is.null(.rnaseq_ui_cache$loaded)) {
    src <- if (exists("rnaseq_pca_meta_path") && file.exists(rnaseq_pca_meta_path)) rnaseq_pca_meta_path
           else if (exists("rnaseq_pca_path") && file.exists(rnaseq_pca_path)) rnaseq_pca_path
           else NULL
    .rnaseq_ui_cache$pca <- if (is.null(src)) NULL else tryCatch(readRDS(src), error = function(e) NULL)
    .rnaseq_ui_cache$loaded <- TRUE
  }
  pca <- .rnaseq_ui_cache$pca

  # Defaults so the panel still renders without data present
  sample_sets   <- c("All samples" = "all")
  dim_choices   <- c("PC1" = "Dim.1", "PC2" = "Dim.2")
  color_choices <- list("Sample attributes" = c("Strain" = "Strain", "Sex" = "Sex"))

  if (!is.null(pca)) {
    have <- pca$set_names %||% names(pca$sets)
    sample_sets <- c(
      if ("all"    %in% have) c("All samples" = "all"),
      if ("male"   %in% have) c("Males only"  = "male"),
      if ("female" %in% have) c("Females only" = "female")
    )
    dim_choices   <- setNames(pca$dims, sub("^Dim\\.", "PC", pca$dims))
    color_choices <- list(
      "Sample attributes" = pca$discrete_vars,
      "Phenotypes"        = pca$pheno_vars
    )
  }

  strain_levels <- c("C57BL/6J", "129S1/SvImJ", "CAST/EiJ", "PWK/PhJ",
                     "B6CASTF1", "129SPWKF1", "B6CAST-129SPWK-F2")
  # Log-scale colour only makes sense for a continuous phenotype, so reveal that
  # checkbox only when Color-by is NOT one of the sample attributes.
  disc_vals <- if (!is.null(pca)) unname(pca$discrete_vars) else c("Strain", "Sex")
  log_cond  <- paste0("[", paste(sprintf("'%s'", disc_vals), collapse = ","),
                      "].indexOf(input.color_by) < 0")

  tagList(
    p("Principal-component analysis of liver transcriptomes. Colour the samples by strain, sex, or any measured phenotype (liver weight, steatosis, fibrosis, ...)."),

    method_panel(
      tags$ul(
        tags$li(tags$b("RNA extraction:"), " livers were pulverised in liquid nitrogen, homogenised in TRIzol and purified with Direct-zol kits; purity (NanoDrop) and integrity (FragmentAnalyzer) were checked."),
        tags$li(tags$b("Sequencing and quantification:"), " strand-specific mRNA-seq on DNBSEQ (at least 50M PE100 reads per sample); reads were mapped to the C57BL/6J reference (GRCm38) with STAR, and gene-level counts and CPM computed with RSEM. Genes were kept if CPM > 1 in at least 50% of samples."),
        tags$li(tags$b("PCA:"), " principal-component analysis of the filtered transcriptomes. In the study, PC1 separated mice by liver-disease severity and PC2 by sex, so disease severity is the dominant axis of transcriptomic variation.")
      )
    ),

    bslib::accordion(
      id = ns("rnaseq_params"),
      bslib::accordion_panel(
        title = "Plot Parameters",
        value = "params",
        open = TRUE,

        # Two aligned rows: four col-3 cells each, so every row-2 field sits directly
        # under the row-1 field above it. Log-color (conditional) tucks under Color by;
        # Height shares its cell with the PDF button.
        fluidRow(
          column(3, selectInput(ns("sample_set"), label_tip("Samples", "Use all samples, or a PCA computed within one sex only"),
                                choices = sample_sets, selected = "all")),
          column(3, selectInput(ns("color_by"), label_tip("Color by", "A sample attribute (strain / sex) or a phenotype such as liver weight or steatosis"),
                                choices = color_choices, selected = "Strain"),
                    conditionalPanel(
                      condition = log_cond, ns = ns,
                      checkboxInput(ns("log_color"), label_tip("Log color", "Apply log1p to the phenotype colour scale (useful for skewed measures like ALT, steatosis, fibrosis)"), value = TRUE))),
          column(3, selectizeInput(ns("strain_filter"), label_tip("Strains", "Show only these strains in the PCA; leave empty for all."),
                                   choices = strain_levels, selected = NULL, multiple = TRUE, width = "100%",
                                   options = list(placeholder = "All strains", plugins = list("remove_button")))),
          column(3, selectInput(ns("shape_by"), label_tip("Shape by", "Optionally use point shape to distinguish sex"),
                                choices = c("None" = "none", "Sex" = "Sex"), selected = "Sex"))
        ),
        fluidRow(
          column(3, selectInput(ns("x_dim"), label_tip("X axis", "Principal component on the horizontal axis"),
                                choices = dim_choices, selected = "Dim.1")),
          column(3, selectInput(ns("y_dim"), label_tip("Y axis", "Principal component on the vertical axis"),
                                choices = dim_choices, selected = "Dim.2")),
          column(3, numericInput(ns("point_size"), label_tip("Point size", "Size of the sample points"),
                                 value = 2, min = 0.5, max = 6, step = 0.5)),
          column(3, fluidRow(
                      column(7, numericInput(ns("plot_height"), label_tip("Height", "Vertical size of the plot in pixels"),
                                             value = 550, min = 300, max = 1200, step = 50)),
                      column(5, div(class = "d-flex align-items-end justify-content-end h-100 pb-1",
                                    downloadButton(ns("download_plot"), "PDF", class = "btn btn-primary")))))
        )
      )
    ),

    div(class = "mt-2 text-center", actionButton(ns("refresh"), "Refresh plot", icon = icon("arrows-rotate"),
                                     class = "btn btn-outline-secondary btn-sm px-4")),
    div(class = "mt-2", tags$small(class = "text-muted", textOutput(ns("pca_caption"), inline = TRUE))),
    div(class = "mt-2", uiOutput(ns("plot_container")))
  )
}
