# Proteomics PCA module UI.
# Control choices (sample sets, PC axes, colour-by options) come from the tiny
# meta sidecar (proteomics_pca_meta.rds) written by preprocess_proteomics_pca.R,
# read once and cached. Falls back to the full bundle, then to a minimal UI.
.proteomics_ui_cache <- new.env(parent = emptyenv())

proteomics_pca_ui <- function(id) {
  ns <- NS(id)

  if (is.null(.proteomics_ui_cache$loaded)) {
    src <- if (exists("proteomics_pca_meta_path") && file.exists(proteomics_pca_meta_path)) proteomics_pca_meta_path
           else if (exists("proteomics_pca_path") && file.exists(proteomics_pca_path)) proteomics_pca_path
           else NULL
    .proteomics_ui_cache$pca <- if (is.null(src)) NULL else tryCatch(readRDS(src), error = function(e) NULL)
    .proteomics_ui_cache$loaded <- TRUE
  }
  pca <- .proteomics_ui_cache$pca

  # Defaults so the panel still renders without data present
  sample_sets   <- c("All samples" = "all")
  dim_choices   <- c("PC1" = "Dim.1", "PC2" = "Dim.2")
  color_choices <- list("Sample attributes" = c("Strain" = "Strain", "Sex" = "Sex", "Generation" = "Generation"))

  if (!is.null(pca)) {
    have <- pca$set_names %||% names(pca$sets)
    sample_sets <- c(
      if ("all"    %in% have) c("All samples"  = "all"),
      if ("male"   %in% have) c("Males only"   = "male"),
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
  disc_vals <- if (!is.null(pca)) unname(pca$discrete_vars) else c("Strain", "Sex", "Generation")
  log_cond  <- paste0("[", paste(sprintf("'%s'", disc_vals), collapse = ","),
                      "].indexOf(input.color_by) < 0")

  tagList(
    p("Principal-component analysis of liver proteomes (log2 protein abundance). Colour the samples by strain, sex, generation, or any measured phenotype (liver weight, steatosis, fibrosis, ...)."),

    method_panel(
      tags$ul(
        tags$li(tags$b("Sample preparation:"), " about 10 mg of cryopulverised liver was extracted in n-butanol/methanol/water; precipitated protein was reduced, alkylated and digested with LysC then trypsin, and the peptides desalted and normalised to 1 mg/mL."),
        tags$li(tags$b("Mass spectrometry:"), " peptides were separated over a 30-min gradient and analysed on an Orbitrap Astral in data-independent acquisition (DIA) mode."),
        tags$li(tags$b("Identification:"), " DIA data were processed in Spectronaut (v18.1) against the reference proteomes of all four F0 founder strains; 9,776 proteins detected in at least 50% of samples were retained."),
        tags$li(tags$b("PCA:"), " principal-component analysis of log2 protein abundance across the retained proteins.")
      )
    ),

    bslib::accordion(
      id = ns("proteomics_params"),
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
          column(3, selectInput(ns("color_by"), label_tip("Color by", "A sample attribute (strain / sex / generation) or a phenotype such as liver weight or steatosis"),
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
