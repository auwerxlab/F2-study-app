# Phenotype choices are read once and reused across tab navigations / sessions
# (qtl_data_dir is defined in app.R).
.qtl_ui_cache <- new.env(parent = emptyenv())

qtl_ui <- function(id) {
  ns <- NS(id)

  # Pre-read phenotype names for immediate population (cached). Prefer the tiny
  # sidecar (qtl_pheno_choices.rds, written by Scripts/preprocess_qtl_data.R) so
  # the UI does not decompress a ~4.6 MB mapping RDS just to list phenotype names;
  # fall back to reading male_005_001.rds if the sidecar is absent (older data).
  if (is.null(.qtl_ui_cache$pheno_choices)) {
    meta_f <- file.path(qtl_data_dir, "qtl_pheno_choices.rds")
    full_f <- file.path(qtl_data_dir, "male_005_001.rds")
    .qtl_ui_cache$pheno_choices <- if (file.exists(meta_f)) {
      tryCatch(readRDS(meta_f), error = function(e) character(0))
    } else if (file.exists(full_f)) {
      sort(unique(readRDS(full_f)$phenotype))
    } else character(0)
  }
  pheno_choices <- .qtl_ui_cache$pheno_choices

  # A sensible default for the single-phenotype view (falls back to the first).
  prof_default <- if ("Liver_g_24" %in% pheno_choices) "Liver_g_24"
                  else if (length(pheno_choices)) unname(pheno_choices[1]) else NULL

  tagList(
    method_panel(
      tags$ul(
        tags$li(tags$b("Genotypes from RNA-seq:"), " variants were called from the RNA-seq alignments with the GATK RNA-seq workflow (HaplotypeCaller) and hard-filtered."),
        tags$li(tags$b("Genotype map:"), " a custom window-based pipeline selected founder-informative SNPs across the four founders (B6, 129, CAST, PWK), assigned each F2 interval to one of four genotype classes (B6-129, B6-PWK, CAST-129, CAST-PWK), and merged concordant windows into genotype blocks."),
        tags$li(tags$b("QTL mapping:"), " performed in R with R/qtl (v1.74) on the RNA-seq-derived autosomal maps, scanning males and females separately. Genome-wide LOD significance thresholds came from 5,000 permutations (alpha = 0.05 and 0.01)."),
        tags$li(tags$b("Confidence intervals:"), " each lead peak carries a 95% Bayesian confidence interval; for the body-composition and liver-disease trait groups, intervals shared by at least two traits were merged into core CI regions used for candidate-gene scoring.")
      )
    ),

  bslib::navset_card_tab(
    id = ns("qtl_tabs"),
    title = tagList(icon("project-diagram"), " QTL"),

    # ---- Tab 1: genome-wide Manhattan plot (all phenotypes) ----
    bslib::nav_panel(
      "By Phenotype",
      p("Interactive Manhattan plot of QTL mapping results across 19 chromosomes."),

      bslib::accordion(
        id = ns("qtl_params"),
        bslib::accordion_panel(
          title = "Plot Parameters",
          value = "params",
          open = TRUE,

          # Two aligned rows: four col-3 cells each, so every row-2 control sits directly
          # under the row-1 field above it.
          fluidRow(
            column(3, selectInput(ns("sex_select"), label_tip("Sex", "Select which sex to display: male, female, or both overlaid"),
                                  choices = c("Male" = "male", "Female" = "female", "Both (overlay)" = "both"),
                                  selected = "male")),
            column(3, selectizeInput(ns("pheno_select"), label_tip("Phenotypes", "Filter to specific phenotypes; leave empty to show all"),
                                     choices = pheno_choices, selected = NULL, multiple = TRUE, width = "100%",
                                     options = list(placeholder = "All phenotypes (leave empty for all)", plugins = list("remove_button")))),
            column(3, numericInput(ns("plot_height"), label_tip("Height", "Adjust the vertical size of the Manhattan plot in pixels"),
                                   value = 500, min = 300, max = 1200, step = 50)),
            column(3, div(class = "d-flex align-items-end justify-content-end h-100 pb-1",
                          downloadButton(ns("download_plot"), "PDF", class = "btn btn-primary")))
          ),
          fluidRow(
            column(3, sliderInput(ns("pval_max"), label_tip("Max P-value", "Only show QTLs with P-value at or below this threshold"),
                                  min = 0.01, max = 0.05, value = 0.025, step = 0.001)),
            column(3, sliderInput(ns("lod_min"), label_tip("Min LOD Score", "Lower Y-axis cutoff: hide points below this LOD score (0 shows all points)"),
                                  min = 0, max = 10, value = 0, step = 0.5)),
            column(3, div(class = "d-flex align-items-end h-100",
                          checkboxInput(ns("show_all_points"), label_tip("Show all points", "Include non-significant points (P > 0.05) in grey"), value = FALSE)))
          )
        )
      ),

      div(class = "mt-3 text-center",
          actionButton(ns("refresh"), "Refresh plot", icon = icon("arrows-rotate"),
                       class = "btn btn-outline-secondary btn-sm px-4")),
      div(class = "mt-3",
          uiOutput(ns("plot_container"))
      ),

      div(class = "mt-4",
          div(class = "d-flex justify-content-between align-items-center",
              h5("QTL Peaks", class = "mb-0"),
              downloadButton(ns("download_peaks"), "Download CSV", class = "btn btn-outline-secondary btn-sm")),
          checkboxInput(ns("show_ns"), label_tip("Include non-significant", "Show all QTL peaks in the table, not just significant ones"), value = FALSE),
          DT::dataTableOutput(ns("peak_table"))
      )
    ),

    # ---- Tab 2: one phenotype, one panel per chromosome (4 per row) ----
    bslib::nav_panel(
      "By Chromosome",
      p("LOD profile of a single phenotype across all 19 chromosomes - one panel per chromosome. The dashed line marks the significance LOD; use the slider to move it, and optionally hide markers below it."),

      bslib::accordion(
        id = ns("prof_params"),
        bslib::accordion_panel(
          title = "Plot Parameters",
          value = "params",
          open = TRUE,

          # Two aligned rows: four col-3 cells each, so the row-2 controls line up
          # under the row-1 fields (download under Height).
          fluidRow(
            column(3, selectInput(ns("prof_pheno"), label_tip("Phenotype", "Phenotype whose LOD profile to display across all chromosomes"),
                                  choices = pheno_choices, selected = prof_default)),
            column(3, selectInput(ns("prof_sex"), label_tip("Sex", "Show male, female, or overlay both sexes"),
                                  choices = c("Male" = "male", "Female" = "female", "Both (overlay)" = "both"),
                                  selected = "male")),
            column(3, sliderInput(ns("prof_lod_min"), label_tip("Significance LOD", "Vertical position of the horizontal dashed significance line. Defaults to this phenotype's own LOD threshold; move it to explore other cutoffs."),
                                  min = 0, max = 10, value = 4, step = 0.1)),
            column(3, numericInput(ns("prof_height"), label_tip("Height", "Adjust the vertical size of the per-chromosome grid"),
                                   value = 800, min = 400, max = 2000, step = 50))
          ),
          fluidRow(
            column(9, div(class = "d-flex align-items-end h-100",
                          checkboxInput(ns("prof_hide_below"), label_tip("Hide points below significance LOD", "Drop markers below the significance line. The x and y axes stay anchored at 0 so the plot does not rescale."), value = FALSE))),
            column(3, div(class = "d-flex align-items-end justify-content-end h-100 pb-1",
                          downloadButton(ns("download_profile"), "Download PDF", class = "btn btn-primary")))
          )
        )
      ),

      div(class = "mt-2 text-center", actionButton(ns("refresh_prof"), "Refresh plot", icon = icon("arrows-rotate"),
                                       class = "btn btn-outline-secondary btn-sm px-4")),
      div(class = "mt-1", tags$small(class = "text-muted", textOutput(ns("profile_caption"), inline = TRUE))),
      div(class = "mt-2",
          uiOutput(ns("profile_container"))
      ),

      div(class = "mt-4",
          div(class = "d-flex justify-content-between align-items-center",
              h5("Chromosome peaks", class = "mb-0"),
              downloadButton(ns("download_profile_table"), "Download CSV", class = "btn btn-outline-secondary btn-sm")),
          checkboxInput(ns("prof_show_ns"), label_tip("Include non-significant", "Also list non-significant (NS) peaks, not just significant ones. Click Refresh to apply."), value = FALSE),
          DT::dataTableOutput(ns("profile_table"))
      )
    )
  )
  )
}
