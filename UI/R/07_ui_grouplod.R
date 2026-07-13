# Group choices are read once and reused across tab navigations / sessions
# (lod_data_dir is defined in app.R).
.grouplod_ui_cache <- new.env(parent = emptyenv())

grouplod_ui <- function(id) {
  ns <- NS(id)

  # Pre-read group names from cached RDS (fast) or fallback to TSV
  if (is.null(.grouplod_ui_cache$group_choices)) {
    rds_file <- file.path(lod_data_dir, "rds", "group_list.rds")
    .grouplod_ui_cache$group_choices <- if (file.exists(rds_file)) {
      readRDS(rds_file)
    } else {
      map_file <- file.path(lod_data_dir, "fix_interst_phenotype_GB")
      if (file.exists(map_file)) {
        mp <- read.table(map_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        sort(unique(mp$Group))
      } else character(0)
    }
  }
  group_choices <- .grouplod_ui_cache$group_choices

  bslib::card(
    header = tagList(icon("chart-bar"), " Group LOD Plot"),
    p("Per-group, per-chromosome LOD visualization with interval bars and LOD lines."),

    method_panel(
      tags$ul(
        tags$li(tags$b("QTL scans:"), " LOD profiles come from R/qtl (v1.74) mapping on the RNA-seq-derived genotype maps, with males and females scanned separately and genome-wide thresholds from 5,000 permutations."),
        tags$li(tags$b("Phenotype groups:"), " related traits are grouped (body-composition traits and liver-disease traits) so that loci shared across a group can be found; other phenotypes are mapped individually."),
        tags$li(tags$b("Shared regions:"), " for each lead locus the app overlaps the 95% Bayesian confidence intervals of the group's phenotypes. The shared-region bar marks where the chosen number of phenotypes overlap (union, all, or at least K) - the same logic used to define the study's core CI regions.")
      )
    ),

    bslib::accordion(
      id = ns("grouplod_params"),
      bslib::accordion_panel(
        title = "Plot Parameters",
        value = "params",
        open = TRUE,

        # Two aligned rows: four col-3 cells each, so every row-2 field sits directly
        # under the row-1 field above it.
        fluidRow(
          column(3, selectInput(ns("sex_select"), label_tip("Sex", "Select male or female QTL data"),
                                choices = c("Male" = "male", "Female" = "female", "Both (overlay)" = "both"), selected = "male")),
          column(3, selectInput(ns("group_select"), label_tip("Phenotype Group", "Group of related phenotypes to visualize together"),
                                choices = group_choices, selected = if (length(group_choices) > 0) group_choices[1] else NULL)),
          column(3, selectInput(ns("chr_select"), label_tip("Chromosome", "Choose a chromosome or auto-select the first with significant QTLs"),
                                choices = c("Auto (significant)" = "auto"), selected = "auto")),
          column(3, selectInput(ns("share_mode"), label_tip("Shared Region Mode", "How to define shared regions: union of any overlap, intersection of all phenotypes, or at least K phenotypes overlapping"),
                                choices = c("Union (any overlap)" = "0", "All phenotypes" = "all", "At least K" = "k"), selected = "k"))
        ),
        fluidRow(
          column(3, conditionalPanel(
                      condition = paste0("input['", ns("share_mode"), "'] == 'k'"),
                      numericInput(ns("share_k"), label_tip("Min. K Phenotypes", "Minimum number of phenotypes that must overlap to define a shared region"),
                                   value = 2, min = 1, max = 20, step = 1),
                      tags$small(class = "text-muted", textOutput(ns("share_k_range"), inline = TRUE)))),
          column(3, numericInput(ns("flank_mb"), label_tip("Flank (Mb)", "Extend the viewing region beyond the significant markers by this many megabases"),
                                 value = 10, min = 0, max = 50, step = 1)),
          column(3, numericInput(ns("plot_height"), label_tip("Height", "Adjust the vertical size of the combined plot in pixels"),
                                 value = 800, min = 400, max = 1500, step = 50)),
          column(3, div(class = "d-flex align-items-end justify-content-end h-100 pb-1",
                        downloadButton(ns("download_plot"), "Download PDF", class = "btn btn-primary")))
        )
      )
    ),

    div(class = "mt-3 text-center",
        actionButton(ns("refresh"), "Refresh plot", icon = icon("arrows-rotate"),
                     class = "btn btn-outline-secondary btn-sm px-4")),
    div(class = "mt-3",
        uiOutput(ns("plot_container"))
    )
  )
}
