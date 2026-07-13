phenotype_ui <- function(id) {
  ns <- NS(id)

  bslib::card(
    header = tagList(icon("dna"), " Phenotype Analysis"),

    # --- Collapsible Methods panel ---
    method_panel(
      h6("Study design"),
      p("A four-way intercross of C57BL/6J, 129S1/SvImJ, CAST/EiJ and PWK/PhJ - strains chosen to span metabolic-disease-resistant to -susceptible backgrounds. Mice were fed a Western diet at thermoneutrality to drive metabolic disease and were phenotyped longitudinally from about 6 to 24 weeks of age."),
      h6("Body composition and energy metabolism"),
      tags$ul(
        tags$li(tags$b("Body composition (EchoMRI):"), " each mouse was weighed and scanned in an EchoMRI 3-in-1; lean mass, fat mass and free-fluid content were reported in grams (three scans per mouse)."),
        tags$li(tags$b("Indirect calorimetry (CLAMS):"), " 48 h in metabolic cages recording movement, oxygen consumption (VO2), carbon dioxide production (VCO2) and food intake every 16 min; the first 24 h were adaptation and the second 24 h were analysed.")
      ),
      h6("Blood, plasma and tissue"),
      tags$ul(
        tags$li(tags$b("Longitudinal sampling:"), " urine and blood were taken at 6, 14 and 22 weeks. Blood glucose (glucometer) and ammonia (PocketChem) were read on fresh blood; plasma was frozen for the other assays."),
        tags$li(tags$b("Terminal collection:"), " mice were sacrificed at 24 weeks after about 17 weeks of diet, following a 4 h morning fast; liver, kidney, heart and spleen were snap-frozen and liver pieces fixed in formalin or OCT for histology."),
        tags$li(tags$b("Plasma biochemistry:"), " creatinine, HDL, LDL, total cholesterol, ASAT, ALT, urea and triglycerides on a Siemens Dimension analyser; TIMP-1 and GDF-15 by Luminex multiplex.")
      ),
      h6("Liver histology"),
      tags$ul(
        tags$li(tags$b("Steatosis (H&E):"), " FFPE liver sections were stained, scanned at 20x, and the area of lipid vacuoles quantified with Visiopharm as a percentage of total liver area."),
        tags$li(tags$b("Fibrosis (Sirius red):"), " OCT cryosections were stained with Sirius Red F3B; whole-slide images were scored in QuPath with a trained pixel classifier, and fibrosis expressed as the percentage of fibrotic area.")
      ),
      h6("Phenotype processing"),
      p("Outliers were removed and each phenotype was transformed toward normality (raw, then log, then Box-Cox or Yeo-Johnson - keeping whichever best passed a normality test). For the genetic mapping, body-composition traits and liver-disease traits were grouped to reveal shared loci, while the remaining phenotypes were analysed individually.")
    ),

    # --- Collapsible parameters box ---
    bslib::accordion(
      id = ns("pheno_params"),
      bslib::accordion_panel(
        title = "Plot Parameters",
        value = "params",
        open = TRUE,

        # Row 1: analysis type + data filters + sizing (always visible)
        fluidRow(
          column(3,
                 selectInput(ns("analysis_type"), label_tip("Analysis Type", "Choose the visualization: heatmaps, body weight plots, variance analysis, or cluster analysis"),
                             choices = list(
                               "Heatmaps" = list(
                                 "Fat Percentage Heatmap" = "fat_heatmap",
                                 "Tissue Weight Z-scores" = "tissue_zscore_heatmap",
                                 "Tissue Relative to Body Weight" = "tissue_rel_heatmap",
                                 "Liver Phenotypes Clustering" = "liver_heatmap",
                                 "Liver F2 Only" = "liver_f2_heatmap"
                               ),
                               "Body Weight Analysis" = list(
                                 "Initial Body Weight" = "bw_initial",
                                 "Final Body Weight" = "bw_final",
                                 "Body Weight Time Course" = "bw_timecourse"
                               ),
                               "Variance Analysis" = list(
                                 "Phenotype Availability" = "pheno_availability",
                                 "Steatosis Variation" = "steatosis_variation",
                                 "Body Weight Variance" = "bw_variance",
                                 "Fat/Lean Variance" = "fatlean_variance",
                                 "Coefficient of Variation" = "cv_analysis"
                               ),
                               "Cluster Analysis" = list(
                                 "Liver Cluster Frequencies" = "cluster_freq",
                                 "Liver Cluster Phenotypes" = "cluster_pheno"
                               )
                             ),
                             selected = "fat_heatmap")
          ),
          column(2, selectInput(ns("sex_filter"), label_tip("Sex", "Filter data by sex"),
                                choices = c("Both" = "both", "Male" = "m", "Female" = "f"), selected = "both")),
          column(2, selectInput(ns("generation_filter"), label_tip("Generation", "Filter by mouse generation"),
                                choices = c("All" = "all", "F0" = "F0", "F1" = "F1", "F2" = "F2"), selected = "all")),
          column(2, selectizeInput(ns("strain_filter"), label_tip("Strain", "Filter by mouse strain; leave empty for all, or pick several"),
                                   choices = c("C57BL/6J" = "C57BL/6J",
                                               "129S1/SvImJ" = "129S1/SvImJ",
                                               "CAST/EiJ" = "CAST/EiJ",
                                               "PWK/PhJ" = "PWK/PhJ",
                                               "B6CASTF1" = "B6CASTF1",
                                               "129SPWKF1" = "129SPWKF1",
                                               "B6CAST-129SPWK-F2" = "B6CAST-129SPWK-F2"),
                                   selected = NULL, multiple = TRUE, width = "100%",
                                   options = list(placeholder = "All strains", plugins = list("remove_button")))),
          column(1, numericInput(ns("plot_width"), label_tip("W (cm)", "Plot width in centimeters"),
                                 value = 25, min = 10, max = 50, step = 1)),
          column(1, numericInput(ns("plot_height"), label_tip("H (cm)", "Plot height in centimeters"),
                                 value = 6, min = 3, max = 20, step = 1)),
          column(1, numericInput(ns("font_size"), label_tip("Font", "Base font size for labels"),
                                 value = 10, min = 6, max = 16, step = 1))
        ),

        # Row 2: analysis-specific options (left) + actions (right).
        # Only one conditionalPanel shows at a time, so this stays a single row.
        fluidRow(
          column(10,
                 # Heatmap options
                 conditionalPanel(
                   condition = paste0("input['", ns("analysis_type"), "'].includes('heatmap')"),
                   fluidRow(
                     column(2,
                            selectInput(ns("color_palette"), label_tip("Color Palette", "Color scale for the heatmap values"),
                                        choices = c("Red-White-Blue" = "rwb",
                                                    "Blue-White-Red" = "bwr",
                                                    "Viridis" = "viridis",
                                                    "Plasma" = "plasma"),
                                        selected = "rwb")
                     ),
                     column(2, div(class = "d-flex align-items-end h-100",
                                   checkboxInput(ns("cluster_columns"), label_tip("Cluster Columns", "Hierarchically cluster and reorder columns by similarity"), value = TRUE))),
                     column(2, div(class = "d-flex align-items-end h-100",
                                   checkboxInput(ns("cluster_rows"), label_tip("Cluster Rows", "Hierarchically cluster and reorder rows by similarity"), value = TRUE))),
                     column(2, div(class = "d-flex align-items-end h-100",
                                   checkboxInput(ns("show_column_names"), label_tip("Show Column Names", "Display individual sample IDs along the heatmap columns"), value = FALSE))),
                     column(2, numericInput(ns("tree_height"), label_tip("Tree ht", "Height of the column dendrogram in points; increase to enlarge the tree relative to the heatmap body"),
                                            value = 150, min = 40, max = 500, step = 20)),
                     column(2,
                            # K-means groups the columns, so it only makes sense when columns
                            # are being clustered - hide it otherwise.
                            div(class = "d-flex align-items-end h-100",
                                conditionalPanel(
                                  condition = paste0("input['", ns("cluster_columns"), "']"),
                                  checkboxInput(ns("enable_kmeans"), label_tip("Enable K-means", "Apply K-means clustering to group columns into K clusters"), value = FALSE),
                                  conditionalPanel(
                                    condition = paste0("input['", ns("enable_kmeans"), "']"),
                                    numericInput(ns("column_km"), label_tip("K Clusters", "Number of clusters for K-means grouping of columns"),
                                                 value = 4, min = 1, max = 10, step = 1)
                                  )
                                )
                            )
                     )
                   )
                 ),
                 # Body-weight time-course options
                 conditionalPanel(
                   condition = paste0("input['", ns("analysis_type"), "'] == 'bw_timecourse'"),
                   fluidRow(
                     column(3,
                            selectInput(ns("bw_metric"), label_tip("Body Weight Metric", "Choose how body weight is expressed: absolute grams, percent gain, or normalized to week 6"),
                                        choices = c("Absolute Weight" = "BW_value",
                                                    "Percent Gain" = "BW_perGain",
                                                    "Normalized to Week 6" = "BW_NormBW6_per"),
                                        selected = "BW_perGain")
                     ),
                     column(2, numericInput(ns("line_alpha"), label_tip("Line Transparency", "Opacity of individual trajectory lines (0.1 = faint, 1 = solid)"),
                                            value = 0.3, min = 0.1, max = 1, step = 0.1)),
                     column(2, numericInput(ns("point_size"), label_tip("Point Size", "Size of data points on the body weight plots"),
                                            value = 0.5, min = 0.1, max = 2, step = 0.1)),
                     column(2, div(class = "d-flex align-items-end h-100",
                                   checkboxInput(ns("show_median_line"), label_tip("Median Line", "Overlay a median trend line per group"), value = FALSE))),
                     column(3,
                            div(class = "d-flex align-items-end h-100",
                                conditionalPanel(
                                  condition = paste0("input['", ns("show_median_line"), "']"),
                                  checkboxInput(ns("show_individual_lines"), label_tip("Individual Lines", "Plot each individual mouse's weight trajectory"), value = TRUE)
                                )
                            )
                     )
                   )
                 ),
                 # Initial / final body-weight options
                 conditionalPanel(
                   condition = paste0(
                     "input['", ns("analysis_type"), "'] == 'bw_initial' || ",
                     "input['", ns("analysis_type"), "'] == 'bw_final'"
                   ),
                   fluidRow(
                     column(4,
                            selectInput(ns("bw_display_type"), label_tip("Display Type", "Bar chart shows group means with error bars; box plot shows distribution"),
                                        choices = c("Bar Chart" = "bar", "Box Plot" = "box"),
                                        selected = "bar")
                     ),
                     column(4, div(class = "d-flex align-items-end h-100",
                                   checkboxInput(ns("separate_by_sex"), label_tip("Separate by Sex", "Split the plot into separate panels for male and female"), value = TRUE)))
                   )
                 )
          ),
          column(2, div(class = "d-flex gap-2 align-items-end justify-content-end h-100 pb-1",
                        actionButton(ns("reset"), "Reset", class = "btn btn-outline-secondary btn-sm"),
                        downloadButton(ns("download_plot"), "Save", class = "btn btn-primary btn-sm")))
        )
      )
    ),

    # --- Refresh: re-render the current plot ---
    div(class = "mt-2 text-center", actionButton(ns("refresh"), "Refresh plot", icon = icon("arrows-rotate"),
                                     class = "btn btn-outline-secondary btn-sm px-4")),

    # --- Plot output ---
    conditionalPanel(
      condition = paste0("input['", ns("analysis_type"), "'].includes('heatmap')"),
      div(class = "mt-3",
          h5("Heatmap Analysis"),
          plotOutput(ns("heatmapPlot"), height = "1000px")
      )
    ),

    conditionalPanel(
      condition = paste0(
        "input['", ns("analysis_type"), "'].includes('bw_') && ",
        "!input['", ns("analysis_type"), "'].includes('variance')"
      ),
      div(class = "mt-3",
          h5("Body Weight Analysis"),
          plotOutput(ns("bodyWeightPlot"), height = "1000px")
      )
    ),

    conditionalPanel(
      condition = paste0(
        "input['", ns("analysis_type"), "'].includes('variance') || ",
        "input['", ns("analysis_type"), "'].includes('cv_') || ",
        "input['", ns("analysis_type"), "'].includes('steatosis') || ",
        "input['", ns("analysis_type"), "'] == 'pheno_availability'"
      ),
      div(class = "mt-3",
          h5("Variance Analysis"),
          plotOutput(ns("variancePlot"), height = "1000px")
      )
    ),

    conditionalPanel(
      condition = paste0("input['", ns("analysis_type"), "'].includes('cluster')"),
      div(class = "mt-3",
          h5("Cluster Analysis"),
          plotOutput(ns("clusterPlot"), height = "1000px")
      )
    ),

    # --- Summary statistics table (conditional) ---
    conditionalPanel(
      condition = paste0("input['", ns("analysis_type"), "'] == 'cv_analysis' || input['", ns("analysis_type"), "'].includes('variance')"),
      div(class = "mt-4",
          div(class = "d-flex justify-content-between align-items-center",
              h5("Summary Statistics", class = "mb-0"),
              downloadButton(ns("download_summary"), "Download CSV", class = "btn btn-outline-secondary btn-sm")),
          DT::dataTableOutput(ns("summaryTable"))
      )
    )
  )
}
