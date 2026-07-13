# Proteomics PCA module server.
# Loads the precomputed PCA (Scripts/preprocess_proteomics_pca.R ->
# Data/Proteomics/proteomics_pca.rds) and renders a PCA scatter coloured by a
# sample attribute (Strain, Sex, Generation) or any phenotype. Reuses the pure
# plot builder rnaseq_build_pca_plot() defined in Server/R/02_srv_rnaseq_pca.R.

# --- Reference data shared across all sessions (proteomics_pca_path from app.R) ---
.proteomics_pca_cache <- new.env(parent = emptyenv())
proteomics_pca_data <- function() {
  if (is.null(.proteomics_pca_cache$loaded)) {
    .proteomics_pca_cache$data <- if (exists("proteomics_pca_path") && file.exists(proteomics_pca_path)) {
      tryCatch(readRDS(proteomics_pca_path), error = function(e) NULL)
    } else NULL
    .proteomics_pca_cache$loaded <- TRUE
  }
  .proteomics_pca_cache$data
}

proteomics_pca_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    library(ggplot2)

    strainColors <- c("#000000", "#FF9300", "#107F40", "#0000FE", "#FF2600", "#942192", "#C0C0C0")
    names(strainColors) <- c("C57BL/6J", "129S1/SvImJ", "CAST/EiJ", "PWK/PhJ",
                             "B6CASTF1", "129SPWKF1", "B6CAST-129SPWK-F2")
    sexColors <- c("m" = "#00A2FF", "f" = "#FF644E")

    # Lazy: load the PCA bundle on first render (when the Proteomics tab opens),
    # not at session startup. proteomics_pca_data() is memoized (loads once).
    delayedAssign("pca", proteomics_pca_data())
    root_input <- session$userData$root_input %||% NULL

    current <- reactive({
      validate(need(!is.null(pca),
        paste0("Proteomics PCA data not found. Run Scripts/preprocess_proteomics_pca.R ",
               "to generate Data/Proteomics/proteomics_pca.rds.")))
      set <- input$sample_set %||% "all"
      s <- pca$sets[[set]]
      if (is.null(s)) s <- pca$sets[["all"]]
      validate(need(!is.null(s), "No data for the selected sample set"))
      s
    })

    pca_plot <- reactive({
      s  <- current()
      df <- s$coord

      # Optional strain subset (empty = all strains).
      sf <- input$strain_filter
      if (length(sf) && "Strain" %in% names(df)) {
        df <- df[df$Strain %in% sf, , drop = FALSE]
        validate(need(nrow(df) > 0, "No samples for the selected strain(s)."))
      }

      color_by    <- input$color_by %||% "Strain"
      is_discrete <- color_by %in% unname(pca$discrete_vars)
      color_label <- if (is_discrete) {
        names(pca$discrete_vars)[match(color_by, pca$discrete_vars)]
      } else {
        nm <- names(pca$pheno_vars)[match(color_by, pca$pheno_vars)]
        if (is.na(nm)) gsub("_", " ", color_by) else nm
      }
      validate(need(color_by %in% names(df),
                    paste0("Column '", color_by, "' is not in the data")))

      is_dark <- !is.null(root_input) && isTRUE(root_input$mode == "dark")

      rnaseq_build_pca_plot(
        df = df, eig = s$eig,
        x_dim = input$x_dim %||% "Dim.1",
        y_dim = input$y_dim %||% "Dim.2",
        color_by = color_by, color_label = color_label,
        is_discrete = is_discrete,
        do_log = isTRUE(input$log_color) && !is_discrete,
        shape_by = input$shape_by %||% "none",
        point_size = input$point_size %||% 2,
        strain_colors = strainColors, sex_colors = sexColors,
        is_dark = is_dark
      )
    })

    output$plot_container <- renderUI({
      h <- input$plot_height %||% 550
      plotOutput(session$ns("pca_plot"), height = paste0(h, "px"))
    })

    # Build only on Refresh click (and once on load); parameter changes wait.
    pca_plot_g <- eventReactive(input$refresh, pca_plot(), ignoreNULL = FALSE)
    output$pca_plot <- renderPlot({ pca_plot_g() }, bg = "transparent") |>
      bindCache(pca_plot_g(), plot_theme_key(session))

    output$pca_caption <- renderText({
      s <- current()
      paste0("n = ", nrow(s$coord), " samples",
             if (!is.null(pca$n_proteins)) paste0("  |  ", pca$n_proteins, " proteins used") else "")
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("proteomics_pca_", input$sample_set %||% "all", "_",
               input$color_by %||% "Strain", "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        p <- tryCatch(pca_plot(), error = function(e) NULL)
        pdf(file, width = 7, height = 5)
        on.exit(dev.off(), add = TRUE)
        if (is.null(p)) {
          grid::grid.newpage()
          grid::grid.text("No plot to render (proteomics PCA data unavailable).")
        } else {
          print(p)
        }
      }
    )
  })
}
