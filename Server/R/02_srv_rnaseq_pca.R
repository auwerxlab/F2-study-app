# RNA-seq PCA module server.
# Loads the precomputed PCA (Scripts/preprocess_rnaseq_pca.R ->
# Data/RNAseq/rnaseq_pca.rds) and renders a PCA scatter coloured by a sample
# attribute (Strain, Sex) or any phenotype (liver weight, steatosis, ...).

# --- Reference data shared across all sessions (rnaseq_pca_path from app.R) ---
.rnaseq_cache <- new.env(parent = emptyenv())
rnaseq_pca_data <- function() {
  if (is.null(.rnaseq_cache$loaded)) {
    .rnaseq_cache$data <- if (exists("rnaseq_pca_path") && file.exists(rnaseq_pca_path)) {
      tryCatch(readRDS(rnaseq_pca_path), error = function(e) NULL)
    } else NULL
    .rnaseq_cache$loaded <- TRUE
  }
  .rnaseq_cache$data
}

# --- Plot builder (pure function: easy to test without Shiny) ---------------
rnaseq_build_pca_plot <- function(df, eig, x_dim, y_dim, color_by, color_label,
                                  is_discrete, do_log, shape_by, point_size,
                                  strain_colors, sex_colors, is_dark = FALSE) {
  text_col  <- if (is_dark) "#E0E0E0" else "black"
  line_col  <- if (is_dark) "#AAAAAA" else "black"
  bg_col    <- if (is_dark) "#222222" else "white"

  x_idx <- suppressWarnings(as.integer(sub("^Dim\\.", "", x_dim)))
  y_idx <- suppressWarnings(as.integer(sub("^Dim\\.", "", y_dim)))
  pct   <- function(i) if (!is.na(i) && i <= length(eig)) paste0(" (", round(eig[i], 2), "%)") else ""
  xlab  <- paste0("PC", x_idx, pct(x_idx))
  ylab  <- paste0("PC", y_idx, pct(y_idx))

  df$pcx <- df[[x_dim]]
  df$pcy <- df[[y_dim]]
  if (is_discrete) {
    df$col_val <- factor(df[[color_by]])
  } else {
    v <- suppressWarnings(as.numeric(df[[color_by]]))
    df$col_val <- if (isTRUE(do_log)) log1p(v) else v
  }
  use_shape <- identical(shape_by, "Sex") && "Sex" %in% names(df)
  if (use_shape) df$shape_val <- factor(df$Sex)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = pcx, y = pcy))
  if (use_shape) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = col_val, shape = shape_val),
                                 size = point_size, na.rm = TRUE)
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = col_val),
                                 size = point_size, na.rm = TRUE)
  }

  if (is_discrete) {
    pal <- if (identical(color_by, "Strain")) strain_colors
           else if (identical(color_by, "Sex")) sex_colors else NULL
    if (!is.null(pal)) {
      p <- p + ggplot2::scale_color_manual(name = color_label, values = pal,
                                           drop = TRUE, na.value = "grey70")
    } else {
      p <- p + ggplot2::scale_color_discrete(name = color_label, na.value = "grey70")
    }
  } else {
    lab <- if (isTRUE(do_log)) paste0(color_label, " (log)") else color_label
    p <- p + ggplot2::scale_color_gradientn(name = lab,
                                            colors = c("blue", "white", "red"),
                                            na.value = "grey70")
  }
  if (use_shape) {
    p <- p + ggplot2::scale_shape_manual(name = "Sex",
                                         values = c("m" = 16, "f" = 17,
                                                    "male" = 16, "female" = 17),
                                         na.value = 4)
  }

  p +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = bg_col, color = NA),
      panel.background  = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.background = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.key        = ggplot2::element_rect(fill = bg_col, color = NA),
      axis.title  = ggplot2::element_text(color = text_col, size = 16),
      axis.text   = ggplot2::element_text(color = text_col, size = 13),
      axis.line   = ggplot2::element_line(color = line_col),
      axis.ticks  = ggplot2::element_line(color = line_col),
      legend.title = ggplot2::element_text(color = text_col),
      legend.text  = ggplot2::element_text(color = text_col),
      panel.grid   = ggplot2::element_blank()
    )
}

rnaseq_pca_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    library(ggplot2)

    strainColors <- c("#000000", "#FF9300", "#107F40", "#0000FE", "#FF2600", "#942192", "#C0C0C0")
    names(strainColors) <- c("C57BL/6J", "129S1/SvImJ", "CAST/EiJ", "PWK/PhJ",
                             "B6CASTF1", "129SPWKF1", "B6CAST-129SPWK-F2")
    sexColors <- c("m" = "#00A2FF", "f" = "#FF644E")

    # Lazy: load the PCA bundle on first render (when the RNA-seq PCA tab opens),
    # not at session startup. rnaseq_pca_data() is memoized, so it loads once.
    delayedAssign("pca", rnaseq_pca_data())
    root_input <- session$userData$root_input %||% NULL

    # Selected sample-set frame (+ eigenvalues)
    current <- reactive({
      validate(need(!is.null(pca),
        paste0("RNA-seq PCA data not found. Run Scripts/preprocess_rnaseq_pca.R ",
               "to generate Data/RNAseq/rnaseq_pca.rds.")))
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
             if (!is.null(pca$n_genes)) paste0("  |  ", pca$n_genes, " genes used") else "")
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("rnaseq_pca_", input$sample_set %||% "all", "_",
               input$color_by %||% "Strain", "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        p <- tryCatch(pca_plot(), error = function(e) NULL)
        pdf(file, width = 7, height = 5)
        on.exit(dev.off(), add = TRUE)
        if (is.null(p)) {
          grid::grid.newpage()
          grid::grid.text("No plot to render (RNA-seq PCA data unavailable).")
        } else {
          print(p)
        }
      }
    )
  })
}
