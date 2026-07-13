# RNA-seq hierarchical clustering module server.
# Renders the liver-transcriptome sample clustering as a ComplexHeatmap: the
# transcriptome dendrogram on top (branches coloured by strain), Sex / Strain
# annotation bars, and a z-scored phenotype heatmap body with thin white tile
# borders. Strains can be subset - the subset is re-clustered from the stored
# sample distance matrix so the dendrogram stays clean.
# (Loads Scripts/preprocess_rnaseq_hc.R -> Data/RNAseq/rnaseq_hc.rds.)

# --- Reference data shared across all sessions (rnaseq_hc_path from app.R) ---
.rnaseq_hc_cache <- new.env(parent = emptyenv())
rnaseq_hc_data <- function() {
  if (is.null(.rnaseq_hc_cache$loaded)) {
    .rnaseq_hc_cache$data <- if (exists("rnaseq_hc_path") && file.exists(rnaseq_hc_path)) {
      tryCatch(readRDS(rnaseq_hc_path), error = function(e) NULL)
    } else NULL
    .rnaseq_hc_cache$loaded <- TRUE
  }
  .rnaseq_hc_cache$data
}

# Z-score a phenotype (optionally log1p first). Returns numeric; NA stays NA.
rnaseq_zscore <- function(x, log_it = FALSE) {
  v <- suppressWarnings(as.numeric(x))
  if (isTRUE(log_it)) v <- log1p(pmax(v, 0))
  z <- (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)
  z[!is.finite(z)] <- NA
  z
}

# Build the clustering + matrices used by the ComplexHeatmap renderer.
# `set` is one entry of bundle$sets: list(hc, meta, dmat?). When `strain_filter`
# is non-empty, the kept samples are re-clustered from `set$dmat` (Euclidean,
# complete - same as the precomputed tree) so the dendrogram matches the subset.
rnaseq_hc_prepare <- function(set, bar_cols, bar_labels = bar_cols,
                              strain_filter = NULL, reorder_by = NULL,
                              show_sex = TRUE, show_strain = TRUE,
                              log_bars = TRUE, strain_colors = NULL, sex_colors = NULL) {

  use_subset <- !is.null(strain_filter) && length(strain_filter) > 0

  if (use_subset) {
    keep_ids <- set$meta$Sample_ID[set$meta$Strain %in% strain_filter]
    avail    <- if (!is.null(set$dmat)) rownames(set$dmat) else set$hc$labels
    keep_ids <- intersect(keep_ids, avail)
    if (length(keep_ids) > 1 && !is.null(set$dmat)) {
      # Preferred: re-cluster the kept samples from the stored distance matrix.
      hc <- stats::hclust(stats::as.dist(set$dmat[keep_ids, keep_ids, drop = FALSE]),
                          method = "complete")
    } else if (length(keep_ids) > 1 && requireNamespace("dendextend", quietly = TRUE) &&
               length(keep_ids) < length(set$hc$labels)) {
      # No stored distance matrix: prune the precomputed tree down to the kept
      # samples so strain filtering still works. as.hclust keeps the same label
      # bookkeeping the rest of this function relies on.
      drop_ids <- setdiff(set$hc$labels, keep_ids)
      hc <- tryCatch(
        stats::as.hclust(dendextend::prune(stats::as.dendrogram(set$hc), drop_ids)),
        error = function(e) set$hc)
    } else {
      hc <- set$hc
    }
  } else {
    hc <- set$hc
  }
  # Align the metadata to the (subset) clustering order so colour bars line up.
  meta <- set$meta[match(hc$labels, set$meta$Sample_ID), , drop = FALSE]
  rownames(meta) <- NULL
  ids <- meta$Sample_ID

  dend <- stats::as.dendrogram(hc)
  if (!is.null(reorder_by) && nzchar(reorder_by) && reorder_by %in% names(meta)) {
    wts <- suppressWarnings(as.numeric(meta[[reorder_by]]))
    if (all(is.na(wts))) wts <- seq_len(nrow(meta))
    wts[is.na(wts)] <- stats::median(wts, na.rm = TRUE)
    dend <- stats::reorder(dend, wts = wts)
  }

  # Colour the dendrogram branches by strain (in leaf order), if dendextend is
  # available; ComplexHeatmap honours the edge colours that color_branches sets.
  if (requireNamespace("dendextend", quietly = TRUE) &&
      !is.null(strain_colors) && "Strain" %in% names(meta)) {
    lo <- stats::order.dendrogram(dend)
    bc <- unname(strain_colors[as.character(meta$Strain[lo])]); bc[is.na(bc)] <- "grey50"
    dend <- tryCatch(dendextend::color_branches(dend, col = bc), error = function(e) dend)
  }

  # z-scored phenotype matrix: rows = phenotype bars, columns = samples (hc order).
  rows <- list(); rlabs <- character(0)
  for (i in seq_along(bar_cols)) {
    if (!bar_cols[i] %in% names(meta)) next
    rows[[length(rows) + 1]] <- rnaseq_zscore(meta[[bar_cols[i]]], log_it = log_bars)
    rlabs <- c(rlabs, bar_labels[i])
  }
  mat <- if (length(rows)) do.call(rbind, rows) else matrix(numeric(0), nrow = 0, ncol = nrow(meta))
  if (nrow(mat)) { rownames(mat) <- rlabs; colnames(mat) <- ids }

  list(
    mat    = mat,
    dend   = dend,
    sex    = if (isTRUE(show_sex)    && "Sex"    %in% names(meta)) as.character(meta$Sex)    else NULL,
    strain = if (isTRUE(show_strain) && "Strain" %in% names(meta)) as.character(meta$Strain) else NULL,
    n      = nrow(meta)
  )
}

# Draw the ComplexHeatmap (dendrogram + Sex/Strain annotations + phenotype tiles
# with white borders). Returns invisibly; falls back to a message if the
# Bioconductor packages are unavailable.
rnaseq_hc_render <- function(prep, title = "", is_dark = FALSE, plot_px = 600,
                             tree_px = 150, strain_colors = NULL, sex_colors = NULL) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "The 'pheatmap' package is required for this plot.", cex = 1.1)
    return(invisible(FALSE))
  }
  if (is.null(prep$mat) || nrow(prep$mat) == 0 || ncol(prep$mat) == 0) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No phenotype bars / samples to display.", cex = 1.1)
    return(invisible(FALSE))
  }

  # Diverging blue-white-red palette centred on 0 via symmetric breaks.
  rng <- suppressWarnings(max(abs(prep$mat), na.rm = TRUE))
  if (!is.finite(rng) || rng == 0) rng <- 1
  pal    <- grDevices::colorRampPalette(c("blue", "white", "red"))(256)
  breaks <- seq(-rng, rng, length.out = length(pal) + 1)

  # Top annotation: Strain + Sex, matched to heatmap columns by sample id.
  ann_df <- data.frame(row.names = colnames(prep$mat))
  ann_cols <- list()
  if (!is.null(prep$strain)) {
    ann_df$Strain <- prep$strain
    sc <- strain_colors[intersect(names(strain_colors), unique(prep$strain))]
    if (length(sc)) ann_cols$Strain <- sc
  }
  if (!is.null(prep$sex)) {
    ann_df$Sex <- prep$sex
    scx <- sex_colors[intersect(names(sex_colors), unique(prep$sex))]
    if (length(scx)) ann_cols$Sex <- scx
  }
  has_ann <- ncol(ann_df) > 0

  # pheatmap clusters columns from an hclust; convert the (possibly reordered /
  # strain-filtered) dendrogram back to one. Branch colours are not carried over.
  hc_cols <- tryCatch(stats::as.hclust(prep$dend), error = function(e) FALSE)

  pheatmap::pheatmap(
    prep$mat,
    color             = pal,
    breaks            = breaks,
    cluster_cols      = hc_cols,
    cluster_rows      = FALSE,
    treeheight_col    = if (isTRUE(inherits(hc_cols, "hclust"))) tree_px else 0,
    annotation_col    = if (has_ann) ann_df else NA,
    annotation_colors = if (length(ann_cols)) ann_cols else NA,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    border_color      = "white",
    na_col            = "grey70",
    main              = title,
    silent            = FALSE
  )
  invisible(TRUE)
}

rnaseq_hc_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    strainColors <- c("#000000", "#FF9300", "#107F40", "#0000FE", "#FF2600", "#942192", "#C0C0C0")
    names(strainColors) <- c("C57BL/6J", "129S1/SvImJ", "CAST/EiJ", "PWK/PhJ",
                             "B6CASTF1", "129SPWKF1", "B6CAST-129SPWK-F2")
    sexColors <- c("m" = "#00A2FF", "f" = "#FF644E")

    # Lazy: load the clustering bundle on first render (when the RNA-seq HC tab
    # opens), not at session startup. rnaseq_hc_data() is memoized (loads once).
    delayedAssign("hc_data", rnaseq_hc_data())
    root_input <- session$userData$root_input %||% NULL

    current <- reactive({
      validate(need(!is.null(hc_data),
        paste0("RNA-seq clustering data not found. Run Scripts/preprocess_rnaseq_hc.R ",
               "to generate Data/RNAseq/rnaseq_hc.rds.")))
      set <- input$sample_set %||% "all"
      s <- hc_data$sets[[set]]
      if (is.null(s)) s <- hc_data$sets[["all"]]
      validate(need(!is.null(s), "No data for the selected sample set"))
      s
    })

    prep <- reactive({
      s <- current()
      sf <- input$strain_filter
      bar_cols <- input$bar_vars %||% hc_data$default_bars
      bar_labels <- vapply(bar_cols, function(col) {
        nm <- names(hc_data$pheno_vars)[match(col, hc_data$pheno_vars)]
        if (is.na(nm)) gsub("_", " ", col) else nm
      }, character(1))
      rnaseq_hc_prepare(
        set = s, bar_cols = bar_cols, bar_labels = bar_labels,
        strain_filter = sf,
        reorder_by = input$reorder_by %||% "",
        show_sex = isTRUE(input$show_sex), show_strain = isTRUE(input$show_strain),
        log_bars = isTRUE(input$log_bars),
        strain_colors = strainColors, sex_colors = sexColors
      )
    })

    title_txt <- reactive({
      set <- input$sample_set %||% "all"
      paste0(switch(set, all = "All", male = "Male", female = "Female", "All"),
             " liver transcriptomes")
    })

    # Gate the (expensive) prep + title AND the plot height behind the Refresh
    # button; only dark mode stays live and just redraws the cached prep. Changing
    # the height slider does nothing until Refresh is clicked.
    hc_ready <- eventReactive(input$refresh, {
      list(prep = prep(), title = title_txt(), height = input$plot_height %||% 600,
           tree = input$tree_height %||% 150,
           filtered = length(input$strain_filter) > 0)
    }, ignoreNULL = FALSE)

    output$plot_container <- renderUI({
      h <- hc_ready()$height %||% 600
      plotOutput(session$ns("hc_plot"), height = paste0(h, "px"))
    })

    output$hc_plot <- renderPlot({
      r <- hc_ready()
      is_dark <- !is.null(root_input) && isTRUE(root_input$mode == "dark")
      rnaseq_hc_render(r$prep, title = r$title, is_dark = is_dark, plot_px = r$height,
                       tree_px = r$tree,
                       strain_colors = strainColors, sex_colors = sexColors)
    }, bg = "transparent", execOnResize = TRUE,
       height = function() hc_ready()$height %||% 600) |>
      bindCache(hc_ready(), plot_theme_key(session))

    output$hc_caption <- renderText({
      # Read the refresh-gated result (NOT prep() directly), so selecting strains does
      # not kick off the expensive re-clustering in the background - the caption, like
      # the plot, only updates when Refresh is clicked.
      r <- hc_ready()
      p <- r$prep
      n <- if (is.null(p)) NA else p$n
      paste0("n = ", n, " samples",
             if (!is.null(hc_data$n_genes)) paste0("  |  ", hc_data$n_genes, " genes") else "",
             "  |  Euclidean distance, complete linkage",
             if (isTRUE(r$filtered)) "  |  strain-filtered (re-clustered)" else "")
    })

    output$download_plot <- downloadHandler(
      filename = function() paste0("rnaseq_hclust_", input$sample_set %||% "all", "_", Sys.Date(), ".pdf"),
      content = function(file) {
        pdf(file, width = 14, height = 8)
        on.exit(dev.off(), add = TRUE)
        rnaseq_hc_render(prep(), title = title_txt(), is_dark = FALSE,
                         tree_px = input$tree_height %||% 150,
                         strain_colors = strainColors, sex_colors = sexColors)
      }
    )
  })
}
