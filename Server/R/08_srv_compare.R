# Protein vs RNA-seq comparison module server (the "Data comparison" tab).
# Loads the precomputed bundle (Scripts/preprocess_prot_rna.R ->
# Data/Comparison/prot_rna.rds) and shows RNA-protein correlation per sample,
# per gene, and a gene-level scatter, as in 20-Prot_RNAseq_comparison.R.

# --- Reference data shared across all sessions (prot_rna_path from app.R) ---
.prot_rna_cache <- new.env(parent = emptyenv())
prot_rna_data <- function() {
  # Retry until the file successfully loads, so a bundle generated *after* the app
  # started is picked up on the next session without a full process restart
  # (instead of caching the "missing" state permanently).
  if (is.null(.prot_rna_cache$data) &&
      exists("prot_rna_path") && file.exists(prot_rna_path)) {
    .prot_rna_cache$data <- tryCatch(readRDS(prot_rna_path), error = function(e) NULL)
  }
  .prot_rna_cache$data
}

# Common theme bits
.pr_theme <- function(is_dark) {
  text_col <- if (is_dark) "#E0E0E0" else "black"
  line_col <- if (is_dark) "#AAAAAA" else "grey40"
  bg_col   <- if (is_dark) "#222222" else "white"
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = bg_col, color = NA),
      panel.background  = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.background = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.key        = ggplot2::element_rect(fill = bg_col, color = NA),
      plot.title   = ggplot2::element_text(color = text_col, size = 13),
      axis.title   = ggplot2::element_text(color = text_col),
      axis.text    = ggplot2::element_text(color = text_col),
      axis.line    = ggplot2::element_line(color = line_col),
      axis.ticks   = ggplot2::element_line(color = line_col),
      legend.title = ggplot2::element_text(color = text_col),
      legend.text  = ggplot2::element_text(color = text_col)
    )
}
pr_corr_col <- function(method) if (identical(method, "spearman")) "corr_spearman" else "corr_pearson"

# Per-sample correlation histogram
protrna_sample_hist <- function(df, method = "pearson", is_dark = FALSE) {
  cc <- pr_corr_col(method)
  df <- df[is.finite(df[[cc]]), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[cc]])) +
    ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "black", linewidth = 0.3) +
    ggplot2::labs(x = paste0("RNA-protein correlation (", method, ")"), y = "Samples (n)") +
    .pr_theme(is_dark)
}

# Per-sample correlation barplot, coloured by strain or cluster
protrna_sample_bar <- function(df, method = "pearson", color_by = "Strain",
                               strain_colors = NULL, cluster_colors = NULL, is_dark = FALSE) {
  cc <- pr_corr_col(method)
  df <- df[is.finite(df[[cc]]), , drop = FALSE]
  pal <- if (color_by == "liver_cluster") cluster_colors else strain_colors
  df$.fill <- factor(as.character(df[[color_by]]))
  df$.id <- factor(df$ID, levels = df$ID[order(df[[cc]])])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .id, y = .data[[cc]], fill = .fill)) +
    ggplot2::geom_col()
  if (!is.null(pal)) p <- p + ggplot2::scale_fill_manual(name = color_by, values = pal, na.value = "grey70")
  p +
    ggplot2::labs(x = "Sample (ranked)", y = paste0("Correlation (", method, ")")) +
    .pr_theme(is_dark) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

# Per-gene correlation histogram
protrna_gene_hist <- function(df, method = "pearson", is_dark = FALSE) {
  cc <- pr_corr_col(method)
  df <- df[is.finite(df[[cc]]), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[cc]])) +
    ggplot2::geom_histogram(binwidth = 0.05, fill = "grey70", color = "black", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = if (is_dark) "#AAAAAA" else "grey40") +
    ggplot2::labs(x = paste0("RNA-protein correlation (", method, ")"), y = "Genes (n)") +
    .pr_theme(is_dark)
}

# Per-gene correlation density by liver cluster
protrna_gene_density <- function(df, method = "pearson", cluster_colors = NULL, is_dark = FALSE) {
  cc <- pr_corr_col(method)
  df <- df[is.finite(df[[cc]]) & !is.na(df$liver_cluster), , drop = FALSE]
  df$Cluster <- factor(as.character(df$liver_cluster))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[cc]], color = Cluster, fill = Cluster)) +
    ggplot2::geom_density(alpha = 0.15)
  if (!is.null(cluster_colors)) {
    p <- p + ggplot2::scale_color_manual(values = cluster_colors, na.value = "grey70") +
      ggplot2::scale_fill_manual(values = cluster_colors, na.value = "grey70")
  }
  p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = if (is_dark) "#AAAAAA" else "grey40") +
    ggplot2::labs(x = paste0("RNA-protein correlation (", method, ")"), y = "Density") +
    .pr_theme(is_dark)
}

# Gene-level scatter: log2 protein vs log2 CPM across mice
protrna_gene_scatter <- function(g, gene_label = "", color_by = "liver_cluster",
                                 strain_colors = NULL, cluster_colors = NULL, is_dark = FALSE) {
  text_col <- if (is_dark) "#E0E0E0" else "black"
  g <- g[is.finite(g$log2_Protein) & is.finite(g$log2CPM_RNA), , drop = FALSE]
  if (nrow(g) < 2) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0, y = 0, label = "Not enough data for this gene",
                               color = text_col, size = 4.5))
  }
  pal <- if (color_by == "Strain") strain_colors else cluster_colors
  g$.col <- factor(as.character(g[[color_by]]))
  r <- suppressWarnings(stats::cor(g$log2CPM_RNA, g$log2_Protein, use = "complete.obs"))
  ttl <- paste0(gene_label, "   (overall r = ", round(r, 2), ", n = ", nrow(g), ")")

  p <- ggplot2::ggplot(g, ggplot2::aes(x = log2CPM_RNA, y = log2_Protein, color = .col)) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = if (is_dark) "#AAAAAA" else "grey50",
                         linewidth = 0.6, formula = y ~ x, na.rm = TRUE) +
    ggplot2::geom_point(size = 2, alpha = 0.85, na.rm = TRUE)
  if (!is.null(pal)) p <- p + ggplot2::scale_color_manual(name = color_by, values = pal, na.value = "grey70")
  p +
    ggplot2::labs(x = "RNA (log2 CPM)", y = "Protein (log2)", title = ttl) +
    .pr_theme(is_dark)
}

compare_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    library(ggplot2)

    # Lazy: prot_rna_data() loads an 82 MB bundle (~4 s to decompress). Binding it
    # as a promise defers that cost until the first "Data comparison" output is
    # rendered (i.e. the user opens this tab) instead of blocking every session at
    # startup. prot_rna_data() is itself memoized, so the load still happens once.
    delayedAssign("pr", prot_rna_data())
    root_input <- session$userData$root_input %||% NULL
    is_dark <- reactive(!is.null(root_input) && isTRUE(root_input$mode == "dark"))

    # Populate the gene selector server-side (efficient for ~7.9k genes: only the
    # matching options are sent as the user types). Server-side selectize does not
    # initialize while its nav_panel is hidden, so defer the one-time population
    # until the "Gene scatter" tab is first shown.
    gene_ready <- reactiveVal(FALSE)
    # Bumped once, when the Gene scatter tab is first opened, so the first gene's
    # scatter draws immediately instead of waiting for a manual Refresh click.
    # Later gene / colour changes still wait for this panel's Refresh button.
    scatter_init <- reactiveVal(0)
    observeEvent(input$tabs, {
      if (identical(input$tabs, "Gene scatter") && !gene_ready() &&
          !is.null(pr) && length(pr$gene_choices)) {
        updateSelectizeInput(session, "gene", choices = pr$gene_choices,
                             selected = unname(pr$gene_choices[1]), server = TRUE)
        gene_ready(TRUE)
        scatter_init(scatter_init() + 1)
      }
    })

    need_data <- function() {
      validate(need(!is.null(pr),
        paste0("Protein/RNA comparison data not found. Run Scripts/preprocess_prot_rna.R ",
               "to generate Data/Comparison/prot_rna.rds.")))
    }

    # ---- Per-sample ----
    output$sample_hist <- renderPlot({
      input$refresh_s   # only re-render on this sub-tab's Refresh (isolate the rest)
      isolate({
        need_data()
        protrna_sample_hist(pr$corr_per_sample, input$smethod %||% "pearson", is_dark())
      })
    }, bg = "transparent")

    output$sample_bar <- renderPlot({
      input$refresh_s   # only re-render on this sub-tab's Refresh (isolate the rest)
      isolate({
        need_data()
        protrna_sample_bar(pr$corr_per_sample, input$smethod %||% "pearson",
                           input$scolor %||% "Strain", pr$strain_colors, pr$cluster_colors, is_dark())
      })
    }, bg = "transparent")

    output$sample_caption <- renderText({
      need_data()
      cc <- pr_corr_col(input$smethod %||% "pearson")
      v <- pr$corr_per_sample[[cc]]; v <- v[is.finite(v)]
      paste0("n = ", length(v), " samples  |  median r = ", round(stats::median(v), 3),
             "  |  range ", round(min(v), 2), " to ", round(max(v), 2))
    })

    # ---- Per-gene ----
    output$gene_hist <- renderPlot({
      input$refresh_g   # only re-render on this sub-tab's Refresh (isolate the rest)
      isolate({
        need_data()
        protrna_gene_hist(pr$corr_per_gene, input$gmethod %||% "pearson", is_dark())
      })
    }, bg = "transparent")

    output$gene_density <- renderPlot({
      input$refresh_g   # only re-render on this sub-tab's Refresh (isolate the rest)
      isolate({
        need_data()
        validate(need(!is.null(pr$corr_per_gene_cluster) && nrow(pr$corr_per_gene_cluster) > 0,
                      "No per-cluster correlations available."))
        protrna_gene_density(pr$corr_per_gene_cluster, input$gmethod %||% "pearson",
                             pr$cluster_colors, is_dark())
      })
    }, bg = "transparent")

    # Per-gene correlation table data, shared by the on-screen table and its CSV download.
    gene_table_data <- reactive({
      need_data()
      df <- pr$corr_per_gene
      df <- df[, intersect(c("Gene.Symbol", "gene_id", "corr_pearson", "corr_spearman", "n_mice"), names(df))]
      for (cc in intersect(c("corr_pearson", "corr_spearman"), names(df))) df[[cc]] <- round(df[[cc]], 3)
      df
    })

    output$gene_table <- DT::renderDataTable({
      df <- gene_table_data()
      DT::datatable(df, rownames = FALSE, filter = "top", selection = "none",
                    options = list(pageLength = 12, scrollX = TRUE,
                                   order = list(list(which(names(df) == "corr_pearson") - 1, "desc"))))
    })

    output$download_gene_table <- downloadHandler(
      filename = function() paste0("prot_rna_per_gene_", Sys.Date(), ".csv"),
      content = function(file) {
        d <- tryCatch(gene_table_data(), error = function(e) data.frame())
        write.csv(d, file, row.names = FALSE)
      }
    )

    # ---- Gene scatter ----
    gene_subset <- reactive({
      need_data()
      # Before the server-side selectize reports its value, input$gene can be NULL
      # OR an empty string; fall back to the first gene in both cases so the initial
      # scatter is a real gene (req("") would otherwise silently blank the plot).
      gid <- input$gene
      if (is.null(gid) || !nzchar(gid))
        gid <- if (length(pr$gene_choices)) unname(pr$gene_choices[1]) else NULL
      req(gid)
      g <- pr$pairs[pr$pairs$gene_id == gid, , drop = FALSE]
      list(g = g, gid = gid,
           label = { nm <- names(pr$gene_choices)[match(gid, pr$gene_choices)]; if (is.na(nm)) gid else nm })
    })

    output$gene_scatter <- renderPlot({
      input$refresh_sc   # re-render on this sub-tab's Refresh...
      scatter_init()     # ...and once, when the tab is first opened (first gene)
      isolate({
        s <- gene_subset()
        protrna_gene_scatter(s$g, gene_label = s$label, color_by = input$pcolor %||% "liver_cluster",
                             strain_colors = pr$strain_colors, cluster_colors = pr$cluster_colors,
                             is_dark = is_dark())
      })
    }, bg = "transparent")

    output$scatter_caption <- renderText({
      s <- gene_subset()
      if (is.null(pr$corr_per_gene_cluster)) return("")
      sub <- pr$corr_per_gene_cluster[pr$corr_per_gene_cluster$gene_id == s$gid, , drop = FALSE]
      if (nrow(sub) == 0) return(paste0(s$label, ": no per-cluster correlations (needs >=5 mice/cluster)."))
      cc <- pr_corr_col(input$gmethod %||% "pearson")
      parts <- vapply(seq_len(nrow(sub)), function(i)
        paste0("cluster ", sub$liver_cluster[i], ": r = ", round(sub[[cc]][i], 2)), character(1))
      paste0(s$label, " - per-cluster correlation  |  ", paste(parts, collapse = "   |   "))
    })

    # ---- Downloads ----
    dl <- function(name, plot_fun) downloadHandler(
      filename = function() paste0("prot_rna_", name, "_", Sys.Date(), ".pdf"),
      content = function(file) {
        p <- tryCatch(plot_fun(), error = function(e) NULL)
        pdf(file, width = 7, height = 5); on.exit(dev.off(), add = TRUE)
        if (is.null(p)) { grid::grid.newpage(); grid::grid.text("No plot.") } else print(p)
      })
    output$download_sample <- dl("per_sample", function()
      protrna_sample_bar(pr$corr_per_sample, input$smethod %||% "pearson", input$scolor %||% "Strain",
                         pr$strain_colors, pr$cluster_colors))
    output$download_gene <- dl("per_gene", function()
      protrna_gene_hist(pr$corr_per_gene, input$gmethod %||% "pearson"))
    output$download_scatter <- dl("gene_scatter", function() {
      s <- gene_subset()
      protrna_gene_scatter(s$g, gene_label = s$label, color_by = input$pcolor %||% "liver_cluster",
                           strain_colors = pr$strain_colors, cluster_colors = pr$cluster_colors)
    })
  })
}
