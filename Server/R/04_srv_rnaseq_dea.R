# RNA-seq differential expression module server.
# Browses the precomputed limma-voom results (Scripts/preprocess_rnaseq_dea.R ->
# Data/RNAseq/rnaseq_dea.rds): per-contrast volcano + results table, plus a
# DEG-count overview across contrasts at the chosen thresholds.

# --- Reference data shared across all sessions (rnaseq_dea_path from app.R) ---
.rnaseq_dea_cache <- new.env(parent = emptyenv())
rnaseq_dea_data <- function() {
  if (is.null(.rnaseq_dea_cache$loaded)) {
    .rnaseq_dea_cache$data <- if (exists("rnaseq_dea_path") && file.exists(rnaseq_dea_path)) {
      tryCatch(readRDS(rnaseq_dea_path), error = function(e) NULL)
    } else NULL
    .rnaseq_dea_cache$loaded <- TRUE
  }
  .rnaseq_dea_cache$data
}

# Classify each gene as Up / Down / NS at the given FDR and |logFC| thresholds.
rnaseq_classify_de <- function(tbl, sig_fdr, lfc) {
  sig <- !is.na(tbl$adj.P.Val) & tbl$adj.P.Val < sig_fdr & abs(tbl$logFC) > lfc
  factor(ifelse(!sig, "NS", ifelse(tbl$logFC > 0, "Up", "Down")),
         levels = c("Up", "Down", "NS"))
}

# Volcano plot for one contrast (pure function; ggplot2 + ggrepel).
rnaseq_build_volcano <- function(tbl, label, sig_fdr = 0.05, lfc = 0, top_n = 15,
                                 show_ns = TRUE, is_dark = FALSE) {
  text_col <- if (is_dark) "#E0E0E0" else "black"
  line_col <- if (is_dark) "#AAAAAA" else "grey40"
  bg_col   <- if (is_dark) "#222222" else "white"

  tbl <- tbl[!is.na(tbl$logFC) & !is.na(tbl$adj.P.Val), , drop = FALSE]
  tbl$de <- rnaseq_classify_de(tbl, sig_fdr, lfc)
  tbl$neglog10 <- -log10(tbl$adj.P.Val)
  fin <- tbl$neglog10[is.finite(tbl$neglog10)]
  cap <- if (length(fin)) max(fin) else 1
  tbl$neglog10[!is.finite(tbl$neglog10)] <- cap        # cap adj.P == 0

  # Hide ALL non-significant (grey) points - those under the FDR threshold OR the
  # |logFC| threshold - leaving only the called Up/Down genes. de == "NS" is
  # exactly (adj.P >= FDR) OR (|logFC| <= lfc), so one filter covers both. The
  # y-axis is anchored at 0 below, so this only removes points, it does not rescale.
  if (!isTRUE(show_ns)) tbl <- tbl[tbl$de != "NS", , drop = FALSE]

  lab <- tbl[tbl$de != "NS" & !is.na(tbl$gene_name) & nzchar(tbl$gene_name), , drop = FALSE]
  lab <- lab[order(lab$adj.P.Val, lab$P.Value), , drop = FALSE]
  if (nrow(lab) > top_n) lab <- lab[seq_len(top_n), , drop = FALSE]

  cols <- c("Up" = "#D55E00", "Down" = "#1f78b4", "NS" = "grey70")

  p <- ggplot2::ggplot(tbl, ggplot2::aes(x = logFC, y = neglog10, color = de)) +
    ggplot2::geom_point(size = 1, alpha = 0.7, na.rm = TRUE) +
    ggplot2::scale_color_manual(name = NULL, values = cols, breaks = c("Up", "Down"),
                                drop = FALSE) +
    ggplot2::geom_hline(yintercept = -log10(sig_fdr), linetype = "dashed", color = line_col)
  if (lfc > 0) {
    p <- p + ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = line_col)
  }
  if (nrow(lab) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = lab, ggplot2::aes(label = gene_name), size = 4.5, color = text_col,
      max.overlaps = 20, box.padding = 0.4, na.rm = TRUE, show.legend = FALSE)
  }
  p +
    ggplot2::expand_limits(y = 0) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::labs(x = "logFC (log2)", y = "-log10(adj. P)", title = label) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = bg_col, color = NA),
      panel.background  = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.background = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.key        = ggplot2::element_rect(fill = bg_col, color = NA),
      plot.title  = ggplot2::element_text(color = text_col, size = 16),
      axis.title  = ggplot2::element_text(color = text_col, size = 15),
      axis.text   = ggplot2::element_text(color = text_col, size = 13),
      axis.line   = ggplot2::element_line(color = line_col),
      axis.ticks  = ggplot2::element_line(color = line_col),
      legend.text = ggplot2::element_text(color = text_col)
    )
}

# DEG counts (Up / Down) per contrast at the chosen thresholds.
rnaseq_deg_counts <- function(contrasts, sig_fdr, lfc) {
  do.call(rbind, lapply(names(contrasts), function(k) {
    de <- rnaseq_classify_de(contrasts[[k]]$table, sig_fdr, lfc)
    data.frame(key = k, label = contrasts[[k]]$label,
               Up = sum(de == "Up"), Down = sum(de == "Down"),
               stringsAsFactors = FALSE)
  }))
}

# Diverging bar of DEG counts across contrasts (Up up, Down down).
rnaseq_build_deg_bar <- function(counts_df, is_dark = FALSE) {
  text_col <- if (is_dark) "#E0E0E0" else "black"
  line_col <- if (is_dark) "#AAAAAA" else "grey40"
  bg_col   <- if (is_dark) "#222222" else "white"
  long <- rbind(
    data.frame(label = counts_df$label, Direction = "Up",   n =  counts_df$Up),
    data.frame(label = counts_df$label, Direction = "Down", n = -counts_df$Down)
  )
  long$label <- factor(long$label, levels = rev(counts_df$label))
  ggplot2::ggplot(long, ggplot2::aes(x = label, y = n, fill = Direction)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = line_col) +
    ggplot2::scale_fill_manual(values = c("Up" = "#D55E00", "Down" = "#1f78b4")) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "DEGs (n)") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = bg_col, color = NA),
      panel.background  = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.background = ggplot2::element_rect(fill = bg_col, color = NA),
      axis.title = ggplot2::element_text(color = text_col, size = 15),
      axis.text  = ggplot2::element_text(color = text_col, size = 13),
      axis.line  = ggplot2::element_line(color = line_col),
      axis.ticks = ggplot2::element_line(color = line_col),
      legend.text  = ggplot2::element_text(color = text_col),
      legend.title = ggplot2::element_text(color = text_col)
    )
}

# --- GSEA (optional layer: Scripts/preprocess_rnaseq_gsea.R -> rnaseq_gsea.rds) ---
.rnaseq_gsea_cache <- new.env(parent = emptyenv())
rnaseq_gsea_data <- function() {
  if (is.null(.rnaseq_gsea_cache$loaded)) {
    .rnaseq_gsea_cache$data <- if (exists("rnaseq_gsea_path") && file.exists(rnaseq_gsea_path)) {
      tryCatch(readRDS(rnaseq_gsea_path), error = function(e) NULL)
    } else NULL
    .rnaseq_gsea_cache$loaded <- TRUE
  }
  .rnaseq_gsea_cache$data
}

# Wrap long GO term descriptions for the axis, capped at 2 lines ("..."-truncated if
# longer) so tall multi-line labels do not overlap on a dense dot plot. Plain ASCII
# "..." (not the U+2026 ellipsis) so it renders cleanly in the PDF downloads too.
rnaseq_wrap_label <- function(x, width = 34, max_lines = 2) {
  vapply(x, function(s) {
    w <- strwrap(s, width = width)
    if (length(w) > max_lines) {
      w <- w[seq_len(max_lines)]
      w[max_lines] <- paste0(sub("\\s+$", "", w[max_lines]), "...")
    }
    paste(w, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# GSEA dot plot: top-N enriched up + down GO terms (NES, coloured by
# -log10(p.adjust), sized by gene ratio), as in 04-RNAseq_DEA.R.
rnaseq_build_gsea_dot <- function(res, label = "", top_n = 10, padj_max = 0.05, is_dark = FALSE) {
  text_col <- if (is_dark) "#E0E0E0" else "black"
  line_col <- if (is_dark) "#AAAAAA" else "grey40"
  bg_col   <- if (is_dark) "#222222" else "white"

  empty_plot <- function(msg) {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg, color = text_col, size = 4.5) +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = bg_col, color = NA))
  }
  if (is.null(res) || nrow(res) == 0) return(empty_plot("No GSEA results for this contrast"))

  res <- res[!is.na(res$p.adjust) & res$p.adjust <= padj_max & !is.na(res$NES), , drop = FALSE]
  if (nrow(res) == 0) return(empty_plot(paste0("No GO terms at p.adjust <= ", padj_max)))

  up <- res[res$NES > 0, , drop = FALSE]; up <- utils::head(up[order(up$p.adjust), , drop = FALSE], top_n)
  dn <- res[res$NES < 0, , drop = FALSE]; dn <- utils::head(dn[order(dn$p.adjust), , drop = FALSE], top_n)
  d <- rbind(dn, up)
  if (nrow(d) == 0) return(empty_plot("No GO terms to display"))

  d <- d[order(d$NES), , drop = FALSE]
  d$Description <- factor(rnaseq_wrap_label(d$Description), levels = rnaseq_wrap_label(d$Description))
  d$neglog_padj <- -log10(d$p.adjust)

  ggplot2::ggplot(d, ggplot2::aes(x = Description, y = NES, color = neglog_padj, size = GeneRatio)) +
    ggplot2::geom_hline(yintercept = 0, color = line_col) +
    ggplot2::geom_point() +
    ggplot2::coord_flip() +
    ggplot2::scale_color_viridis_c(name = "-log10(p.adj)", option = "magma", direction = -1,
                                   limits = c(0, NA)) +
    ggplot2::scale_size_continuous(name = "Gene ratio", range = c(2, 8)) +
    ggplot2::labs(x = NULL, y = "Normalized Enrichment Score (NES)", title = label) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = bg_col, color = NA),
      panel.background  = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.background = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.key        = ggplot2::element_rect(fill = bg_col, color = NA),
      plot.title  = ggplot2::element_text(color = text_col, size = 14),
      axis.title  = ggplot2::element_text(color = text_col, size = 14),
      axis.text   = ggplot2::element_text(color = text_col, size = 12),
      axis.line   = ggplot2::element_line(color = line_col),
      axis.ticks  = ggplot2::element_line(color = line_col),
      legend.text  = ggplot2::element_text(color = text_col),
      legend.title = ggplot2::element_text(color = text_col)
    )
}

# --- Set-specific ORA (optional: Scripts/preprocess_rnaseq_ora.R -> rnaseq_ora.rds) ---
.rnaseq_ora_cache <- new.env(parent = emptyenv())
rnaseq_ora_data <- function() {
  if (is.null(.rnaseq_ora_cache$loaded)) {
    .rnaseq_ora_cache$data <- if (exists("rnaseq_ora_path") && file.exists(rnaseq_ora_path)) {
      tryCatch(readRDS(rnaseq_ora_path), error = function(e) NULL)
    } else NULL
    .rnaseq_ora_cache$loaded <- TRUE
  }
  .rnaseq_ora_cache$data
}

# Friendly label for a coefficient key from the ORA bundle's choice map.
ora_label <- function(ora, key) {
  nm <- names(ora$coef_choices)[match(key, ora$coef_choices)]
  if (is.na(nm)) key else nm
}

# Diverging up/down GO dot plot for one combination (ORA), as in 04-RNAseq_DEA.R:
# up terms get +(-log10 p.adjust), down terms get -(-log10 p.adjust).
rnaseq_build_ora_dot <- function(up_df, down_df, title = "", top_n = 10,
                                 padj_max = 0.05, is_dark = FALSE) {
  text_col <- if (is_dark) "#E0E0E0" else "black"
  line_col <- if (is_dark) "#AAAAAA" else "grey40"
  bg_col   <- if (is_dark) "#222222" else "white"
  empty_plot <- function(msg) {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg, color = text_col, size = 4.5) +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = bg_col, color = NA))
  }
  pick <- function(df, dir, strict) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- df[!is.na(df$p.adjust), , drop = FALSE]
    if (isTRUE(strict)) df <- df[df$p.adjust <= padj_max, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df <- utils::head(df[order(df$p.adjust), , drop = FALSE], top_n)
    df$Direction <- dir
    df$logP <- if (dir == "up") -log10(df$p.adjust) else log10(df$p.adjust)
    df[, c("Description", "Direction", "logP", "GeneRatioNum", "p.adjust")]
  }
  # Strict only: show exactly the terms at/below the p.adjust cutoff (top_n per
  # direction), so the plot always matches the table. No "best available" fallback -
  # drawing terms the table has filtered out is what made the two disagree.
  dd <- rbind(pick(up_df, "up", TRUE), pick(down_df, "down", TRUE))
  if (is.null(dd) || nrow(dd) == 0)
    return(empty_plot(paste0("No GO terms at p.adjust <= ", padj_max, " for this set")))
  dd$Description <- rnaseq_wrap_label(dd$Description)

  # Faint 0 reference (up right / down left) rather than a heavy divider splitting
  # the panel in two.
  zero_col <- if (is_dark) "#555555" else "grey85"
  ggplot2::ggplot(dd, ggplot2::aes(x = stats::reorder(Description, logP), y = logP,
                                   color = Direction, size = GeneRatioNum)) +
    ggplot2::geom_hline(yintercept = 0, color = zero_col, linewidth = 0.25) +
    ggplot2::geom_point() +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = c("up" = "#D55E00", "down" = "#1f78b4")) +
    ggplot2::scale_size_continuous(name = "Gene ratio", range = c(2, 8)) +
    ggplot2::labs(x = NULL, y = "Signed -log10(p.adjust)", title = title) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = bg_col, color = NA),
      panel.background  = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.background = ggplot2::element_rect(fill = bg_col, color = NA),
      legend.key        = ggplot2::element_rect(fill = bg_col, color = NA),
      plot.title  = ggplot2::element_text(color = text_col, size = 14),
      axis.title  = ggplot2::element_text(color = text_col, size = 14),
      axis.text   = ggplot2::element_text(color = text_col, size = 12),
      axis.line   = ggplot2::element_line(color = line_col),
      axis.ticks  = ggplot2::element_line(color = line_col),
      legend.text  = ggplot2::element_text(color = text_col),
      legend.title = ggplot2::element_text(color = text_col)
    )
}

rnaseq_dea_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    library(ggplot2)
    library(ggrepel)

    # Lazy: load the DEA / GSEA / ORA bundles on first render (when the RNA-seq
    # differential-expression tab opens), not at session startup. Each accessor is
    # memoized, so each bundle loads at most once per process.
    delayedAssign("dea",  rnaseq_dea_data())
    delayedAssign("gsea", rnaseq_gsea_data())
    delayedAssign("ora",  rnaseq_ora_data())
    root_input <- session$userData$root_input %||% NULL

    current_tbl <- reactive({
      validate(need(!is.null(dea),
        paste0("RNA-seq DEA data not found. Run Scripts/preprocess_rnaseq_dea.R ",
               "to generate Data/RNAseq/rnaseq_dea.rds.")))
      key <- input$contrast %||% names(dea$contrasts)[1]
      ct <- dea$contrasts[[key]]
      validate(need(!is.null(ct), "Unknown contrast"))
      ct
    })

    is_dark <- reactive(!is.null(root_input) && isTRUE(root_input$mode == "dark"))

    volcano <- reactive({
      ct <- current_tbl()
      rnaseq_build_volcano(
        tbl = ct$table, label = ct$label,
        sig_fdr = input$sig_fdr %||% 0.05,
        lfc = input$lfc %||% 0,
        top_n = input$top_n %||% 15,
        show_ns = isTRUE(input$show_ns),
        is_dark = is_dark()
      )
    })

    output$plot_container <- renderUI({
      h <- input$plot_height %||% 520
      plotOutput(session$ns("volcano"), height = paste0(h, "px"))
    })
    # Build only when this tab's Refresh is clicked (and once on load): changing
    # parameters no longer auto-redraws - the user clicks Refresh to apply them.
    # eventReactive caches the built plot, so a height resize redraws the cached
    # plot rather than silently picking up un-refreshed parameter changes.
    volcano_g <- eventReactive(input$refresh, volcano(), ignoreNULL = FALSE)
    output$volcano <- renderPlot({ volcano_g() }, bg = "transparent") |>
      bindCache(volcano_g(), plot_theme_key(session))

    output$dea_caption <- renderText({
      ct <- current_tbl()
      de <- rnaseq_classify_de(ct$table, input$sig_fdr %||% 0.05, input$lfc %||% 0)
      paste0(ct$description, "  |  ", sum(de == "Up"), " up / ", sum(de == "Down"),
             " down at FDR < ", input$sig_fdr %||% 0.05,
             if (!is.null(ct$n_samples)) paste0("  |  n = ", ct$n_samples, " samples") else "")
    })

    # Results table
    output$dea_table <- DT::renderDataTable({
      ct <- current_tbl()
      df <- ct$table
      num <- intersect(c("logFC", "AveExpr", "t", "B"), names(df))
      for (cc in num) df[[cc]] <- round(df[[cc]], 3)
      for (cc in intersect(c("P.Value", "adj.P.Val"), names(df))) df[[cc]] <- signif(df[[cc]], 3)
      DT::datatable(df, rownames = FALSE, filter = "top", selection = "none",
                    options = list(pageLength = 15, scrollX = TRUE,
                                   order = list(list(which(names(df) == "adj.P.Val") - 1, "asc"))))
    })

    # DEG overview across all contrasts
    output$deg_bar <- renderPlot({
      input$refresh   # only re-render on this tab's Refresh (isolate the rest)
      isolate({
        validate(need(!is.null(dea), ""))
        rnaseq_build_deg_bar(
          rnaseq_deg_counts(dea$contrasts, input$sig_fdr %||% 0.05, input$lfc %||% 0),
          is_dark = is_dark())
      })
    }, bg = "transparent")

    # --- GO GSEA for the selected contrast (optional layer) ---
    gsea_res <- reactive({
      validate(need(!is.null(gsea),
        paste0("GSEA data not found. Run Scripts/preprocess_rnaseq_gsea.R ",
               "to generate Data/RNAseq/rnaseq_gsea.rds.")))
      key <- input$contrast %||% names(dea$contrasts)[1]
      g <- gsea$gsea[[key]]
      validate(need(!is.null(g), "No GSEA results available for this contrast."))
      g
    })

    output$gsea_dot <- renderPlot({
      input$refresh_gsea   # only re-render on this panel's Refresh (isolate the rest)
      isolate({
        g <- gsea_res()
        rnaseq_build_gsea_dot(g$result, label = g$label,
                              top_n = input$gsea_top_n %||% 10,
                              padj_max = input$gsea_padj %||% 0.05,
                              is_dark = is_dark())
      })
    }, bg = "transparent")

    output$gsea_table <- DT::renderDataTable({
      g <- gsea_res()
      df <- g$result
      df <- df[is.na(df$p.adjust) | df$p.adjust <= (input$gsea_padj %||% 0.05), , drop = FALSE]
      if ("core_enrichment" %in% names(df)) df$core_enrichment <- NULL  # too wide for display
      for (cc in intersect(c("NES", "GeneRatio"), names(df))) df[[cc]] <- round(df[[cc]], 3)
      for (cc in intersect(c("pvalue", "p.adjust"), names(df))) df[[cc]] <- signif(df[[cc]], 3)
      DT::datatable(df, rownames = FALSE, filter = "top", selection = "none",
                    options = list(pageLength = 10, scrollX = TRUE,
                                   order = list(list(which(names(df) == "p.adjust") - 1, "asc"))))
    })

    # --- Set-specific ORA (pick two model coefficients) ---
    ora_pair <- reactive({
      validate(need(!is.null(ora),
        paste0("Set-specific ORA data not found. Run Scripts/preprocess_rnaseq_ora.R ",
               "to generate Data/RNAseq/rnaseq_ora.rds.")))
      a <- input$ora_coef_a %||% names(ora$coef_choices)[1]
      b <- input$ora_coef_b %||% names(ora$coef_choices)[2]
      validate(need(a != b, "Pick two different coefficients."))
      key <- paste(sort(c(a, b)), collapse = "|")
      pr <- ora$pairs[[key]]
      validate(need(!is.null(pr), "No ORA available for this pair."))
      list(pair = pr, a = a, b = b)
    })

    # Map the chosen segment to the stored (canonical) segment.
    ora_seg <- reactive({
      pp <- ora_pair(); pr <- pp$pair
      seg <- input$ora_segment %||% "a"
      stored <- if (seg == "common") "common"
                else if (seg == "a") (if (pp$a == pr$a) "a_only" else "b_only")
                else (if (pp$b == pr$a) "a_only" else "b_only")
      lab_a <- pr$a_label; lab_b <- pr$b_label
      focus <- switch(seg,
                      a = paste0(ora_label(ora, pp$a), " specific"),
                      b = paste0(ora_label(ora, pp$b), " specific"),
                      common = "Shared (both)")
      list(up = pr$seg[[stored]]$up, down = pr$seg[[stored]]$down, focus = focus)
    })

    # Auto-tune the ORA p.adjust cutoff so a plot always appears for the chosen
    # coefficients / segment: when they change, relax the threshold to include the
    # most significant terms (expanding the slider's max if even the best term sits
    # above the default 0.25 ceiling). Fires only on coefficient/segment change, so
    # the user can still fine-tune the slider afterwards.
    observeEvent(list(input$ora_coef_a, input$ora_coef_b, input$ora_segment), {
      s <- tryCatch(ora_seg(), error = function(e) NULL)
      if (is.null(s)) return()
      padj <- c(s$up$p.adjust, s$down$p.adjust)
      padj <- padj[is.finite(padj)]
      if (!length(padj)) return()                      # nothing to show at any cutoff
      target <- max(sort(padj)[min(length(padj), 10L)], 0.01)  # show up to ~10 top terms
      slider_max <- max(0.25, ceiling(target * 100) / 100)
      updateSliderInput(session, "ora_padj", max = slider_max, value = target)
    })

    # Number of terms the plot will draw (top_n per direction, at/below the cutoff),
    # refresh-gated so it tracks the drawn plot. Drives the plot height below.
    ora_n_shown <- eventReactive(input$refresh_ora, {
      s <- tryCatch(ora_seg(), error = function(e) NULL)
      if (is.null(s)) return(0L)
      pm <- input$ora_padj %||% 0.05; tn <- input$ora_top_n %||% 10
      nb <- function(df) if (!is.null(df) && nrow(df)) min(sum(!is.na(df$p.adjust) & df$p.adjust <= pm), tn) else 0L
      nb(s$up) + nb(s$down)
    }, ignoreNULL = FALSE)

    # Height grows with the term count so each 2-line GO label has room (no overlap).
    output$ora_plot_container <- renderUI({
      h <- max(320, min(1600, 90 + ora_n_shown() * 34))
      plotOutput(session$ns("ora_dot"), height = paste0(h, "px"))
    })

    output$ora_dot <- renderPlot({
      input$refresh_ora   # only re-render on this panel's Refresh (isolate the rest)
      isolate({
        s <- ora_seg()
        rnaseq_build_ora_dot(s$up, s$down, title = s$focus,
                             top_n = input$ora_top_n %||% 10,
                             padj_max = input$ora_padj %||% 0.05,
                             is_dark = is_dark())
      })
    }, bg = "transparent")

    output$ora_caption <- renderText({
      pp <- ora_pair(); pr <- pp$pair
      cu <- pr$counts$up; cd <- pr$counts$down
      paste0(pr$a_label, " specific: ", cu["a_only"], " up / ", cd["a_only"], " down   |   ",
             pr$b_label, " specific: ", cu["b_only"], " up / ", cd["b_only"], " down   |   ",
             "Shared: ", cu["common"], " up / ", cd["common"], " down",
             "   (FDR < ", ora$fdr, ")")
    })

    # Table shares the plot's Refresh + p.adjust cutoff, so the two never disagree
    # (previously the table updated live while the refresh-gated plot lagged, and the
    # plot's fallback showed terms the table had filtered out).
    ora_table_data <- eventReactive(input$refresh_ora, {
      s <- ora_seg()
      up <- s$up; dn <- s$down
      if (nrow(up)) up$Direction <- "up"
      if (nrow(dn)) dn$Direction <- "down"
      df <- rbind(up, dn)
      df <- df[!is.na(df$p.adjust) & df$p.adjust <= (input$ora_padj %||% 0.05), , drop = FALSE]
      if ("geneID" %in% names(df)) df$geneID <- NULL
      for (cc in intersect(c("GeneRatioNum"), names(df))) df[[cc]] <- round(df[[cc]], 3)
      for (cc in intersect(c("pvalue", "p.adjust"), names(df))) df[[cc]] <- signif(df[[cc]], 3)
      df[order(df$p.adjust), , drop = FALSE]
    }, ignoreNULL = FALSE)

    output$ora_table <- DT::renderDataTable({
      df <- ora_table_data()
      validate(need(nrow(df) > 0, "No GO terms at or below this p.adjust. Raise the threshold, then click Refresh."))
      DT::datatable(df, rownames = FALSE, filter = "top", selection = "none",
                    options = list(pageLength = 10, scrollX = TRUE,
                                   order = list(list(which(names(df) == "p.adjust") - 1, "asc"))))
    })

    # Downloads
    output$download_plot <- downloadHandler(
      filename = function() paste0("rnaseq_volcano_", input$contrast %||% "contrast", "_", Sys.Date(), ".pdf"),
      content = function(file) {
        p <- tryCatch(volcano(), error = function(e) NULL)
        pdf(file, width = 7, height = 6); on.exit(dev.off(), add = TRUE)
        if (is.null(p)) { grid::grid.newpage(); grid::grid.text("No plot to render.") } else print(p)
      }
    )
    output$download_table <- downloadHandler(
      filename = function() paste0("rnaseq_DEA_", input$contrast %||% "contrast", "_", Sys.Date(), ".csv"),
      content = function(file) {
        ct <- tryCatch(current_tbl(), error = function(e) NULL)
        if (is.null(ct)) write.csv(data.frame(), file, row.names = FALSE)
        else write.csv(ct$table, file, row.names = FALSE)
      }
    )
    output$download_gsea_plot <- downloadHandler(
      filename = function() paste0("rnaseq_GSEA_", input$contrast %||% "contrast", "_", Sys.Date(), ".pdf"),
      content = function(file) {
        g  <- tryCatch(gsea_res(), error = function(e) NULL)
        tn <- input$gsea_top_n %||% 10
        n  <- if (is.null(g) || is.null(g$result)) 0 else {
          r <- g$result
          r <- r[!is.na(r$p.adjust) & r$p.adjust <= (input$gsea_padj %||% 0.05) & !is.na(r$NES), , drop = FALSE]
          min(sum(r$NES > 0), tn) + min(sum(r$NES < 0), tn)
        }
        # Larger page, and height that grows with the term count (matches the screen).
        pdf(file, width = 11, height = max(7, min(24, 2.5 + n * 0.45))); on.exit(dev.off(), add = TRUE)
        if (is.null(g)) { grid::grid.newpage(); grid::grid.text("No GSEA results.") }
        else print(rnaseq_build_gsea_dot(g$result, label = g$label,
                                         top_n = tn,
                                         padj_max = input$gsea_padj %||% 0.05))
      }
    )
    output$download_gsea_table <- downloadHandler(
      filename = function() paste0("rnaseq_GSEA_", input$contrast %||% "contrast", "_", Sys.Date(), ".csv"),
      content = function(file) {
        g <- tryCatch(gsea_res(), error = function(e) NULL)
        if (is.null(g)) write.csv(data.frame(), file, row.names = FALSE)
        else write.csv(g$result, file, row.names = FALSE)
      }
    )
    output$download_ora_plot <- downloadHandler(
      filename = function() paste0("rnaseq_ORA_", input$ora_coef_a %||% "a", "_vs_", input$ora_coef_b %||% "b", "_", Sys.Date(), ".pdf"),
      content = function(file) {
        s  <- tryCatch(ora_seg(), error = function(e) NULL)
        tn <- input$ora_top_n %||% 10; pm <- input$ora_padj %||% 0.05
        nb <- function(df) if (!is.null(df) && nrow(df)) min(sum(!is.na(df$p.adjust) & df$p.adjust <= pm), tn) else 0
        n  <- if (is.null(s)) 0 else nb(s$up) + nb(s$down)
        # Larger page, and height that grows with the term count (matches the screen).
        pdf(file, width = 11, height = max(7, min(24, 2.5 + n * 0.45))); on.exit(dev.off(), add = TRUE)
        if (is.null(s)) { grid::grid.newpage(); grid::grid.text("No ORA results.") }
        else print(rnaseq_build_ora_dot(s$up, s$down, title = s$focus,
                                        top_n = tn,
                                        padj_max = pm))
      }
    )
    output$download_ora_table <- downloadHandler(
      filename = function() paste0("rnaseq_ORA_", input$ora_coef_a %||% "a", "_vs_", input$ora_coef_b %||% "b", "_", Sys.Date(), ".csv"),
      content = function(file) {
        s <- tryCatch(ora_seg(), error = function(e) NULL)
        if (is.null(s)) { write.csv(data.frame(), file, row.names = FALSE); return(invisible()) }
        up <- s$up; dn <- s$down
        if (nrow(up)) up$Direction <- "up"
        if (nrow(dn)) dn$Direction <- "down"
        write.csv(rbind(up, dn), file, row.names = FALSE)
      }
    )
  })
}
