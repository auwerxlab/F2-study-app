# Reference QTL data is read once and shared across all user sessions
# (qtl_data_dir is defined in app.R).
.qtl_ref <- new.env(parent = emptyenv())
qtl_reference <- function() {
  if (is.null(.qtl_ref$chrom_df)) {
    .qtl_ref$chrom_df   <- readRDS(file.path(qtl_data_dir, "chr_GRCm38.rds"))
    .qtl_ref$male_raw   <- readRDS(file.path(qtl_data_dir, "male_005_001.rds"))
    .qtl_ref$female_raw <- readRDS(file.path(qtl_data_dir, "female_005_001.rds"))
  }
  .qtl_ref
}

qtl_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    library(ggplot2)
    library(dplyr)
    library(ggrepel)
    library(DT)

    # Reference data, shared across sessions. Bound as promises (delayedAssign) so
    # the ~9 MB of QTL RDS that qtl_reference() decompresses is loaded on the first
    # Manhattan render - when the user opens the QTL tab - rather than blocking
    # every session at startup. The derived frames are promises too: reading any of
    # them (male_raw, rect_df, x_min, ...) inside a reactive forces the load once.
    delayedAssign("ref",        qtl_reference())
    delayedAssign("chrom_df",   ref$chrom_df)
    delayedAssign("male_raw",   ref$male_raw)
    delayedAssign("female_raw", ref$female_raw)
    delayedAssign("rect_df", {
      rd <- chrom_df
      rd$shade <- seq_len(nrow(chrom_df)) %% 2 == 0
      rd$xmin  <- chrom_df$cum_start
      rd$xmax  <- chrom_df$cum_start + as.numeric(chrom_df$chr_length)
      rd
    })
    delayedAssign("x_min", min(rect_df$xmin))
    delayedAssign("x_max", max(rect_df$xmax))

    # --- Reactive: filtered data ---
    filtered_data <- reactive({
      sex <- input$sex_select
      req(sex)
      
      if (sex == "male") {
        df <- male_raw
        df$sex_label <- "Male"
      } else if (sex == "female") {
        df <- female_raw
        df$sex_label <- "Female"
      } else {
        m <- male_raw;   m$sex_label <- "Male"
        f <- female_raw; f$sex_label <- "Female"
        df <- rbind(m, f)
      }
      
      # LOD filter (biggest performance win - drops most rows)
      lod_min <- input$lod_min %||% 2
      df <- df[df$lod >= lod_min, ]
      
      # Phenotype filter
      sel_phenos <- input$pheno_select
      if (!is.null(sel_phenos) && length(sel_phenos) > 0) {
        df <- df[df$phenotype %in% sel_phenos, ]
      }
      
      # P-value filter
      show_all <- isTRUE(input$show_all_points)
      if (!show_all) {
        pval_max <- input$pval_max %||% 0.025
        if (pval_max <= 0.01) {
          df <- df[df$sig_cat == "P<=0.01", ]
        } else {
          df <- df[df$sig_cat %in% c("P<=0.01", "0.01<=P<=0.05"), ]
        }
      }
      
      df
    })
    
    # --- Significant peaks table data ---
    peak_data <- reactive({
      sex <- input$sex_select
      req(sex)
      
      if (sex == "male") {
        df <- male_raw
        df$sex_label <- "Male"
      } else if (sex == "female") {
        df <- female_raw
        df$sex_label <- "Female"
      } else {
        m <- male_raw;   m$sex_label <- "Male"
        f <- female_raw; f$sex_label <- "Female"
        df <- rbind(m, f)
      }
      
      # Significance filter for table
      show_ns  <- isTRUE(input$show_ns)
      pval_max <- input$pval_max %||% 0.025
      keep <- character(0)
      if (pval_max <= 0.01) {
        keep <- "P<=0.01"
      } else {
        keep <- c("P<=0.01", "0.01<=P<=0.05")
      }
      if (show_ns) {
        keep <- c(keep, "NS")
      }
      df <- df[df$sig_cat %in% keep, ]
      
      # Phenotype filter
      sel_phenos <- input$pheno_select
      if (!is.null(sel_phenos) && length(sel_phenos) > 0) {
        df <- df[df$phenotype %in% sel_phenos, ]
      }
      
      if (nrow(df) == 0) return(data.frame())
      
      # Top peak per phenotype per chromosome - use base R aggregate
      split_key <- paste(df$phenotype, df$chr, df$sig_cat, sep = "|||")
      peaks_list <- lapply(split(df, split_key), function(chunk) {
        idx <- which.max(chunk$lod)
        data.frame(
          phenotype     = chunk$phenotype[1],
          chr           = chunk$chr[1],
          peak_pos      = chunk$pos[idx],
          peak_lod      = chunk$lod[idx],
          lod_threshold = chunk$lod_threshold[1],
          significance  = chunk$sig_cat[1],
          n_markers     = nrow(chunk),
          stringsAsFactors = FALSE
        )
      })
      peaks <- do.call(rbind, peaks_list)
      rownames(peaks) <- NULL
      
      if (sex == "both") {
        # Add sex column
        sex_list <- lapply(split(df, split_key), function(chunk) {
          chunk$sex_label[which.max(chunk$lod)]
        })
        peaks$sex <- unlist(sex_list)
        peaks <- peaks[, c("sex", "phenotype", "chr", "peak_pos", "peak_lod",
                           "lod_threshold", "significance", "n_markers")]
      }
      
      peaks <- peaks[order(-peaks$peak_lod), ]
      peaks
    })
    
    # --- Reactive: highlighted points from table selection ---
    highlight_data <- reactive({
      sel_rows <- input$peak_table_rows_selected
      if (is.null(sel_rows) || length(sel_rows) == 0) return(NULL)
      
      pd <- peak_data()
      if (nrow(pd) == 0) return(NULL)
      
      sel <- pd[sel_rows, , drop = FALSE]
      df  <- filtered_data()
      
      # Match all points belonging to selected phenotype+chr combos
      sel_keys <- paste(sel$phenotype, sel$chr, sep = "|||")
      df_keys  <- paste(df$phenotype, df$chr, sep = "|||")
      df[df_keys %in% sel_keys, ]
    })
    
    # --- Build the Manhattan plot ---
    # Access dark mode from root session
    root_input <- session$userData$root_input %||% NULL

    manhattan_plot <- reactive({
      df <- filtered_data()
      req(nrow(df) > 0)

      # Detect dark mode
      is_dark <- FALSE
      if (!is.null(root_input)) {
        is_dark <- isTRUE(root_input$mode == "dark")
      }

      bg_col   <- if (is_dark) "#222222" else "white"
      text_col <- if (is_dark) "#E0E0E0" else "black"
      line_col <- if (is_dark) "#AAAAAA" else "black"
      shade_col <- if (is_dark) "white" else "grey"

      sex <- input$sex_select
      
      # Peak labels: top significant hit per phenotype per chromosome
      sig_df <- df[df$sig_cat != "NS", ]
      max_df_red <- data.frame()
      if (nrow(sig_df) > 0) {
        max_df <- sig_df %>%
          group_by(phenotype, chr) %>%
          filter(lod == max(lod)) %>%
          slice(1) %>%
          ungroup()
        max_df_red <- max_df[max_df$sig_cat == "P<=0.01", ]
      }
      
      dummy_df <- data.frame(x = c(x_min, x_max), y = c(0, 0))
      
      p <- ggplot() +
        geom_rect(data = subset(rect_df, shade),
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                  fill = shade_col, alpha = 0.15, inherit.aes = FALSE) +
        geom_blank(data = dummy_df, aes(x = x, y = y))
      
      if (sex == "both") {
        p <- p +
          geom_point(data = df,
                     aes(x = global_pos, y = lod, color = sig_cat, shape = sex_label),
                     size = 1, alpha = 0.6, na.rm = TRUE) +
          scale_shape_manual(name = "Sex",
                             values = c("Male" = 16, "Female" = 17))
      } else {
        p <- p +
          geom_point(data = df,
                     aes(x = global_pos, y = lod, color = sig_cat),
                     size = 1, alpha = 0.7, na.rm = TRUE)
      }
      
      p <- p +
        scale_color_manual(name = "P value",
                           values = c("NS" = "#7F7F7F",
                                      "0.01<=P<=0.05" = "#1f78b4",
                                      "P<=0.01" = "#D55E00"),
                           breaks = c("0.01<=P<=0.05", "P<=0.01")) +
        guides(color = guide_legend(override.aes = list(size = 4)))
      
      if (nrow(max_df_red) > 0) {
        p <- p +
          geom_text_repel(data = max_df_red,
                          aes(x = global_pos, y = lod, label = phenotype),
                          size = 3.5, color = text_col,
                          nudge_y = 1, force = 3, box.padding = 0.5,
                          na.rm = TRUE, max.overlaps = 20,
                          show.legend = FALSE)
      }
      
      # Highlight selected table rows
      hl <- highlight_data()
      if (!is.null(hl) && nrow(hl) > 0) {
        p <- p +
          geom_point(data = hl,
                     aes(x = global_pos, y = lod),
                     color = "#FFD700", size = 2.5, alpha = 0.8,
                     shape = 21, fill = "#FFD700", stroke = 0.5,
                     na.rm = TRUE, show.legend = FALSE) +
          geom_text_repel(data = hl[!duplicated(paste(hl$phenotype, hl$chr)), ],
                          aes(x = global_pos, y = lod, label = paste0(phenotype, " (chr", chr, ")")),
                          size = 3, color = "#B8860B", fontface = "bold",
                          nudge_y = 0.5, force = 2, box.padding = 0.4,
                          na.rm = TRUE, max.overlaps = 30,
                          show.legend = FALSE)
      }
      
      lod_min <- input$lod_min %||% 2
      
      p <- p +
        scale_x_continuous(breaks = chrom_df$cum_mid,
                           labels = chrom_df$chr,
                           limits = c(x_min, x_max),
                           expand = c(0, 0)) +
        scale_y_continuous(limits = c(lod_min, NA),
                           expand = expansion(mult = c(0, 0.05))) +
        xlab("Chromosome") +
        ylab("LOD Score") +
        theme_classic() +
        theme(
          plot.background   = element_rect(fill = bg_col, color = NA),
          panel.background  = element_rect(fill = bg_col, color = NA),
          axis.title        = element_text(size = 14, color = text_col),
          axis.text         = element_text(size = 12, color = text_col),
          axis.line         = element_line(color = line_col),
          axis.ticks        = element_line(color = line_col),
          legend.title      = element_text(size = 12, color = text_col),
          legend.text       = element_text(size = 10, color = text_col),
          legend.background = element_rect(fill = bg_col, color = NA),
          legend.key        = element_rect(fill = bg_col, color = NA)
        )
      
      p
    })
    
    # --- Dynamic plot height container ---
    output$plot_container <- renderUI({
      h <- input$plot_height %||% 500
      plotOutput(session$ns("manhattan_plot"), height = paste0(h, "px"))
    })
    
    # --- Render the plot ---
    manhattan_g <- eventReactive(input$refresh, manhattan_plot(), ignoreNULL = FALSE)
    output$manhattan_plot <- renderPlot({ manhattan_g() }, bg = "transparent") |>
      bindCache(manhattan_g(), plot_theme_key(session))
    
    # --- Render the peaks table ---
    output$peak_table <- DT::renderDataTable({
      pd <- peak_data()
      req(nrow(pd) > 0)
      
      # Format numbers for display
      pd$peak_lod      <- round(pd$peak_lod, 3)
      pd$lod_threshold <- round(pd$lod_threshold, 3)
      pd$peak_pos      <- format(pd$peak_pos, big.mark = ",")
      
      DT::datatable(pd,
                    rownames = FALSE,
                    filter = "top",
                    selection = "multiple",
                    options = list(
                      pageLength = 15,
                      scrollX = TRUE,
                      order = list(list(which(names(pd) == "peak_lod") - 1, "desc"))
                    ))
    })

    output$download_peaks <- downloadHandler(
      filename = function() {
        paste0("qtl_peaks_", input$sex_select %||% "male", "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        d <- tryCatch(peak_data(), error = function(e) data.frame())
        write.csv(d, file, row.names = FALSE)
      }
    )

    # --- Download handler ---
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("manhattan_qtl_", input$sex_select, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        p <- manhattan_plot()
        if (is.null(p)) {
          pdf(file, width = 27, height = 10)
          plot.new()
          text(0.5, 0.5, "No data to plot", cex = 2)
          dev.off()
        } else {
          pdf(file, width = 27, height = 10)
          print(p)
          dev.off()
        }
      }
    )

    # ===================================================================
    # --- "By Chromosome" tab: one phenotype, one panel per chromosome ---
    # ===================================================================

    # Full LOD profile for the selected phenotype (the whole curve is shown - no
    # LOD / p-value filtering). Reading male_raw / female_raw here forces the lazy
    # QTL reference load on the first render of this view.
    profile_data <- reactive({
      ph <- input$prof_pheno
      req(ph)
      sex <- input$prof_sex %||% "male"

      pick <- function(df, lab) {
        d <- df[df$phenotype == ph, , drop = FALSE]
        d$sex_label <- rep(lab, nrow(d))   # length-safe: 0 rows -> empty column (no rbind mismatch)
        d
      }
      df <- switch(sex,
                   male   = pick(male_raw,   "Male"),
                   female = pick(female_raw, "Female"),
                   both   = rbind(pick(male_raw, "Male"), pick(female_raw, "Female")),
                   pick(male_raw, "Male"))

      req(nrow(df) > 0)
      df$chr    <- factor(df$chr, levels = as.character(1:19))
      df$pos_mb <- df$pos / 1e6
      df[order(df$chr, df$pos_mb), , drop = FALSE]
    })

    # Default the significance-LOD slider to the selected phenotype's own genome-wide
    # threshold whenever the phenotype / sex changes, so the dashed line starts in the
    # right place (the user can still move it). Reading *_raw here forces the lazy load.
    observeEvent(list(input$prof_pheno, input$prof_sex), {
      ph <- input$prof_pheno
      if (is.null(ph)) return()
      sex <- input$prof_sex %||% "male"
      get_thr <- function(df) {
        t <- unique(df$lod_threshold[df$phenotype == ph]); t <- t[is.finite(t)]
        if (length(t)) t else numeric(0)
      }
      thr <- tryCatch(switch(sex,
                             male   = get_thr(male_raw),
                             female = get_thr(female_raw),
                             both   = c(get_thr(male_raw), get_thr(female_raw)),
                             get_thr(male_raw)),
                      error = function(e) numeric(0))
      thr <- thr[is.finite(thr)]
      if (length(thr)) updateSliderInput(session, "prof_lod_min", value = round(mean(thr), 2))
    })

    # Faceted LOD profile: one panel per chromosome, 4 panels per row.
    profile_plot <- reactive({
      df  <- profile_data()
      sex <- input$prof_sex %||% "male"

      # Significance-line position (slider). Do NOT drop rows here: the full set
      # defines the facets and per-chromosome axis ranges. When "hide below" is on we
      # only omit below-line markers from the DRAWN points, so every chromosome panel
      # still shows (dropping rows collapsed the grid to a single panel before).
      sig  <- input$prof_lod_min %||% 0
      hide <- isTRUE(input$prof_hide_below) && is.finite(sig) && sig > 0

      is_dark  <- !is.null(root_input) && isTRUE(root_input$mode == "dark")
      bg_col   <- if (is_dark) "#222222" else "white"
      text_col <- if (is_dark) "#E0E0E0" else "black"
      line_col <- if (is_dark) "#AAAAAA" else "grey40"

      # Colour relative to the movable line: markers below it are grey (NS); above it
      # they keep the P-value tiers (blue = 0.01<=P<=0.05, orange = P<=0.01). So moving
      # the slider recolours points while the blue tier stays visible.
      df$col_cat <- ifelse(is.na(df$lod) | df$lod < sig, "NS",
                           ifelse(df$sig_cat == "P<=0.01", "P<=0.01", "0.01<=P<=0.05"))
      df$col_cat <- factor(df$col_cat, levels = c("NS", "0.01<=P<=0.05", "P<=0.01"))

      # Points actually drawn (all, or only at/above the line when hiding).
      df_pts <- if (hide) df[!is.na(df$lod) & df$lod >= sig, , drop = FALSE] else df

      # geom_blank() on the full data fixes the facets + axis ranges for every
      # chromosome, so hiding points never collapses the grid to a single panel.
      p <- ggplot(df, aes(x = pos_mb, y = lod)) + geom_blank()

      # Colour always encodes the P-value tier. When overlaying both sexes, tell them
      # apart by point SHAPE (circle = male, triangle = female) rather than colour, so
      # the P-value colours are kept in both single-sex and overlay modes.
      if (sex == "both") {
        p <- p +
          geom_point(data = df_pts, aes(color = col_cat, shape = sex_label), size = 1.1, na.rm = TRUE) +
          scale_shape_manual(name = "Sex", values = c("Male" = 16, "Female" = 17))
      } else {
        p <- p +
          geom_point(data = df_pts, aes(color = col_cat), size = 0.9, na.rm = TRUE)
      }
      p <- p +
        scale_color_manual(name = "P value",
                           values = c("NS" = "#B0B0B0", "0.01<=P<=0.05" = "#1f78b4", "P<=0.01" = "#D55E00"),
                           breaks = c("NS", "0.01<=P<=0.05", "P<=0.01"), drop = FALSE)

      # Significance line at the slider position (defaults to the phenotype's own LOD
      # threshold; see the observeEvent that syncs the slider on phenotype change).
      if (is.finite(sig) && sig > 0) {
        p <- p + geom_hline(yintercept = sig, linetype = "dashed",
                            color = line_col, linewidth = 0.5)
      }

      p +
        facet_wrap(~ chr, ncol = 4, scales = "free_x", drop = FALSE) +
        expand_limits(x = 0, y = 0) +
        labs(x = "Position (Mb)", y = "LOD score",
             title = paste0(input$prof_pheno, " - LOD profile by chromosome")) +
        theme_bw(base_size = 11) +
        theme(
          plot.background   = element_rect(fill = bg_col, color = NA),
          panel.background  = element_rect(fill = bg_col, color = NA),
          legend.background = element_rect(fill = bg_col, color = NA),
          legend.key        = element_rect(fill = bg_col, color = NA),
          legend.position   = "bottom",
          strip.background  = element_rect(fill = if (is_dark) "#333333" else "grey90", color = NA),
          strip.text        = element_text(color = text_col, face = "bold", size = 13),
          plot.title        = element_text(color = text_col, size = 14),
          axis.title        = element_text(color = text_col),
          axis.text         = element_text(color = text_col, size = 8),
          axis.line         = element_line(color = line_col),
          axis.ticks        = element_line(color = line_col),
          panel.grid.minor  = element_blank(),
          legend.title      = element_text(color = text_col),
          legend.text       = element_text(color = text_col)
        )
    })

    output$profile_caption <- renderText({
      df <- tryCatch(profile_data(), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) return("")
      thr <- unique(df$lod_threshold); thr <- thr[is.finite(thr)]
      paste0(input$prof_pheno, "  |  ", nrow(df), " markers across ",
             length(unique(df$chr[!is.na(df$chr)])), " chromosomes  |  significance LOD ",
             if (length(thr)) paste(round(thr, 2), collapse = " / ") else "NA",
             "  |  ", sum(df$sig_cat != "NS"), " markers above P<=0.05")
    })

    output$profile_container <- renderUI({
      h <- input$prof_height %||% 800
      plotOutput(session$ns("profile_plot"), height = paste0(h, "px"))
    })

    profile_g <- eventReactive(input$refresh_prof, profile_plot(), ignoreNULL = FALSE)
    output$profile_plot <- renderPlot({ profile_g() }, bg = "transparent") |>
      bindCache(profile_g(), plot_theme_key(session))

    # Per-chromosome peak table for the selected phenotype, gated behind Refresh like
    # the plot. Uses the SAME aggregation, columns and formatting as the "By Phenotype"
    # QTL Peaks table (peak_data) - top marker per chromosome x significance tier, with
    # n_markers = markers in that chr/tier group - so both tables report identical
    # peak_pos / peak_lod / n_markers for a given phenotype and sex. Both significant
    # tiers are kept (P<=0.01 and 0.01<=P<=0.05); ticking "Include non-significant" also
    # lists the NS tier. This keys off sig_cat (the genome-wide threshold), NOT the
    # movable slider line, so dragging the plot's dashed line never changes the table.
    profile_table_data <- eventReactive(input$refresh_prof, {
      df  <- profile_data()
      if (!nrow(df)) return(data.frame())
      sex  <- input$prof_sex %||% "male"
      both <- (sex == "both") && "sex_label" %in% names(df)

      keep <- c("P<=0.01", "0.01<=P<=0.05")
      if (isTRUE(input$prof_show_ns)) keep <- c(keep, "NS")
      df <- df[df$sig_cat %in% keep, , drop = FALSE]
      if (!nrow(df)) return(data.frame())

      df$chr <- as.character(df$chr)
      split_key <- paste(df$chr, df$sig_cat, sep = "|||")
      rows <- lapply(split(df, split_key), function(chunk) {
        idx <- which.max(chunk$lod)
        data.frame(
          phenotype     = chunk$phenotype[1],
          chr           = chunk$chr[1],
          peak_pos      = chunk$pos[idx],
          peak_lod      = chunk$lod[idx],
          lod_threshold = chunk$lod_threshold[1],
          significance  = chunk$sig_cat[1],
          n_markers     = nrow(chunk),
          stringsAsFactors = FALSE
        )
      })
      peaks <- do.call(rbind, rows); rownames(peaks) <- NULL

      if (both) {
        sex_list <- lapply(split(df, split_key), function(chunk) chunk$sex_label[which.max(chunk$lod)])
        peaks$sex <- unlist(sex_list)
        peaks <- peaks[, c("sex", "phenotype", "chr", "peak_pos", "peak_lod",
                           "lod_threshold", "significance", "n_markers")]
      }
      peaks[order(-peaks$peak_lod), , drop = FALSE]
    }, ignoreNULL = FALSE)

    output$profile_table <- DT::renderDataTable({
      d <- profile_table_data()
      validate(need(nrow(d) > 0, "No significant peaks for this phenotype (tick 'Include non-significant' to see all)."))
      # Same display formatting as the By Phenotype table.
      d$peak_lod      <- round(d$peak_lod, 3)
      d$lod_threshold <- round(d$lod_threshold, 3)
      d$peak_pos      <- format(d$peak_pos, big.mark = ",")
      DT::datatable(d, rownames = FALSE, filter = "top", selection = "none",
                    options = list(pageLength = 10, scrollX = TRUE, deferRender = TRUE,
                                   order = list(list(which(names(d) == "peak_lod") - 1, "desc"))))
    })

    output$download_profile_table <- downloadHandler(
      filename = function() {
        paste0("qtl_chr_peaks_", input$prof_pheno %||% "phenotype", "_",
               input$prof_sex %||% "male", "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        d <- tryCatch(profile_table_data(), error = function(e) data.frame())
        write.csv(d, file, row.names = FALSE)
      }
    )

    output$download_profile <- downloadHandler(
      filename = function() {
        paste0("qtl_profile_", input$prof_pheno %||% "phenotype", "_",
               input$prof_sex %||% "male", "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        p <- tryCatch(profile_plot(), error = function(e) NULL)
        pdf(file, width = 14, height = 10)
        on.exit(dev.off(), add = TRUE)
        if (is.null(p)) { plot.new(); text(0.5, 0.5, "No data to plot", cex = 2) }
        else print(p)
      }
    )
  })
}
