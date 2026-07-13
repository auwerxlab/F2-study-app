# Group LOD reference data is shared across all user sessions
# (lod_data_dir is defined in app.R).
.lod_group_cache <- new.env(parent = emptyenv())   # per-(sex, group) data
.lod_ref         <- new.env(parent = emptyenv())   # phenotype map
lod_phenotype_map <- function(rds_dir) {
  if (is.null(.lod_ref$mp)) .lod_ref$mp <- readRDS(file.path(rds_dir, "phenotype_map.rds"))
  .lod_ref$mp
}

grouplod_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    library(ggplot2)
    library(RColorBrewer)
    library(grid)

    rds_dir <- file.path(lod_data_dir, "rds")

    # --- Phenotype map (small, loaded once and shared across sessions) ---
    # Lazy: read on first use (when a group is plotted), not at session startup.
    delayedAssign("mp", lod_phenotype_map(rds_dir))

    # --- Per-group RDS cache (shared across sessions) ---
    group_cache <- .lod_group_cache

    load_group_data <- function(sex, group) {
      cache_key <- paste0(sex, "_", group)
      if (!is.null(group_cache[[cache_key]])) return(group_cache[[cache_key]])

      safe_group <- gsub("[^A-Za-z0-9_]", "_", group)
      rds_file <- file.path(rds_dir, paste0(sex, "_", safe_group, ".rds"))
      if (!file.exists(rds_file)) return(NULL)

      df <- readRDS(rds_file)
      group_cache[[cache_key]] <- df
      df
    }

    # --- Helpers ---
    is_sig <- function(status) status %in% c("Significant_005", "Significant_001")

    chr_key <- function(x) {
      suppressWarnings(as.numeric(gsub("^chr", "", x, ignore.case = TRUE)))
    }

    make_palette_no_yellow <- function(levels_vec) {
      n <- length(levels_vec)
      if (n <= 0) return(setNames(character(0), character(0)))
      if (n <= 8) {
        cols <- brewer.pal(max(3, n), "Dark2")[seq_len(n)]
        return(setNames(cols, levels_vec))
      }
      hues <- seq(15, 375, length.out = n + 1)[seq_len(n)]
      cols <- grDevices::hcl(h = hues, c = 115, l = 55)
      setNames(cols, levels_vec)
    }

    merge_intervals <- function(starts, ends) {
      ok <- is.finite(starts) & is.finite(ends) & (ends > starts)
      starts <- starts[ok]; ends <- ends[ok]
      if (length(starts) == 0) return(data.frame(start = numeric(0), end = numeric(0)))
      o <- order(starts, ends); starts <- starts[o]; ends <- ends[o]
      out_s <- c(); out_e <- c()
      cur_s <- starts[1]; cur_e <- ends[1]
      if (length(starts) >= 2) {
        for (i in 2:length(starts)) {
          if (starts[i] <= cur_e) { cur_e <- max(cur_e, ends[i]) }
          else { out_s <- c(out_s, cur_s); out_e <- c(out_e, cur_e); cur_s <- starts[i]; cur_e <- ends[i] }
        }
      }
      out_s <- c(out_s, cur_s); out_e <- c(out_e, cur_e)
      data.frame(start = out_s, end = out_e)
    }

    get_sig_runs <- function(dat_ph, region_start, region_end, min_bar_mb = 0.05) {
      dat_ph <- dat_ph[order(dat_ph$x_mb), ]
      if (nrow(dat_ph) == 0) return(data.frame(start = numeric(0), end = numeric(0)))
      sig <- is_sig(dat_ph$Status)
      if (!any(sig)) return(data.frame(start = numeric(0), end = numeric(0)))
      run_id <- cumsum(sig & !c(FALSE, sig[-length(sig)]))
      sub <- dat_ph[sig, , drop = FALSE]
      sub$run_id <- run_id[sig]
      runs <- do.call(rbind, lapply(split(sub, sub$run_id), function(chunk) {
        data.frame(start = min(chunk$x_mb), end = max(chunk$x_mb))
      }))
      runs$mid <- (runs$start + runs$end) / 2
      runs$len <- runs$end - runs$start
      runs$start <- ifelse(runs$len < min_bar_mb, runs$mid - min_bar_mb / 2, runs$start)
      runs$end   <- ifelse(runs$len < min_bar_mb, runs$mid + min_bar_mb / 2, runs$end)
      runs$start <- pmax(runs$start, region_start)
      runs$end   <- pmin(runs$end, region_end)
      runs <- runs[runs$end > runs$start, c("start", "end")]
      runs
    }

    get_overlap_regions <- function(interval_df, k) {
      if (nrow(interval_df) == 0 || k <= 1) return(merge_intervals(interval_df$start, interval_df$end))
      events <- rbind(
        data.frame(pos = interval_df$start, type = "start", pheno = interval_df$display_pheno, stringsAsFactors = FALSE),
        data.frame(pos = interval_df$end,   type = "end",   pheno = interval_df$display_pheno, stringsAsFactors = FALSE)
      )
      events <- events[order(events$pos), ]
      positions <- sort(unique(events$pos))
      if (length(positions) < 2) return(data.frame(start = numeric(0), end = numeric(0)))
      counts <- character(0); out_s <- c(); out_e <- c()
      for (i in seq_len(length(positions) - 1)) {
        pos_now <- positions[i]; pos_next <- positions[i + 1]
        ev_here <- events[events$pos == pos_now, , drop = FALSE]
        for (p in ev_here$pheno[ev_here$type == "end"]) counts <- counts[counts != p]
        for (p in ev_here$pheno[ev_here$type == "start"]) counts <- c(counts, p)
        if (length(unique(counts)) >= k) { out_s <- c(out_s, pos_now); out_e <- c(out_e, pos_next) }
      }
      if (length(out_s) == 0) return(data.frame(start = numeric(0), end = numeric(0)))
      merge_intervals(out_s, out_e)
    }

    build_segments <- function(dat) {
      dat <- dat[order(dat$Chr, dat$display_pheno, dat$x_mb), ]
      do.call(rbind, lapply(split(dat, paste(dat$Chr, dat$display_pheno)), function(chunk) {
        if (nrow(chunk) < 2) return(NULL)
        chunk$x_next  <- c(chunk$x_mb[-1], NA)
        chunk$y_next  <- c(chunk$Lod[-1], NA)
        chunk$sig_now <- is_sig(chunk$Status)
        chunk$line_group <- ifelse(chunk$sig_now, "Colored", "NS")
        chunk[is.finite(chunk$x_next) & is.finite(chunk$y_next), ]
      }))
    }

    # --- Reactive: load only the needed group data ---
    plot_data <- reactive({
      sex <- input$sex_select
      group <- input$group_select
      req(sex, group)

      if (sex == "both") {
        df_m <- load_group_data("male", group)
        df_f <- load_group_data("female", group)
        if (!is.null(df_m) && nrow(df_m) > 0) df_m$sex_label <- "Male"
        if (!is.null(df_f) && nrow(df_f) > 0) df_f$sex_label <- "Female"
        df <- rbind(df_m, df_f)
      } else {
        df <- load_group_data(sex, group)
        if (!is.null(df) && nrow(df) > 0) df$sex_label <- ifelse(sex == "male", "Male", "Female")
      }
      req(!is.null(df), nrow(df) > 0)
      df
    })

    # --- Update chromosome selector based on significant chromosomes ---
    observe({
      df <- plot_data()
      req(nrow(df) > 0)

      sig_chrs <- unique(df$Chr[is_sig(df$Status)])
      sig_chrs <- sig_chrs[order(chr_key(sig_chrs))]

      chr_choices <- c("Auto (first significant)" = "auto")
      if (length(sig_chrs) > 0) {
        named <- setNames(sig_chrs, paste("Chr", sig_chrs))
        chr_choices <- c(chr_choices, named)
      }
      all_chrs <- sort(unique(df$Chr), decreasing = FALSE)
      all_chrs <- all_chrs[order(chr_key(all_chrs))]
      all_chrs <- setdiff(all_chrs, sig_chrs)
      if (length(all_chrs) > 0) {
        named_all <- setNames(all_chrs, paste("Chr", all_chrs, "(NS)"))
        chr_choices <- c(chr_choices, named_all)
      }

      updateSelectInput(session, "chr_select", choices = chr_choices, selected = "auto")
    })

    # --- Update share_k range based on number of significant phenotypes ---
    observe({
      df <- plot_data()
      req(nrow(df) > 0)
      n_sig_pheno <- length(unique(df$display_pheno[is_sig(df$Status)]))
      max_k <- max(n_sig_pheno, 1)
      cur_k <- min(input$share_k %||% 2, max_k)
      updateNumericInput(session, "share_k", min = 1, max = max_k, value = cur_k)
      output$share_k_range <- renderText(paste0("Range: 1 - ", max_k, " phenotypes"))
    })

    # Access dark mode from root session
    root_input <- session$userData$root_input %||% NULL

    # --- Build the two-panel plot ---
    group_plot <- reactive({
      df <- plot_data()
      req(nrow(df) > 0)

      # Detect dark mode
      is_dark <- FALSE
      if (!is.null(root_input)) {
        is_dark <- isTRUE(root_input$mode == "dark")
      }
      bg_col     <- if (is_dark) "#222222" else "white"
      text_col   <- if (is_dark) "#E0E0E0" else "black"
      line_col   <- if (is_dark) "#AAAAAA" else "black"
      shared_fill <- if (is_dark) "grey50" else "grey80"

      sex_sel   <- input$sex_select
      is_both   <- (sex_sel == "both")
      group     <- input$group_select
      flank_mb  <- input$flank_mb %||% 10
      chr_sel   <- input$chr_select %||% "auto"

      share_mode_in <- input$share_mode %||% "k"
      share_k_val   <- input$share_k %||% 2

      # Determine which chromosome to plot
      df_sig <- df[is_sig(df$Status), ]
      if (nrow(df_sig) == 0) {
        return(ggplot() + theme_void() +
          geom_text(aes(x = 0, y = 0, label = "No significant QTL found for this group")))
      }

      sig_chrs <- unique(df_sig$Chr)
      sig_chrs <- sig_chrs[order(chr_key(sig_chrs))]

      if (chr_sel == "auto") {
        ch <- sig_chrs[1]
      } else {
        ch <- chr_sel
      }

      # Group-level phenotype list
      pheno_all_group_display <- ifelse(
        is.na(mp$Simply_name[mp$Group == group]) | mp$Simply_name[mp$Group == group] == "",
        mp$Raw_name[mp$Group == group],
        mp$Simply_name[mp$Group == group]
      )
      pheno_all_group_display <- unique(pheno_all_group_display)

      # Filter to selected chromosome only
      sig_chr <- df_sig[df_sig$Chr == ch, ]
      all_chr <- df[df$Chr == ch, ]
      if (nrow(sig_chr) == 0 || nrow(all_chr) == 0) {
        return(ggplot() + theme_void() +
          geom_text(aes(x = 0, y = 0, label = paste("No significant QTL on Chr", ch))))
      }

      chr_min <- min(all_chr$x_mb); chr_max <- max(all_chr$x_mb)
      region_start <- max(min(sig_chr$x_mb) - flank_mb, chr_min)
      region_end   <- min(max(sig_chr$x_mb) + flank_mb, chr_max)
      if (region_end - region_start < 0.2) {
        mid <- (region_start + region_end) / 2
        region_start <- mid - 0.1; region_end <- mid + 0.1
      }

      sig_pheno_chr <- unique(sig_chr$display_pheno)
      df_win <- all_chr[all_chr$x_mb >= region_start & all_chr$x_mb <= region_end, ]
      if (nrow(df_win) == 0) {
        return(ggplot() + theme_void() +
          geom_text(aes(x = 0, y = 0, label = "No data in region")))
      }

      # === TOP PANEL: interval bars ===
      if (is_both) {
        # Build bars per sex
        sexes_present <- unique(df_win$sex_label)
        bars_list <- list()
        for (sx in sexes_present) {
          for (p in sig_pheno_chr) {
            dat_p <- df_win[df_win$display_pheno == p & df_win$sex_label == sx, ]
            if (nrow(dat_p) == 0) next
            runs <- get_sig_runs(dat_p, region_start, region_end)
            if (nrow(runs) == 0) next
            runs$display_pheno <- p
            runs$sex_label <- sx
            bars_list[[paste(sx, p)]] <- runs
          }
        }
        bars_df <- do.call(rbind, bars_list)
      } else {
        bars_list <- lapply(sig_pheno_chr, function(p) {
          dat_p <- df_win[df_win$display_pheno == p, ]
          if (nrow(dat_p) == 0) return(NULL)
          runs <- get_sig_runs(dat_p, region_start, region_end)
          if (nrow(runs) == 0) return(NULL)
          runs$display_pheno <- p
          runs
        })
        bars_df <- do.call(rbind, bars_list)
      }
      if (is.null(bars_df) || nrow(bars_df) == 0) {
        return(ggplot() + theme_void() +
          geom_text(aes(x = 0, y = 0, label = "No significant intervals found")))
      }

      # Order: significant phenotypes on top, then rest
      sig_ordered <- pheno_all_group_display[pheno_all_group_display %in% sig_pheno_chr]
      nonsig_ordered <- pheno_all_group_display[!pheno_all_group_display %in% sig_ordered]
      y_order <- c(sig_ordered, nonsig_ordered)

      y_map <- data.frame(display_pheno = y_order, y = seq_along(y_order),
                          stringsAsFactors = FALSE)

      bars_df <- merge(bars_df, y_map, by = "display_pheno")
      bar_height <- 0.65
      if (is_both) {
        # Offset: male upper half, female lower half
        half <- bar_height / 2
        bars_df$ymin <- ifelse(bars_df$sex_label == "Male",
                               bars_df$y - half, bars_df$y)
        bars_df$ymax <- ifelse(bars_df$sex_label == "Male",
                               bars_df$y, bars_df$y + half)
      } else {
        bars_df$ymin <- bars_df$y - bar_height / 2
        bars_df$ymax <- bars_df$y + bar_height / 2
      }

      pal_all <- make_palette_no_yellow(sig_ordered)

      # Shared regions
      merged_by_pheno <- do.call(rbind, lapply(split(bars_df, bars_df$display_pheno), function(sub) {
        m <- merge_intervals(sub$start, sub$end)
        if (nrow(m) == 0) return(NULL)
        data.frame(display_pheno = sub$display_pheno[1], start = m$start, end = m$end,
                   stringsAsFactors = FALSE)
      }))

      shared_regions <- data.frame(start = numeric(0), end = numeric(0))
      if (!is.null(merged_by_pheno) && nrow(merged_by_pheno) > 0) {
        n_ph <- length(unique(merged_by_pheno$display_pheno))
        if (share_mode_in == "0") {
          shared_regions <- merge_intervals(merged_by_pheno$start, merged_by_pheno$end)
        } else if (share_mode_in == "all") {
          shared_regions <- get_overlap_regions(merged_by_pheno, k = n_ph)
        } else {
          shared_regions <- get_overlap_regions(merged_by_pheno, k = min(share_k_val, n_ph))
        }
        if (nrow(shared_regions) > 0) {
          shared_regions$start <- pmax(shared_regions$start, region_start)
          shared_regions$end   <- pmin(shared_regions$end, region_end)
          shared_regions <- shared_regions[shared_regions$end > shared_regions$start, ]
        }
      }

      subtitle_txt <- ""
      if (nrow(shared_regions) > 0) {
        labs_vec <- paste0(sprintf("%.2f", shared_regions$start), "M-",
                          sprintf("%.2f", shared_regions$end), "M")
        subtitle_txt <- paste0("Shared region: ", paste(labs_vec, collapse = "; "))
      }

      pretty_group <- gsub("_", " ", group)
      chr_label <- paste0("Chr ", gsub("^chr", "", ch, ignore.case = TRUE))
      if (is_both && nchar(subtitle_txt) > 0) {
        subtitle_txt <- paste0(subtitle_txt, "  |  solid = Male, faded = Female")
      } else if (is_both) {
        subtitle_txt <- "solid = Male, faded = Female"
      }
      br <- pretty(c(region_start, region_end), n = 5)
      br <- br[br >= region_start & br <= region_end]

      # Top plot
      p_top <- ggplot()
      if (nrow(shared_regions) > 0) {
        p_top <- p_top +
          geom_rect(data = shared_regions,
                    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                    fill = shared_fill, alpha = 0.28, inherit.aes = FALSE)
      }
      if (is_both) {
        bars_male <- bars_df[bars_df$sex_label == "Male", ]
        bars_female <- bars_df[bars_df$sex_label == "Female", ]
        if (nrow(bars_male) > 0) {
          p_top <- p_top +
            geom_rect(data = bars_male,
                      aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = display_pheno),
                      color = NA, alpha = 1)
        }
        if (nrow(bars_female) > 0) {
          p_top <- p_top +
            geom_rect(data = bars_female,
                      aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = display_pheno),
                      color = NA, alpha = 0.5, linetype = "dashed")
        }
        # Sex indicator in subtitle instead of annotation
      } else {
        p_top <- p_top +
          geom_rect(data = bars_df,
                    aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = display_pheno),
                    color = NA, alpha = 1)
      }
      p_top <- p_top +
        scale_fill_manual(values = pal_all, guide = "none") +
        scale_x_continuous(name = "Position (Mb)", position = "top",
                           limits = c(region_start, region_end), breaks = br,
                           labels = function(x) paste0(x, "M"), expand = c(0, 0)) +
        scale_y_reverse(name = NULL, breaks = y_map$y, labels = y_map$display_pheno,
                        limits = c(length(y_order) + 0.5, 0.5)) +
        labs(title = paste0(pretty_group, " QTL summary | ", chr_label),
             subtitle = subtitle_txt) +
        theme_classic() +
        theme(
          plot.title      = element_text(size = 16, hjust = 0.5, color = text_col),
          plot.subtitle   = element_text(size = 11, hjust = 0.5, color = if (is_dark) "#AAAAAA" else "grey40"),
          axis.title.x    = element_text(size = 12, color = text_col),
          axis.text.x     = element_text(size = 10, color = text_col),
          axis.text.y     = element_text(size = 9, color = text_col),
          axis.line        = element_line(color = line_col),
          axis.ticks       = element_line(color = line_col),
          panel.background = element_rect(fill = bg_col, color = NA),
          plot.background  = element_rect(fill = bg_col, color = NA),
          plot.margin      = margin(t = 5, r = 8, b = 2, l = 5)
        )

      # === BOTTOM PANEL: LOD lines ===
      df_bottom <- all_chr[all_chr$display_pheno %in% sig_pheno_chr &
                           all_chr$x_mb >= region_start & all_chr$x_mb <= region_end, ]
      df_bottom$display_pheno <- factor(df_bottom$display_pheno, levels = pheno_all_group_display)

      if (nrow(df_bottom) == 0) return(p_top)

      # Build segments - split by sex too if both mode
      if (is_both) {
        seg_list <- lapply(split(df_bottom, df_bottom$sex_label), function(sub) {
          s <- build_segments(sub)
          if (!is.null(s) && nrow(s) > 0) s$sex_label <- sub$sex_label[1]
          s
        })
        seg_df <- do.call(rbind, seg_list)
      } else {
        seg_df <- build_segments(df_bottom)
        if (!is.null(seg_df) && nrow(seg_df) > 0) seg_df$sex_label <- df_bottom$sex_label[1]
      }
      if (!is.null(seg_df)) {
        seg_df <- seg_df[seg_df$x_mb >= region_start & seg_df$x_mb <= region_end &
                         seg_df$x_next >= region_start & seg_df$x_next <= region_end, ]
      }

      pheno_present <- sig_ordered[sig_ordered %in% unique(as.character(df_bottom$display_pheno))]
      pheno_present <- pheno_present[pheno_present %in% names(pal_all)]
      if (length(pheno_present) == 0) return(p_top)

      pal_bottom <- pal_all[pheno_present]
      color_values <- c(pal_bottom, "Non-significant" = "#7F7F7F")

      p_bottom <- ggplot()
      if (nrow(shared_regions) > 0) {
        p_bottom <- p_bottom +
          geom_rect(data = shared_regions,
                    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                    fill = shared_fill, alpha = 0.28, inherit.aes = FALSE)
      }
      if (!is.null(seg_df) && nrow(seg_df) > 0) {
        seg_col <- seg_df[seg_df$line_group == "Colored", ]
        seg_ns  <- seg_df[seg_df$line_group == "NS", ]
        if (is_both) {
          # Use linetype to distinguish sex
          if (nrow(seg_col) > 0) {
            p_bottom <- p_bottom +
              geom_segment(data = seg_col,
                           aes(x = x_mb, y = Lod, xend = x_next, yend = y_next,
                               group = interaction(display_pheno, sex_label),
                               color = as.character(display_pheno),
                               linetype = sex_label),
                           linewidth = 1, na.rm = TRUE)
          }
          if (nrow(seg_ns) > 0) {
            p_bottom <- p_bottom +
              geom_segment(data = seg_ns,
                           aes(x = x_mb, y = Lod, xend = x_next, yend = y_next,
                               group = interaction(display_pheno, sex_label),
                               color = "Non-significant",
                               linetype = sex_label),
                           linewidth = 1, na.rm = TRUE)
          }
        } else {
          if (nrow(seg_col) > 0) {
            p_bottom <- p_bottom +
              geom_segment(data = seg_col,
                           aes(x = x_mb, y = Lod, xend = x_next, yend = y_next,
                               group = display_pheno, color = as.character(display_pheno)),
                           linewidth = 1, na.rm = TRUE)
          }
          if (nrow(seg_ns) > 0) {
            p_bottom <- p_bottom +
              geom_segment(data = seg_ns,
                           aes(x = x_mb, y = Lod, xend = x_next, yend = y_next,
                               group = display_pheno, color = "Non-significant"),
                           linewidth = 1, na.rm = TRUE)
          }
        }
      }
      # Significant points
      df_sig_pts <- df_bottom[is_sig(df_bottom$Status), ]
      if (nrow(df_sig_pts) > 0) {
        if (is_both) {
          p_bottom <- p_bottom +
            geom_point(data = df_sig_pts,
                       aes(x = x_mb, y = Lod, color = as.character(display_pheno),
                           shape = sex_label),
                       size = 1.2, na.rm = TRUE) +
            scale_shape_manual(name = "Sex", values = c("Male" = 16, "Female" = 17))
        } else {
          p_bottom <- p_bottom +
            geom_point(data = df_sig_pts,
                       aes(x = x_mb, y = Lod, color = as.character(display_pheno)),
                       size = 0.9, na.rm = TRUE)
        }
      }

      p_bottom <- p_bottom +
        scale_color_manual(name = "Phenotype", values = color_values,
                           breaks = c(names(pal_bottom), "Non-significant"), drop = TRUE)
      if (is_both) {
        p_bottom <- p_bottom +
          scale_linetype_manual(name = "Sex", values = c("Male" = "solid", "Female" = "dashed"))
      }
      p_bottom <- p_bottom +
        scale_x_continuous(name = "Position (Mb)", limits = c(region_start, region_end),
                           breaks = br, labels = function(x) paste0(x, "M"),
                           expand = c(0, 0)) +
        scale_y_continuous(name = "LOD Score", expand = expansion(mult = c(0, 0.05))) +
        theme_classic() +
        theme(
          axis.title       = element_text(size = 12, color = text_col),
          axis.text        = element_text(size = 10, color = text_col),
          axis.line        = element_line(color = line_col),
          axis.ticks       = element_line(color = line_col),
          legend.title     = element_text(size = 11, color = text_col),
          legend.text      = element_text(size = 9, color = text_col),
          legend.background = element_rect(fill = bg_col, color = NA),
          legend.key       = element_rect(fill = bg_col, color = NA),
          panel.background = element_rect(fill = bg_col, color = NA),
          plot.background  = element_rect(fill = bg_col, color = NA),
          plot.margin      = margin(t = 2, r = 8, b = 5, l = 5)
        ) +
        guides(color = guide_legend(ncol = 1, byrow = TRUE))

      list(top = p_top, bottom = p_bottom)
    })

    # --- Dynamic plot height ---
    output$plot_container <- renderUI({
      h <- input$plot_height %||% 800
      plotOutput(session$ns("grouplod_plot"), height = paste0(h, "px"))
    })

    # --- Render the combined plot ---
    group_plot_g <- eventReactive(input$refresh, group_plot(), ignoreNULL = FALSE)
    output$grouplod_plot <- renderPlot({
      result <- group_plot_g()
      if (is.list(result) && !is.null(result$top) && !is.null(result$bottom)) {
        g_top    <- ggplotGrob(result$top)
        g_bottom <- ggplotGrob(result$bottom)
        max_widths <- grid::unit.pmax(g_top$widths, g_bottom$widths)
        g_top$widths    <- max_widths
        g_bottom$widths <- max_widths
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1,
                           heights = grid::unit(c(0.4, 0.6), "npc"))))
        grid::pushViewport(grid::viewport(layout.pos.row = 1))
        grid::grid.draw(g_top)
        grid::popViewport()
        grid::pushViewport(grid::viewport(layout.pos.row = 2))
        grid::grid.draw(g_bottom)
        grid::popViewport(2)
      } else {
        print(result)
      }
    }, bg = "transparent") |>
      bindCache(group_plot_g(), plot_theme_key(session))

    # --- Download handler ---
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("grouplod_", input$sex_select, "_", input$group_select, "_",
               input$chr_select, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        result <- group_plot()
        if (is.list(result) && !is.null(result$top) && !is.null(result$bottom)) {
          pdf(file, width = 14, height = 10)
          g_top    <- ggplotGrob(result$top)
          g_bottom <- ggplotGrob(result$bottom)
          max_widths <- grid::unit.pmax(g_top$widths, g_bottom$widths)
          g_top$widths    <- max_widths
          g_bottom$widths <- max_widths
          grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1,
                             heights = grid::unit(c(0.4, 0.6), "npc"))))
          grid::pushViewport(grid::viewport(layout.pos.row = 1))
          grid::grid.draw(g_top)
          grid::popViewport()
          grid::pushViewport(grid::viewport(layout.pos.row = 2))
          grid::grid.draw(g_bottom)
          grid::popViewport(2)
          dev.off()
        } else {
          ggsave(file, result, width = 14, height = 10, device = "pdf")
        }
      }
    )
  })
}
