phenotype_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Load required libraries (library() fails fast if a package is missing).
    # Note: plyr is intentionally NOT loaded - it masks dplyr's summarise/mutate,
    # which are used throughout this module, and it is not needed here.
    # ComplexHeatmap (a heavy Bioconductor package, ~2-4 s to attach) is NOT loaded
    # here: every call to it below is ComplexHeatmap::-qualified, so its namespace
    # loads lazily on the first heatmap render instead of at every session startup.
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(DT)
    library(readxl)
    library(reshape2)
    
    # Define color palettes and themes
    strainColors <- c("#000000","#FF9300","#107F40","#0000FE","#FF2600", "#942192","#C0C0C0")
    names(strainColors) <- c("C57BL/6J","129S1/SvImJ", "CAST/EiJ","PWK/PhJ", "B6CASTF1","129SPWKF1","B6CAST-129SPWK-F2")
    sexColors <- c("m" = "#00A2FF", "f" = "#FF644E")

    # Root input - lets the header "Refresh" button re-render this module's plots.
    root_input <- session$userData$root_input %||% NULL

    # Reactive values for data storage
    values <- reactiveValues(
      data_wide = NULL,
      data_melt = NULL,
      bw_time_course = NULL,
      current_plot = NULL,
      load_error = NULL
    )
    
    safe_numeric <- function(x) {
      if (inherits(x, "Date") || inherits(x, "POSIXt")) return(rep(NA_real_, length(x)))
      if (is.numeric(x)) return(x)
      xx <- as.character(x)
      xx <- gsub(",", ".", xx)                    # fix decimal commas
      keep <- grepl("^\\s*-?\\d+(\\.\\d+)?\\s*$", xx)
      out <- rep(NA_real_, length(xx))
      out[keep] <- as.numeric(xx[keep])
      out
    }
    
    read_bw_file <- function(path) {
      # read header only to decide types
      hdr <- suppressWarnings(readxl::read_excel(path, n_max = 0))
      ct  <- rep("guess", ncol(hdr))
      wk  <- grep("^Week\\b", names(hdr))
      if (length(wk) > 0) ct[wk] <- "text"       # avoid readxl date warnings
      
      bw <- suppressWarnings(readxl::read_excel(path, col_types = ct))
      colnames(bw)[colnames(bw) == "Animal ID"] <- "Pyrat_ID"
      bw
    }
    
    # Load data on initialization
    observe({
      values$load_error <- NULL
      tryCatch({
        # Load phenotype data (data_dir defined in app.R)
        if(file.exists(paste0(data_dir, '/allData_wide_clean_noOutlier.csv'))) {
          values$data_wide <- read.csv(paste0(data_dir, '/allData_wide_clean_noOutlier.csv'))
        }
        
        if(file.exists(paste0(data_dir, '/allData_melt_clean_noOutlier.csv'))) {
          values$data_melt <- read.csv(paste0(data_dir, '/allData_melt_clean_noOutlier.csv'))
        }
        
        # Load body weight time course data (phenotype_dir defined in app.R)
        if (file.exists(paste0(phenotype_dir, "/A_BW/Aim2Step2_BodyWeight.xlsx"))) {
          bw <- read_bw_file(paste0(phenotype_dir, "/A_BW/Aim2Step2_BodyWeight.xlsx"))
          bw$Strain[which((bw$Strain == "PWK/PhJ") & (bw$Sex == "m"))] <- "129SPWKF1"
          
          if (file.exists(paste0(phenotype_dir, "/A_BW/Aim2a_BodyWeight_pwk.xlsx"))) {
            bw_aim2a <- read_bw_file(paste0(phenotype_dir, "/A_BW/Aim2a_BodyWeight_pwk.xlsx"))
            bw <- dplyr::bind_rows(bw, bw_aim2a)
          }
          
          bw_time_course <- bw[bw$Parameter == "BW", ]
          bw_time_course <- bw_time_course[!is.na(bw_time_course$Pyrat_ID),]
          bw_time_course <- bw_time_course[, !grepl("Vet|Comments", names(bw_time_course))]
          bw_time_course <- as.data.frame(bw_time_course)
          
          # robust numeric conversion of all Week columns
          week_cols <- grep("^Week\\b", names(bw_time_course))
          if (length(week_cols) > 0) {
            bw_time_course[week_cols] <- lapply(bw_time_course[week_cols], safe_numeric)
          }
          
          # cleanups
          bw_time_course$Generation <- "F0"
          bw_time_course[bw_time_course$Strain == "B6CAST-129SPWK-F2", "Generation"] <- "F2"
          bw_time_course[bw_time_course$Strain %in% c("B6CASTF1", "129SPWKF1"), "Generation"] <- "F1"
          
          # reshape
          id_vars  <- names(bw_time_course)[!grepl("^Week\\b", names(bw_time_course))]
          meas_vars <- names(bw_time_course)[grepl("^Week\\b", names(bw_time_course))]
          if (length(meas_vars) > 0) {
            bw_melt <- reshape2::melt(bw_time_course, id.vars = id_vars,
                                      measure.vars = meas_vars, variable.name = "Week",
                                      value.name = "BW_value")
            bw_melt$Week <- as.numeric(gsub("^Week[ .]", "", as.character(bw_melt$Week)))
            bw_melt <- bw_melt[!is.na(bw_melt$BW_value), ]
            # clean odd strain names
            bw_melt$Strain <- gsub("^\\*", "", bw_melt$Strain)
            
            # derive initial/final BW per animal
            bw_init <- bw_melt %>%
              dplyr::group_by(Pyrat_ID) %>%
              dplyr::summarise(bw_initial = if (any(Week == 6)) BW_value[Week == 6][1] else NA_real_, .groups = "drop")
            
            bw_fin <- bw_melt %>%
              dplyr::group_by(Pyrat_ID) %>%
              dplyr::summarise(bw_final = if (any(Week == 24)) BW_value[Week == 24][1] else NA_real_, .groups = "drop")
            
            bw_melt <- bw_melt %>%
              dplyr::left_join(bw_init, by = "Pyrat_ID") %>%
              dplyr::left_join(bw_fin,  by = "Pyrat_ID") %>%
              dplyr::mutate(
                BW_NormBW6_per = 100 * BW_value / bw_initial,
                BW_perGain     = 100 * (BW_value - bw_initial) / bw_initial
              )
            
            # factors for plotting
            bw_melt$Sex        <- factor(bw_melt$Sex, levels = c("m","f"))
            bw_melt$Generation <- factor(bw_melt$Generation, levels = c("F0","F1","F2"))
            bw_melt$Strain     <- factor(bw_melt$Strain, levels = names(strainColors))
            
            # store for plotting
            values$bw_time_course <- bw_melt
          }
        }

        # If the primary table never loaded, the data directory is likely missing
        if (is.null(values$data_wide)) {
          values$load_error <- paste0(
            "Phenotype data not found under '", data_dir,
            "'. Check that the data directory is available (mounted)."
          )
        }
      }, error = function(e) {
        values$load_error <- paste("Error loading phenotype data:", conditionMessage(e))
        message(values$load_error)
      })
    })
    
    # Filter data based on user inputs
    filtered_data <- reactive({
      req(values$data_wide)
      
      data <- values$data_wide
      
      # Apply sex filter
      if(input$sex_filter != "both") {
        data <- data[data$Sex == input$sex_filter, ]
      }
      
      # Apply generation filter
      if(input$generation_filter != "all") {
        data <- data[data$Generation == input$generation_filter, ]
      }
      
      # Apply strain filter
      if(length(input$strain_filter) > 0) {
        data <- data[data$Strain %in% input$strain_filter, ]
      }
      
      data
    })
    
    # Build heatmap object (reused by renderPlot and the download handler)
    heatmap_obj <- reactive({
        validate(need(is.null(values$load_error), values$load_error))
        req(filtered_data())

        dfPlot <- NULL
        if (input$analysis_type == "fat_heatmap") {
          cols_needed <- c("ID", "Strain", "Sex")
          fat_cols <- grep("Fat_perc", names(filtered_data()), value = TRUE)
          if (length(fat_cols) > 0) dfPlot <- filtered_data()[, c(cols_needed, fat_cols)]
        } else if (input$analysis_type == "tissue_zscore_heatmap") {
          cols_needed <- c("ID", "Strain", "Sex", "Heart_g_24", "Liver_g_24", "Kidney_g_24", "Spleen_g_24")
          cols_exist <- cols_needed[cols_needed %in% names(filtered_data())]
          if (length(cols_exist) >= 4) dfPlot <- filtered_data()[, cols_exist]
        } else if (input$analysis_type == "tissue_rel_heatmap") {
          cols_needed <- c("ID", "Strain", "Sex", "Heart_g_RelBWSac_24", "Liver_g_RelBWSac_24", "Kidney_g_RelBWSac_24", "Spleen_g_RelBWSac_24")
          cols_exist <- cols_needed[cols_needed %in% names(filtered_data())]
          if (length(cols_exist) >= 4) dfPlot <- filtered_data()[, cols_exist]
        } else if (input$analysis_type == "liver_heatmap") {
          cols_needed <- c("ID","Strain","Sex","BW_perGain_22","Fat_perc_22","Liver_g_RelBWSac_24","Liver_g_24",
                           "ALAT_24","TIMP1_24","sum_all_vacuoles_percentage_24","Fibrosis_perc_24")
          cols_exist <- cols_needed[cols_needed %in% names(filtered_data())]
          if (length(cols_exist) >= 4) dfPlot <- filtered_data()[, cols_exist]
        } else if (input$analysis_type == "liver_f2_heatmap") {
          data_f2 <- filtered_data()
          if ("Generation" %in% names(data_f2)) data_f2 <- data_f2[data_f2$Generation == "F2", ]
          cols_needed <- c("ID","Strain","Sex","Fat_perc_22","Liver_g_RelBWSac_24","ALAT_24","TIMP1_24",
                           "sum_all_vacuoles_percentage_24","Fibrosis_perc_24")
          cols_exist <- cols_needed[cols_needed %in% names(data_f2)]
          if (length(cols_exist) >= 4) dfPlot <- data_f2[, cols_exist]
        }
        
        validate(need(!is.null(dfPlot) && nrow(dfPlot) > 0, "No data available for this heatmap"))
        dfPlot <- stats::na.omit(dfPlot)
        validate(need(nrow(dfPlot) > 0, "No complete cases to plot"))
        
        hm.identity <- dfPlot[, c("ID","Strain","Sex")]
        rownames(hm.identity) <- hm.identity$ID
        hm.identity <- hm.identity[, -1, drop = FALSE]
        
        dfPlot_melt <- reshape2::melt(dfPlot, id.vars = c("ID","Strain","Sex"))
        
        if (input$analysis_type != "fat_heatmap") {
          dfPlot_processed <- dfPlot_melt %>%
            dplyr::group_by(variable) %>%
            dplyr::mutate(zscore = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
            dplyr::select(variable, ID, zscore) %>%
            tidyr::spread(ID, zscore)
          mat <- as.matrix(dfPlot_processed[, -1, drop = FALSE])
          rownames(mat) <- dfPlot_processed$variable
          legend_name <- "z-score"
        } else {
          dfPlot_processed <- dfPlot_melt %>%
            dplyr::select(variable, ID, value) %>%
            tidyr::spread(ID, value)
          mat <- as.matrix(dfPlot_processed[, -1, drop = FALSE])
          rownames(mat) <- dfPlot_processed$variable
          legend_name <- "Fat %"
        }
        
        # palette
        col <- switch(
          input$color_palette,
          rwb = colorRampPalette(c("red","white","blue"))(256),
          bwr = colorRampPalette(c("blue","white","red"))(256),
          viridis = viridis::viridis(256),
          plasma = viridis::plasma(256)
        )
        
        # annotation data.frame for pheatmap (rownames must match the heatmap columns)
        ann <- hm.identity[colnames(mat), c("Strain", "Sex"), drop = FALSE]
        ann_colors <- list(
          Strain = strainColors[intersect(names(strainColors), unique(as.character(ann$Strain)))],
          Sex    = sexColors[intersect(names(sexColors),       unique(as.character(ann$Sex)))]
        )

        km <- if (isTRUE(input$cluster_columns) &&
                  isTRUE(input$enable_kmeans) &&
                  !is.null(input$column_km) &&
                  !is.na(input$column_km) &&
                  input$column_km > 1) as.integer(input$column_km) else 1L
        km <- min(km, ncol(mat))   # do not exceed number of columns

        # Return the drawing ingredients; draw_phenotype_heatmap() renders them.
        list(
          mat = mat, ann = ann, ann_colors = ann_colors, col = col,
          cluster_columns   = isTRUE(input$cluster_columns),
          cluster_rows      = isTRUE(input$cluster_rows),
          km                = km,
          tree_height       = input$tree_height %||% 150,
          show_column_names = isTRUE(input$show_column_names)
        )
    })

    # Draw the phenotype heatmap with pheatmap (CRAN). Replaces ComplexHeatmap so the
    # app carries no Bioconductor dependency and deploys cleanly to Connect Cloud.
    draw_phenotype_heatmap <- function(obj) {
      if (is.null(obj) || is.null(obj$mat) || !nrow(obj$mat) || !ncol(obj$mat)) {
        grid::grid.newpage(); grid::grid.text("No data to plot"); return(invisible())
      }
      pheatmap::pheatmap(
        obj$mat,
        color             = obj$col,
        cluster_cols      = obj$cluster_columns,
        cluster_rows      = obj$cluster_rows,
        treeheight_col    = if (isTRUE(obj$cluster_columns)) obj$tree_height else 0,
        treeheight_row    = if (isTRUE(obj$cluster_rows))    obj$tree_height else 0,
        cutree_cols       = if (obj$cluster_columns && obj$km > 1) obj$km else NA,
        annotation_col    = obj$ann,
        annotation_colors = obj$ann_colors,
        show_colnames     = obj$show_column_names,
        show_rownames     = TRUE,
        border_color      = "white",
        na_col            = "black",
        silent            = FALSE
      )
    }

    heatmap_g <- eventReactive(input$refresh, heatmap_obj(), ignoreNULL = FALSE)
    output$heatmapPlot <- renderPlot({
      draw_phenotype_heatmap(heatmap_g())
    }, bg = "transparent") |>
      bindCache(heatmap_g(), plot_theme_key(session))

    # Build body weight plot (reused by renderPlot and the download handler)
    bodyWeightPlot_obj <- reactive({
      if (!isTRUE(input$analysis_type %in% c("bw_initial","bw_final","bw_timecourse"))) {
        return(NULL)
      }
      validate(need(is.null(values$load_error), values$load_error))
      req(values$bw_time_course)
      validate(need(nrow(values$bw_time_course) > 0, "Body weight data not loaded yet"))
      
      bw_data <- values$bw_time_course
      
      # Apply filters
      if(input$sex_filter != "both") {
        bw_data <- bw_data[bw_data$Sex == input$sex_filter, ]
      }
      if(input$generation_filter != "all") {
        bw_data <- bw_data[bw_data$Generation == input$generation_filter, ]
      }
      if(length(input$strain_filter) > 0) {
        bw_data <- bw_data[bw_data$Strain %in% input$strain_filter, ]
      }
      
      if(nrow(bw_data) == 0) return(NULL)
      
      if(input$analysis_type == "bw_initial") {
        # Initial body weight
        bw_week6 <- bw_data[bw_data$Week == 6, ]
        
        if(input$bw_display_type == "bar") {
          p <- ggplot(bw_week6, aes(x = reorder(Pyrat_ID, bw_initial), 
                                    y = bw_initial, fill = Strain)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = strainColors) +
            theme_minimal() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()) +
            labs(x = "Animals", y = "Initial Body Weight (g)")
        } else {
          p <- ggplot(bw_week6, aes(x = Strain, y = bw_initial, fill = Strain)) +
            geom_boxplot() +
            geom_jitter(width = 0.2, alpha = 0.5) +
            scale_fill_manual(values = strainColors) +
            theme_minimal() +
            labs(x = "Strain", y = "Initial Body Weight (g)")
        }
        
        if(input$separate_by_sex && input$sex_filter == "both") {
          p <- p + facet_wrap(~Sex)
        }
        
      } else if(input$analysis_type == "bw_final") {
        # Final body weight
        bw_week24 <- bw_data[bw_data$Week == 24, ]
        
        if(input$bw_display_type == "bar") {
          p <- ggplot(bw_week24, aes(x = reorder(Pyrat_ID, bw_final), 
                                     y = bw_final, fill = Strain)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = strainColors) +
            theme_minimal() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()) +
            labs(x = "Animals", y = "Final Body Weight (g)")
        } else {
          p <- ggplot(bw_week24, aes(x = Strain, y = bw_final, fill = Strain)) +
            geom_boxplot() +
            geom_jitter(width = 0.2, alpha = 0.5) +
            scale_fill_manual(values = strainColors) +
            theme_minimal() +
            labs(x = "Strain", y = "Final Body Weight (g)")
        }
        
        if(input$separate_by_sex && input$sex_filter == "both") {
          p <- p + facet_wrap(~Sex)
        }
        
      } else if(input$analysis_type == "bw_timecourse") {
        # Body weight time course
        metric_col <- input$bw_metric
        
        p <- ggplot(bw_data, aes(x = Week, y = .data[[metric_col]],
                                 color = Strain, group = Pyrat_ID))
        
        if (isTRUE(input$show_individual_lines)) {
          # outline line behind colored line
          p <- p +
            geom_line(linewidth = 1.2, colour = "black", alpha = input$line_alpha, show.legend = FALSE) +
            geom_line(linewidth = 0.8, alpha = input$line_alpha) +
            geom_point(size = input$point_size + 0.6, colour = "black", alpha = input$line_alpha, show.legend = FALSE) +
            geom_point(size = input$point_size, alpha = input$line_alpha)
        }
        
        if (isTRUE(input$show_median_line)) {
          p <- p +
            stat_summary(aes(group = Strain), fun = median,
                         geom = "line", linewidth = 2.0, colour = "black", show.legend = FALSE) +
            stat_summary(aes(group = Strain), fun = median,
                         geom = "line", linewidth = 1.2)
        }
        
        p <- p +
          scale_color_manual(values = strainColors) +
          scale_x_continuous(breaks = seq(6, 24, by = 2)) +
          theme_minimal() +
          labs(x = "Age (weeks)",
               y = switch(metric_col,
                          "BW_value" = "Body Weight (g)",
                          "BW_perGain" = "Body Weight Gain (% from start)",
                          "BW_NormBW6_per" = "Body Weight (% of week 6)"))
        
        if (input$separate_by_sex && input$sex_filter == "both") {
          p <- p + facet_grid(Sex ~ Strain)
        } else {
          p <- p + facet_wrap(~Strain, nrow = 1)
        }
      }
      
      p + theme(text = element_text(size = input$font_size))
    })

    bw_plot_g <- eventReactive(input$refresh, bodyWeightPlot_obj(), ignoreNULL = FALSE)
    output$bodyWeightPlot <- renderPlot({
      bw_plot_g()
    }, res = 144, bg = "white") |>
      bindCache(bw_plot_g(), plot_theme_key(session))

    # Build variance plot (reused by renderPlot and the download handler)
    variancePlot_obj <- reactive({
      validate(need(is.null(values$load_error), values$load_error))
      req(values$data_melt)

      p <- NULL
      if(input$analysis_type == "pheno_availability") {
        # Phenotype availability plot
        counts <- values$data_melt %>%
          dplyr::group_by(Week, Variable) %>%
          dplyr::summarise(n = n(), .groups = "drop") %>%
          dplyr::filter(n > 0)
        
        p <- ggplot(counts, aes(x = as.factor(Week), y = Variable)) +
          geom_point(size = 3, color = "darkgrey") +
          scale_y_discrete(limits = rev) +
          theme_minimal() +
          labs(x = "Week", y = "Phenotype")
        
      } else if(input$analysis_type == "steatosis_variation") {
        # Steatosis variation
        steatosis_data <- values$data_melt %>%
          dplyr::filter(Variable == "sum_all_vacuoles_percentage" & Week == 24)
        
        if(nrow(steatosis_data) > 0) {
          p <- ggplot(steatosis_data, aes(x = reorder(ID, Value), y = Value, 
                                          fill = Sex)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = sexColors) +
            theme_minimal() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()) +
            labs(x = "Animals", y = "Steatosis (%)")
        } else {
          p <- ggplot() + theme_void() + 
            geom_text(aes(x = 0, y = 0, label = "No steatosis data available"))
        }
        
      } else if(input$analysis_type == "bw_variance") {
        # Body weight variance
        data_w  <- filtered_data()
        bw_cols <- grep("^BW_[0-9]+$", names(data_w), value = TRUE)
        
        if (length(bw_cols) > 0) {
          df <- data_w %>%
            dplyr::select(Sex, dplyr::all_of(bw_cols)) %>%
            tidyr::gather(Variable, Value, -Sex) %>%
            dplyr::mutate(WeekNum = as.integer(sub("^BW_(\\d+)$", "\\1", Variable)))
          
          if (input$sex_filter != "both" || isTRUE(input$separate_by_sex)) {
            grp <- df %>%
              dplyr::group_by(Variable, Sex, WeekNum) %>%
              dplyr::summarise(variance = var(Value, na.rm = TRUE), .groups = "drop")
            lvl <- grp %>% dplyr::distinct(Variable, WeekNum) %>% dplyr::arrange(WeekNum) %>% dplyr::pull(Variable)
            grp$Variable <- factor(grp$Variable, levels = lvl)
            
            p <- ggplot(grp, aes(x = Variable, y = variance, fill = Sex)) +
              geom_bar(stat = "identity", position = "dodge") +
              scale_fill_manual(values = sexColors, drop = FALSE) +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(x = "Week", y = "Variance")
          } else {
            grp <- df %>%
              dplyr::group_by(Variable, WeekNum) %>%
              dplyr::summarise(variance = var(Value, na.rm = TRUE), .groups = "drop")
            lvl <- grp %>% dplyr::distinct(Variable, WeekNum) %>% dplyr::arrange(WeekNum) %>% dplyr::pull(Variable)
            grp$Variable <- factor(grp$Variable, levels = lvl)
            
            p <- ggplot(grp, aes(x = Variable, y = variance)) +
              geom_bar(stat = "identity", fill = "grey50") +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(x = "Week", y = "Variance")
          }
        } else {
          p <- ggplot() + theme_void() +
            geom_text(aes(x = 0, y = 0, label = "No body weight data available"))
        }
        
      } else if(input$analysis_type == "fatlean_variance") {
        # Fat/Lean variance
        data_w <- filtered_data()
        
        fat_cols  <- grep("Fat_g_",  names(data_w), value = TRUE)
        lean_cols <- grep("Lean_g_", names(data_w), value = TRUE)
        fatlean_cols <- c(fat_cols, lean_cols)
        fatlean_cols <- fatlean_cols[!grepl("24$", fatlean_cols)]
        
        if (length(fatlean_cols) > 0) {
          df <- data_w %>%
            dplyr::select(Sex, dplyr::all_of(fatlean_cols)) %>%
            tidyr::gather(Variable, Value, -Sex) %>%
            dplyr::mutate(
              WeekNum = as.integer(sub("^.*_(\\d+)$", "\\1", Variable)),
              Measure = factor(ifelse(grepl("^Fat_", Variable), "Fat", "Lean"), levels = c("Fat","Lean"))
            )
          
          if (input$sex_filter != "both" || isTRUE(input$separate_by_sex)) {
            grp <- df %>%
              dplyr::group_by(Variable, Sex, WeekNum, Measure) %>%
              dplyr::summarise(variance = var(Value, na.rm = TRUE), .groups = "drop")
            lvl <- grp %>% dplyr::distinct(Variable, WeekNum, Measure) %>% dplyr::arrange(WeekNum, Measure) %>% dplyr::pull(Variable)
            grp$Variable <- factor(grp$Variable, levels = lvl)
            
            p <- ggplot(grp, aes(x = Variable, y = variance, fill = Sex)) +
              geom_bar(stat = "identity", position = "dodge") +
              scale_fill_manual(values = sexColors, drop = FALSE) +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(x = "Week (Fat / Lean)", y = "Variance")
          } else {
            grp <- df %>%
              dplyr::group_by(Variable, WeekNum, Measure) %>%
              dplyr::summarise(variance = var(Value, na.rm = TRUE), .groups = "drop")
            lvl <- grp %>% dplyr::distinct(Variable, WeekNum, Measure) %>% dplyr::arrange(WeekNum, Measure) %>% dplyr::pull(Variable)
            grp$Variable <- factor(grp$Variable, levels = lvl)
            
            p <- ggplot(grp, aes(x = Variable, y = variance)) +
              geom_bar(stat = "identity", fill = "grey50") +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(x = "Week (Fat / Lean)", y = "Variance")
          }
        } else {
          p <- ggplot() + theme_void() +
            geom_text(aes(x = 0, y = 0, label = "No fat/lean data available"))
        }
        
      } else if(input$analysis_type == "cv_analysis") {
        # Coefficient of variation
        numeric_cols <- names(values$data_wide)[sapply(values$data_wide, is.numeric)]
        numeric_cols <- numeric_cols[!numeric_cols %in% c("ID")]
        
        if(length(numeric_cols) > 0) {
          cv_data <- values$data_wide %>%
            dplyr::select(all_of(numeric_cols)) %>%
            dplyr::summarise_all(list(
              mean = ~mean(., na.rm = TRUE),
              sd = ~sd(., na.rm = TRUE)
            )) %>%
            tidyr::gather(key, value) %>%
            tidyr::separate(key, c("Variable", "stat"), sep = "_(?=[^_]+$)", fill = "right") %>%
            tidyr::spread(stat, value) %>%
            dplyr::mutate(cv = sd/mean) %>%
            dplyr::filter(!is.na(cv) & !is.infinite(cv))
          
          p <- ggplot(cv_data, aes(x = cv, y = reorder(Variable, cv))) +
            geom_point(size = 2, shape = 18, color = "darkblue") +
            theme_minimal() +
            labs(x = "Coefficient of Variation", y = "")
        } else {
          p <- ggplot() + theme_void() + 
            geom_text(aes(x = 0, y = 0, label = "No numeric data available"))
        }
      }
      
      validate(need(!is.null(p), "No variance plot for this selection"))
      p + theme(text = element_text(size = input$font_size))
    })

    variance_g <- eventReactive(input$refresh, variancePlot_obj(), ignoreNULL = FALSE)
    output$variancePlot <- renderPlot({
      variance_g()
    }, bg = "white", res = 144) |>
      bindCache(variance_g(), plot_theme_key(session))

    # Cluster plot - not yet implemented; show an explicit placeholder
    # instead of a blank panel.
    output$clusterPlot <- renderPlot({
      grid::grid.newpage()
      grid::grid.text("Cluster analysis is not yet implemented.",
                      gp = grid::gpar(col = "grey50", fontsize = 14))
    }, res = 144)
    
    # Summary statistics (shared by the on-screen table and its CSV download).
    summary_stats_data <- reactive({
      req(values$data_wide)
      numeric_cols <- names(values$data_wide)[sapply(values$data_wide, is.numeric)]
      numeric_cols <- numeric_cols[!numeric_cols %in% c("ID")]
      if (length(numeric_cols) == 0) return(NULL)
      values$data_wide %>%
        dplyr::select(all_of(numeric_cols)) %>%
        dplyr::summarise_all(list(
          n = ~sum(!is.na(.)),
          mean = ~mean(., na.rm = TRUE),
          sd = ~sd(., na.rm = TRUE),
          min = ~min(., na.rm = TRUE),
          max = ~max(., na.rm = TRUE),
          cv = ~sd(., na.rm = TRUE)/mean(., na.rm = TRUE)
        )) %>%
        tidyr::gather(key, value) %>%
        tidyr::separate(key, c("Variable", "Statistic"), sep = "_(?=[^_]+$)", fill = "right") %>%
        tidyr::spread(Statistic, value) %>%
        dplyr::mutate(across(where(is.numeric), ~round(., 3)))
    })

    # Summary statistics table
    output$summaryTable <- DT::renderDT({
      summary_stats <- summary_stats_data()
      if (!is.null(summary_stats)) {
        datatable(summary_stats,
                  options = list(
                    pageLength = 20,
                    scrollX = TRUE,
                    scrollY = "400px"
                  ),
                  rownames = FALSE)
      } else {
        # Return empty datatable with message
        datatable(data.frame(Message = "No numeric data available"),
                  options = list(dom = 't'),
                  rownames = FALSE)
      }
    })

    output$download_summary <- downloadHandler(
      filename = function() paste0("phenotype_summary_stats_", Sys.Date(), ".csv"),
      content = function(file) {
        d <- tryCatch(summary_stats_data(), error = function(e) NULL)
        if (is.null(d)) d <- data.frame(Message = "No numeric data available")
        write.csv(d, file, row.names = FALSE)
      }
    )
    
    # Reset to defaults
    observeEvent(input$reset, {
      updateSelectInput(session, "analysis_type", selected = "fat_heatmap")
      updateCheckboxInput(session, "cluster_columns", value = TRUE)
      updateCheckboxInput(session, "cluster_rows", value = TRUE)
      updateNumericInput(session, "column_km", value = 4)
      updateSelectInput(session, "color_palette", selected = "rwb")
      updateCheckboxInput(session, "show_column_names", value = FALSE)
      updateCheckboxInput(session, "show_individual_lines", value = TRUE)
      updateCheckboxInput(session, "show_median_line", value = FALSE)
      updateSelectInput(session, "bw_metric", selected = "BW_perGain")
      updateNumericInput(session, "line_alpha", value = 0.3)
      updateNumericInput(session, "point_size", value = 0.5)
      updateSelectInput(session, "bw_display_type", selected = "bar")
      updateCheckboxInput(session, "separate_by_sex", value = TRUE)
      updateSelectInput(session, "sex_filter", selected = "both")
      updateSelectInput(session, "generation_filter", selected = "all")
      updateSelectizeInput(session, "strain_filter", selected = character(0))
      updateNumericInput(session, "plot_width", value = 25)
      updateNumericInput(session, "plot_height", value = 6)
      updateNumericInput(session, "font_size", value = 10)
      updateCheckboxInput(session, "enable_kmeans", value = FALSE)
    })
    
    # Download plot - render the same object the on-screen plot uses.
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("phenotype_", input$analysis_type, ".pdf")
      },
      content = function(file) {
        atype <- input$analysis_type
        pdf(file, width = input$plot_width / 2.54, height = input$plot_height / 2.54)
        on.exit(dev.off(), add = TRUE)
        tryCatch({
          if (grepl("heatmap", atype)) {
            draw_phenotype_heatmap(heatmap_obj())
          } else if (atype %in% c("bw_initial", "bw_final", "bw_timecourse")) {
            print(bodyWeightPlot_obj())
          } else if (grepl("variance|cv_analysis|steatosis|pheno_availability", atype)) {
            print(variancePlot_obj())
          } else {
            grid::grid.newpage()
            grid::grid.text("No plot available for this analysis type.")
          }
        }, error = function(e) {
          grid::grid.newpage()
          grid::grid.text(paste("Could not render plot:", conditionMessage(e)))
        })
      }
    )
  })
}