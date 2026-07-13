# Libraries
library(shiny)
library(bslib)
library(shinydashboardPlus)
thematic::thematic_shiny(font = NA)


# ---- Render cache: reuse the instance RAM to skip recomputing identical plots ----

shinyOptions(cache = cachem::cache_mem(max_size = 1024 * 1024^2))   # ~1 GB / worker

# Light/dark signature for plot cache keys, so a cached image matches the current bslib colour mode
plot_theme_key <- function(session) {
  ri <- session$userData$root_input
  if (!is.null(ri) && isTRUE(ri$mode == "dark")) "dark" else "light"
}


# ---------- CSS / JS -----------
navbar_background <- source("./UI/CSS/navbar_background.R")$value
dark_sidebar_fix <- source("./UI/CSS/dark_sidebar_fix.R")$value
menu_styles      <- source("./UI/CSS/menu_styles.R")$value   # accordion + active-tab highlight
active_tab_js    <- source("./UI/JS/active_tab.R")$value      # toggles .active-tab on the selected menu button


# ----- Starting theme -----
cerulean_theme <- bs_theme(version = 5, bootswatch = "cerulean")


# ------ Data directories (single source of truth, used by all modules) ------
data_dir      <- "Aim2_Step2_manuscript/Data"            # phenotype CSVs (external mount)
phenotype_dir <- "Aim2_Step2_manuscript/phenotype_data"  # body-weight Excel (external mount)
qtl_data_dir  <- "Data/Figure1_all_strain"               # QTL Manhattan RDS
lod_data_dir  <- "Data/LOD"                               # Group LOD RDS (rds/ subfolder)
rnaseq_data_dir <- "Data/RNAseq"                          # RNA-seq RDS (PCA, clustering, DEA)
rnaseq_pca_path <- file.path(rnaseq_data_dir, "rnaseq_pca.rds")
rnaseq_hc_path  <- file.path(rnaseq_data_dir, "rnaseq_hc.rds")
rnaseq_dea_path <- file.path(rnaseq_data_dir, "rnaseq_dea.rds")
rnaseq_gsea_path <- file.path(rnaseq_data_dir, "rnaseq_gsea.rds")
rnaseq_ora_path  <- file.path(rnaseq_data_dir, "rnaseq_ora.rds")
# Tiny UI-choices sidecars (written next to the bundles by the preprocess scripts)
rnaseq_pca_meta_path <- file.path(rnaseq_data_dir, "rnaseq_pca_meta.rds")
rnaseq_hc_meta_path  <- file.path(rnaseq_data_dir, "rnaseq_hc_meta.rds")
rnaseq_dea_meta_path <- file.path(rnaseq_data_dir, "rnaseq_dea_meta.rds")
rnaseq_ora_meta_path <- file.path(rnaseq_data_dir, "rnaseq_ora_meta.rds")
compare_data_dir <- "Data/Comparison"                     # protein vs RNA-seq RDS
prot_rna_path    <- file.path(compare_data_dir, "prot_rna.rds")
prot_rna_meta_path <- file.path(compare_data_dir, "prot_rna_meta.rds")  # tiny UI-choices sidecar
proteomics_data_dir <- "Data/Proteomics"                  # liver proteomics RDS (PCA)
proteomics_pca_path <- file.path(proteomics_data_dir, "proteomics_pca.rds")
proteomics_pca_meta_path <- file.path(proteomics_data_dir, "proteomics_pca_meta.rds")


# ------ Functions ------
source("./UI/Functions/tooltip_label.R")      # Helper for parameter tooltips
source("./UI/Functions/method_panel.R")   # Collapsible per-tab "Methods" panel
source("./Server/Functions/menu_btn.R")  # Helper to build a sidebar button with icon
source("./Server/Functions/get_app_version.R")  # Version discovery helpers
APP_VERSION     <- get_app_version()
APP_COMMIT_MSG  <- {
  out <- if (has_git()) safe_sys("git", c("log", "-1", "--pretty=%B")) else character()
  if (length(out) == 0) "Unknown" else paste(out, collapse = "\n")
}


# ------ UI Components ------
source("./UI/R/01_ui_phenotype.R", local = TRUE)
source("./UI/R/02_ui_rnaseq_pca.R", local = TRUE)
source("./UI/R/03_ui_rnaseq_hc.R", local = TRUE)
source("./UI/R/04_ui_rnaseq_dea.R", local = TRUE)
source("./UI/R/05_ui_qtl.R", local = TRUE)
source("./UI/R/06_ui_info.R", local = TRUE)
source("./UI/R/07_ui_grouplod.R", local = TRUE)
source("./UI/R/08_ui_compare.R", local = TRUE)
source("./UI/R/09_ui_proteomics_pca.R", local = TRUE)


# ------ Server Components ------
source("./Server/R/01_srv_phenotype.R", local = TRUE)
source("./Server/R/02_srv_rnaseq_pca.R", local = TRUE)
source("./Server/R/03_srv_rnaseq_hc.R", local = TRUE)
source("./Server/R/04_srv_rnaseq_dea.R", local = TRUE)
source("./Server/R/05_srv_qtl.R", local = TRUE)
source("./Server/R/06_srv_info.R", local = TRUE)
source("./Server/R/07_srv_grouplod.R", local = TRUE)
source("./Server/R/08_srv_compare.R", local = TRUE)
source("./Server/R/09_srv_proteomics_pca.R", local = TRUE)

ui <- function(request) {
  page_sidebar(
    theme = cerulean_theme,
    title = "A novel mouse cross uncovers candidate therapeutic targets for hepatic steatosis, adiposity, and dyslipidemia",
    # Spinner shown over any plot while it is recomputing
    tags$style(HTML("
      .shiny-plot-output.recalculating,
      .shiny-image-output.recalculating {
        opacity: 1 !important;
        position: relative;
      }
      .shiny-plot-output.recalculating::before,
      .shiny-image-output.recalculating::before {
        content: ''; position: absolute; inset: 0;
        background: rgba(127,127,127,0.12); z-index: 4; pointer-events: none;
      }
      .shiny-plot-output.recalculating::after,
      .shiny-image-output.recalculating::after {
        content: ''; position: absolute; top: 50%; left: 50%;
        width: 44px; height: 44px; margin: -22px 0 0 -22px;
        border: 4px solid rgba(127,127,127,0.3); border-top-color: #764ba2;
        border-radius: 50%; animation: cc-spin 0.8s linear infinite;
        z-index: 5; pointer-events: none;
      }
      @keyframes cc-spin { to { transform: rotate(360deg); } }
      /* Wider hover tooltips so multi-line help (e.g. Pearson vs Spearman) reads
         across the screen instead of wrapping into a narrow column. */
      .tooltip-inner { max-width: 540px; }
    ")),
    navbar_background,
    dark_sidebar_fix,
    menu_styles,
    active_tab_js,
    # navset_card_tab: keep the sub-tabs next to the card title (left)
    tags$style(HTML(
      ".bslib-navs-card-title { justify-content: flex-start !important; gap: 1.25rem; }
       .bslib-navs-card-title > .nav { margin-left: 0 !important; }"
    )),
    sidebar = sidebar(
      title = "Menu",
      open = "desktop",
      
      # Menu buttons
      menu_btn("tab_intro",       "Introduction",  "home"),
      menu_btn("tab_pheno",       "Phenotypes",    "dna"),
      menu_btn("tab_rnaseq",      "RNA-seq",       "dna"),
      menu_btn("tab_proteomics",  "Proteomics",    "flask"),
      menu_btn("tab_qtl",         "QTL",           "project-diagram"),
      menu_btn("tab_grouplod",    "Group LOD",     "chart-bar"),
      menu_btn("tab_compare",     "Data comparison","exchange-alt"),
      tags$hr(),
      menu_btn("tab_info",        "Info",          "info-circle"),
      
      # Version footer
      tags$div(
        style = "position: fixed; bottom: 10px; left: 0; width: var(--bslib-sidebar-width, 250px);
           text-align: center; font-size: 12px; color: #6c757d;
           border-top: 1px solid #dee2e6; padding-top: 8px;",
        actionLink(
          inputId = "go_info",
          label = if (grepl("^[0-9]", APP_VERSION)) paste0("Version ", APP_VERSION) else APP_VERSION,
          class = "btn btn-link p-0",
          style = "text-decoration: underline; color: inherit;"
        ),
        tags$div(
          style = "margin-top: 4px; font-size: 11px; color: #6c757d; letter-spacing: 0.5px;",
          tags$style(HTML("
            .coffee-cup {
              display: inline-block;
              font-size: 14px;
              position: relative;
              animation: coffee-tilt 3s ease-in-out infinite;
            }
            @keyframes coffee-tilt {
              0%, 100% { transform: rotate(0deg); }
              25% { transform: rotate(-10deg); }
              50% { transform: rotate(0deg); }
              75% { transform: rotate(5deg); }
            }
            .steam {
              display: inline-block;
              position: relative;
              top: -2px;
              font-size: 11px;
              font-weight: bold;
              color: #4A3020;            /* dark coffee brown - solid, no fade */
              opacity: 1;
            }
            .steam:nth-child(2) { animation-delay: 0.4s; }
            .steam:nth-child(3) { animation-delay: 0.8s; }
            @keyframes steam-rise {
              0% { opacity: 0; transform: translateY(0px); }
              30% { opacity: 0.7; }
              100% { opacity: 0; transform: translateY(-10px); }
            }
            .footer-text {
              opacity: 1;
              transition: opacity 0.3s ease;
            }
            .footer-text:hover {
              opacity: 1;
            }
            .footer-text:hover .coffee-cup {
              animation: coffee-sip 0.6s ease forwards;
            }
            @keyframes coffee-sip {
              0% { transform: rotate(0deg); }
              50% { transform: rotate(-30deg); }
              100% { transform: rotate(0deg); }
            }
          ")),
          tags$span(class = "footer-text",
            HTML("brewed with "),
            tags$img(
              src = "https://media.tenor.com/Uq_-tDUQlJkAAAAj/hot-beverage-joypixels.gif",
              alt = "hot beverage",
              style = "height: 22px; vertical-align: middle; background: transparent;"
            ),
            HTML(" at EPFL")
          )
        )
      )
    ),
    
    # Header with time, date, dark mode
    div(
      style = "display: flex; justify-content: space-between; align-items: center; width: 100%; padding: 0 0 15px 0; border-bottom: 1px solid #dee2e6; margin-bottom: 20px;",
      uiOutput("main_title"),
      div(
        style = "display: flex; align-items: center; gap: 20px;",
        div(
          style = "text-align: right;",
          div(textOutput("currentDateTime", inline = TRUE), style = "font-weight: 500; font-size: 16px;")
        ),
        input_dark_mode(id = "mode", mode = "light")
      )
    ),
    
    uiOutput("main_content")
  )
}



server <- function(input, output, session) {
  # Share root input with modules (for dark mode detection)
  session$userData$root_input <- input

  # Clock
  output$currentDateTime <- renderText({
    invalidateLater(30000, session)
    paste(
      format(Sys.time(), "%A, %B %d, %Y"),
      "|",
      format(Sys.time(), "%H:%M")
    )
  })
  
  # Navigation state
  current_tab <- reactiveVal("intro")
  
  observeEvent(input$tab_intro,      { current_tab("intro") })
  observeEvent(input$tab_pheno,      { current_tab("pheno") })
  observeEvent(input$tab_rnaseq,     { current_tab("rnaseq") })
  observeEvent(input$tab_proteomics, { current_tab("proteomics") })
  observeEvent(input$tab_qtl,        { current_tab("qtl") })
  observeEvent(input$tab_grouplod,   { current_tab("grouplod") })
  observeEvent(input$tab_compare,    { current_tab("compare") })
  observeEvent(input$tab_info,       { current_tab("info") })

  # Highlight the selected menu button in the sidebar. Driven off current_tab()
  # so it tracks both button clicks and programmatic jumps (e.g. the Info link).
  observe({
    session$sendCustomMessage("setActiveTab", paste0("tab_", current_tab()))
  })

  # Persistent page title in the top header strip, so it is always clear which
  # tab is open. Mirrors the sidebar menu labels / icons.
  output$main_title <- renderUI({
    meta <- switch(
      current_tab(),
      intro      = list(icon = "home",            label = "Introduction"),
      pheno      = list(icon = "dna",             label = "Phenotypes"),
      rnaseq     = list(icon = "dna",             label = "RNA-seq"),
      proteomics = list(icon = "flask",           label = "Proteomics"),
      qtl        = list(icon = "project-diagram", label = "QTL"),
      grouplod   = list(icon = "chart-bar",       label = "Group LOD"),
      compare    = list(icon = "exchange-alt",    label = "Data comparison"),
      info       = list(icon = "info-circle",     label = "Info"),
      list(icon = "circle", label = "")
    )
    tags$h4(class = "m-0 fw-bold d-flex align-items-center gap-2",
            icon(meta$icon), tags$span(meta$label))
  })

  # Main content per tab
  output$main_content <- renderUI({
    switch(
      current_tab(),
      intro = card(
        header = tagList(icon("home"), " Introduction"),
        p("Interactive app to the F2 multi-strain liver study. Use the menu on the left to move between analyses."),

        # --- Study abstract ---
        div(
          class = "mb-3 p-3",
          style = "border-left: 4px solid #764ba2; background: rgba(118, 75, 162, 0.05); border-radius: 6px;",
          tags$p(class = "text-muted fst-italic mb-2", "Abstract"),
          tags$p(
            class = "mb-0", style = "text-align: justify;",
            "Metabolic dysfunction\u2013associated steatotic liver disease (MASLD) is a major component of cardiometabolic disease, yet the genetic factors that determine susceptibility to hepatic steatosis and progression to advanced liver injury remain incompletely understood. We generated a four-way mouse intercross spanning resistant and susceptible genetic backgrounds and profiled longitudinal metabolic, histological, transcriptomic, and proteomic phenotypes during dietary induction of metabolic disease. The resulting population of 572 mice captured the full spectrum of liver disease severity and revealed substantial genetic separation between adiposity and liver injury. Integration of mouse and human genetic evidence prioritized candidate genes for steatosis, adiposity, and plasma lipid traits, including ",
            tags$i("Atp8b1"), ", ", tags$i("Lgr4"), ", ", tags$i("Hsd17b12"), ", ",
            tags$i("C1qtnf4"), ", ", tags$i("Aqp9"), ", and ", tags$i("Cd163"),
            ". These findings show that studies in mouse populations provide a framework for identifying therapeutic targets underlying cardiometabolic disease."
          )
        ),

        bslib::accordion(
          id = "intro_info",
          bslib::accordion_panel(
            title = "Analysis Information",
            value = "info",
            open = TRUE,
            h6("Study design"),
            p("A multi-generational mouse cross profiled for liver phenotypes, transcriptome and proteome:"),
            tags$ul(
              tags$li(tags$b("F0"), " parental strains: C57BL/6J, 129S1/SvImJ, CAST/EiJ, PWK/PhJ"),
              tags$li(tags$b("F1"), " crosses: B6CASTF1, 129SPWKF1"),
              tags$li(tags$b("F2"), " intercross: B6CAST-129SPWK-F2"),
              tags$li("Measurements from weeks 6-24: body composition, metabolic markers and liver histology")
            ),
            h6("What each tab does"),
            tags$ul(
              tags$li(tags$b("Phenotypes"), " - heatmaps, body-weight, variance and cluster views of the physiological measurements."),
              tags$li(tags$b("RNA-seq"), " - liver transcriptome PCA, hierarchical clustering, and differential expression (limma-voom) with GO GSEA / ORA enrichment."),
              tags$li(tags$b("Proteomics"), " - liver proteome PCA."),
              tags$li(tags$b("QTL"), " - genome-wide Manhattan mapping across the 19 chromosomes."),
              tags$li(tags$b("Group LOD"), " - per-group, per-chromosome LOD profiles with interval bars."),
              tags$li(tags$b("Data comparison"), " - protein vs RNA-seq abundance correlation.")
            ),
            h6("Colour coding"),
            p("Strain colours follow a consistent scheme across every plot. Sex is coded blue (male) and red (female).")
          )
        )
      ),
      pheno = phenotype_ui("pheno"),
      rnaseq = bslib::navset_card_tab(
        title = tagList(icon("dna"), " RNA-seq"),
        bslib::nav_panel("PCA", rnaseq_pca_ui("rnaseq_pca")),
        bslib::nav_panel("Hierarchical clustering", rnaseq_hc_ui("rnaseq_hc")),
        bslib::nav_panel("Differential expression", rnaseq_dea_ui("rnaseq_dea"))
      ),
      proteomics = bslib::navset_card_tab(
        title = tagList(icon("flask"), " Proteomics"),
        bslib::nav_panel("PCA", proteomics_pca_ui("proteomics_pca"))
      ),
      qtl = qtl_ui("qtl"),
      grouplod = grouplod_ui("grouplod"),
      compare = compare_ui("compare"),
      info = info_ui("info")
    )
  })
  
  # Phenotype
  phenotype_server("pheno")

  # RNA-seq (PCA + hierarchical clustering + differential expression)
  rnaseq_pca_server("rnaseq_pca")
  rnaseq_hc_server("rnaseq_hc")
  rnaseq_dea_server("rnaseq_dea")

  # QTL
  qtl_server("qtl")

  # Group LOD
  grouplod_server("grouplod")

  # Data comparison (protein vs RNA-seq)
  compare_server("compare")

  # Proteomics (PCA)
  proteomics_pca_server("proteomics_pca")

  # INFO
  observeEvent(input$go_info, { current_tab("info") })
  info_server("info")
}



shinyApp(ui = ui, server = server)