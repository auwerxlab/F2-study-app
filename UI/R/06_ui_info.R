# ./UI/R/06_ui_info.R
info_ui <- function(id) {
  ns <- NS(id)
  
  bslib::card(
    header = tagList(icon("info-circle"), " Info"),
    
    # Current version
    div(
      style = "margin-bottom: 12px;",
      tags$b("Current version: "), textOutput(ns("app_version"), inline = TRUE)
    ),
    
    # Latest commit message
    div(
      style = "white-space: pre-wrap; border: 1px solid #eee; 
               padding: 8px; border-radius: 6px; margin-bottom: 16px;",
      tags$b("Latest commit message:"), tags$br(),
      textOutput(ns("latest_commit_msg"))
    ),
    
    # Contact
    tags$h5("Contact"),
    div(
      style = "margin-bottom: 16px; padding: 12px; border: 1px solid #eee; border-radius: 6px;",
      tags$p(
        icon("laptop-code"), " For questions related to the ",
        tags$b("app"), ", contact: ",
        tags$a(href = "mailto:Alaa.Badreddine@epfl.ch", "Alaa Badreddine"),
        " - ", tags$code("Alaa.Badreddine@epfl.ch")
      ),
      tags$p(
        style = "margin-bottom: 0;",
        icon("flask"), " For questions related to the ",
        tags$b("research"), ", contact: ",
        tags$a(href = "mailto:giorgia.benegiamo@epfl.ch", "Giorgia Benegiamo"),
        " - ", tags$code("giorgia.benegiamo@epfl.ch")
      )
    ),

    # Release history
    tags$h5("Release history"),
    uiOutput(ns("tagHistory"))
  )
}
