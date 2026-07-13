# ./Server/R/06_srv_info.R
info_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    output$app_version <- renderText({ APP_VERSION })
    output$latest_commit_msg <- renderText({ APP_COMMIT_MSG })
    
    output$tagHistory <- renderUI({
      hist <- get_tag_history(limit = 50)
      if (nrow(hist) == 0) return(div("No tags found."))
      
      tagList(lapply(seq_len(nrow(hist)), function(i) {
        tags$div(
          style = "padding:10px; border:1px solid #eee; border-radius:8px; margin-bottom:8px;",
          tags$div(
            tags$b(paste0("Version ", hist$tag[i])),
            HTML(sprintf("&nbsp;&nbsp;(<code>%s</code> | %s)", hist$commit[i], hist$date[i]))
          ),
          tags$pre(style = "white-space: pre-wrap; margin:6px 0 0 0;", hist$message[i])
        )
      }))
    })
  })
}
