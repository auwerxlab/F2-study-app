menu_btn <- function(id, label, icon_name) {
  tags$div(
    style = "margin-bottom:8px;",
    actionLink(
      inputId = id,
      label = tagList(icon(icon_name), span(label)),
      class = "btn btn-outline-secondary w-100 text-start"
    )
  )
}