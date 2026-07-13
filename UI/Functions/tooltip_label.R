# Helper: create a label with hover tooltip icon
# Usage: label_tip("Sex", "Select male, female, or both sexes for analysis")
label_tip <- function(label_text, tip_text) {
  tags$span(
    label_text, " ",
    bslib::tooltip(
      tags$i(
        class = "fas fa-question-circle",
        style = "color: #6c757d; font-size: 11px; cursor: help;"
      ),
      tip_text,
      placement = "top"
    )
  )
}
