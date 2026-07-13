# ./UI/Functions/method_accordion.R
# Collapsible "Methods" text shown at the top of each analysis tab.
#
# Deliberately NOT a bslib accordion - just a native <details>/<summary> with a
# plain grey ">" that rotates to "v" when opened. Text is condensed from the
# study's Materials and Methods (manuscript F2_MASH_manuscript_v15). Collapsed by
# default so it never pushes the plot parameters or figures down.
#
# Styling lives in UI/CSS/menu_styles.R (.method-details / .method-summary).
method_panel <- function(..., title = "Methods") {
  args <- list(...)
  # Back-compat: ignore an optional leading intro string.
  if (length(args) && is.character(args[[1]]) && length(args[[1]]) == 1) args <- args[-1]

  tags$details(
    class = "method-details mb-2",
    tags$summary(class = "method-summary", title),
    tags$div(class = "method-body", args)
  )
}
