# Help page ---------------------------------------------------------------
# Renders the HELP.md shipped with the app: navigation, page-by-page
# description, endpoint families, statistical caveats and runnable R examples
# that use the package's own functions and real cohorts.
ui_mod_help <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Help", class = "section-title"),
    div(class = "cpas-help", .md_to_html(.help_md_path()))
  )
}

server_mod_help <- function(id) {
  moduleServer(id, function(input, output, session) {})
}
