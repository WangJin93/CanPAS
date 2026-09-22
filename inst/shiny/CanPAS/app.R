# ================================================================
# CanPAS Shiny app
#
#   library(CanPAS); run_cpas_app()
#
# Top-level navigation (two groups, no deeper nesting):
#   Dashboard | Datasets
#   Single dataset analysis : KM analysis | COX analysis | COX by genes
#   Multi-datasets analysis : COX by datasets | Pooled KM | Meta-analysis
#   Methods | Help (HELP.md rendered in the app)
#
# A single shared selection (dataset + endpoint family) is synchronised across
# all pages: choosing a dataset anywhere updates the others, and Datasets rows
# can be pushed straight into the analysis pages. Every result page shows the
# figure above the table and prints 4 decimals.
# ================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bs4Dash)
  library(shinyWidgets)
  library(shinycssloaders)
  library(DT)
  library(CanPAS)
})

# ---- source helpers & modules (works from any cwd) ----------------
.appdir <- function() {
  if (dir.exists("apps")) return(getwd())
  system.file("shiny", "CanPAS", package = "CanPAS")
}
ad <- .appdir()
for (f in c("core.R", "mod_welcome.R", "mod_datasets.R", "mod_km.R", "mod_cox.R",
            "mod_cox_by_genes.R", "mod_cox_by_datasets.R", "mod_pooled_km.R",
            "mod_meta.R", "mod_methods.R", "mod_help.R")) {
  source(file.path(ad, "apps", f))
}

# ---- brand logo (top-left of the app) ----------------------------
# www/logo.jpg is served through an explicit resource path so the logo also
# resolves when the app is started from an arbitrary working directory.
.logo_path <- file.path(ad, "www", "logo.jpg")
if (file.exists(.logo_path) && !("cpas-assets" %in% names(shiny::resourcePaths())))
  shiny::addResourcePath("cpas-assets", file.path(ad, "www"))
.brand_logo <- if (file.exists(.logo_path))
  tags$img(src = "cpas-assets/logo.jpg", class = "cpas-logo", alt = "CanPAS logo")

dataset_info <- .norm_catalog(.pkg_data2("dataset_info"))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# ---- UI ----------------------------------------------------------
ui <- bs4Dash::dashboardPage(
  dark = FALSE,
  bs4Dash::dashboardHeader(
    title = tags$div(
      class = "cpas-brand",
      .brand_logo,
      tags$span(class = "cpas-brand-text", "Cancer Prognosis Analysis Suite")
    ),
    skin = "dark"
  ),
  bs4Dash::dashboardSidebar(
    skin = "light",
    bs4Dash::sidebarMenu(
      id = "sidebar",
      bs4Dash::menuItem("Dashboard", tabName = "welcome",
                        icon = icon("info-circle")),
      bs4Dash::menuItem("Datasets", tabName = "ds",
                        icon = icon("database")),
      bs4Dash::menuItem("Single dataset analysis",
                        icon = icon("chart-line"), startExpanded = TRUE,
        bs4Dash::menuSubItem("KM analysis", tabName = "sd_km",
                             icon = icon("chart-area")),
        bs4Dash::menuSubItem("COX analysis", tabName = "sd_cox",
                             icon = icon("cogs")),
        bs4Dash::menuSubItem("COX by genes", tabName = "sd_genes",
                             icon = icon("dna"))
      ),
      bs4Dash::menuItem("Multi-datasets analysis",
                        icon = icon("layer-group"), startExpanded = TRUE,
        bs4Dash::menuSubItem("COX by datasets", tabName = "md_cox",
                             icon = icon("table-columns")),
        bs4Dash::menuSubItem("Pooled KM", tabName = "md_km",
                             icon = icon("chart-line")),
        bs4Dash::menuSubItem("Meta-analysis", tabName = "md_meta",
                             icon = icon("chart-bar"))
      ),
      bs4Dash::menuItem("Methods", tabName = "methods",
                        icon = icon("book")),
      bs4Dash::menuItem("Help", tabName = "help",
                        icon = icon("question-circle"))
    )
  ),
  bs4Dash::dashboardBody(
    tags$head(tags$style(HTML(
      "/* ---- brand logo (top-left) ----------------------------------------
         logo.jpg has a light (#f3f5f5) background, so the brand block is
         painted in the same colour: the logo blends into the header instead
         of showing a white rectangle on the coloured bar. */
       .main-sidebar .cpas-brand{
         display:flex; align-items:center; flex-wrap:nowrap; gap:8px;
         padding:8px 10px 8px 12px; background:#f3f5f5;
         border-bottom:1px solid #dde3e8;
       }
       .main-sidebar .cpas-brand .cpas-logo{
         display:block; flex:0 0 auto; height:46px; width:auto; max-width:100%;
       }
       .main-sidebar .cpas-brand .cpas-brand-text{
         font-weight:600; color:#0C4B72; font-size:.78rem; line-height:1.2;
         letter-spacing:.01em;
       }
       /* collapsed (mini) sidebar: keep the logo inside the narrow rail */
       body.sidebar-mini.sidebar-collapse .main-sidebar .cpas-brand{
         justify-content:center; padding:6px 4px;
       }
       body.sidebar-mini.sidebar-collapse .main-sidebar:not(:hover) .cpas-brand-text{display:none;}
       body.sidebar-mini.sidebar-collapse .main-sidebar:not(:hover) .cpas-logo{height:30px;}
       @media (max-width: 991.98px){
         .main-sidebar .cpas-brand{justify-content:center;}
         .main-sidebar .cpas-brand .cpas-brand-text{display:none;}
       }
       /* short labels with a hover explanation (native title tooltip) */
       .main-sidebar .cpas-tip{font-size:.85em; opacity:.85;}
       .main-sidebar .cpas-tip:hover{opacity:1;}
       .cpas-hint{font-size:.78rem; color:#5b6b7b; line-height:1.25;
                  margin:-6px 0 6px 0;}"
    ))),
    tags$head(tags$script(HTML(
      "Shiny.addCustomMessageHandler('open_url', function(u){ window.open(u,'_blank'); });"
    ))),
    bs4Dash::tabItems(
      bs4Dash::tabItem("welcome", ui_mod_welcome("welcome")),
      bs4Dash::tabItem("ds", ui_mod_datasets("ds")),
      bs4Dash::tabItem("sd_km", ui_mod_km("km")),
      bs4Dash::tabItem("sd_cox", ui_mod_cox("cox")),
      bs4Dash::tabItem("sd_genes", ui_mod_cox_by_genes("genes")),
      bs4Dash::tabItem("md_cox", ui_mod_cox_by_datasets("bycohort")),
      bs4Dash::tabItem("md_km", ui_mod_pooled_km("pooled")),
      bs4Dash::tabItem("md_meta", ui_mod_meta("meta")),
      bs4Dash::tabItem("methods", ui_mod_methods("methods")),
      bs4Dash::tabItem("help", ui_mod_help("help"))
    )
  )
)

# ---- server ------------------------------------------------------
server <- function(input, output, session) {
  rv <- reactiveValues(
    sel = list(acc = NULL, type = NULL, family = NULL, ep = NULL),  # shared selection
    marker = NULL,                                                  # marker spec from KM
    km_last = NULL
  )

  jump <- function(tab) bs4Dash::updateTabItems(session, "sidebar", selected = tab)

  server_mod_welcome("welcome")
  server_mod_datasets("ds", rv, jump, dataset_info)
  server_mod_km("km", rv, dataset_info)
  server_mod_cox("cox", rv, dataset_info)
  server_mod_cox_by_genes("genes", rv, dataset_info)
  server_mod_cox_by_datasets("bycohort", rv, dataset_info)
  server_mod_pooled_km("pooled", rv, dataset_info)
  server_mod_meta("meta", rv, dataset_info)
  server_mod_methods("methods", dataset_info)
  server_mod_help("help")
}

shinyApp(ui = ui, server = server)
