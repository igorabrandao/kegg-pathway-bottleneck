#***************************************************#
# Pipeline to generate the network visualization app #
#***************************************************#

# ---- IMPORT SECTION ----

# app.R #

#' This is the pipeline script to generate the shiny app
#' for networks visualizations
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(shiny)
library(shinydashboard)

#*******************************************************************************************#

# ---- SETTINGS SECTION ----

#*************************#
# Pipeline basic settings #
#*************************#

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

# Get the list of files
folder = paste0("./output/network/")
network_list <- list.files(path=folder, pattern='*.html')
menu_idx <- 0
file_idx <- 0

# Define UI object for app
ui <- dashboardPage(
  dashboardHeader(title = "Networks menu"),
  dashboardSidebar(sidebarMenu(id="mytabs", sidebarMenuOutput("menu"))),

  # Set the networks html files dynamically
  dashboardBody(
    tabItems <- lapply(network_list, function(file) {
      file_idx <<- file_idx + 1
      return(tabItem(tabName = paste0("dashboard", file_idx), includeHTML(paste0(folder, network_list[file_idx]))))
    })
  )
)

# Define server logic required to show the graphs
server <- function(input, output, session) {
  output$menu <- renderMenu({

    # Set the menu text dynamically
    menu_list <- lapply(network_list, function(file) {
      menu_idx <<- menu_idx + 1
      menuText <- onlyNumber(file)
      return(menuItem(menuText[2], tabName = paste0("dashboard", menu_idx), icon = icon("project-diagram")))
    })

    sidebarMenu(menu_list)
  })
  isolate({updateTabItems(session, "mytabs", "dashboard1")})
}

# Run the APP!
shinyApp(ui, server)
