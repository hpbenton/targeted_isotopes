#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

require(shiny) || stop("Cannot find shiny package please install")
require(markdown) || stop("Cannot find all needed requirements please check you have everything installed")

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # App title ----
    titlePanel("Uploading Files"),
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        # Sidebar panel for inputs ----
        sidebarPanel(
            # Input: Select a file ----
            fileInput("files", "Choose mzML file(s)",
                      multiple = TRUE,
                      accept = "mzML"),
            
            # # Horizontal line ----
            # tags$hr(),
            
        ),
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Data file ----
            tableOutput("contents")
            
        )
    ),
    tabPanel(
        titlePanel("Peak Picking"),
        sidebarLayout(
            radioButtons("plotType", "Plot type",
                         c("Scatter"="p", "Line"="l")
            )
        )
    ),
))
