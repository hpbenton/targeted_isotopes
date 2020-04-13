#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
require(MSnbase) || stop("Cannot open app, please make sure all dependencies are met")
require(xcms)|| stop("Cannot open app, please make sure all dependencies are met")

options(shiny.maxRequestSize=5000*1024^2)
# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    output$contents <- renderTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        # tryCatch(
        #     {
        #         df <- read.csv(input$file1$datapath,
        #                        header = input$header,
        #                        sep = input$sep,
        #                        quote = input$quote)
        #     },
        #     error = function(e) {
        #         # return a safeError if a parsing error occurs
        #         stop(safeError(e))
        #     }
        # )
        
        
        return(dir())
        
    })

})
