#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(lubridate)
library(rsconnect)
library(shinyBS)
library(matrixStats)
library(ggplot2)
library(scales)
library(reshape2)
library(nlme)
library(mgcv)
library(BCEA)
library(inlabru)
library(devtools)
library(tibble)
library(tidyverse)
library(ggpubr)
library(rms)
library(Rmisc)
library(fmsb)
library(remotes)
library(readxl)
library(plyr)
library(stats)
library(triangle)
library(EnvStats)
library(e1071)
library(meta)
library(metafor)
library(gridExtra)
library(here)
library(dplyr)
library(ellipse)
library(ggplot2)
library(lazyeval)
library(igraph)
library(ggraph)
library(knitr)
library(plyr)
library(stats)
library(diagram)
library(blscrapeR)
library(foreach)
library(mondate)
library(dampack)

library(darthtools) # a R package that contains tools frequently used by the DARTH workgroup



# Define UI for application that draws a histogram
ui <- fluidPage(
    
    theme = shinytheme(theme = "sandstone"),
    useShinyjs(),
    withMathJax(),
    
    tags$head(
        includeCSS("www/style_EK.css")
    ),
    

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 10,
                        value = 5)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           plotOutput("distPlotE")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
        
    })
    
    output$distPlotE <- renderPlot({
        # generate bins based on input$bins from ui.R
        # draw the histogram with the specified number of bins
        w    <- faithful[, 2]
        x <- rate_to_prob(w/100)
        bins <- seq(min(x), max(x), length.out = input$bins + 1)
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
