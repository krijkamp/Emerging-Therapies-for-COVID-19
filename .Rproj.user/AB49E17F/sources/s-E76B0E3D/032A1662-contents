#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


# Load the packages 
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(lubridate)
library(rsconnect)
library(shinyBS)
if (!require('pacman')) install.packages('pacman'); library(pacman) 

# use this package to conveniently install other packages
p_load("matrixStats", "ggplot2", 
       "scales", "reshape2",
       "nlme", "mgcv", "BCEA", 
       "inlabru", "devtools",
       "tibble", "tidyverse", "ggpubr",
       "rms", "europepmc", "Rmisc",
       "fmsb", "remotes", "readxl", "plyr", "stats", "triangle",
       "EnvStats", "e1071", "meta","metafor", "gridExtra", 
       "here", "dplyr","ellipse", "ggplot2", "lazyeval", 
       "igraph", "ggraph","knitr", "plyr", "stats", "diagram",
       "triangle", "HMDHFDplus", "blscrapeR", "here", "gridExtra", "foreach", "mondate",  "parallel")


p_load_gh("DARTH-git/dampack") # coding framework to construct model-based cost-effectiveness analysis in R
#p_load_gh("DARTH-git/darthpack") # package for analyzing and visualizing the health economic outputs of mathematical models
p_load_gh("DARTH-git/darthtools") # a R package that contains tools frequently used by the DARTH workgroup


if (!require("remotes"))
  install.packages("remotes")
remotes::install_github("jcrodriguez1989/shinyParallel")
library("shinyParallel")

rm(list=ls())

# Load the data = output from the main model
load("output/l_param_trt.rda")  # 
load("output/l_names.rda")


# Select for each drug one list -> using the drug name
l_param_shiny <- l_param_trt[c(l_names$v_names_trt_report)]
names(l_param_shiny) <- l_names$v_names_trt_report_ful


jscode <- "shinyjs.refresh = function() { history.go(0); }"




# Define UI for application that draws a histogram
ui <- function(req) {
  Lang = "EN"
  source(paste0("layout ", Lang, ".R"), local= T)
  source(paste0("values ", Lang, ".R"), local= T)
  source(paste0("documentation ", Lang, ".R"), local= T)
  #source("documentation.R") #@EK changed for language hope that works
  
fluidPage(
    tags$head(
      includeCSS("www/style_EK.css")
    ),
    theme = shinytheme(theme = "sandstone"),
    useShinyjs(),
    withMathJax(),
          

          h1(paste(title, lng, sep =", "), .noWs ="outside"),
                
                # Application title
                Credits2,
                Credits3,
                Credits4,
          br(),
                Subtitle1,
                Subtitle1.1,
          hr(),
                Instructions0,
                Instructions1,
                # Sidebar with a slider input for number of bins 
                sidebarLayout(
                  sidebarPanel(
                    tabsetPanel(id = "inputs",
                                type = "tabs",
                                tabPanel(title = Tab0_0,
                                         Tab1_1,
                                         Tab1_2,
                                         numericInput("hr_D_Trt_timespan1", 
                                                      Slider1, 
                                                      value = n_eff_size,
                                                      step = 0.001),
                                         
                                         sliderInput("Uncertainty around effect size",
                                                     Slider1_1,
                                                     min = 0,
                                                     max = 2,
                                                     value = c(ci_eff_size),
                                                     step = 0.001),
                                         Tab1_3,
                                         
                                         numericInput("c_Trt",
                                                     Slider2,
                                                     value = n_costs,
                                                     step = 1),
                                         
                                         sliderInput("ci_c_Trt",
                                                     Slider2_1,
                                                     min = 0,
                                                     max = 1e6,
                                                     value = c(ci_costs),
                                                     step = 1),
                                         
                                         numericInput("n_Trt",
                                                      Slider2_2,
                                                      value = n_Trt,
                                                      step = 1),
                                  
                                         
                                         Tab_instr1,
                                         
                                         sliderInput("LOS_Trt",
                                                     Slider3,
                                                     min = 0,
                                                     max = 100,
                                                     value = LOS_Trt,
                                                     step = 1),
                                         
                                         sliderInput("LOS_noTrt",
                                                     Slider4,
                                                     min = 0,
                                                     max = 100,
                                                     value = LOS_noTrt,
                                                     step = 1)
                                ),
                                tabPanel(title = Tab0_1,
                                         Tab_instr2,
                                         sliderInput("r_discount",
                                                     Slider5,
                                                     min = 0,
                                                     max = 10,
                                                     value = r_discount,
                                                     step = 0.5),
                                         
                                         sliderInput("n_wtp",
                                                      Slider7,
                                                     min = 0,
                                                     max = 3*1e5,
                                                     step = 1e3,
                                                     value = n_wtp),

                                         Tab_instr3,
                                         sliderInput("n_iter",
                                                     Slider6,
                                                     min = 0,
                                                     max = 10000,
                                                     value = n_iter,
                                                     step = 50)
                                         
                                ),
                                tabPanel(title = Tab0_2,
                                         Tab_instr4,

                                         numericInput(inputId = "n_H_current",
                                                      label = Population1,
                                                      value = n_H_current),
                                         
                                         numericInput(inputId = "n_H_future",
                                                      label = Population0,
                                                      value = n_H_future),
                                         
                                         actionButton("RunVOI", "Run VOI")

                                )
                                
                                
                    ),
                    
                    actionButton("RunPSA", "Run PSA")
                    ),
               mainPanel (
                 
                tabsetPanel(id = "outputs", type = "tabs",
                  # Display outputs
                             Main_Output1,
                  tabPanel(title = Tab0_3,
                           Documentation0,
                           Documentation1,
                           Documentation2,
                           Documentation2_1,
                           Documentation2_2,
                           Documentation3,
                           Documentation3_1,
                           Documentation4,
                           Documentation4_1,
                           Documentation5,
                           Documentation5_1,
                           Documentation6,
                           Documentation6_1,
                           Documentation7,
                           Documentation7_1,
                           Documentation8,
                           Documentation8_1,
                           Documentation9,
                           Documentation10,
                           Documentation11,
                           Documentation12,
                           Documentation13,
                           Documentation14,
                           Documentation14_1,
                           Documentation15,
                           Documentation15_1,
                           Documentation16,
                           Documentation16_1,
                           Documentation17,
                           Documentation17_1,
                           Documentation18,
                           Documentation18_1
                           ), 
                  tabPanel(title = "Disclaimer",
                           DocumentationDisclaimer1,
                           DocumentationDisclaimer2,
                           DocumentationDisclaimer3,
                           DocumentationDisclaimer4,
                           DocumentationDisclaimer5,
                           DocumentationDisclaimer6),
                  tabPanel(title = "Acknowledgement",
                           DocumentationAck1,
                           DocumentationAck2,
                           DocumentationAck3)
                  )
                )
                  
                ),
                hr(),
                License
)
}
# Define server logic required to draw a histogram
server <- function(input, output) {
  source("functions/00_general_mycolour.R")
  source("functions/00_functions_darthtools.R")
  source("functions/01_model_input_data_hospitalization_functions.R") # functions for the model
  source("functions/01_model_inputs_functions.R")  # functions for the data
  source("functions/02_decision_model_functions.R") # functions for model structure & to run the model 
  source("functions/02_decision_model_plot_functions.R") # plot functions
  source("functions/02_decision_model_calcout_functions.R") # the entire model
  Lang <- read.csv("changes.csv")$x
  if (is.na(Lang)) Lang = "EN"

  source(paste0("layout ", Lang, ".R"), local= T)
  source(paste0("values ", Lang, ".R"), local= T)
  source(paste0("documentation ", Lang, ".R"), local= T)
  
  write.csv("EN", "changes.csv")
  
  
  # Load the data = output from the main model
  load("output/l_param_trt_basecase.rda")  # 
  load("output/l_names.rda")
  load("output/l_m_Parameters.rda")  
  load("output/l_input_general.rda")
  
  
  data <- reactive({
    validate(
      need(input$c_Trt > input$ci_c_Trt[1] & input$c_Trt < input$ci_c_Trt[2] , "Please select a data set")
    )
  })
  
  
  # Select for each drug one list -> using the drug name
  l_param_shiny <- l_param_trt_basecase[c(l_names$v_names_trt_report)]
  names(l_param_shiny) <- l_names$v_names_trt_report_ful
  
  l_m_Parameters_shiny <- l_m_Parameters[c(l_names$v_names_trt_report)]
  names(l_m_Parameters_shiny) <- l_names$v_names_trt_report_ful
  l_param <- l_param_shiny$Dexamethasone
  
  # run the basecase function
  output$df_basecase <- renderTable({
    #update the values 
    l_param$hr_D_Trt_timespan1_novent<- input$hr_D_Trt_timespan1
    l_param$hr_D_Trt_timespan1_vent  <- input$hr_D_Trt_timespan1
    l_param$c_Trt_private            <- input$c_Trt
    l_param$c_Trt_public             <- input$c_Trt
    l_param$n_Trt                    <- input$n_Trt
    l_param$LOS_Trt                  <- input$LOS_Trt
    l_param$LOS_noTrt                <- input$LOS_noTrt 
    l_param$d_c                      <- input$r_discount/100
    l_param$d_e                      <- input$r_discount/100
    
    df_ce   <- calculate_cea_output_VOI_COVID(l_param, 
                                            n_wtp = input$n_wtp)
    df_ce %>% arrange(desc(Strategy))
    
  })
  
  # run the calculate ICER function
  output$df_CEA_basecase<- renderTable({
    l_param$hr_D_Trt_timespan1_novent<- input$hr_D_Trt_timespan1
    l_param$hr_D_Trt_timespan1_vent  <- input$hr_D_Trt_timespan1
    l_param$c_Trt_private            <- input$c_Trt
    l_param$c_Trt_public             <- input$c_Trt
    l_param$n_Trt                    <- input$n_Trt
    l_param$LOS_Trt                  <- input$LOS_Trt
    l_param$LOS_noTrt                <- input$LOS_noTrt
    l_param$d_c                      <- input$r_discount/100
    l_param$d_e                      <- input$r_discount/100
    
    df_ce <- calculate_cea_output_VOI_COVID(l_param, 
                                            n_wtp = input$n_wtp)
    
    df_cea_QALY <- calculate_icers(cost       = df_ce$Cost,
                                   effect     = df_ce$Effect,
                                   strategies = df_ce$Strategy)
    df_cea_QALY %>% 
      arrange(desc(Strategy))
  })
  
 

  

  
  

  

  
  # Trick file date creation update
  onStop(function() {
    
    # File name
    p <- paste0(getwd(), "/app.R")
    
    # Update file 'date creation'
    Sys.setFileTime(p, now())
    print("test")
  
  }) # onStop
}

# Run the application 
shinyApp(ui = ui, server = server)
