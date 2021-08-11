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
  source(paste0("layout EN.R"), local= T)
  source(paste0("values EN.R"), local= T)
  source(paste0("documentation EN.R"), local= T)
  #source("documentation.R") #@EK changed for language hope that works
  

  
fluidPage(
    tags$head(
      includeCSS("www/style_EK.css")
    ),
    theme = shinytheme(theme = "sandstone"),
    useShinyjs(),
    withMathJax(),
          

          h1(paste(title, sep =", "), .noWs ="outside"),
                
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
                                         numericInput("hr_D_Trt_timespan1_novent", 
                                                      Slider1, 
                                                      value = n_eff_size,
                                                      step = 0.001,
                                                      min = 0),
                                         
                                         sliderInput("Uncertainty around effect size",
                                                     Slider1_1,
                                                     min = 0,
                                                     max = 3,
                                                     value = c(ci_eff_size),
                                                     step = 0.001),
                                         Tab1_2_1,
                                         numericInput("hr_D_Trt_timespan1_vent", 
                                                      Slider1_2, 
                                                      value = n_eff_size_vent,
                                                      step = 0.001,
                                                      min = 0),
                                         
                                         sliderInput("Uncertainty around effect size",
                                                     Slider1_3,
                                                     min = 0,
                                                     max = 3,
                                                     value = c(ci_eff_size_vent),
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
                                                     max = 73,
                                                     value = LOS_Trt,
                                                     step = 1),
                                         
                                         sliderInput("LOS_noTrt",
                                                     Slider4,
                                                     min = 0,
                                                     max = 73,
                                                     value = LOS_noTrt,
                                                     step = 1),
                                         Tab_instr1_1,
                                         
                                         sliderInput("p_IC",
                                                     Slider4_1,
                                                     min = 0,
                                                     max = 1,
                                                     value = 0.165,
                                                     step = 0.01)
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
                                                     max = 1e3,
                                                     value = n_iter,
                                                     step = 100)
                                         
                                ),
                                tabPanel(titl = Tab0_1_1,
                                         h3("Triangular distributions"),
                                         h5("Make sure the mean fits in the CI range"),
                                         
                                         h4("Being hospitalized "),
                                         sliderInput("u_H",
                                                     "Mean ",
                                                     min = 0,
                                                     max = 1, 
                                                     value = u_H,
                                                     step = 0.01),
                                         
                                         sliderInput("ci_u_H",
                                                     "upper and lower limit",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(ci_u_H),
                                                     step = 0.01),
                                         
                                         h4("Being in the intensive care unit "),
                                         sliderInput("u_I",
                                                     "Mean",
                                                     min = 0,
                                                     max = 1, 
                                                     value = u_I,
                                                     step = 0.01),
                                         
                                         sliderInput("ci_u_I",
                                                     "upper and lower limit",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(ci_u_I),
                                                     step = 0.01),
                                         
                                         h3("Beta distributions"),
                                         h4("After recovery from hospitalization"),
                                         sliderInput("u_R_H",
                                                     "mean",
                                                     min = 0,
                                                     max = 1, 
                                                     value = u_R_H,
                                                     step = 0.01),
                                         
                                         numericInput("sd_u_R_H",
                                                      "Standard deviatie",
                                                      value = sd_u_R_H,
                                                      step = 0.01),
                                         
                                         h4("After recovery from the ICU"),
                                         
                                         sliderInput("u_R_IC",
                                                     "Mean",
                                                     min = 0,
                                                     max = 1, 
                                                     value = u_R_IC,
                                                     step = 0.01),
                                         
                                         numericInput("sd_u_R_IC",
                                                      "Standard deviatie",
                                                      value = sd_u_R_IC,
                                                      step = 0.01),
                                         Tab_instr3.1,
                                         ), 
                                
                                tabPanel(title = Tab0_2,
                                         Tab_instr4,

                                         numericInput(inputId = "n_H_current",
                                                      label = Population1,
                                                      value = n_H_current),
                                         Population1_1,
                                         
                                         numericInput(inputId = "n_H_future",
                                                      label = Population0,
                                                      value = n_H_future),
                                         Population0_1,
                                         
                                         actionButton("RunVOI", "Run VOI")

                                )
                                
                                
                    ),
                    
                    actionButton("RunPA", "Run PA")
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
  source("functions/06_VOI_functions.R") # load VOI functions
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
    l_param$hr_D_Trt_timespan1_novent<- input$hr_D_Trt_timespan1_novent
    l_param$hr_D_Trt_timespan1_vent  <- input$hr_D_Trt_timespan1_vent
    l_param$c_Trt_private            <- input$c_Trt
    l_param$c_Trt_public             <- input$c_Trt
    l_param$n_Trt                    <- input$n_Trt
    l_param$LOS_Trt                  <- input$LOS_Trt
    l_param$LOS_noTrt                <- input$LOS_noTrt 
    l_param$p_IC                     <- input$p_IC
    l_param$d_c                      <- input$r_discount/100
    l_param$d_e                      <- input$r_discount/100
    l_param$u_H                      <- input$u_H
    l_param$u_I                      <- input$u_I
    l_param$u_R_H                    <- input$u_R_H
    l_param$u_R_IC                   <- input$u_R_IC
    
    df_ce   <- calculate_cea_output_VOI_COVID(l_param, 
                                            n_wtp = input$n_wtp)
    df_ce %>% arrange(desc(Strategy))
    
  })
  
  # run the calculate ICER function
  output$df_CEA_basecase<- renderTable({
    l_param$hr_D_Trt_timespan1_novent<- input$hr_D_Trt_timespan1_novent
    l_param$hr_D_Trt_timespan1_vent  <- input$hr_D_Trt_timespan1_vent
    l_param$c_Trt_private            <- input$c_Trt
    l_param$c_Trt_public             <- input$c_Trt
    l_param$n_Trt                    <- input$n_Trt
    l_param$LOS_Trt                  <- input$LOS_Trt
    l_param$LOS_noTrt                <- input$LOS_noTrt
    l_param$p_IC                     <- input$p_IC
    l_param$d_c                      <- input$r_discount/100
    l_param$d_e                      <- input$r_discount/100
    l_param$u_H                      <- input$u_H
    l_param$u_I                      <- input$u_I
    l_param$u_R_H                    <- input$u_R_H
    l_param$u_R_IC                   <- input$u_R_IC
    
    df_ce <- calculate_cea_output_VOI_COVID(l_param, 
                                            n_wtp = input$n_wtp)
    
    df_cea_QALY <- calculate_icers(cost       = df_ce$Cost,
                                   effect     = df_ce$Effect,
                                   strategies = df_ce$Strategy)
    df_cea_QALY %>% 
      arrange(desc(Strategy))
  })
  
 

  

  
  observeEvent(input$RunPA, {  #
    # add all the action that are needed before the model needs to start running
    
    withProgress(message = 'Performing probabilistic analysis', 
                 detail = "After the PA, making the figures takes some time as well", value = 0, {
      
      # select the treatment of interest
      m_Parameters <<- l_m_Parameters_shiny$Dexamethasone[1:input$n_iter, ]
      v_output <- c("LY", "QALY", "Costs")            # Vector of output names
      n_str <- 2
      m_C_psa <- m_E_psa <- matrix(NA, ncol = 2, nrow = input$n_iter)
      colnames(m_C_psa) <- colnames(m_E_psa) <- c("notrt", "trt")
      
      #Initiate list
      l_df_par <- l_m_output_par <-  list()
      
      m_output_par <- matrix(NA, 
                             ncol = as.numeric(length(v_output)) * as.numeric(n_str), 
                             nrow = input$n_iter, 
                             dimnames = list(paste("Iteration", 1:input$n_iter, sep = " "), 
                                             (paste(rep(v_output, n_str), rep(c("notrt", "trt"), 
                                                                              each = length(v_output)), sep = " "))))
     
      
        # Parameters with a distribution - sample fist 
          v_c_trt <- runif(n = input$n_iter, 
                           min =  input$ci_c_Trt[1], 
                           max =  input$ci_c_Trt[2])
          
          # Utility values 
          v_u_H    <- rtriangle(n = input$n_iter,
                                a = input$ci_u_H[1],
                                b = input$ci_u_H[2],
                                c = input$u_H)
          
          v_u_I     <- rtriangle(n = input$n_iter,
                                 a = input$ci_u_I[1],
                                 b = input$ci_u_I[2],
                                 c = input$u_I)
          
          v_u_R_H  <- rbeta(n = input$n_iter,
                            shape1 = beta_params(input$u_R_H, input$sd_u_R_H)$alpha,
                            shape2 = beta_params(input$u_R_H, input$sd_u_R_H)$beta)
         
          v_u_R_IC <- rbeta(n = input$n_iter,
                             shape1 = beta_params(input$u_R_IC, input$sd_u_R_IC)$alpha,
                             shape2 = beta_params(input$u_R_IC, input$sd_u_R_IC)$beta)
          
        
        # Replace baseline items in a list
        m_Parameters[, "hr_D_Trt_timespan1_novent"]<- input$hr_D_Trt_timespan1_novent
        m_Parameters[, "hr_D_Trt_timespan1_vent"]  <- input$hr_D_Trt_timespan1_vent
        m_Parameters[, "c_Trt_private"]            <- v_c_trt
        m_Parameters[, "c_Trt_public"]             <- v_c_trt
        m_Parameters[, "n_Trt"]                    <- input$n_Trt
        m_Parameters[, "LOS_Trt"]                  <- input$LOS_Trt
        m_Parameters[, "LOS_noTrt"]                <- input$LOS_noTrt
        m_Parameters[, "p_IC"]                     <- input$p_IC
        m_Parameters[, "d_c"]                      <- input$r_discount/100
        m_Parameters[, "d_e"]                      <- input$r_discount/100
        m_Parameters[, "u_H"]                      <- v_u_H
        m_Parameters[, "u_I"]                      <- v_u_I
        m_Parameters[, "u_R_H"]                    <- v_u_R_H
        m_Parameters[, "u_R_IC"]                   <- v_u_R_IC
        
        

        for (g in 1:input$n_iter){
         
    
          
          l_param_psa <- as.list(m_Parameters[g, ])
          l_param_psa <- c(l_param_psa, l_input_general) # Combine with the general input

          
          # Run the model
          l_out_temp <- calculate_cea_output_VOI_COVID(l_param_psa, n_wtp = l_param_psa$wtp, verbose = FALSE)
          
          df_ce <- c(l_out_temp$Cost, 
                     l_out_temp$Effect, 
                     l_out_temp$LY,
                     l_out_temp$NMB,
                     l_out_temp$NHB)
          
        
        
        # Save the output
        m_output_par[g, "Costs notrt"]   <- df_ce[ 1]
        m_output_par[g, "Costs trt"]     <- df_ce[ 2]
        
        m_output_par[g, "QALY notrt"]    <- df_ce[ 3]
        m_output_par[g, "QALY trt"]      <- df_ce[ 4]
        
        m_output_par[g, "LY notrt"]      <- df_ce[ 5]
        m_output_par[g, "LY trt"]        <- df_ce[ 6]
        
        # Add storing the values of each iterations 
        m_C_psa[g, "notrt"] <- m_output_par[g, "Costs notrt"]
        m_C_psa[g, "trt"]   <- m_output_par[g, "Costs trt"]
        
        m_E_psa[g, "notrt"] <- m_output_par[g, "QALY notrt"]
        m_E_psa[g, "trt"]   <- m_output_par[g, "QALY trt"]
        
     
        incProgress(amount = 1/input$n_iter) # Update the progress bar   
        Sys.sleep(0.25)  # pretending to execute code
        
      
      }
    
      
  })
    
    # Do the CEA calculation 
    output$df_CEA_PSA <- renderTable({

      
      df_cea_QALY_PSA <- calculate_icers(cost   =     colMeans(m_C_psa),
                                         effect     = colMeans(m_E_psa),
                                         strategies = l_out_temp$Strategy)
      df_cea_QALY_PSA
      
      df_cea_QALY_PSA %>% arrange(desc(Strategy))  # Have trt first row

    })
    
    output$plot_CE_PSA <- renderPlot({
      
      plot_CE_PSA <- make_psa_obj(cost  = as.data.frame(m_C_psa), 
                                  effectiveness  = as.data.frame(m_E_psa),
                                  parameters = m_Parameters[1:input$n_iter, ],
                                  strategies = l_out_temp$Strategy)
      plot(plot_CE_PSA) + 
        labs(title     = paste("Cost-effectiveness plane"), 
             subtitle   = paste("intervention vs control", sep = " ")) + 
        xlab("Effectiveness (in QALY)") +
        ylab("Costs ($)") +
        geom_vline(xintercept = 0, color = my_darkgray, size = 0.6) +
        geom_hline(yintercept = 0, color = my_darkgray, size = 0.6) +
        scale_color_manual(values = c(my_lightred, my_green, my_black)) + 
        scale_fill_manual(values = c(my_lightred, my_green, my_black))

      
    })
    
    
  })
  
  observeEvent(input$RunVOI, {
    withProgress(message = 'Performing value of information analysis', value = 0, {
      
      v_wtp <- seq(0, 3*1e5, 1e3)
      l_pop <- input$n_H_year * mean(m_Parameters[, "p_IC_notrt"])

      # calculate population EVPI for QALY
      l_out_evpi_QALY <- evpi(v.wtp = v_wtp,
                              m.e   = m_E_psa, 
                              m.c   = m_C_psa, 
                              pop   = l_pop) 
      
  
      l_out_evpi_QALY$EVPI[l_out_evpi_QALY$WTP == 100000] # print QALY result for a specific WTP
      
      # store in a list
      l_evpi_value <- cbind(WTP       = l_out_evpi$WTP, 
                          EVPI_LY   = l_out_evpi_LY$EVPI, 
                          EVPI_QALY = l_out_evpi_QALY$EVPI) 
    }) 
    
    #output$plotVOI <- renderPlot({})  
  
    
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
