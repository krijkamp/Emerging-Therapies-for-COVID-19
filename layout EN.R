
# layout for Shiny
title = "Emerging Therapies for COVID-19: the value of performing clinical trials" 

# subtitle/explanations
Subtitle1  = h3("This tool interactively shows the most updated result of the manuscript Emerging Therapies for COVID-19: the value of performing clinical trials.") 

Subtitle1.1 = h4("In the panel on the left, you can adjust parameter values for a COVID-19 therapy of interest. The results of the model will be shown in the results panel on the right. For details about the structure and assumptions of the model, see the panel called [About the tool] and peer-reviewed publication [ADD CITATION].")

# Explanation 
Credits2  = h5("Research collaboration of Stijntje Dijk, Eline Krijkamp, Natalia Kunst, Cary Gross, John Wong, Myriam Hunink. Corresponding author: m.hunink@erasmusmc.nl") 

Credits3  = h5("Affiliations: Department of Epidemiology, Erasmus University Medical Center, Rotterdam, The Netherlands. Department of Health Management and Health Economics, Institute of Health and Society, Faculty of Medicine, Oslo, Norway. Section of General Internal Medicine, Yale University School of Medicine, New Haven, CT, USA. Division of Clinical Decision Making, Tufts Medical Center, Boston, USA. Center for Health Decision Science, Harvard University T.H. Chan School of Public Health, Boston, USA. ")

#### Add others as appropriate ###
Credits4 = h5("Shiny application created by Eline Krijkamp. Please report any bugs and inconsistencies to e.krijkamp@erasmusmc.nl") 

Instructions0 =  h3("Instructions:") 
Instructions1 =  p("Adjust the parameters of interest. The basecase code will update automatically. To run the probabilistic analysis, click on the [run PA] botton.") 

Tab0_0 = "Drug inputs" 
Tab0_1 = "Model inputs"
Tab0_1_1 = "Utility inputs"
Tab0_2 = "VOI inputs"
Tab0_3 = "About the Tool"


Tab1_1 = h3("Modify to see how the results vary due to changes in the parameter values.") 
Tab1_2 = h4("Specify effect size")
Tab1_2_1 = h5("when on mechanical ventilation")
Tab1_3 = h4("Specify treatment costs")


# Order of importance 
# 1 Treatment effect = 1. effect therapy, effect type + 95% CI
# Cost therapy + uncertainty
# Length of stay
# Hoeveel patienten in de toemot COVID therapy geven
# Mortality control group 

# Effect size
Slider1    = "Hazard ratio of mortality without ventilation(treatment compared to no treatment)" 
Slider1_1  = "95%-CI Hazard ratio" 
Slider1_2    = "Hazard ratio of mortality WITH ventilation (treatment compared to no treatment)" 
Slider1_3  = "95%-CI Hazard ratio WITH  ventilation"

# Cost of the therapy
Slider2    = "Cost therapy per dose ($)" 
Slider2_1  = "95%-CI costs therapy per dose:" 
Slider2_2  = "Number of doses per patient"

# Length of stay
Slider3    = "Length of stay with therapy (days):" 
Slider3_1  = "CI length of stay"
Slider4    = "Length of stay without therapy (days):" 
Slider4_1  = "Probability of being admitted to the ICU"
Slider5    = "Discount rate (%)" 
Slider6    = "Number of PA iteration:" 

Slider7   = "Willingness-to-pay ($/QALY)" 

Population0  = "Future patients"
Population0_1 = h5("The number of patients expected to be hospitalized after the trial results are available")
Population1  = "Current patients"
Population1_1 = h5("Expected number of patients to be hospitalized in the USA while awaiting trial results and their implementation over a 3 months period")


Tab_instr1   = h4("Specify the length of stay") 
Tab_instr1_1 = h4("Disease severity: ICU admission")
Tab_instr2   = h4("Change these values if you would like to consider different model assumptions in the calculations.")
Tab_instr3   = h5("Click [run PA] after changing the number of PA iterations. Note: the simulation takes a couple of minutes")
Tab_instr3.1 = h5("NOTE: these utiliy values are only evaluated in the basecase analysis. For the PA analysis, the distribution as specified in the manuscript is used")
Tab_instr4   = h4("Change these values for the VOI analysis.")
Tab_instr5   = h5("Click [run VOI] to run the analysis, but make sure you first run the PA analysis")


Main_Output1 =  tabPanel( title = "Results",
  br(),
  h3(strong("With this therapy:")),
  h4(strong("Basecase values:")),
  h5("Cost = US dollar, Effect = QALY, LY = Life years, NMB: Net monetary benefit for QALYs, NHB: Net health benefit for QALYs, NMB_LY: Net monetary benefit for LYs, NHB_LY: Net health benefit for LYs"),
  tableOutput("df_basecase"),
  h4(strong("Basecase cost-effectiveness results:")),
  h5("Effect = QALY and Cost = US dollar"),
  tableOutput("df_CEA_basecase"),
  h5("Lengend: trt: treatment, notrt: no treatment, Inc_Cost: Incremental costs, Inc:Effect: Incremental effects in QALY, QALY:Quality adjusted life years, ICER: Incremental cost effectivenes ratio, ND: no-dominated, D:dominated, ED: extended dominated"),
  h4(strong("PA cost-effectiveness results:")),
  p("Click the [Run PA] button to see the (changes in the) probabilistic analysis results. Please note that this simulation takes some time."),
  plotOutput("plot_CE_PSA"),
  p(""),
  tableOutput("df_CEA_PSA"),
  h4(strong("Expected net benefit results:")),
  textOutput("test"),
  p("Click the [Run VOI] button to see (changes in the) VOI results. Please note that this simulation takes some time.")
  
)




License = h6("This software is licensed via the GNU General Public License 3.0.") 



# IDEA: Add error messages





