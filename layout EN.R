lng = "USA"
# layout for Shiny
title = "Emerging Therapies for COVID-19: the Value of More Clinical Research vs Treatment Implementation" 

# subtitle/explanations
Subtitle1  = h3("This tool interactively shows the most updated result of the manuscript Emerging Therapies for COVID-19: the Value of More Clinical Research vs Treatment Implementation.") 

Subtitle1.1 = h4("In the panel on the left, you can adjust parameter values for a COVID-19 therapy of interest. The results of the model will be shown in the results panel on the right. For details about the structure and assumptions of the model, see the panel called [About the tool].")

# Explanation 
Credits2  = h5("Research collaboration of Stijntje Dijk, Eline Krijkamp, Natalia Kunst, Cary Gross, John Wong, Myriam Hunink. Corresponding author: m.hunink@erasmusmc.nl") 

Credits3  = h5("Affiliations: Department of Epidemiology, Erasmus University Medical Center, Rotterdam, The Netherlands. Department of Health Management and Health Economics, Institute of Health and Society, Faculty of Medicine, Oslo, Norway. Section of General Internal Medicine, Yale University School of Medicine, New Haven, CT, USA. Division of Clinical Decision Making, Tufts Medical Center, Boston, USA. Center for Health Decision Science, Harvard University T.H. Chan School of Public Health, Boston, USA. ")

#### Add others as appropriate ###
Credits4 = h5("Shiny application created by Eline Krijkamp. Please report any bugs and inconsistencies to e.krijkamp@erasmusmc.nl") 

Instructions0 =  h3("Instructions:") 
Instructions1 =  p("Adjust the parameters of interest. The basecase code will update automatically. For the probabilistic sensitivity analysis click on the run botton.") 

Tab0_0 = "Drug inputs" 
Tab0_1 = "Model inputs"
Tab0_2 = "VOI inputs"
Tab0_3 = "About the Tool"


Tab1_1 = h3("Modify to see how the results vary due to changes in the parameter values.") 
Tab1_2 = h4("Specifyeffect size")
Tab1_3 = h4("Specify treatment costs")


# Order of importance 
# 1 Treatment effect = 1. effect therapy, effect type + 95% CI
# Cost therapy + uncertainty
# Length of stay
# Hoeveel patienten in de toemot COVID therapy geven
# Mortality control group 

# Effect size
Slider1    = "Hazard ratio of mortality" 
Slider1_1  = "95%-CI Hazard ratio" 
# Cost of the therapy
Slider2    = "Cost therapy per dose ($)" 
Slider2_1  = "95%-CI costs therapy per dose:" 
Slider2_2  = "Number of doses per patient"

# Length of stay
Slider3    = "Length of stay with therapy (days):" 
Slider3_1  = "CI length of stay"
Slider4    = "Length of stay without therapy (days):" 
Slider5    = "Discount rate (%)" 
Slider6    = "Number of PSA iteration:" 

Slider7   = "Willingness-to-pay ($/QALY)" 

Population0  = "Future patients"
Population1  = "Current patients"


Tab_instr1 = h4("Specify the length of stay") 
Tab_instr2 = h4("Change these values if you would like to consider different model assumptions in the calculations.")
Tab_instr3 = h5("Click run model after changing the number of PSA iterations. Note, the simulation takes a couple of minutes")
Tab_instr4 = h4("Change these values for the VOI analysis.")
Tab_instr5 = h5("Click run VOI to run the analysis, but make sure you first run the PSA analysis")


Main_Output1 =  tabPanel( title = "Results",
  br(),
  h3(strong("With this therapy:")),
  h4(strong("Basecase values:")),
  tableOutput("df_basecase"),
  h4(strong("Basecase cost-effectiveness results:")),
  tableOutput("df_CEA_basecase"),
  h5("Effect = QALY and Cost = US dollar"),
  h4(strong("PSA cost-effectiveness results:")),
  p("Click the [Run PSA] button to see the (changes in the) PSA results. Please note that this simulation takes some time."),
  plotOutput("plot_CE_PSA"),
  tableOutput("df_CEA_PSA"),
  h4(strong("Expected net benefit results:")),
  p("Click the [Run VOI] button to see (changes in the) VOI results. Please note that this simulation takes some time.")
  
)




License = h6("This software is licensed via the GNU General Public License 3.0.") 



# IDEA: Add error messages





