Documentation0 <- h1("Documentation")

Documentation1 <- p("The goal of this tool is to clearly present the results of the performed study and to transparently communicate the data used and assumptions made. We also would like to provide an easy way to see how changes in parameters values might influence the model results and consequently maybe even the conclusion.")

Documentation2 <- p("Below, we elaborate on the drug inputs, model parameters, and calculations.")

Documentation2_1 <-  p("Note: This tool is for illustrative purposes. Use for decision-making or policy planning is at your own risk. Please also read the Disclaimer before using it.")

Documentation2_2 <- h2("Drug Inputs")

Documentation3 <- h4("Hazard rate ratio")
Documentation3_1 <- p("In general, a hazard rate ratio (HRR) is a measure of an effect of an intervention on an outcome of interest over time. In our case it is the treatment effect of the drug on the survival of hospitalized covid patients. An HRR of 1 means that at any particular point in time, the event rate are the same in both the control and the intervention group. In other words, the treatment is not different from the control group of the outcome dying. A HRR of 0.5 indicates that at any particular point in time, half as many patient in the treatment group experienced the event (dying) compared to the control group. This means the treatment has a protective effect on dying compared to the control. While, when the treatment is causing more deaths compared the control group the HRR is above 1. For example, an HRR of 2 indicated that at any particual point in time, twice as many patients in the treatment group experience the event (dying) compared to the control.")

Documentation4   <- h4("95%-CI hazard rate ratio")
Documentation4_1 <- p("This is the 95%-confidence interval of the hazard rate ratio of the treatment effect (treatment vs control). These parameter values are used in the sensitivity analysis of the model to incorporate the uncertainty around the estimated HRR based on the current treatment effect from a trial or meta-analysis of trials.")

Documentation5   <- h4("Cost therapy per dose")
Documentation5_1 <- p("The mean expected costs of a dose of the drug in US dollar. Please note, users might like to calculate the mean price of the drug based on the % of users having a private vs a public insurance. In our simulations we estimate the mean drug costs using data suggesting that 66,2% of the patient is privately insured. (source: https://www.census.gov/library/publications/2019/demo/p60-267.html.) ")

Documentation6   <- h4("Number of doses per patient")
Documentation6_1 <- p("Mean expected number of required vials, estimated as total used over period of 6 days")

Documentation7   <- h4("Length of stay with therapy")
Documentation7_1 <-p("Mean length of stay in the hospital in days for a patient under treatment. This parameter is used to estimate the utility and cost of patients in the hospital. The maximum length of stay in the hospital is 73 days. This hospitalization period is modeled in the first cycle of the model. ")
.
Documentation8   <- h4("Length of stay without therapy")
Documentation8_1 <-p("Mean length of stay in the hospital in days for a COVID patient under usual care. This parameter is used to estimate the utility and cost of patients in the hospital. The maximum length of stay in the hospital is 73 days. This hospitalization period is modeled in the first cycle of the model.")

Documentation9<- h2("Model Inputs")

Documentation10 <- h3("General input - fixed")
Documentation11 <- p("The general input parameters define the settings for the model. Part of these parameters values are fixed in this shiny application. However, usefull to know what they mean. For our strategies, we use the terms trt and notrt to indicate the treatment strategy or care as usual respectively. Our health states are indicated using the names Hospitalized, Recovered from Hospital ward, Recovered from ICU, Dead and in short they are indicated using H, R_H, R_IC and R, respectively. We make use of a cycle length of 73 days and assume a year had 365 days. In our simulation we start with hospitalized patient, since this is also the place were patients are treated with the drugs for COVID-19.") 

Documentation12 <- p("We simulate a life time horizon which is defined as a maximum age of 120. Users are not able to change the time horizon of interest.") 


Documentation13 <- h3("General input - flexible")

Documentation14   <- h4("Discount rate")
Documentation14_1 <- p("Since we are interested in the USA setting, we used the recommended discounting factor of 3% for both costs and effects.")

Documentation15   <- h4("Number of PSA iterations")
Documentation15_1 <- p("The PSA is technique used to propagate uncertainty from model inputs to model outcomes. This is done by running the model using different sets of model parameters based on their distribution. The number of PSA iterations indicated how many different parameter sets are generate to run the model with. More iterations give a more stable result, but comes with the costs of computational time. ")

Documentation16   <- h4("Willingness-to-pay ($/QALY)")
Documentation16_1 <- p("A threshold that represents what the decision maker or society is willing-to-pay for a unit of health outcome. The threshold is expressed in monetary units per health outcome.")

Documentation17 <- h4("Current patients")
Documentation17_1 <- p("Number of patients expected to be hospitalized overall while another RCT would be conducted. Assumed to be 3 months in our analysis. Source: IHME projections. ")

Documentation18 <- h4("Future patients")
Documentation18_1 <- p("Number of patients expected to be hospitalized excluding current patients. Source: IHME projections.")

Documentation19 = h4("Limitations")

Documentation20 = p("Limitation description")




# Disclaimer section

DocumentationDisclaimer1 = h3("Disclaimer")

DocumentationDisclaimer2 = h4("The presented model on this website is the result of carefully conducted research. However, the corresponding code and manuscript have not been through peer-review yet.")

DocumentationDisclaimer3 = h5("Before you use this app make sure that you read the model description as well as the peer reviewed paper [link will be added once ready](ADD LINK ), which will help you to interpret its output. For questions, please don't hesitate to contact the authors (see title page) and please report any inconsistencies or doubts directly. Thank you.")

DocumentationDisclaimer4 = p("As a model can never replace clinical judgment, it can only be used as a decision-support tool. This decision tool should be used exclusively by health care professionals as a complementary tool to estimate urgency of procedures in order to triage effectively and ethically. Any responsibility for using this model and its results will rest solely by the health care professional using the model. Using it you should understand and agree that this site is not responsible or liable for any claim, loss, or damage resulting from its use. While we try to keep the information on the site as accurate as possible, we disclaim any warranty concerning its accuracy, timeliness, and completeness, and any other warranty, express or implied, including warranties of merchantability or fitness for a particular purpose. ")

DocumentationDisclaimer5 = p("This site is not an attempt to practice medicine or provide specific medical advice, nor does the use of the site establish a doctor-patient relationship. For medical treatment or answers to personal questions, we strongly encourage you to consult with a qualified health care provider. For advice about your own care, please ask your doctor.")

DocumentationDisclaimer6 = p("You assume full responsibility for using the information on this site, and you understand and agree that this site is not responsible or liable for any claim, loss, or damage resulting from its use. While we try to keep the information on the site as accurate as possible, we disclaim any warranty concerning its accuracy, timeliness, and completeness, and any other warranty, express or implied, including warranties of merchantability or fitness for a particular purpose. We do not warrant that access to the site will be error- or virus-free.")


# Documentation Acknowledgment 
DocumentationAck1 = h3("Acknowledgment")
DocumentationAck2 = h4("The authors would like to thank the Society for Medical Decision Making COVID-19 Decision and the Gordon and Betty Moore Foundation for their support through the COVID-19 Decision Modeling Initiative")

DocumentationAck3 = p("We like to thank the DARTH workgroup for providing the foundations for our R code and Shiny structure template. And a special thanks goes to Merlijn Mac Gillavry for helping out with the technical issue around placing Shiny on the server. ")











