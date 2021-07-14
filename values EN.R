# The values for the shiny
# This depends on the selected drug -> how can we make that?

# Save the values in the csv -> with all the possible names and input
# Load the values in shiny app applications

n_eff_size    <- 0.82             # effect size treatment
ci_eff_size   <- c(0.47, 1.13)    # confidence interval of the effect size 


n_costs     <- 1.9                 # cost of the drug 
ci_costs    <- c(1.425, 2.564)     # confidence interval costs drug 
n_Trt       <- 7                   # number of required doses of the treatment

LOS_Trt      <- 11            # Length of stay with the therapy
LOS_noTrt    <- 12            # Length of stay without the therapy

r_discount    <- 3             # the discount rate
n_iter       <- 100            # then number of PSA iteration 

n_wtp        <- 100000        # the willingness-to-pay in $/QALY
n_H_future   <- 88000         # the future number of hospitalized patients 
n_H_current  <- 741000        # the number of COVID patients during the trail 

## Not possible to adjust right now 

n_age_max <- 120       # max. age for life time horizon analysis
n_DpC     <- 73        # days per cycle
n_DpY     <- 365       # days per year
d_c       <- 0.03      # discount rate 

u_H    <- 0.5
u_I    <- 0.05
u_R_H  <- 0.88
u_R_IC <- 0.677
