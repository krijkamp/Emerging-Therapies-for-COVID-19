# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)
# Petros Pechlivanoglou, PhD (6) 

# In collaboration of: 		
# 1 Center for Research and Teaching in Economics (CIDE), Drug Policy Program, Mexico
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada



# Installation and loading packages
install_and_load <- function(packages) {
  # Modified from https://www.listendata.com/2018/12/install-load-multiple-r-packages.html
  # The function below performs the following operations -
  #  - First it finds all the already installed R packages
  #  - Check packages which we want to install are already installed or not.
  #  - If package is already installed, it does not install it again.
  #  - If package is missing (not installed), it installs the package.
  #  - Loop through steps 2, 3 and 4 for multiple packages we want to install
  #  - Load all the packages (both already available and new ones).
  k <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(k)){
    install.packages(k, 
                     repos = "https://cran.rstudio.com/", 
                     dependencies = TRUE)
  }
  
  for(package_name in packages){
    library(package_name,
            character.only = TRUE, 
            quietly = TRUE)
  }
}

# used for sensitivity analysis

calculate_ce_out <- function (l_params_all, n_wtp = 100000) {
  with(as.list(l_params_all), {
    
    # Static characteristics
    v_x      <- runif(n_i, min = 0.95, max = 1.05) # treatment effect modifier at baseline                                         
    v_age0   <- sample(x = dist_Age$age, prob = dist_Age$prop, size = n_i, replace = TRUE) # sample from age distribution an initial age for every individual
    df_X     <- data.frame(ID = 1:n_i, x = v_x, Age = v_age0)
    
    
    #### 05.1 Probability function ####
    # The Probs function that updates the transition probabilities of every cycle is shown below.
    
    Probs <- function(M_t, df_X, v_Ts, t) { 
      # Arguments:
      # M_t: health state occupied  at cycle t (character vector)
      # df_X: dataframe with individual caracteristics
      # v_Ts: time an individual is sick
      # t:     current cycle 
      # Returns: 
      #   transition probabilities for that cycle
      
      m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  # create matrix of state transition probabilities
      rownames(m_p_t) <-  v_n                               # give the state names to the rows
      
      # lookup baseline probability and rate of dying based on individual characteristics
      p_HD_all <- inner_join(df_X, p_mort, by = c("Age"))
      p_HD     <- p_HD_all[M_t == "H","p_HD"]
      
      # update the v_p with the appropriate probabilities   
      m_p_t[, M_t == "H"]  <- rbind(1 - p_HS1 - p_HD, p_HS1, 0, p_HD)                             # transition probabilities when healthy
      m_p_t[, M_t == "S1"] <- rbind(p_S1H, 1 - p_S1H - p_S1S2 - p_S1D[v_Ts], p_S1S2, p_S1D[v_Ts]) # transition probabilities when sick
      m_p_t[, M_t == "S2"] <- rbind(0, 0, 1 - p_S2D, p_S2D)                                       # transition probabilities when sicker
      m_p_t[, M_t == "D"]  <- rbind(0, 0, 0, 1)                                                   # transition probabilities when dead   
      return(t(m_p_t))
    }       
    
    #### 05.2 Cost function ####
    # The Costs function estimates the costs at every cycle.
    
    Costs <- function (M_t, Trt = FALSE) {
      # M_t: health state occupied by individual i at cycle t (character variable)
      # Trt:  is the individual being treated? (default is FALSE) 
      
      c_t <- 0                                 # by default the cost for everyone is zero 
      c_t[M_t == "H"]  <- c_H                  # update the cost if healthy
      c_t[M_t == "S1"] <- c_S1 + c_Trt * Trt   # update the cost if sick conditional on treatment
      c_t[M_t == "S2"] <- c_S2 + c_Trt * Trt   # update the cost if sicker conditional on treatment
      c_t[M_t == "D"]  <- c_D                  # update the cost if dead
      
      return(c_t)        		                   # return the costs
    }
    
    #### 05.3 Health outcome function ####
    # The Effs function to update the utilities at every cycle.
    
    Effs <- function (M_t, df_X, Trt = FALSE, cl = 1) {
      # M_t: health state occupied by individual i at cycle t (character variable)
      # df_Pop: inidividual characteristics inclusing Age, Sex and the effect mofifier of the treatment effect
      # Trt:  is the individual treated? (default is FALSE) 
      # cl:   cycle length (default is 1)
      
      e_t <- 0                                       # by default the utility for everyone is zero
      e_t[M_t == "H"]  <- e_H                        # update the utility if healthy
      e_t[M_t == "S1" & Trt == FALSE] <- e_S1        # update the utility if sick
      e_t[M_t == "S1" & Trt == TRUE]  <- e_Trt * df_X$x[M_t == "S1"]  # update the utility if sick but on treatment (adjust for individual effect modifier) 
      e_t[M_t == "S2"] <- e_S2                       # update the utility if sicker
      e_t[M_t == "D"]  <- e_D                        # update the utility if dead
      
      QALYs <-  e_t * cl            # calculate the QALYs during cycle t
      return(QALYs)                 # return the QALYs
    }
    
    #### 06 Run Microsimulation ####
    MicroSim <- function(n_i, df_X , Trt = FALSE, seed = 1) {
      # Arguments:  
      # n_i:     number of individuals
      # df_X     data frame with individual data 
      ## Age      age of the individuals
      ## Sex      sex of the indivuduals 
      ## x        effect modifier  
      # Trt:     is this the individual receiving treatment? (default is FALSE)
      # seed:    defauls is 1
      
      set.seed(seed) # set the seed
      
      # create three matrices called m_M, m_C and m_E
      # number of rows is equal to the n_i, the number of columns is equal to n_t  (the initial state and all the n_t cycles)
      # m_M is used to store the health state information over time for every individual
      # m_C is used to store the costs information over time for evey individual
      # m_E is used to store the effects information over time for every individual
      
      m_M <- m_C <- m_E <- m_Ts <-  matrix(nrow = n_i, ncol = n_t + 1, 
                                           dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                           paste("cycle", 0:n_t, sep = " ")))  
      
      m_M [, 1] <- v_M_init    # initial health state at cycle 0 for individual i
      v_Ts      <- v_Ts_init   # initialize time since illnes onset for individual i
      
      m_C[, 1]  <- Costs(m_M[, 1], Trt)         # calculate costs per individual during cycle 0
      m_E[, 1]  <- Effs (m_M[, 1], df_X, Trt)   # calculate QALYs per individual during cycle 0
      
      # open a loop for time running cycles 1 to n_t 
      for (t in 1:n_t) {
        v_p <- Probs(m_M[, t], df_X, v_Ts, t)             # calculate the transition probabilities for the cycle based on  health state t
        m_M[, t + 1]  <- samplev(v_p, 1)                  # sample the current health state and store that state in matrix m_M 
        m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt)         # calculate costs per individual during cycle t + 1
        m_E[, t + 1]  <- Effs(m_M[, t + 1], df_X, Trt)    # calculate QALYs per individual during cycle t + 1
        
        v_Ts <- if_else(m_M[, t + 1] == "S1", v_Ts + 1, 0) # update time since illness onset for t + 1 
        df_X$Age[m_M[,t + 1] != "D"]  <- df_X$Age[m_M[, t + 1] != "D"] + 1
        
        
      } # close the loop for the time points 
      
      # calculate  
      tc <- m_C %*% v_dwc    # total (discounted) cost per individual
      te <- m_E %*% v_dwe    # total (discounted) QALYs per individual 
      tc_hat <- mean(tc)     # average (discounted) cost 
      te_hat <- mean(te)     # average (discounted) QALYs
      
      # store the results from the simulation in a list
      results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, tc_hat = tc_hat, te_hat = te_hat)   
      return(results)  # return the results
    } # end of the MicroSim function  
    
    ### Run the simulation for both no treatment and treatment options
    outcomes_no_trt  <- MicroSim(n_i, df_X, Trt = FALSE, seed = 1)
    outcomes_trt     <- MicroSim(n_i, df_X, Trt = TRUE, seed = 1)
    
    ## Vector with total discounted mean Costs and QALYs
    v_tc_d      <- c(outcomes_no_trt$tc_hat, outcomes_trt$tc_hat)
    v_te_d      <- c(outcomes_no_trt$te_hat, outcomes_trt$te_hat)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d     <- v_te_d * n_wtp - v_tc_d
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_te_d,
                        NMB      = v_nmb_d)
    return(df_ce)
  }
  )
}



# plot health state trace
plot_m_TR <- function(m_M, title = "Health state trace" ) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  matplot(m_TR, type = "l", main = title, col= 1:n_states,
          ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  legend("topright", v_n, col = 1:n_states,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)
  
}

#-----------------------------------------------------------------------------------------------#
#### R function to extract the parameters of a beta distribution from mean and st. deviation ####
#-----------------------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
betaPar <- function(m, s) 
{
  a <- m * ((m * (1 - m) / s ^ 2) - 1)
  b <- (1 - m) * ((m * (1 - m) / s ^ 2) - 1)
  list(a = a, b = b)
}

beta_params <- function(mean, sigma) {
  alpha <- ((1 - mean) / sigma ^ 2 - 1 / mean) * mean ^ 2
  beta  <- alpha * (1 / mean - 1)
  params <- list(alpha = alpha, beta = beta)
  return(params)
}

#-------------------------------------------------------------------------------------------------#
#### R function to extract the parameters of a gamma distribution from mean and st. deviation  ####
#-------------------------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
gammaPar <- function(m, s) {   
  # m: mean  
  # s: standard deviation 
  shape <- m ^ 2 / s ^ 2
  scale <- s ^ 2 / m
  list(shape = shape, scale = scale)
}

#-----------------------------------------------------------------------------------------------#
#### R function to create Markov Trace ####
#-----------------------------------------------------------------------------------------------#
CalculateMarkovTrace <- function(m_Trans, v_MT1, n_s, v_n, n_t){
  # Arguments
  # m_Trans Transition probability matrix (control or intervention)
  # v_MT1   Starting cohort allocation
  # n_s    No. of health states
  # v_n    Vector with health state names
  # n_t:    No. of cycles
  # Return
  # m_Trace The cohort trace matrix 
  
  m_Trace <- matrix(NA, nrow = n_t + 1, ncol = n_s,
                    dimnames = list(paste("Cycle", 0:n_t, sep = " "), v_n))
  # Creates a trace matrix for the allocation of the cohort in each cycle for each health state
  m_Trace[1, ] <- v_MT1
  # First cycle is the starting allocation of the cohort
  for(i in 2:(n_t + 1)) {
    m_Trace[i, ] <- t(m_Trace[i - 1, ]) %*% m_Trans
  }
  # Fills the rows of the trace matrix by multiplying the first row 
  # of the trace matrix with the transition probability matrix
  m_Trace
  # Function returns the trace matrix
}

#-----------------------------------------------------------------------------------------------#
#### R function to plot Markov Trace ####
#-----------------------------------------------------------------------------------------------#
PlotTrace <- function(trace, xlab, title, txtsize = 12) {
  # Plots the Markov trace
  # Args:
  #  trace:   Markov trace generated by `CalculateMarkovTrace` function of Micro trace generated by 'CalculateMicroTrace' 
  #  xlab:    x-axis label (e.g. "years", "days" etc.)
  #  title:   Title of the plot, (e.g. "Markov Trace" or "Microsimulation Trace")
  #  txtsize: Text size for plot, default = 12
  #
  # Return
  #  plot_trace: ggplot of Markov trace
  require(reshape2)
  require(ggplot2)
  
  trace <- data.frame(time = seq(1, (nrow(trace))), trace)
  trace <- melt(trace, id.vars = "time")
  plot_trace <- ggplot(trace, aes(x = time, y = value, colour = variable)) +
    geom_line() +
    scale_colour_hue("States", l = 50) +
    ggtitle(title) +
    xlab(xlab) +
    ylab("Proportion") +
    theme_bw() +
    theme(title = element_text(face = "bold", size = 14),
          axis.title.x = element_text(face = "bold", size = txtsize),
          axis.title.y = element_text(face = "bold", size = txtsize),
          axis.text.y  = element_text(size = txtsize),
          axis.text.x  = element_text(size = txtsize))
  
  return(plot_trace)
}

#-----------------------------------------------------------------------------------------------#
#### R function to plot Markov Trace with set limits x/y axis  ####
#-----------------------------------------------------------------------------------------------#
PlotTrace2 <- function(trace, xlab, title, txtsize = 12) {
  # Plots the Markov trace
  # Args:
  #  trace:   Markov trace generated by `CalculateMarkovTrace` function of Micro trace generated by 'CalculateMicroTrace' 
  #  xlab:    x-axis label (e.g. "years", "days" etc.)
  #  title:   Title of the plot, (e.g. "Markov Trace" or "Microsimulation Trace")
  #  txtsize: Text size for plot, default = 12
  #
  # Return
  #  plot_trace: ggplot of Markov trace
  require(reshape2)
  require(ggplot2)
  v_names <- colnames(trace)
  trace <- data.frame(time = seq(1, (nrow(trace))), trace)
  trace <- melt(trace, id.vars = "time")
 
  
  plot_trace <- ggplot(trace, aes(x = time, y = value, group = variable)) +
    geom_line(aes(linetype = variable, color = variable)) +
    scale_colour_hue("States", l = 50, labels = v_names) +
    scale_linetype_discrete("States", labels = v_names)+
    scale_y_continuous(limits = c(0, 1)) +
    ggtitle(title) +
    xlab(xlab) +
    ylab("Proportion") +
    theme_bw() +
    theme(title = element_text(face = "bold", size = 14),
          axis.title.x = element_text(face = "bold", size = txtsize),
          axis.title.y = element_text(face = "bold", size = txtsize),
          axis.text.y  = element_text(size = txtsize),
          axis.text.x  = element_text(size = txtsize))
  
  return(plot_trace)
}


#-----------------------------------------------------------------------------------------------#
#### R function to create distribution  ####
#-----------------------------------------------------------------------------------------------#

# Creating beta distribution function
beta_mom <- function(mean, var){
  # Arguments
  # mean: mean of the probability
  # var: variance of the probability
  # Returns
  # beta distribution for the probability
  term <- mean * (1 - mean) / var - 1
  alpha <- mean * term
  beta <- (1 - mean) * term
  if (var >= mean * (1 - mean)) stop("var must be less than mean * (1 - mean)")
  return(list(alpha = alpha, beta = beta))
}

beta_mom(0.8, 0.1)

# Creating a lognormal distributions
lnorm_mom <- function(mean, sd){
  # Arguments
  # mean : mean of the relative risk (in our case)
  # sd : standard deviation of the relative risk
  # Returns
  # lognormal distribution for the probability
  if (mean > 0){
    sigma2 <- log((sd ^ 2 + mean ^ 2) /mean ^ 2)
    mu <- log(mean) - 1/2 * sigma2
  } else{
    stop("Mean must be positive")
  }
  return(list(mu = mu, sigma2 = sigma2))
}

#-----------------------------------------------------------------------------------------------#
#### R function to change probabilities to rates and rates to probabilities  ####
#-----------------------------------------------------------------------------------------------#
ProbToRate <- function(p , t){
  # argument
  # p: the probability of the event 
  # t:  time in which the event took place 
  # Retunrs:
  # r : rate
  r <- -(1/t) * log(1 - p)
  r
}

RateToProb <- function(r, t){
  # ARgument 
  p <- 1 - exp(-r * t)
  return(p)
}


#-----------------------------------------------------------------------------------------------#
#### R function to generate net monetary and net heatlh benefit  ####
#-----------------------------------------------------------------------------------------------#

# Function for net health benefit (NBM)
calculateNHB <- function(effectiveness, costs, WTP){
  NHB <- effectiveness - costs / WTP
  return(NHB)
}

# Function for net monetary benefit
calculateNMB <- function(effectiveness, costs, WTP) {
  NMB <- effectiveness * WTP - costs
  return(NMB)
}

#-----------------------------------------------------------------------------------------------#
#### R function that converts a VAS score into a standard gamble equivalent  ####
#-----------------------------------------------------------------------------------------------#
convertVAStoUtility <- function(vas_score, r){
  # Arguments
  ## vas_score: as reported by individual 
  ## r: conversion factor ranging between 1.6 and 2.3
  # Returns 
  ## The utility values 
  utility <- 1 - (1 - vas_score/100) ^ r
  return(utility)
}

#-----------------------------------------------------------------------------------------------#
#### R function for Value of Information Analysis  ####
#-----------------------------------------------------------------------------------------------#
#### Formatting functions ####
# Run them all before continuing!
# Function for number of axis ticks in ggplot2 graphs
number_ticks <- function(n) {function(limits) pretty(limits, n)} 
# Total population affected by the decision
TotPop <- function(time, prev, incid, disc = 0){
  # Computes total population afected by technology
  #
  # Args:
  #   time:  vector with time points defining technology lifetime
  #   prev:  present prevalence
  #   incid: incidence
  #   disc:  discount factor; deafult = 0.
  #
  # Returns:
  #   tot.pop: total population afected by technology over technology lifetime
  #  
  # Technology Life Time, the last entry of vector `time`
  LT            <- time[length(time)]
  # Vector with population afected by the technolgy at each time point
  pop_time      <- c(prev, rep(incid, (length(time)-1))) 
  # Vector with present value of population afected by the technolgy at each time point
  disc_pop_time <- pop_time/(1+disc)^time
  # Total population afected by the technology
  tot_pop <-sum(disc_pop_time)
}
# Cost of Research
CostRes <- function(fixed.cost = 0, 
                    samp_size, 
                    cost_per_patient, 
                    INMB, 
                    clin_trial = TRUE, n_arms = 2){
  # Computes the cost of collecting information (i.e., through a research study)
  #
  # Args:
  #   fixed_cost:       fixed cost of collecting information
  #                     (e.g., fixed cost of a clinical trial); default = 0
  #   samp_size:               vector with sample sizes
  #   cost_per_patient: cost per patient in research study
  #   INMB:             Incremental Net Monetary Benefit
  #   clin_trial:       indicator whether calculation is for a clinical trial;
  #                     default = TRUE
  #   n_arms:           Number of arms in research study design; default = 2
  #
  # Returns:
  #   cost_res: vector with the total cost of collecting information for each simple size
  #
  if (clin_trial){
    Cost_Res <- fixed_cost + n_arms*samp_size*cost_per_patient + samp_size*INMB
  } else { # E.g., cohort study
    Cost_Res <- fixed_cost + samp_size*cost_per_patient
  }
  return(Cost_Res)
}



beta_params <- function(mean, sigma) {
  alpha <- ((1 - mean) / sigma ^ 2 - 1 / mean) * mean ^ 2
  beta  <- alpha * (1 / mean - 1)
  params <- list(alpha = alpha, beta = beta)
  return(params)
}


#------------------------------------------------------------------------------#
####                         Decision Model                                 ####
#------------------------------------------------------------------------------#
#' Decision Model
#'
#' \code{decision_model} implements the decision model used.
#'
#' @param l_params_all List with all parameters of decision model
#' @param verbose Logical variable to indicate print out of messages
#' @return The transition probability array and the cohort trace matrix.
#' 
decision_model <- function(l_params_all, verbose = FALSE) {
  with(as.list(l_params_all), {
    # compute internal paramters as a function of external parameters
    r_HD    = - log(1 - p_HD) # rate of death in healthy
    r_S1D   = hr_S1 * r_HD 	  # rate of death in sick
    r_S2D   = hr_S2 * r_HD  	# rate of death in sicker
    p_S1D   = 1 - exp(-r_S1D) # probability to die in sick
    p_S2D   = 1 - exp(-r_S2D) # probability to die in sicker
    
    ####### INITIALIZATION ##########################################
    # create the cohort trace
    m_M <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n))     # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    
    m_M[1, ] <- c(1, 0, 0, 0)                      # initialize Markov trace
    
    # create transition probability matrix for NO treatment
    m_P <- matrix(0,
                  nrow = n_s, 
                  ncol = n_s,
                  dimnames = list(v_n, v_n))
    # fill in the transition probability array
    ### From Healthy
    m_P["H", "H"]   <- 1 - (p_HS1 + p_HD)
    m_P["H", "S1"]  <- p_HS1
    m_P["H", "D"]   <- p_HD
    ### From Sick
    m_P["S1", "H"]  <- p_S1H
    m_P["S1", "S1"] <- 1 - (p_S1H + p_S1S2 + p_S1D)
    m_P["S1", "S2"] <- p_S1S2
    m_P["S1", "D"]  <- p_S1D
    ### From Sicker
    m_P["S2", "S2"] <- 1 - p_S2D
    m_P["S2", "D"]  <- p_S2D
    ### From Dead
    m_P["D", "D"]   <- 1
    
    # check rows add up to 1
    if (!isTRUE(all.equal(as.numeric(rowSums(m_P)), as.numeric(rep(1, n_s))))) {
      stop("This is not a valid transition Matrix")
    }
    
    ############# PROCESS ###########################################
    
    for (t in 1:n_t){                     # throughout the number of cycles
      m_M[t + 1, ] <- m_M[t, ] %*% m_P    # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    
    ####### EPIDEMIOLOGICAL OUTPUT  ###########################################
    #### Overall Survival (OS) ####
    v_os      <- 1 - m_M[, "D"]           # calculate the overall survival (OS) probability for no treatment
    
    #### Disease prevalence #####
    v_prev    <- rowSums(m_M[, c("S1", "S2")])/v_os
    
    #### Proportion of sick in S1 state #####
    v_prop_S1 <- m_M[, "S1"] / v_prev
    
    ####### RETURN OUTPUT  ###########################################
    out <- list(m_M      = m_M,
                m_P      = m_P,
                Surv     = v_os[-1],
                Prev     = v_prev[-1],
                PropSick = v_prop_S1[c(11, 21, 31)])
    
    return(out)
  }
  )
}

#------------------------------------------------------------------------------#
####              Calculate cost-effectiveness outcomes                     ####
#------------------------------------------------------------------------------#
#' Calculate cost-effectiveness outcomes
#'
#' \code{calculate_ce_out} calculates costs and effects for a given vector of parameters using a simulation model.
#' @param l_params_all List with all parameters of decision model
#' @param n_wtp Willingness-to-pay threshold to compute net benefits
#' @return A data frame with discounted costs, effectiveness and NMB.
#' 
calculate_ce_out <- function(l_params_all, n_wtp = 100000){ # User defined
  with(as.list(l_params_all), {
    ## Create discounting vectors
    v_dwc <- 1 / ((1 + d_e) ^ (0:(n_t))) # vector with discount weights for costs
    v_dwe <- 1 / ((1 + d_c) ^ (0:(n_t))) # vector with discount weights for QALYs
    
    ## Run STM model at a parameter set for each intervention
    l_model_out_no_trt <- decision_model(l_params_all = l_params_all)
    l_model_out_trt    <- decision_model(l_params_all = l_params_all)
    
    ## Cohort trace by treatment
    m_M_no_trt  <- l_model_out_no_trt$m_M # No treatment
    m_M_trt     <- l_model_out_trt$m_M    # Treatment
    
    ## Vectors with costs and utilities by treatment
    v_e_no_trt  <- c(e_H, e_S1, e_S2, e_D)
    v_e_trt     <- c(e_H, e_trt, e_S2, e_D)
    
    v_c_no_trt  <- c(c_H, c_S1, c_S2, c_D)
    v_c_trt     <- c(c_H, c_S1 + c_trt, c_S2 + c_trt, c_D)
    
    ## Mean Costs and QALYs for Treatment and NO Treatment
    v_te_no_trt <- m_M_no_trt %*% v_e_no_trt
    v_te_trt    <- m_M_trt %*% v_e_trt
    
    v_tc_no_trt <- m_M_no_trt %*% v_c_no_trt
    v_tc_trt    <- m_M_trt %*% v_c_trt
    
    ## Total discounted mean Costs and QALYs
    te_d_no_trt <- t(v_te_no_trt) %*% v_dwe 
    te_d_trt    <- t(v_te_trt) %*% v_dwe
    
    tc_d_no_trt <- t(v_tc_no_trt) %*% v_dwc
    tc_d_trt    <- t(v_tc_trt)    %*% v_dwc
    
    ## Vector with total discounted mean Costs and QALYs
    v_tc_d      <- c(tc_d_no_trt, tc_d_trt)
    v_te_d      <- c(te_d_no_trt, te_d_trt)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d     <- v_te_d * n_wtp - v_tc_d
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_te_d,
                        NMB      = v_nmb_d)
    
    return(df_ce)
  }
  )
}

#------------------------------------------------------------------------------#
#### R function to extract the parameters of a beta distribution            ####
####                from mean and st. deviation                             ####
#------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
betaPar <- function(m, s) 
{
  a <- m * ((m * (1 - m) / s ^ 2) - 1)
  b <- (1 - m) * ((m * (1 - m) / s ^ 2) - 1)
  list(a = a, b = b)
}

#------------------------------------------------------------------------------#
#### R function to extract the parameters of a gamma distribution           ####
####                   from mean and st. deviation                          ####
#------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
gammaPar <- function(m, s) {   
  # m: mean  
  # s: standard deviation 
  shape <- m ^ 2 / s ^ 2
  scale <- s ^ 2 / m
  list(shape = shape, scale = scale)
}

# Rowmaxs function (from https://rdrr.io/github/ejanalysis/analyze.stuff/src/R/rowMaxs.R)
#' @title Returns the max value of each row of a data.frame or matrix
#'
#' @description
#' Returns maximum value of each row of a data.frame or matrix.
#' @param df Data.frame or matrix, required.
#' @param na.rm Logical value, optional, TRUE by default. Defines whether NA values should be removed first. Otherwise result will be NA when any NA is in the given vector.
#' @return Returns a vector of numbers of length equal to number of rows in df.
#' @template maxmin
#' @export
rowMaxs <- function(df, na.rm=TRUE) {
  
  if (is.matrix(df)) {df <- data.frame(df, stringsAsFactors=FALSE, drop = FALSE)}
  
  valid.cols <- sapply(df, function(x) { is.numeric(x) || is.logical(x) || is.character(x)})
  stopifnot(any(valid.cols))
  # or could just return NA?:
  # if (!any(valid.cols)) {return(NA)}
  if (any(!valid.cols) ) {warning('using only numeric (double or integer) or logical or character columns -- ignoring other columns ')}
  
  result <- do.call(pmax, c(df[ , valid.cols, drop = FALSE], na.rm=na.rm))
  
  result[nononmissing <- rowSums(!is.na(df[ , valid.cols, drop = FALSE]))==0] <- -Inf
  if (any(nononmissing)) {warning('where no non-missing arguments, returning -Inf')}
  return(result)
  
  # df = data.frame of numeric values, i.e. a list of vectors passed to pmax
  # Value returned is vector, each element is max of a row of df
}

#### Formatting functions ####
# Run them all before continuing!
# Function for number of axis ticks in ggplot2 graphs
number_ticks <- function(n) {function(limits) pretty(limits, n)} 
# Total population affected by the decision
TotPop <- function(time, prev, incid, disc = 0){
  # Computes total population afected by technology
  #
  # Args:
  #   time:  vector with time points defining technology lifetime
  #   prev:  present prevalence
  #   incid: incidence
  #   disc:  discount factor; deafult = 0.
  #
  # Returns:
  #   tot.pop: total population afected by technology over technology lifetime
  #  
  # Technology Life Time, the last entry of vector `time`
  LT            <- time[length(time)]
  # Vector with population afected by the technolgy at each time point
  pop.time      <- c(prev, rep(incid, (length(time)-1))) 
  # Vector with present value of population afected by the technolgy at each time point
  disc.pop.time <- pop.time/(1+disc)^time
  # Total population afected by the technology
  tot.pop <-sum(disc.pop.time)
}
# Cost of Research
CostRes <- function(fixed.cost = 0, 
                    samp.size, 
                    cost.per.patient, 
                    INMB, 
                    clin.trial = TRUE, n.arms = 2){
  # Computes the cost of collecting information (i.e., through a research study)
  #
  # Args:
  #   fixed.cost:       fixed cost of collecting information
  #                     (e.g., fixed cost of a clinical trial); default = 0
  #   samp.size:               vector with sample sizes
  #   cost.per.patient: cost per patient in research study
  #   INMB:             Incremental Net Monetary Benefit
  #   clin.trial:       indicator whether calculation is for a clinical trial;
  #                     default = TRUE
  #   n.arms:           Number of arms in research study design; default = 2
  #
  # Returns:
  #   cost.res: vector with the total cost of collecting information for each simple size
  #
  if (clin.trial){
    Cost.Res <- fixed.cost + n.arms*samp.size*cost.per.patient + samp.size*INMB
  } else { # E.g., cohort study
    Cost.Res <- fixed.cost + samp.size*cost.per.patient
  }
  return(Cost.Res)
  
  
  ## Function to fit multiple functional forms to survival data
  fit.fun <- function(time, status, data = data , add = FALSE, extrapolate = FALSE, times)  
  {
    #Extact the right data columns 
    data$time   <-   data[,   time]  
    data$status <-   data[, status]  
    
    if (extrapolate == TRUE)  {
      plot.times <- max(times)
    } else if  (extrapolate == FALSE) {
      plot.times <- max(data$time)
    }
    
    # Progression free survival  
    KM.fit     <-     survfit(Surv(time, status) ~ 1, data = data)                         # fit Kaplan-Meier curve 
    fit.llogis <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "llogis" )       # fit model with loglogistic distribution
    fit.weib   <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "weibull")       # fit model with Weibull distribution
    fit.lnorm  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "lnorm"  )       # fit model with lognormal distribution
    fit.gamma  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "gamma"  )       # fit model with gamma distribution 
    fit.exp    <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "exp"    )       # fit model with exponential distribution
    fit.gengamma  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "gengamma"  ) # fit model with gamma distribution  
    
    
    # extarapolate all models beyond the KM curve
    if(add){ lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
    if(!add){ plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F, mark.time= T)}
    lines(fit.llogis,   t = times, col = 2, ci = F)
    lines(fit.weib,     t = times, col = 3, ci = F)
    lines(fit.lnorm,    t = times, col = 4, ci = F)
    lines(fit.gamma,    t = times, col = 5, ci = F)
    lines(fit.gengamma,    t = times, col = 6, ci = F)
    lines(fit.exp,      t = times, col = 7, ci = F)
    legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal", "Gamma","GenGamma", "Exponential"), col = 1:7, lty = rep(1, 7), bty="n")
    
    # compare AIC values
    AIC <- c(    Loglogistic = AIC(fit.llogis),                                         
                 Weibull     = AIC(fit.weib), 
                 Lognormal   = AIC(fit.lnorm), 
                 Gamma       = AIC(fit.gamma),
                 GenGamma       = AIC(fit.gengamma),
                 Exponentail = AIC(fit.exp))
    AIC= round(AIC,3)
    
    # compare BIC values
    BIC <- c(    Loglogistic = BIC(fit.llogis),                                         
                 Weibull     = BIC(fit.weib), 
                 Lognormal   = BIC(fit.lnorm), 
                 Gamma       = BIC(fit.gamma),
                 GenGamma    = BIC(fit.gengamma),
                 Exponential = BIC(fit.exp))
    
    BIC <- round(BIC,3)
    
    res <- list(Loglogistic = fit.llogis,
                Weibull     = fit.weib,
                Lognormal   = fit.lnorm, 
                Gamma       = fit.gamma,
                GenGamma       = fit.gengamma,
                Exponential = fit.exp, 
                AIC         = AIC,
                BIC         = BIC)
    res
  }
  
  
  fit.mstate <- function(time, status, trans,  data = data , add = FALSE, extrapolate = FALSE, times)  
  {
    data$time  <- data[, time  ]
    data$tatus <- data[, status]
    
    if (extrapolate == TRUE)  {
      plot.times <- max(times)
    } else if  (extrapolate == FALSE) {
      plot.times <- max(data$time)
    }
    
    # Progression free survival  
    KM.fit     <-     survfit(Surv(time, status) ~ trans , data = data)                                 # fit Kaplan-Meier curve 
    fit.llogis <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "llogis" ) # fit model with loglogistic distribution
    fit.weib   <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "weibull") # fit model with Weibull distribution
    fit.lnorm  <- flexsurvreg(Surv(time, status) ~ trans + sdlog(trans), data = data, dist = "lnorm"  ) # fit model with lognormal distribution
    fit.gamma  <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "gamma"  ) # fit model with gamma distribution 
    fit.gengamma  <- flexsurvreg(Surv(time, status) ~ trans + Q(trans) + sigma(trans), data = data, dist = "gengamma"  ) # fit model with gamma distribution 
    
    
    # extarapolate all models beyond the KM curve
    if(add){ lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
    if(!add){ plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
    lines(fit.llogis,   t = times, col = 2, ci = F)
    lines(fit.weib,     t = times, col = 3, ci = F)
    lines(fit.lnorm,    t = times, col = 4, ci = F)
    lines(fit.gamma,    t = times, col = 5, ci = F)
    lines(fit.gengamma,    t = times, col = 6, ci = F)
    legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal", "Gen.Gamma"), col = 1:5, lty = rep(1, 5), bty="n")
    
    # compare AIC values
    AIC <- c(    Loglogistic = AIC(fit.llogis),                                         
                 Weibull     = AIC(fit.weib), 
                 Lognormal   = AIC(fit.lnorm), 
                 Gamma       = AIC(fit.gamma),
                 GenGamma    = AIC(fit.gengamma))
    AIC= round(AIC,3)
    
    # compare BIC values
    BIC <- c(    Loglogistic = BIC(fit.llogis),                                         
                 Weibull     = BIC(fit.weib), 
                 Lognormal   = BIC(fit.lnorm), 
                 Gamma       = BIC(fit.gamma),
                 GenGamma    = BIC(fit.gengamma))
    
    BIC <- round(BIC,3)
    
    res <- list(Loglogistic = fit.llogis,
                Weibull     = fit.weib,
                Lognormal   = fit.lnorm, 
                Gamma       = fit.gamma,
                GenGamma    = fit.gengamma,
                AIC         = AIC,
                BIC         = BIC)
    res
  }
  
  
  trace.DES = function(msm_sim = des_sim, tmat, n_i, times )
  {
    # Restructure the data to extract markov trace
    data.mstate.sim <- data.frame(cbind(matrix(t(msm_sim$st), ncol=1),
                                        matrix(t(msm_sim$t) , ncol=1)))
    colnames(data.mstate.sim) <- c("state","time")
    data.mstate.sim$subject <- rep(1:n_i, each = ncol(msm_sim$st))
    
    data.mstate.sim = na.omit(data.mstate.sim)
    data.mstate.sim = data.mstate.sim[!duplicated(data.mstate.sim), ] # remove duplicate entries in the dataset
    
    # create transition intensitiy matrix with initial values based on the structure of tmat
    twoway7.q               <- tmat
    twoway7.q[!is.na(tmat)] <- 0.5
    twoway7.q[is.na(tmat)]  <- 0
    # fit msm model only so that we can extract the prevalence (i.e. trace) thrrough the prevalence.msm function
    
    fit.msm.sim <- msm(state ~ time,subject = subject, data = data.mstate.sim, qmatrix = twoway7.q, 
                       exacttimes = T, use.deriv = TRUE, analyticp = FALSE, fixedpars = TRUE, hessian = F)
    
    M.tr.des <- prevalence.msm(fit.msm.sim, times = times) # Markov trace when DES model is used
    
    
    return(M.tr.des[[3]]/100)
  }
  
  
  partsurv <- function(fit.pfs, fit.os, time = times){
    # Input
    # fit.pfs: flexsurv obj fitting pfs
    # fit.os: flexsurv obj fitting os
    # title:
    # time = numeric vector of time to estimate probabilities
    # output:
    #  res a list w/ one entry of a data frame w/ probabilities associated w/ stable ,prog and dead.
    
    pfs.surv <- summary(fit.pfs, t = time, ci = F)[[1]]$est
    os.surv  <- summary(fit.os,  t = time, ci = F)[[1]]$est
    sick                 <- os.surv - pfs.surv      # estimate the probability of remaining in the progressed state
    sick[sick < 0]       <- 0                       # in cases where the probability is negative replace with zero
    healthy              <- pfs.surv                # probability of remaining stable
    dead                 <- 1 - os.surv             # probability of being dead
    trace <- cbind(healthy, sick, dead)
    res   <- list(trace = trace)
    
    return(res)
  }
  
  
  flexsurvreg_prob <- function(object, newparams = NULL, times){
    
    if(is.null(newparams) == T ){
      params <- object$res[,1]  
      params <- as.matrix(t(params))
    }else {
      params <- newparams 
      params <- as.matrix(params)
    }
    
    if (ncol(params)== 1){
      surv <- object$dfns$p(times, params[,1], lower.tail = F)
    }else if (ncol(params)== 2){
      surv <- object$dfns$p(times,params[,1],params[,2], lower.tail = F)
    }else if (ncol(params)== 2){
      surv <- object$dfns$p(times,params[,1],params[,2],params[,3], lower.tail = F)
    } else{
      surv <- object$dfns$p(times,params[,1],params[,2],params[,3], lower.tail = F)
    }
    
    t.p <- 1- surv[-1]/(surv[-length(surv)])
    return(t.p = t.p)
  }
  
  
  gen_data <- function(n_pat, n_years)
  {
    # specification of hazard functions to generate data from
    hazardf <- gems::generateHazardMatrix(n_s)
    colnames(hazardf@list.matrix) <- 
      rownames(hazardf@list.matrix) <- v_n
    
    # specifying the transition hazard from healthy -> sick
    hazardf[["healthy","sick"]] <- function (t, r1, r2){
      hweibull(t,r1, r2)
    }
    
    # specifying the transition hazard from healthy -> dead 
    hazardf[["healthy","dead"]] <- function (t, r1, r2){
      flexsurv::hgompertz(t,r1, r2)
    }
    
    # specifying the transition hazard from sick -> dead 
    hazardf[["sick","dead"]] <- function (t, r1, r2){
      hweibull(t,r1, r2)
    }
    
    
    # list of parameters for the hazard functions defined above
    mu        <- gems::generateParameterMatrix(hazardf) 
    rownames(mu@list.matrix) <- 
      colnames(mu@list.matrix) <- v_n
    
    mu[["healthy", "sick"]] <- list(1.5, 6)      #  the Weibull parameters for H -> S
    mu[["healthy", "dead"]] <- list(0.25, 0.08)  # the Gompertz params for H -> D
    mu[["sick",    "dead"]] <- list(0.5,4)       #  the Weibull parameters for S -> D
    
    
    
    # simulate the cohort
    cohort <- gems::simulateCohort(
      transitionFunctions = hazardf,
      parameters = mu,
      cohortSize = n_pat,
      to = n_years)
    
    # extract the simulated true data 
    true_data <- cohort@time.to.state
    colnames(true_data) <- v_n
    
    true_data$dead[is.na(true_data$dead)] <- n_years
    true_data$sick[is.na(true_data$sick)] <- true_data$dead[is.na(true_data$sick)]
    
    
    # create a status variable that will capture the transition events
    true_status         <- matrix(NA, nrow = n_pat, ncol = n_s, dimnames = list(1:n_pat,v_n))
    true_status         <- as.data.frame(true_status)
    true_status$healthy <- ifelse(is.na(true_data$healthy),0,1)
    true_status$dead    <- ifelse(true_data$dead == n_years, 0, 1)
    true_status$sick    <- ifelse(true_data$dead == true_data$sick, 0, 1)
    
    
    censtime <- runif(n = n_pat, 0, n_years)
    
    censored_sick <- ifelse(censtime      <= true_data$sick |
                              true_data$sick >  5, 1, 0)
    censored_dead <- ifelse(censtime <= true_data$dead|
                              true_data$dead >5, 1, 0)
    
    sim_data <- true_data
    
    sim_data$sick[censored_sick == 1] <-  censtime[censored_sick == 1]
    sim_data$sick[sim_data$sick >5 ]  <-  5
    
    sim_data$dead[censored_dead == 1] <-  censtime[censored_dead == 1]
    sim_data$dead[sim_data$dead >5] <-  5
    
    status <- true_status
    
    status$sick[censored_sick == 1] = 0
    status$dead[censored_dead == 1] = 0
    
    # Usually trials report OS/PFS outcomes so we will recreate those
    
    OS_PFS_data <- data.frame(row.names = 1:n_pat)
    
    OS_PFS_data$PFS_time        <- apply(sim_data[, c("sick","dead")], 1, min) 
    OS_PFS_data$PFS_status      <- ifelse(status$dead == 1 | status$sick == 1, 1, 0 )
    
    OS_PFS_data$OS_time         <- sim_data$dead
    OS_PFS_data$OS_status       <- status$dead 
    list(cohort = cohort, true_data = true_data, true_status = true_status, 
         sim_data =  sim_data,      status = status, OS_PFS_data = OS_PFS_data)
  }
  
  samplev <- function(m.Probs, m) {
    # Arguments
    # m.Probs: matrix with probabilities (n.i * n.s)
    # m:       number of states than need to be sampled per individual  
    # Return
    # ran:    n.i x m matrix filled with sampled health state(s) per individual
    
    d <- dim(m.Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
    n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
    k <- d[2]          # second dimension - number of columns (number of health states considered)
    lev <- dimnames(m.Probs)[[2]]  # extract the names of the health states considered for sampling
    if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
      lev <- 1:k       # create a sequence from 1:k (number of health states considered)
    # create a matrix 
    ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
    U <- t(m.Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
    
    for(i in 2:k) {    # start loop, from the 2nd health states
      U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
    }
    if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
      stop("error in multinom: probabilities do not sum to 1")
    
    for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
      un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
      ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
    }
    ran # return the new health state per individual n.i x m
  } # close the function 
  
  
  #plot health state trace
  plot_m_TR <- function(m_M) {
    # plot the distribution of the population across health states over time (trace)
    # count the number of individuals in each health state at each cycle
    m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
    m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
    colnames(m_TR) <- v_n                                    # name the rows of the matrix
    rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
    # Plot trace of first health state
    matplot(m_TR, type = "l", main = "Health state trace", col= 1:n_s,
            ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
    legend("topright", v_n, col = 1:n_s,    # add a legend to current plot
           lty = rep(1, 3), bty = "n", cex = 0.65)
    
  }
  
  
  
}




#----------------------------------------------------------------------------#
####   Function to check if transition probability array/matrix  is valid ####
#----------------------------------------------------------------------------#
#' Check if transition array is valid
#'
#' \code{check_transition_probability} checks if transition probabilities are in \[0, 1\].
#'
#' @param a_P A transition probability array or matrix.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = FALSE
#'
#' @return
#' This function stops if transition probability array is not valid and shows 
#' what are the entries that are not valid
#' @import utils
#' @export
check_transition_probability <- function(a_P,
                                         err_stop = FALSE, 
                                         verbose = FALSE) {
  
  a_P <- as.array(a_P)
  m_indices_notvalid <- arrayInd(which(a_P < 0 | a_P > 1), 
                                 dim(a_P))
  
  if(dim(m_indices_notvalid)[1] != 0){
    v_rows_notval   <- rownames(a_P)[m_indices_notvalid[, 1]]
    v_cols_notval   <- colnames(a_P)[m_indices_notvalid[, 2]]
    v_cycles_notval <- dimnames(a_P)[[3]][m_indices_notvalid[, 3]]
    
    df_notvalid <- data.frame(`Transition probabilities not valid:` = 
                                matrix(paste0(paste(v_rows_notval, v_cols_notval, sep = "->"),
                                              "; at cycle ",
                                              v_cycles_notval), ncol = 1), 
                              check.names = FALSE)
    
    if(err_stop) {
      stop("Not valid transition probabilities\n",
           paste(capture.output(df_notvalid), collapse = "\n"))
    }
    
    if(verbose){
      warning("Not valid transition probabilities\n",
              paste(capture.output(df_notvalid), collapse = "\n"))
    } 
  }
}

#----------------------------------------------------------------------------#
####   Function to check if sum of transition probabilities equal to one  ####
#----------------------------------------------------------------------------#
#' Check if the sum of transition probabilities equal to one. 
#'
#' \code{check_sum_of_transition_array} checks if each of the rows of the 
#' transition matrices sum to one. 
#' 
#' @param a_P A transition probability array.
#' @param n_states Number of health states.
#' @param n_t Number of cycles.
#' @param err_stop Logical variable to stop model run if set up as TRUE. 
#' Default = TRUE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = TRUE
#' @return 
#' The transition probability array and the cohort trace matrix.
#' @import dplyr
#' @export
check_sum_of_transition_array <- function(a_P,
                                          n_states,
                                          n_t,  
                                          err_stop = TRUE, 
                                          verbose  = TRUE) {
  
  a_P <- as.array(a_P)
  d <- length(dim(a_P))
  # For matrix
  if (d == 2) {
    valid <- sum(rowSums(a_P))
    if (abs(valid - n_states)> 1e-04 ) {
      if(err_stop) {
        browser()
        stop("This is not a valid transition Matrix")
      }
      
      if(verbose){
        warning("This is not a valid transition Matrix")
      } 
    }
  } else {
    # For array
    valid <- (apply(a_P, d, function(x) sum(rowSums(x))) == n_states)
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
      if(err_stop) {
        stop("This is not a valid transition Matrix")
      }
      
      if(verbose){
        warning("This is not a valid transition Matrix")
      } 
    }
  }
}



#----------------------------------------------------------------------------#
####   Function to check if all the items in the list contains information  ####
#----------------------------------------------------------------------------#
#' Check if the each item in the list contains information. 
#'
#' \code{check_list_items} checks if item in the list contains a value 
#' 
#' @param l_param A list with parameter values.
#' @param err_stop Logical variable to stop model run if set up as TRUE. 
#' Default = TRUE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = TRUE
#' @return 
#' The transition probability array and the cohort trace matrix.
#' @import dplyr
#' @export
check_list_items <- function(l_param,
                             err_stop = TRUE, 
                             verbose  = TRUE) {
  
    # check if each compontent of the list is valie, not NULL, not NA etc.
    valid <- !is.null(l_param) & class(l_param) != "NULL" & class(l_param) != "logical" & length(l_param) != 0 & !is.na(l_param)
    
    m_valid <- count(valid) # sum the logic vale for each parameters
    
  if(m_valid$freq[m_valid$x == TRUE] != length(names(l_param))) {
    if(err_stop) {
#      browser()
      stop("This is not a valid list. At least one parameter does not contain information.")
    }
    
    if(verbose){
      warning("This is not a valid list. At least one parameter does not contain information.")
    }
    
  }
    else {
      print("This is a valid list")}
  } # close the function
  
  
#---------#
# Odds Ratio to Relative Risk 
# Example https://clincalc.com/Stats/ConvertOR.aspx?example
# RR = risk ratio, OR = Odds ratio, P_ref = prevalence outcome in reference group
# RR = prob exposed/prob non-exposed, OR = odds exposed/odds non-exposed

# OR <- 0.51
# P_ref <- 0.93

# RR <- OR / ( (1-P_ref) + (P_ref* OR) ) 
# RR # 0.94

ORtoRR <- function(OR, P_ref){
  # Argument 
  RR <- OR / ( (1-P_ref) + (P_ref* OR) )
  return(RR)
}
#---------#


