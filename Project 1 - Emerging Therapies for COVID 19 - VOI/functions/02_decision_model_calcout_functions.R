# Note to myself:
# Started with ## 04.1 Cohort trace



#' The decision model for the COVID therapies 
#'
#' \code{calculate_cea_output_VOI_COVID} implements the decision model used.
#'
#' @param l_params List with all parameters of decision model
#' @param n_wtp The willingness to pay threshold, default 100000 
#' @param return_all Logical variable to indicate to return all model structures or the CE dataframe only. Default is FALSE; ce data frame only. 
#' @param verbose Logical variable to indicate print out of messages
#' @return A cost effectiveness dataframe 
#' 
calculate_cea_output_VOI_COVID <- function(l_list, n_wtp = 100000, return_all = FALSE, verbose = FALSE){
  
  l_param_trt_f <- c()
  l_param_trt_f <- l_list
 
   with(as.list(l_list), {
    
     ##### Data Generation ####
    # Estimate general model input
    l_input_general <- estimate_input_general(l_param_trt_f)
    
    # Add the general to the drug specific list 
    l_param_trt_f <- c(l_param_trt_f, l_input_general)
    
    ## Estimate the treatment cost of the drug
    l_param_trt_f <- estimate_cost_drug(l_list = l_param_trt_f)
    
    # Estimate the age specific probability of dying 
    l_param_trt_f$v_p_RD <- estimate_v_p_RD(df_r_HD = l_param_trt_f$df_r_HD, 
                                              n_age = l_param_trt_f$n_age, 
                                              p_men = l_param_trt_f$p_men,
                                              n_YpC = l_param_trt_f$n_YpC)

    # Calculate the treatment specific cost
    l_param_trt_f <- calc_costs_model(l_list = l_param_trt_f, verbose = verbose) 
    
    # split the first cycle for each treatment
    l_param_trt_f <- split_firstcycle(l_list = l_param_trt_f) 
    
    l_param_trt_f <- estimate_trt_parameters(l_list = l_param_trt_f)
  
    #### Estimate and store effects ####
    l_effects <- c()
    l_effects <- estimate_effects(l_list = l_param_trt_f)
    l_param_trt_f <- c(l_param_trt_f, l_effects)
      ###############################
      # Extract essential parameters
      n_t            <- l_param_trt_f$n_t
      n_YpC          <- l_param_trt_f$n_YpC
      v_names_states <- l_param_trt_f$v_names_states
      n_states       <- l_param_trt_f$n_states
      v_dwe          <- l_param_trt_f$v_dwe
      v_dwc          <- l_param_trt_f$v_dwc
      
      with(l_param_trt_f, {
      ####### INITIALIZATION ##########################################
      # Create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
      m_M_notrt <- m_M_trt <- matrix(NA, 
                                     nrow     = n_t + 1, ncol = n_states,
                                     dimnames = list(paste("cycle", 0:n_t, sep = " "), v_names_states))
      
      # Initialize multidimensional array for both strategies
      a_A_notrt <- a_A_trt <- a_A <- array(0, 
                                           dim = c(n_states, n_states, n_t + 1), 
                                           dimnames = list(v_names_states, v_names_states, 0:n_t)) 
      
      # The cohort starts as hospitalized, initiate first cycle of the cohort trace therefore as all patients in state H
      m_M_notrt[1, ] <- m_M_trt[1, ] <- diag(a_A_notrt[, , 1]) <- diag(a_A_trt[, , 1]) <- diag(a_A[, , 1]) <- v_s_init  # initiate first cycle of cohort trace 
      
      # Create the transition probability matrix
      m_P_trt <- m_P_notrt  <- matrix(0,
                                      nrow = n_states,
                                      ncol = n_states,
                                      dimnames = list(v_names_states, v_names_states)) # name the columns and rows of the matrix
      
      })
      ##########################################
      
      
      #### Fill m_P ####
      l_param_trt_f <- fill_m_P(l_list = l_param_trt_f)
    
      #### Run Markov model ####
      l_markov_dyn <- list() # create a list to store all the dynamics in the dynamics array (l_a_A) as well as the cohort trace matrix (l_m_M)
      
      l_trace <- l_array <- list()
      
      m_P_trt   <- l_param_trt_f$m_P_trt
      m_P_notrt <- l_param_trt_f$m_P_notrt
      
      results_trt <- run_markov(m_P = m_P_trt, df_r_HD = df_r_HD, l_list = l_param_trt_f) # use the function to run the model, take the n_age from the list (l_param) which is the same for all treatments
      
      results_notrt <- run_markov(m_P = m_P_notrt, df_r_HD = df_r_HD, l_list = l_param_trt_f) # use the function to run the model, take the n_age from the list (l_param) which is the same for all treatments
      
      l_trace$trt   <- l_markov_dyn$m_M_trt   <- results_trt$m_M # store the matrices in a list
      l_trace$notrt <- l_markov_dyn$m_M_notrt <- results_notrt$m_M # store the matrices in a list
      
      l_array$trt   <- l_markov_dyn$a_A_trt   <- results_trt$a_A # store the arrays in a list
      l_array$notrt <- l_markov_dyn$a_A_notrt <- results_notrt$a_A # store the arrays in a list
    
      ### Cost - Effectiveness Analysis ###
      l_m_R_costs <- l_a_R_costs <- l_a_Y_costs <- l_v_u <- l_v_e <- l_v_tc <- l_v_tu <- l_v_te <- list() # create a list to store all the dynamics
      
      l_a_rewards <- list()
      
      l_a_rewards <- fill_a_R(l_list = l_param_trt_f)
      
      ## Create multidimensional arrays to store expected outcomes 
      a_Y_costs_trt <- a_Y_costs_notrt <- with(l_param_trt_f, 
                                               array(0, 
                                                     dim = c(n_states, n_states, n_t + 1), 
                                                     dimnames = list(v_names_states, 
                                                                     v_names_states, 
                                                                     0:n_t))
      )  # close the with function
      
      # Initialize arrays
      a_Y_costs_trt   [, , 1] <- l_array$trt[, , 1]   * l_a_rewards$a_R_costs_trt[, , 1]
      a_Y_costs_notrt [, , 1] <- l_array$notrt[, , 1] * l_a_rewards$a_R_costs_notrt[, , 1]
        
      
      v_e_trt   <- c(l_param_trt_f$e_H_mix_Trt,   e_R_H, e_R_IC, e_D)
      v_e_notrt <- c(l_param_trt_f$e_H_mix_noTrt, e_R_H, e_R_IC, e_D)
      
      
      # Vectors with utility 
      v_u_trt   <- c(l_param_trt_f$u_H_mix_Trt,   u_R_H, u_R_IC, u_D) 
      v_u_notrt <- c(l_param_trt_f$u_H_mix_noTrt, u_R_H, u_R_IC, u_D)
      

      
      v_te_trt    <- (l_markov_dyn$m_M_trt     %*%  v_e_trt)   * n_YpC # Expressed in LY
      v_te_notrt  <- (l_markov_dyn$m_M_notrt   %*%  v_e_notrt) * n_YpC # Expressed in LY
      
      v_tu_trt    <- (l_markov_dyn$m_M_trt     %*%  v_u_trt)   * n_YpC # Expressed in QALY
      v_tu_notrt  <- (l_markov_dyn$m_M_notrt   %*%  v_u_notrt) * n_YpC # Expressed in QALY
      
      

      for (t in 1:n_t){
        # element-wise-multiplication of array A with the rewards matrices for costs
        
        a_Y_costs_trt[, , t + 1]   <- l_markov_dyn$a_A_trt[, , t + 1] * l_a_rewards$a_R_costs_trt[, , t + 1] 
        a_Y_costs_notrt[, , t + 1] <- l_markov_dyn$a_A_notrt[, , t + 1] * l_a_rewards$a_R_costs_notrt[, , t + 1]  
      }
      
      
      # sum the costs to get the cost per cycle for the strategy
      v_tc_trt             <- rowSums(t(colSums(a_Y_costs_trt)))
      v_tc_notrt           <- rowSums(t(colSums(a_Y_costs_notrt)))

      
      te_d_trt    <- t(v_te_trt)   %*%  v_dwe # LY
      te_d_notrt  <- t(v_te_notrt) %*%  v_dwe # LY   
      
      tu_d_trt    <- t(v_tu_trt)   %*%  v_dwe  # QALY
      tu_d_notrt  <- t(v_tu_notrt) %*%  v_dwe  # QALY 
      
      tc_d_trt    <- t(v_tc_trt)   %*%  v_dwc  # costs
      tc_d_notrt  <- t(v_tc_notrt) %*%  v_dwc  # costs
      
      
      # store them into a vector, of cost per year and effect expressed in LY
      v_tc      <- c(tc_d_notrt, tc_d_trt)         # total costs in a year
      v_te      <- c(te_d_notrt, te_d_trt)         # express in LY per years
      v_tu      <- c(tu_d_notrt, tu_d_trt)         # express in QALY
      names(v_tc) <- names(v_te) <- names(v_tu) <- v_names_str
      
      # Dataframe with discounted costs and effectiveness in LY
      df_ce       <- data.frame(Strategy = v_names_str,
                                Cost     = v_tc,
                                Effect   = v_te) # effects are lifeyears (NOT QALYs)
      
      # Dataframe with discounted costs and effectiveness in QALYs
      df_cu       <- data.frame(Strategy = v_names_str,
                                Cost     = v_tc,
                                Effect   = v_tu) # effects are in utility (QALYs)
      
      # combine results of utility and effects in one dataframe 
      df_results   <- data.frame(Strategy = v_names_str,
                                   Cost   = v_tc,
                                   LY     = v_te, 
                                   QALYs  = v_tu) 
  
      ## Vector with discounted net monetary benefits (NMB) - QALY
      v_nmb    <- (v_tu * n_wtp) - v_tc
      
      ## Vector with discounted net health benefits (NHB) - QALY
      v_nhb    <- v_tu - (v_tc/n_wtp)
      
      ## Vector with discounted net monetary benefits (NMB) - LY
      v_nmb_LY <- (v_te * n_wtp) - v_tc
      
      ## Vector with discounted net health benefits (NHB) - LY
      v_nhb_LY <- v_te - (v_tc/n_wtp)
      
      ## Dataframe with discounted costs, effectiveness and NMB
      df_ce_combined <- data.frame(Strategy = v_names_str,
                                   Cost     = v_tc,
                                   Effect   = v_tu,
                                   LY       = v_te,
                                   NMB      = v_nmb,
                                   NHB      = v_nhb,
                                   NMB_LY   = v_nmb_LY,
                                   NHB_LY   = v_nhb_LY)
      
      #### Create a summary data frame  #####
    
      l_out_all <- list(l_param_trt_f = l_param_trt_f,
                       results_trt    = results_trt,
                       results_notrt  = results_notrt,
                       l_trace        = l_trace, 
                       l_array        = l_array,
                       v_tc_notrt     = v_tc_notrt,
                       v_tc_trt       = v_tc_trt,
                       df_ce          = df_ce, 
                       df_cu          = df_cu,
                       df_results     = df_results,
                       df_ce_combined = df_ce_combined)
       
      # dataframe with combined CE results
       l_out_ce <- df_ce_combined
       
       # if the Return all is true, report the full information
       if(return_all == TRUE){
         l_out_ce <-  l_out_all
       }
       
       #### RETURN OUTPUT #####
    return(l_out_ce) 
  })
}


