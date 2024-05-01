# This R file contains the model specific function for the VOI for COVID-therapy







#-------------------------------------------------------------------#
#### Function that splits the first cycle in separate time points  ####
#-------------------------------------------------------------------#
#'Split the first cycle (n_DpC) in multiple time terms to allow adjustment for early mortality and different treatment effects 
#'
#' \code{split_firstcycle} splits the first cycle of the model in separate time points based on the duration in which the treatment effect is applicable 
#' 
#' @param l_list  The list with all the model parameters
#' @return 
#' A list where the time points t1, t2 and t3 are added 
#' @import 
#' @export
#' 
split_firstcycle <- function(l_list){
  # The names of all the lists you have in the environment to get data from
  l_list_f <- l_list
  l_out <- c()
  with(as.list(l_list), {

    t1 <- t_Die_mj  # The first period is from the cohort when the majority 
  
    if (n_days_timespan1 > n_DpC){
      n_days_timespan1 <- n_DpC }
    
    t2 <- max(n_days_timespan1  - t1, 0)  # duration of the RCT - period 
    t3 <- max(n_DpC - t1 -  t2, 0)        # cycle - other time period 
    
    # Test if adds up to total days in a cycle
    if (sum(t1 + t2 + t3) != n_DpC) {
      warning('Splitting the first cycle in terms did not work')
    }
    
    l_out <-  c(l_list,
          t1 = t1, 
          t2 = t2, 
          t3 = t3)
    return(l_out)
})

}



# Make a sub-function for each part 
# If we have a treatment effect that is a relative risk RR
apply_trt_RR <- function(l_list, vent = FALSE){
  with(as.list(l_list), {
    
    if(vent == FALSE) {
      # if there is no ventilation information - temporary create values
      rr_D_Trt_timespan1_vent <- rr_D_Trt_timespan1
      venttype <- "no distinction in applied RR to patients with or without ventilation" #overwrite
      rr_D_Trt_timespan1_novent <- rr_D_Trt_timespan1
      # after trial RR         - mechanical vent
      rr_D_Trt_timespan2_vent   <- rr_D_Trt_timespan2
      rr_D_Trt_timespan2_novent <- rr_D_Trt_timespan2
      
      p_R_IC_D_vent   <- p_R_IC_D  
      p_R_IC_D_novent <- p_R_IC_D
    } 
    
    # Calculate transition probabilities for effect type RR   
    # Probability of dying in ICU in treat group for RR
    p_R_IC_D_Rx <- 
      (p_vent  *  p_R_IC_D_vent    * rr_D_Trt_timespan1_vent     *  p_Die_1 + 
       p_vent  *  p_R_IC_D_vent    * rr_D_Trt_timespan1_vent     * (p_Die_2 * (t2/(t2+t3)))  + 
       p_vent  *  p_R_IC_D_vent    * rr_D_Trt_timespan2_vent     * (p_Die_2 * (t3/(t2+t3))) ) +  
 ((1 - p_vent) *  p_R_IC_D_novent  * rr_D_Trt_timespan1_novent   *  p_Die_1 + 
  (1 - p_vent) *  p_R_IC_D_novent  * rr_D_Trt_timespan1_novent   * (p_Die_2 * (t2/(t2+t3))) + 
  (1 - p_vent) *  p_R_IC_D_novent  * rr_D_Trt_timespan2_novent   * (p_Die_2 * (t3/(t2+t3))) )
    

    # Probability dying in the hospital ward in treat group for RR
    p_R_H_D_Rx <- p_R_H_D * rr_D_Trt_timespan1 *  p_Die_1 + 
                  p_R_H_D * rr_D_Trt_timespan1 * (p_Die_2 * (t2/(t2+t3))) +
                  p_R_H_D * rr_D_Trt_timespan2 * (p_Die_2 * (t3/(t2+t3)))   
    
    ### REMOVE  after code is ocrrect ### 
    # @EK @SD
    if(p_R_IC_D_Rx>1){
      p_R_IC_D_Rx <- 1
    }
    
    # @EK -> solve bug when p>1 but need to be done better
    if(p_R_H_D_Rx>1){
      p_R_H_D_Rx <- 1
    }
    
    
    ###### RETURN OUTPUT  
    l_list$p_R_H_D_Rx      <- p_R_H_D_Rx
    l_list$p_R_IC_D_Rx     <- p_R_IC_D_Rx
    l_list$p_R_IC_D_vent   <- p_R_IC_D_vent
    l_list$p_R_IC_D_novent <- p_R_IC_D_novent
    out <- as.list(l_list)
    return(out) 
  })
}

# If we have a treatment effect that is a risk difference RD
apply_trt_RD <- function(l_list, vent = FALSE){
  with(as.list(l_list), {
    
    if(vent == FALSE) {
      # if there is no ventilation information - temporary create values
      rd_D_Trt_timespan1_vent <- rd_D_Trt_timespan1
      venttype <- "no distinction in applied RD to patients with or without ventilation" #overwrite
      rd_D_Trt_timespan1_novent <- rd_D_Trt_timespan1
      # after trial RR         - mechanical vent
      rd_D_Trt_timespan2_vent   <- rd_D_Trt_timespan2
      rd_D_Trt_timespan2_novent <- rd_D_Trt_timespan2
      
      p_R_IC_D_vent   <- p_R_IC_D  
      p_R_IC_D_novent <- p_R_IC_D
    } 
  
    # Calculate transition probabilities for effect type RD   
    # Probability of dying in ICU in treat group for RD
    p_R_IC_D_Rx <- 
        (p_vent  *  (p_R_IC_D_vent   + rd_D_Trt_timespan1_vent  ) *  p_Die_1 + 
         p_vent  *  (p_R_IC_D_vent   + rd_D_Trt_timespan1_vent  ) * (p_Die_2 * (t2/(t2+t3)))  + 
         p_vent  *  (p_R_IC_D_vent   + rd_D_Trt_timespan2_vent  ) * (p_Die_2 * (t3/(t2+t3))) ) +  
   ((1 - p_vent) *  (p_R_IC_D_novent + rd_D_Trt_timespan1_novent) *  p_Die_1 + 
    (1 - p_vent) *  (p_R_IC_D_novent + rd_D_Trt_timespan1_novent) * (p_Die_2 * (t2/(t2+t3))) + 
    (1 - p_vent) *  (p_R_IC_D_novent + rd_D_Trt_timespan2_novent) * (p_Die_2 * (t3/(t2+t3))) )
    
    # Probability dying in the hospital ward in treat group for RD
    p_R_H_D_Rx <- (p_R_H_D + rd_D_Trt_timespan1) *  p_Die_1 + 
                  (p_R_H_D + rd_D_Trt_timespan1) * (p_Die_2 * (t2/(t2+t3))) +
                  (p_R_H_D + rd_D_Trt_timespan2) * (p_Die_2 * (t3/(t2+t3)))  
    
    
    ###### RETURN OUTPUT  
    l_list$p_R_H_D_Rx      <- p_R_H_D_Rx
    l_list$p_R_IC_D_Rx     <- p_R_IC_D_Rx
    l_list$p_R_IC_D_vent   <- p_R_IC_D_vent
    l_list$p_R_IC_D_novent <- p_R_IC_D_novent
    out <- as.list(l_list)
    
    return(out) 
  })
}

# If we have a treatment effect that is a hazard rate HR
apply_trt_HR <- function(l_list, vent = FALSE){
  with(as.list(l_list), {
    # Else part of function
    if(vent == FALSE) {
      
      hr_D_Trt_timespan1_vent   <- hr_D_Trt_timespan1 #  mechanical vent
      hr_D_Trt_timespan1_novent <- hr_D_Trt_timespan1 #  oxygen only
      hr_D_Trt_timespan2_vent   <- hr_D_Trt_timespan2 # mechanical vent
      hr_D_Trt_timespan2_novent <- hr_D_Trt_timespan2 ##  oxygen only
      #null - do not use the weighted probs to calculated the transition probability
      p_R_IC_D_vent  <- p_R_IC_D  
      p_R_IC_D_novent <- p_R_IC_D
    } 
    
    # Calculate transition probabilities for effecttype HR     
    ## Number of people dying per time interval - with ventilation
    n_D_t1_vent      <- p_R_IC_D_vent * p_Die_1
    n_D_t2_vent      <- p_R_IC_D_vent * (1 - p_Die_1) * (t2/(t2+t3))
    
    ## Proportion of people dead in each time interval
    n_D_t12_vent <- n_D_t1_vent + n_D_t2_vent
    n_D_t3_vent  <- p_R_IC_D_vent * (1 - p_Die_1) * (t3/(t2+t3))
    
    ## Probability to die in each time interval
    p_D_t12_vent <- n_D_t12_vent
    p_D_t3_vent  <- n_D_t3_vent / (1 - p_D_t12_vent)
    
    ## Rate to die in each time interval
    r_D_t12_vent <- prob_to_rate(p = p_D_t12_vent, t = (t1 + t2)) 
    r_D_t3_vent  <- prob_to_rate(p = p_D_t3_vent,  t = t3)
    ## Apply HR as treatment effect & calculate those that die
    n_D_t12_trt_vent <- (1 - exp(-r_D_t12_vent * hr_D_Trt_timespan1_vent * (t1+t2))) * 1 
    n_D_t3_trt_vent  <- (1 - exp(-r_D_t3_vent  * hr_D_Trt_timespan2_vent * t3))      * (1 - n_D_t12_trt_vent)
    ## Probability to die total period vent and Rx (p_R_IC_D_vent_1_Rx)
    p_R_IC_D_vent_1_Rx <- n_D_t12_trt_vent + n_D_t3_trt_vent
    
    # Calculate p_R_IC_D_novent_Rx
    ## Number of people dying per time interval
    n_D_t1_novent      <- p_R_IC_D_novent * p_Die_1
    n_D_t2_novent      <- p_R_IC_D_novent * (1 - p_Die_1) * (t2/(t2+t3))
    ## Proportion of people dead in each time interval
    n_D_t12_novent <- n_D_t1_novent + n_D_t2_novent
    n_D_t3_novent  <- p_R_IC_D_novent * (1 -p_Die_1) * (t3/(t2+t3))
    ## Probability to die in each time interval
    p_D_t12_novent <- n_D_t12_novent
    p_D_t3_novent  <- n_D_t3_novent / (1 - p_D_t12_novent)
    ## Rate to die in each time interval
    r_D_t12_novent <- prob_to_rate(p = p_D_t12_novent, t = (t1 + t2)) 
    r_D_t3_novent  <- prob_to_rate(p = p_D_t3_novent,  t = t3)
    ## Apply HR as treatment effect  & calculate # of dying
    n_D_t12_trt_novent <- (1 - exp(-r_D_t12_novent * hr_D_Trt_timespan1_novent * (t1+t2))) * 1 
    n_D_t3_trt_novent  <- (1 - exp(-r_D_t3_novent  * hr_D_Trt_timespan2_novent * t3))      * (1 - n_D_t12_trt_novent)
    ## Probability to die total period vent and Rx (p_R_IC_D_vent_1_Rx)
    p_R_IC_D_novent_1_Rx <- n_D_t12_trt_novent + n_D_t3_trt_novent
    
    # Calculate probability to die ICU Rx
    p_R_IC_D_Rx <- (p_vent * p_R_IC_D_vent_1_Rx) + 
      ( (1 - p_vent) * p_R_IC_D_novent_1_Rx )       
    
    ## Number of people dying per time interval
    n_D_t1      <- p_R_H_D * p_Die_1
    n_D_t2      <- p_R_H_D * (1 - p_Die_1) * (t2/(t2+t3))
    
    ## Proportion of people dead in each time interval
    n_D_t12 <- n_D_t1 + n_D_t2
    n_D_t3  <- p_R_H_D * (1 - p_Die_1) * (t3/(t2+t3))
    ## Probability to die in each time interval
    p_D_t12 <- n_D_t12
    p_D_t3  <- n_D_t3 / (1 - p_D_t12)
    ## Rate to die in each time interval
    r_D_t12 <- prob_to_rate(p = p_D_t12, t = (t1 + t2)) 
    r_D_t3  <- prob_to_rate(p = p_D_t3,  t = t3)
    ## Apply HR as treatment effect
    n_D_t12_trt <- (1 - exp(-r_D_t12 * hr_D_Trt_timespan1 * (t1+t2))) * 1 
    n_D_t3_trt  <- (1 - exp(-r_D_t3  * hr_D_Trt_timespan2 * t3))      * (1 - n_D_t12_trt)
    ## Probability to die total period vent and Rx (p_R_H_D_1_Rx)
    p_R_H_D_Rx <- n_D_t12_trt + n_D_t3_trt         
    
    # Calculate probability to die ICU Rx
    p_R_H_D_Rx <- p_R_H_D_Rx  
    
    #### RETURN OUTPUT  
    l_list$p_R_H_D_Rx      <- p_R_H_D_Rx
    l_list$p_R_IC_D_Rx     <- p_R_IC_D_Rx
    l_list$p_R_IC_D_vent   <- p_R_IC_D_vent
    l_list$p_R_IC_D_novent <- p_R_IC_D_novent
    out <- as.list(l_list)
    
    return(out) 
  })
}

# If we have a treatment effect that is a odds ratio OD 
apply_trt_OR <- function(l_list, vent = FALSE){
  with(as.list(l_list), {
    
    if(vent == FALSE) {
      # if there is no ventilation information - temporary create values
      or_D_Trt_timespan1_vent <- or_D_Trt_timespan1
      venttype <- "no distinction in applied RR to patients with or without ventilation" #overwrite
      or_D_Trt_timespan1_novent <- or_D_Trt_timespan1
      # after trial RR - mechanical vent
      or_D_Trt_timespan2_vent   <- or_D_Trt_timespan2
      or_D_Trt_timespan2_novent <- or_D_Trt_timespan2
      
      p_R_IC_D_vent   <- p_R_IC_D  
      p_R_IC_D_novent <- p_R_IC_D
    } 
    
    # Calculate transition probabilities for effect type RD   
    # Probability of dying in ICU in treat group for RD
    
    
    p_R_IC_D_Rx <- 
      (p_vent  *  odds_to_prob(prob_to_odds(p_R_IC_D_vent) * or_D_Trt_timespan1_vent) *  p_Die_1 + 
       p_vent  *  odds_to_prob(prob_to_odds(p_R_IC_D_vent) * or_D_Trt_timespan1_vent) * (p_Die_2 * (t2/(t2+t3)))  + 
       p_vent  *  odds_to_prob(prob_to_odds(p_R_IC_D_vent) * or_D_Trt_timespan2_vent) * (p_Die_2 * (t3/(t2+t3))) ) +  
 ((1 - p_vent) *  odds_to_prob(prob_to_odds(p_R_IC_D_novent) * or_D_Trt_timespan1_novent) *  p_Die_1 + 
  (1 - p_vent) *  odds_to_prob(prob_to_odds(p_R_IC_D_novent) * or_D_Trt_timespan1_novent) * (p_Die_2 * (t2/(t2+t3))) + 
  (1 - p_vent) *  odds_to_prob(prob_to_odds(p_R_IC_D_novent) * or_D_Trt_timespan2_novent) * (p_Die_2 * (t3/(t2+t3))) )
    
    # Probability dying in the hospital ward in treat group for RD
    p_R_H_D_Rx <- odds_to_prob(prob_to_odds(p_R_H_D) * or_D_Trt_timespan1) *  p_Die_1 + 
                  odds_to_prob(prob_to_odds(p_R_H_D) * or_D_Trt_timespan1) * (p_Die_2 * (t2/(t2+t3))) +
                  odds_to_prob(prob_to_odds(p_R_H_D) * or_D_Trt_timespan2) * (p_Die_2 * (t3/(t2+t3)))  
    
  
    ###### RETURN OUTPUT  
    l_list$p_R_H_D_Rx      <- p_R_H_D_Rx
    l_list$p_R_IC_D_Rx     <- p_R_IC_D_Rx
    l_list$p_R_IC_D_vent   <- p_R_IC_D_vent
    l_list$p_R_IC_D_novent <- p_R_IC_D_novent
    
    out <- as.list(l_list)
    
    return(out) 
  })
}




# function to run for list
# Applying the appropriate treatment effect measures
# Step 1: Check whether analysis concerns a hazard rate or a relative risk for treatment effect
# Step 2: If it concerns a Hazard Rate, Then check if there is a vent/novent distinction
# Step 3: Once the model knows the relevant RR or HR parameters (conditional on vent/novent or no subdivision), then it needs to continue the calculation with the appropriate measure 




estimate_trt_parameters <- function(l_list){
  l_list_f <- c()
  l_list_f <- l_list
  
  # Determine the treatment effect
  
  v_names_param <- names(l_list_f)
  
  v_trt_type_name <- c("rr_D_Trt_timespan1", 
                       "rd_D_Trt_timespan1", 
                       "or_D_Trt_timespan1", 
                       "hr_D_Trt_timespan1")
  v_trt_type <- c("RR",
                  "RD",
                  "OR",
                  "HR")
  # Make a dateframe
  m_trt_type <- as.data.frame(cbind(v_trt_type_name, v_trt_type))
  # select the full name of the treatment type
  name <- v_names_param[ which(v_names_param %in% v_trt_type_name)]
  # match it with the corresponding treatment type
  trt_effect <- m_trt_type$v_trt_type[m_trt_type$v_trt_type_name == name] #
  
  
  vent <- FALSE  # make a variable that by default is FALSE
  if(trt_effect == "RR" & !is.null(l_list_f$rr_D_Trt_timespan1_vent)){
    vent <- TRUE # if  we have information about ventilation, overwrite
  }
  if(trt_effect == "RD" & !is.null(l_list_f$rd_D_Trt_timespan1_vent)){
    vent <- TRUE
  }
  if(trt_effect == "OR" & !is.null(l_list_f$or_D_Trt_timespan1_vent)){
    vent <- TRUE
  }
  if(trt_effect == "HR" & !is.null(l_list_f$hr_D_Trt_timespan1_vent)){
    vent <- TRUE
  }
  
  # Use the correct treatment effect function
  if (trt_effect == "RR"){
    l_list_f <- apply_trt_RR(l_list = l_list_f, vent = vent)
  } else if (trt_effect == "HR"){
    l_list_f <- apply_trt_HR(l_list = l_list_f, vent = vent)
  } else if (trt_effect == "RD"){
    l_list_f <- apply_trt_RD(l_list = l_list_f, vent = vent)
  } else if (trt_effect == "OR"){
    l_list_f <- apply_trt_OR(l_list = l_list_f, vent = vent)
  }
  
  with(as.list(l_list_f), {
    # Probability of dying in ICU in no treat group; 
    # the same for both RR and HR effect groups as no treatment effect is applied, therefore, outside of the loop
    p_R_IC_D <- 
      (p_vent  * p_R_IC_D_vent   *  p_Die_1 + 
         p_vent  * p_R_IC_D_vent   * (p_Die_2 * (t2/(t2 + t3)))  + 
         p_vent  * p_R_IC_D_vent   * (p_Die_2 * (t3/(t2 + t3))) ) +  
      ((1 - p_vent) * p_R_IC_D_novent *  p_Die_1 + 
         (1 - p_vent) * p_R_IC_D_novent * (p_Die_2 * (t2/(t2 + t3))) + 
         (1 - p_vent) * p_R_IC_D_novent * (p_Die_2 * (t3/(t2 + t3))) )  
    # Probability of dying in H in no treat group;
    # the same for both RR and HR effect groups as no treatment effect is applied, therefore, outside of the loop
    p_R_H_D <- p_R_H_D * p_Die_1 + 
      p_R_H_D * (p_Die_2 * (t2/(t2 + t3))) +
      p_R_H_D * (p_Die_2 * (t3/(t2 + t3)))   
    
    
    # Probability of being discharged alive
    # From ICU
    p_R_IC       <-      p_IC  * (1 - p_R_IC_D)    # Without Rx
    p_R_IC_Rx    <-      p_IC  * (1 - p_R_IC_D_Rx) # With Rx
    # From Hospital ward (1 - probability to be in ICU)
    p_R_H        <- (1 - p_IC) * (1 - p_R_H_D)    # without Rx
    p_R_H_Rx     <- (1 - p_IC) * (1 - p_R_H_D_Rx) # With Rx, Note: dexamethasone has no patient in IC
    # Proportion of the cohort that dies in cycle one; 
    # Compared to the proportion of people who die in the cohort study ~ 5165 / (20133-6769) 
    p_HD    <- 1 - p_R_H    - p_R_IC 
    p_HD_Rx <- 1 - p_R_H_Rx - p_R_IC_Rx
    p_RD    <- v_p_RD[n_age] # age specific prob to die 
    
    #### Return output
    l_list_f$p_R_IC     <- p_R_IC
    l_list_f$p_R_IC_Rx  <- p_R_IC_Rx
    l_list_f$p_R_H      <- p_R_H
    l_list_f$p_R_H_Rx   <- p_R_H_Rx
    l_list_f$p_HD       <- p_HD
    l_list_f$p_HD_Rx    <- p_HD_Rx
    l_list_f$p_RD       <- p_RD
    l_list_f$p_R_IC_D   <- p_R_IC_D
    l_list_f$p_R_H_D    <- p_R_H_D
    l_list_f$trt_effect <- trt_effect
    l_list_f$vent       <- vent 
    
    # ISSUE: some parameters are in the list - have to be replaces 
    
    out <- as.list(l_list_f)
    
    return(out) 
  })
}



old_estimate_trt_parameters <- function(l_list){
    l_list_f <- c()
    l_list_f <- l_list
  
  # Determine the treatment effect
  if(is.null(l_list_f$hr_D_Trt_timespan1)){ # if we don't have hr info = RR
    trt_effect <- c("RR")
  } else { #else if the list does contain infor about hr = HR
    trt_effect <- c("HR")
  }
  
  vent <- FALSE  # make a variable that by default is FALSE
  if(trt_effect == "RR" & !is.null(l_list_f$rr_D_Trt_timespan1_vent)){
    vent <- TRUE # if  we have information about ventilation, overwrite
  }
  if(trt_effect == "HR" & !is.null(l_list_f$hr_D_Trt_timespan1_vent)){
    vent <- TRUE
  }
  
  if(trt_effect == "RR"){
    l_list_f <- apply_trt_RR(l_list = l_list_f, vent = vent)
  } else if (trt_effect == "HR"){
    l_list_f <- apply_trt_HR(l_list = l_list_f, vent = vent)
  }
  
  with(as.list(l_list_f), {
    # Probability of dying in ICU in no treat group; 
    # the same for both RR and HR effect groups as no treatment effect is applied, therefore, outside of the loop
    p_R_IC_D <- 
      (p_vent  * p_R_IC_D_vent   *  p_Die_1 + 
         p_vent  * p_R_IC_D_vent   * (p_Die_2 * (t2/(t2 + t3)))  + 
         p_vent  * p_R_IC_D_vent   * (p_Die_2 * (t3/(t2 + t3))) ) +  
      ((1 - p_vent) * p_R_IC_D_novent *  p_Die_1 + 
         (1 - p_vent) * p_R_IC_D_novent * (p_Die_2 * (t2/(t2 + t3))) + 
         (1 - p_vent) * p_R_IC_D_novent * (p_Die_2 * (t3/(t2 + t3))) )  
    # Probability of dying in H in no treat group;
    # the same for both RR and HR effect groups as no treatment effect is applied, therefore, outside of the loop
    p_R_H_D <- p_R_H_D * p_Die_1 + 
      p_R_H_D * (p_Die_2 * (t2/(t2 + t3))) +
      p_R_H_D * (p_Die_2 * (t3/(t2 + t3)))   
    
    
    # Probability of being discharged alive
    # From ICU
    p_R_IC       <-      p_IC  * (1 - p_R_IC_D)    # Without Rx
    p_R_IC_Rx    <-      p_IC  * (1 - p_R_IC_D_Rx) # With Rx
    # From Hospital ward (1 - probability to be in ICU)
    p_R_H        <- (1 - p_IC) * (1 - p_R_H_D)    # without Rx
    p_R_H_Rx     <- (1 - p_IC) * (1 - p_R_H_D_Rx) # With Rx, Note: dexamethasone has no patient in IC
    # Proportion of the cohort that dies in cycle one; 
    # Compared to the proportion of people who die in the cohort study ~ 5165 / (20133-6769) 
    p_HD    <- 1 - p_R_H    - p_R_IC 
    p_HD_Rx <- 1 - p_R_H_Rx - p_R_IC_Rx
    p_RD    <- v_p_RD[n_age] # age specific prob to die 
    
    #### Return output
    l_list_f$p_R_IC    <- p_R_IC
    l_list_f$p_R_IC_Rx <- p_R_IC_Rx
    l_list_f$p_R_H     <- p_R_H
    l_list_f$p_R_H_Rx  <- p_R_H_Rx
    l_list_f$p_HD      <- p_HD
    l_list_f$p_HD_Rx   <- p_HD_Rx
    l_list_f$p_RD      <- p_RD
    l_list_f$p_R_IC_D  <- p_R_IC_D
    l_list_f$p_R_H_D   <- p_R_H_D
    
    # ISSUE: some parameters are in the list - have to be replaces 
    
    out <- as.list(l_list_f)
    
    return(out) 
  })
}




#-------------------------------------------------------------------#
#### Function to fill the structures with rewards  ####
#-------------------------------------------------------------------#
#' This code creates and populates the transition probability matrix
#'
#' \code{fill_m_P} function that populate the transition probability matrix
#' 
#' @param l_list  The list with all the model parameters
#' @return 
#' A list with all the transition probability matrix added to the input list
#' @import 
#' @export


fill_m_P <- function(l_list){
  
  l_list_f <- l_list
  

  
  with(as.list(l_list), {
    
    n_states <- length(v_names_states)
    
    # Create the transition probability matrix
    m_P_trt <- m_P_notrt  <- matrix(0,
                                    nrow = n_states,
                                    ncol = n_states,
                                    dimnames = list(v_names_states, v_names_states)) # name the columns and rows of the matrix
    
    # Create vectors to store the names of the matrices that will be generated as well as make a sequece to index.This are the odd numbers (in case of 2 strategies) as we have to index for 1 for the no_trt and 2 for the trt
    #v_names_matrix   <- vector(length = n_str * length(names(l_list))) # 
    #v_sequence_index <- seq(from = 1, to = n_str * length(names(l_list)), by = n_str)

    
    # from Hospitalized
    m_P_notrt["H", "H"]    <- 0
    m_P_notrt["H", "R_IC"] <- p_R_IC
    m_P_notrt["H", "R_H"]  <- p_R_H
    m_P_notrt["H", "D"]    <- p_HD
    
    # from Recovered from ICU
    m_P_notrt["R_IC", "R_IC"] <- 1 - p_RD
    m_P_notrt["R_IC", "D"]    <-     p_RD 
    
    # from Recovered from Hospital
    m_P_notrt["R_H", "R_H"] <- 1 - p_RD
    m_P_notrt["R_H", "D"]   <-     p_RD 
    
    # from Dead
    m_P_notrt["D", "D"] <- 1
    
    ### Check if transition matrix is valid (i.e., each row should add up to 1)
   darthtools::check_transition_probability( a_P = m_P_notrt, verbose = TRUE)
     darthtools::check_sum_of_transition_array(a_P = m_P_notrt,  
                                              n_states = n_states,
                                              n_cycles  = n_t, 
                                             verbose = TRUE)
    
    ## TREAT GROUP ##
    # from Hospitalized
    m_P_trt["H", "H"]    <- 0
    m_P_trt["H", "R_IC"] <- p_R_IC_Rx
    m_P_trt["H", "R_H"]  <- p_R_H_Rx
    m_P_trt["H", "D"]    <- p_HD_Rx
    
    # from Recovered from ICU
    m_P_trt["R_IC", "R_IC"] <- 1 - p_RD
    m_P_trt["R_IC", "D"]    <-     p_RD 
    
    # from Recovered from Hospital
    m_P_trt["R_H", "R_H"] <- 1 - p_RD
    m_P_trt["R_H", "D"]   <-     p_RD 
    
    # from Dead
    m_P_trt["D", "D"] <- 1
    
    
    ### Check if transition matrix is valid (i.e., each row should add up to 1)
    darthtools::check_transition_probability( a_P = m_P_trt, verbose = TRUE)
    darthtools::check_sum_of_transition_array(a_P = m_P_trt, 
                                              n_states = n_states,
                                              n_cycles = n_t, 
                                              verbose = TRUE)
    
    l_list_f$m_P_trt   <- m_P_trt   # store the matrix TP for trt
    l_list_f$m_P_notrt <- m_P_notrt # store the matrix TP no trt
  
    out <- as.list(l_list_f)
    
    return(out) 
  })
}


#-------------------------------------------------------------------#
#### Function to fill the structures with rewards  ####
#-------------------------------------------------------------------#
#' This code created all the model structures needed to apply the rewards
#'
#' \code{fill_a_R} populate the array (matrix and vector) structures for the rewards
#' 
#' @param l_list  The list with all the model parameters
#' @return 
#' A list with all the reward structures 
#' @import 
#' @export


fill_a_R <- function(l_list){

  with(as.list(l_list), {
    
    n_states <- length(v_names_states)
    
    ## Create matrices to store rewards
    m_R_costs_trt <- m_R_costs_notrt  <- matrix(NA, 
                                                nrow = n_states, 
                                                ncol = n_states,
                                                dimnames = list(v_names_states, v_names_states))
      # NO TREATMENT
      # All Hospitalization-related costs are written as a transition cost
      # when moving from H to R_IC or H to R_H
      # from Hospitalized
      m_R_costs_notrt["H", "H"]        <- 0 # Cost of being in the hospital during LOS
      m_R_costs_notrt["H", "R_IC"]     <- 0 + c_Hospital_mix  + c_Alive_post + c_recovery        # Costs of having been in the ICU + costs of rehabilitation post ICU stay
      m_R_costs_notrt["H", "R_H"]      <- 0 + c_Hospital_mix  + c_Alive_post                # Costs of having been in the hospital ward
      m_R_costs_notrt["H", "D"]        <- 0 + c_Hospital_mix  
      
      # from Recovered ICU
      # The recovery costs have been added in the transition from hospital to post-ICU and are therefore not added to this section of the matrix, which only has state rewards and not transition rewards 
      
      m_R_costs_notrt["R_IC", "H"]     <- 0 + c_Alive_post
      m_R_costs_notrt["R_IC", "R_IC"]  <- 0 + c_Alive_post
      m_R_costs_notrt["R_IC", "R_H"]   <- 0 + c_Alive_post
      m_R_costs_notrt["R_IC", "D"]     <- 0 
      
      # from Recovered Hospital
      m_R_costs_notrt["R_H", "H"]      <- 0 + c_Alive_post
      m_R_costs_notrt["R_H", "R_IC"]   <- 0 + c_Alive_post
      m_R_costs_notrt["R_H", "R_H"]    <- 0 + c_Alive_post
      m_R_costs_notrt["R_H", "D"]      <- 0 
      
      # from Dead
      m_R_costs_notrt["D", "H"]        <- 0
      m_R_costs_notrt["D", "R_IC"]     <- 0
      m_R_costs_notrt["D", "R_H"]      <- 0
      m_R_costs_notrt["D", "D"]        <- c_D
      
      
      # TREAT
      # from Hospitalized
      m_R_costs_trt["H", "H"]        <- 0            
      m_R_costs_trt["H", "R_IC"]     <- 0 + c_Hospital_mix_trt + c_recovery + c_Trt  + c_Alive_post     # Costs of having been in the ICU + costs of rehabilitation post ICU stay
      m_R_costs_trt["H", "R_H"]      <- 0 + c_Hospital_mix_trt + c_Trt + c_Alive_post    # Costs of having been in the hospital ward
      m_R_costs_trt["H", "D"]        <- 0 + c_Hospital_mix_trt + c_Trt    # Weighted average ICU/H cost which is given to those that die
      
      # from Recovered ICU
      m_R_costs_trt["R_IC", "H"]     <- 0 + c_Alive_post
      m_R_costs_trt["R_IC", "R_IC"]  <- 0 + c_Alive_post
      m_R_costs_trt["R_IC", "R_H"]   <- 0 + c_Alive_post
      m_R_costs_trt["R_IC", "D"]     <- 0 
      
      # from Recovered Hospital
      m_R_costs_trt["R_H", "H"]      <- 0 + c_Alive_post
      m_R_costs_trt["R_H", "R_IC"]   <- 0 + c_Alive_post
      m_R_costs_trt["R_H", "R_H"]    <- 0 + c_Alive_post
      m_R_costs_trt["R_H", "D"]      <- 0 
      
      # from Dead
      m_R_costs_trt["D", "H"]        <- 0     
      m_R_costs_trt["D", "R_IC"]     <- 0
      m_R_costs_trt["D", "R_H"]      <- 0
      m_R_costs_trt["D", "D"]        <- c_D
      
      
      # Create array
      a_R_costs_trt <- a_R_costs_notrt <-  array(0, 
                                                 dim = c(n_states, n_states, n_t + 1), 
                                                 dimnames = list(v_names_states, v_names_states, 0:n_t))
      
      # store the costs matrix in the array. 
      # Note, now: all dimensions of time are the same matrix
      a_R_costs_trt  [, , 1:(n_t + 1)] <- m_R_costs_trt    
      a_R_costs_notrt[, , 1:(n_t + 1)] <- m_R_costs_notrt
      

      # The effects only have state rewards, therefore the vectors are enough
      # Vectors with effect 
      #@EK, right now we have the effects that apply to all from the l_params. But we like to make sure we have them all in the l_temporary in order to have them copies to the list where also all the parameters for dexametason ar etc. That makes it easier to run the correct value
      
      v_e_trt   <- c(e_H_mix_Trt,   e_R_H, e_R_IC, e_D)
      v_e_notrt <- c(e_H_mix_noTrt, e_R_H, e_R_IC, e_D)

  
      # Vectors with utility 
      v_u_trt   <- c(u_H_mix_Trt,   u_R_H, u_R_IC, u_D) 
      v_u_notrt <- c(u_H_mix_noTrt, u_R_H, u_R_IC, u_D)
      #  NB: adjustment for cycle length is done when calculating mean costs and effects
      

      l_out <- list(m_R_costs_notrt = m_R_costs_notrt, 
                 a_R_costs_notrt = a_R_costs_notrt,
                 v_e_notrt       = v_e_notrt,
                 v_u_notrt       = v_u_notrt,
                 m_R_costs_trt   = m_R_costs_trt,
                 a_R_costs_trt   = a_R_costs_trt,
                 v_e_trt         = v_e_trt,
                 v_u_trt         = v_u_trt )
      
    
      out <- l_out
    
    return(out) 
  })
}
# 


#-------------------------------------------------------------------#
#### Function to simulate the cohort in the markov model  ####
#-------------------------------------------------------------------#
#' This code runs the simulation of the cohort for the markov model 
#'
#' \code{run_markov} runs the simulation for the markov model 
#' 
#' @param m_P     Transition probability matrix 
#' @param df_r_HD  data frame with the age and gender specific mortality rates 
#' @param l_list  list with parameters 
#' @return 
#' A list with the cohort markov trace matrix (m_M) and the cohort transition dynamics-array (a_A)
#' @import 
#' @export


run_markov <- function(m_P, df_r_HD, l_list){

  with(as.list(l_list), {
    
    m_M <- matrix(NA, 
                  nrow     = n_t + 1, ncol = n_states,
                  dimnames = list(paste("cycle", 0:n_t, sep = " "), v_names_states))
    
    a_A <- array(0, 
                 dim = c(n_states, n_states, n_t + 1), 
                 dimnames = list(v_names_states, v_names_states, 0:n_t)) 
    
    m_M[1, ]  <- diag(a_A[, , 1]) <- v_s_init  # initiate first cycle of cohort trace 
    
    
    for (t in 1:n_t){     # Loop through the number of cycles
      
      # Update the transition probability matrix for age
      # Select the age specific mortality rate
      
      n_age_index <- if_else(floor(n_age + (t * n_DpC)/n_DpY) <110,  floor(n_age + (t * n_DpC)/n_DpY), 110) # update the age with the cycles, round to the lower bound, after being >110, keep age with 110 otherwise you can not index the rate anymore
      
      p_RD <- v_p_RD[n_age_index] # age specific prob to die 
      
      # update the transition probability matrix with the age specific prob of dying 
      # from Recovered from ICU, both trt and notrt
      m_P["R_IC", "R_IC"] <- 1 - p_RD
      m_P["R_IC", "D"]    <- p_RD 
      
      # from Recovered from Hospital, both trt and notrt
      m_P["R_H", "R_H"] <- 1 - p_RD
      m_P["R_H", "D"]   <- p_RD 
      
      ### Check if transition matrix is valid (i.e., each row should add up to 1)
      darthtools::check_transition_probability(a_P = m_P, verbose = TRUE)
      darthtools::check_sum_of_transition_array(a_P = m_P, 
                                                n_states = n_states,
                                                n_cycles = n_t, 
                                                verbose = TRUE)
      
      m_M[t + 1, ]   <- m_M[t, ]  %*% m_P      # estimate the Markov trace 
      # for the next cycle (t + 1)
      a_A[, , t + 1] <- diag(m_M[t, ]) %*% m_P # estimate the transition dynamics at t + 1
    } # close the loop for time
    
    return(list(m_M = m_M, a_A = a_A))
  })
} # close the function




