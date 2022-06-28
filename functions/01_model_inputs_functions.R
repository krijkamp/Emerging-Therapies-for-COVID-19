# This file contains all the functions that are needed in the data generation process of the model for treatment of COVID-19 



#-------------------------------------------------------------------#
#### Function remove date from a name  ####
#-------------------------------------------------------------------#
#' Remove characters at the end of a name  
#'
#' \code{remove_date_from_name} remove if unit is missing 
#' 
#' @param name A character of the name
#' @param date The character of the data structure
#' @return 
#' The character without the part that had to be removed 

remove_date_from_name <- function(name, date = "_01_01_2020"){
  
  n_char_name_trt <- nchar(name)
  n_char_data     <- nchar(date)
  n_new_name_trt  <- substr(name, 1, (n_char_name_trt - n_char_data))
  return(n_new_name_trt)

  }


#-------------------------------------------------------------------#
#### Function remove values from the data that do not have a unit ####
#-------------------------------------------------------------------#
#' Remove lines of code where there is no unit value in the data. 
#'
#' \code{remove_param_without_unit} remove if unit is missing 
#' 
#' @param df_param A dataframe with parameter values for the model
#' @return 
#' The dataframe without the parameter line for which the unit was missing
#' @import 
#' @export

remove_param_without_unit <- function(df_param){
  
  # Document which parameters do not have a unit
  param_NA <- df_param[which(is.na(df_param$Unit)), ]
  
  # Remove parameters without unit
  df_param <- df_param[!is.na(df_param$Unit), ]
  
  # Make variables numeric
  numvar <- c("Med", "Lo_alpha", "Hi_beta") # give it the names of the Excel columns
  for (i in numvar) {df_param[, i] <- as.numeric(df_param[, i])} # make the values numbers

  return(df_param) # Return the dataframe
}


#-------------------------------------------------------------------#
#### Function to split the dataframe in treatment specific lists with only the mean parameter value ####
#-------------------------------------------------------------------#
#' Make a list with contains for each treatment a list with the mean values for all parameters based on the information in the dataframe 
#'
#' \code{split_df_in_trt_lists}  
#' 
#' @param df_param A dataframe with parameter values for the model.
#' @return 
#' A list with a list for each treatment a list with the mean values of the parameters
#' @import 
#' @export

split_df_to_lists_mean_param <- function(df_param){
  l_param <- c() # initiate a list
  
  for (n in sort(unique(df_param$Treatment))){
    sub_param   <- filter(df_param, Treatment == n) # select the values corresponding to each treatment
    l_temporary <- as.list(sub_param$Med) # Subtract the estimated value of parameter values, save as list
    names(l_temporary) <- sub_param$Param # Name the parameters
    l_param[[n]] <- l_temporary # remove the temporary files 
  }
  return(l_param)  # return the list 
}


#-------------------------------------------------------------------#
#### Function to combine lists ####
#-------------------------------------------------------------------#
#' Make a list where all the parameters of one list are merged into all other lists 
#'
#' \code{combine_lists}  
#' 
#' @param list    A list with parameter values for the model
#' @param combine The data you like to add to each list 
#' @return 
#' A list with a list for each treatment with the general information merged into the lists
#' @import 
#' @export

combine_lists <- function(list, combine = "All"){
  
  names <- sort(names(list))  # select the names of the lists
  names_treatment <- names[-(names == combine)] # remove the name that needs to be added to all lists
  l_param <- c() # initiate the new list
  
  for (n in names_treatment){
    l_param[[n]] <- c(list[[n]], list[[combine]]) # combine the treatment specific list with the list specified via the combine variable
  }
  return(l_param) # return the list 
}


#----------------------------------------------------------------------------#
####    Function to estimate general input parameters ####
#----------------------------------------------------------------------------#
#' Estimate general input parameters for the model
#'
#' \code{estimate_input_general} Estimate general input parameters
#' 
#' @param l_list A list with all the model parameters
#' @return 
#' The a list with the general input. These elements have to bee added to the list of treatment specific parameters
#' @import 
#' @export
#'

estimate_input_general <- function(l_list){
  l_list_f <- c()
  with(as.list(l_list), {
    
    n_y     <- 120 - n_age        # time horizon in years, 120 years - mean age of 73
    
    # Number of strategies
    n_str     <- length(v_names_str)
    n_CpY     <- n_DpY / n_DpC          # number of cycles per year
    n_YpC     <- n_DpC / n_DpY          # number of years per cycle
    n_t       <- n_y * n_CpY            # number of cycles in model
    n_states  <- length(v_names_states) # number of health states 
    
    times <- seq(0, n_t/n_CpY, by = n_YpC) # factor to adjust discount weights for cycle length.
    
    # Calculate discount weights for effects for each cycle based on discount rate
    v_dwe      <- 1 / ((1 + d_e) ^ (times))    
    v_dwc      <- 1 / ((1 + d_c) ^ (times))    
    
    # Add the parameters to the lists
    l_list_f$n_str    <- n_str
    l_list_f$n_CpY    <- n_CpY
    l_list_f$n_YpC    <- n_YpC
    l_list_f$n_t      <- n_t
    l_list_f$n_states <- n_states
    l_list_f$times    <- times
    l_list_f$v_dwe    <- v_dwe
    l_list_f$v_dwc    <- v_dwc
    
    return(l_list_f) 
  })
}



#----------------------------------------------------------------------------#
####   Function to estimate treatment cost  ####
#----------------------------------------------------------------------------#
#' Estimate the cost of the treatment based on the drug price and the percentage private insurance and required doses  
#'
#' \code{estimate_cost_drug} Estimate the costs of the drug
#' 
#' @param l_list A list with all the parameter values for the model.
#' @return 
#' The new list with the new estimated drug price added to the input list
#' @import 
#' @export

estimate_cost_drug <- function(l_list){
  l_list_f <- c()
  l_list_f <- l_list

   with(as.list(l_list), {
    
    # Calculate the cost of the treatment and add it to the list 
    c_Trt <- (c_Trt_private *     p_Private_insurance + 
              c_Trt_public * (1 - p_Private_insurance)) * n_Trt
    
    l_list_f$c_Trt <- c_Trt
    out <- as.list(l_list_f)
    return(out) 
  })
}

#----------------------------------------------------------------------------#
####    Function calculate the costs parameters in the model  ####
#----------------------------------------------------------------------------#
#' Estimate the cost parameters in the model  
#'
#' \code{calc_costs_model} Calculates all costs 
#' 
#' @param l_list  A list with all the parameter values for the model.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @return 
#' The new list with new estimated costs added
#' @import 
#' @export
#' 

calc_costs_model <- function(l_list, verbose = TRUE){
  l_temporary <- c()    # initiate list 
  l_temporary <- l_list # store the input list 
  
  with(as.list(l_list), {
    
    # Calculate the cost of the treatment
    # Costs post hospital period = = in the cost per cycle in the recovery health state
    c_Alive_post <- c_Healthcare * n_YpC # Cost per cycle for being alive after the first cycle of being in the study; Calculated as healthcare mean expected healthcare costs US citizen 65+ 
    
    # Total cost of hospital stay noTrt group
    c_H_noTrt <- LOS_noTrt * c_Hospital   # Calculate the total cost of  Hospital ward for the entire LOS
    c_I_noTrt <- LOS_noTrt * (c_I_vent * p_vent + 
                              c_I_noVent * (1 - p_vent))   # Calculate the total cost of ICU for the entire LOS, adjust cost of ICU ward per day based on the probability of requiring mechanical ventilation
    
    # Total cost of hospital stay Trt group
    c_H_Trt <-  LOS_Trt * c_Hospital    # Calculate the total cost of  Hospital ward for the entire LOS
    c_I_Trt <-  LOS_Trt * (c_I_vent * p_vent + 
                          c_I_noVent * (1 - p_vent))  # # Calculate the total cost of ICU for the entire LOS, adjust cost of ICU ward per day based on the probability of requiring mechanical ventilation
    
    # Cost of first cycle (in the hospital), composed of a proportion in the hospital ward and a proportion in the ICU
    c_Hospital_mix     <- (     p_IC   *  c_I_noTrt) +  
                          ((1 - p_IC)  *  c_H_noTrt)
    
    c_Hospital_mix_trt <- (     p_IC   *  c_I_Trt)  +  
                          ((1 - p_IC)  *  c_H_Trt)
    
    
    # Add the parameters to the lists
    l_out <-  c(l_temporary,
                c_Alive_post       = c_Alive_post,
                c_H_noTrt          = c_H_noTrt,
                c_I_noTrt          = c_I_noTrt,
                c_H_Trt            = c_H_Trt,
                c_I_Trt            = c_I_Trt,
                c_Hospital_mix     = c_Hospital_mix,
                c_Hospital_mix_trt = c_Hospital_mix_trt)
   
    l_out <- as.list(l_out)
    # before storing the information - check if all items in the list contain information
    check_list_elements(l_list = l_out, err_stop = TRUE, verbose = verbose)
    return(l_out) 
  })

}


#----------------------------------------------------------------------------#
####    Function to estimate the age and gender specific probability to die after recovery ####
#----------------------------------------------------------------------------#
#' Estimate the age and gender specific mortality probability after recovery
#'
#' \code{estimate_p_RD} Estimate the age specific mortality probability
#' 
#' @param df_r_HD A dataframe with the age and gender specific mortality
#' @param n_age The age of the population
#' @param p_man The proportion of men in the population
#' @param n_YpC The fraction of years per cycle 
#' @return 
#' The age and gender specific mortality of the cohort 
#' @import 
#' @export
#' 

estimate_p_RD <- function(df_r_HD, n_age, p_men, n_YpC){
  # Mortality related transition probability (post hospitalization)
  # Calculate weighed average of the age specific mortality rate based on 2017 numbers
  # Get the vector of weighted average based on male/female for all ages
  v_r_D_age        <- df_r_HD$Female * (1 - p_men) + df_r_HD$Male * p_men
  names(v_r_D_age) <- df_r_HD$Age   # name the vector, 
  v_r_D_age        <- v_r_D_age[-1] # to avoid issues with indexing remove the age of 0 as we don't need it for this model
  
  # Make a vector of probs to die per cycle from the mortality rates 
  v_p_D_age_cycle <- darthtools::rate_to_prob(v_r_D_age, t = n_YpC) 
  p_RD            <- v_p_D_age_cycle[floor(n_age)] # The probability to die in recover at the starting age of the cohort 
  
  return(p_RD) # return the probability to die in recovery at the start fo the model
}

#----------------------------------------------------------------------------#
####    Function to estimate the age and gender specific probability to die after recovery ####
#----------------------------------------------------------------------------#
#' Estimate the age and gender specific mortality probability after recovery
#'
#' \code{estimate_v_p_RD} Estimate the age specific mortality probability
#' 
#' @param df_r_HD A dataframe with the age and gender specific mortality
#' @param n_age The age of the population
#' @param p_man The proportion of men in the population
#' @param n_YpC The fraction of years per cycle 
#' @return 
#' A vector with the age and gender specific mortality
#' @import 
#' @export
#' 

estimate_v_p_RD <- function(df_r_HD, n_age, p_men, n_YpC){
  # Mortality related transition probability (post hospitalization)
  # Calculate weighed average of the age specific mortality rate based on 2017 numbers
  # Get the vector of weighted average based on male/female for all ages
  v_r_D_age        <- df_r_HD$Female * (1 - p_men) + df_r_HD$Male * p_men
  names(v_r_D_age) <- df_r_HD$Age   # name the vector, 
  v_r_D_age        <- v_r_D_age[-1] # to avoid issues with indexing remove the age of 0 as we don't need it for this model
  
  # Make a vector of probs to die per cycle from the mortality rates 
  v_p_D_age_cycle <- darthtools::rate_to_prob(v_r_D_age, t = n_YpC) 
  p_RD            <- v_p_D_age_cycle[floor(n_age)] # The probability to die in recover at the starting age of the cohort 
  
  return(v_p_D_age_cycle) # return the probability to die in recovery at the start fo the model
}



#----------------------------------------------------------------------------#
####    Function to remove the list items that are non-numeric ####
#----------------------------------------------------------------------------#
#' Remove non numeric list items from a list 
#'
#' \code{remove_non_numeric_list_items} removes all the non numeric list items from a list. This is needed with the current version of dampack, although that is a bug and might be solved soon which makes this function less relevent 
#' 
#' @param l_list A list with all the model parameters
#' @return 
#' A list with only numeric values 
#' @import 
#' @export
#' 

remove_non_numeric_list_items <- function(l_list){
  
  m_summary <- as.matrix(summary(l_list))  # investigate the caracteristics of the list
  which(m_summary[, "Mode"] != "numeric")  # find out which are non-numeric
  
  l_param_trt_numeric_only <- l_list[-c(which(m_summary[, "Mode"] != "numeric"))] # Remove the non-numeric list items
  return(l_param_trt_numeric_only) 
  
}

#----------------------------------------------------------------------------#
####    Function to estimate effects LY and utilities ####
#----------------------------------------------------------------------------#
#' Estimate the LY and utilities in each health state
#'
#' \code{estimate_effects} Estimate effects (LY and utilities)
#' 
#' @param l_list A list with all the model parameters
#' @return 
#' The list with the effects added to the list
#' @import 
#' @export
#' 

estimate_effects <- function(l_list){
  
  with(as.list(l_list), {
    
      # weighed LY for no Trt based =  on LOS + those that survive * remaining time alive (n_DpC - LOS)
      e_H_mix_noTrt <- e_H * (LOS_noTrt/n_DpC +
                              (1 - p_HD) * (n_DpC - LOS_noTrt)/n_DpC) 
      
      # weighed LY for Trt based =  on LOS + those that survive * remaining time alive (n_DpC - LOS)
      e_H_mix_Trt   <- e_H * (LOS_Trt/n_DpC +
                                            (1 - p_HD_Rx) * (n_DpC - LOS_Trt)/n_DpC)
      

      
      # Calculation of the weighted average of utility in hospital 
      # weighted by time spend in a condition based on mean utility - treat group and time spend in that condition. 
      u_pH_Trt <- (1 - p_IC)  * u_H  * (LOS_Trt/n_DpC) # Utility of in hospital during length of stay
      u_pI_Trt <- p_IC * u_I         * (LOS_Trt/n_DpC) # Utility of in ICU during length of stay
      u_pD_Trt <- p_HD * u_D         * ((n_DpC - LOS_Trt)/n_DpC) # Utility of dying remainder of first cycle
      u_pRIC_Trt <- p_R_IC * u_R_IC  * ((n_DpC - LOS_Trt)/n_DpC) # Utility of post ICU remainder of first cycle
      u_pRH_Trt <- p_R_H * u_R_H     * ((n_DpC - LOS_Trt)/n_DpC) # Utility of post ICU remainder of first cycle
      
      u_H_mix_Trt <- u_pH_Trt + u_pI_Trt + u_pD_Trt + u_pRIC_Trt + u_pRH_Trt
      

      u_pH_noTrt <- (1 - p_IC)  * u_H  * (LOS_noTrt/n_DpC) # Utility of in hospital during length of stay
      u_pI_noTrt <- p_IC * u_I         * (LOS_noTrt/n_DpC) # Utility of in ICU during length of stay
      u_pD_noTrt <- p_HD * u_D         * ((n_DpC - LOS_noTrt)/n_DpC) # Utility of dying remainder of first cycle
      u_pRIC_noTrt <- p_R_IC * u_R_IC  * ((n_DpC - LOS_noTrt)/n_DpC) # Utility of post ICU remainder of first cycle
      u_pRH_noTrt <- p_R_H * u_R_H     * ((n_DpC - LOS_noTrt)/n_DpC) # Utility of post ICU remainder of first cycle
      
      u_H_mix_noTrt <- u_pH_noTrt + u_pI_noTrt + u_pD_noTrt + u_pRIC_noTrt + u_pRH_noTrt

    
    l_out <- c(e_H_mix_noTrt  = e_H_mix_noTrt,
               e_H_mix_Trt    = e_H_mix_Trt,
               u_H_mix_noTrt  = u_H_mix_noTrt,
               u_H_mix_Trt    = u_H_mix_Trt)  # return the list 
    l_out <- as.list(l_out)
    return(l_out) 
  })
}




#-------------------------------------------------------------------#
####   FUNCTIONS FOR DARTHTOOLS ####
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
####   Function to check if all the items in the list contains information  ####
#-------------------------------------------------------------------#
#' Check if the each item in the list contains information. 
#'
#' \code{check_list_items} checks if item in the list contains a value 
#' 
#' @param l_param A list with parameter values.
#' @param err_stop Logical variable to stop model run if set up as TRUE. 
#' Default = TRUE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' Default = TRUE
#' @return 
#' Information about the list
#' @import 
#' @export
check_list_items <- function(l_param,
                             err_stop = TRUE, 
                             verbose  = TRUE) {
  
  # check if each component of the list is valid, not NULL, not NA etc.
  valid <- !is.null(l_param) & class(l_param) != "NULL" & class(l_param) != "logical" & length(l_param) != 0 & !is.na(l_param)
  
  m_valid <- count(valid) # sum the logic vale for each parameters
  
  if(m_valid$freq[m_valid$x == TRUE] != length(names(l_param))) {
    if(err_stop) {
      stop("This is not a valid list. At least one parameter does not contain information.")
    }
    
    if(verbose){
      warning("This is not a valid list. At least one parameter does not contain information.")
    }
    
  }
  else if (verbose){
    print("This is a valid list")}
} # close the function


#' Check if the each item in the list contains information.
#'
#' \code{check_list_elements} checks if item in the list contains a value
#'
#' @param l_list A list with parameter values.
#' @param err_stop Logical variable to stop model run if set up as TRUE.
#' Default = TRUE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @return
#' Information about the validity of the list
#' @import
#' @export
check_list_elements <- function(l_list,
                                err_stop = TRUE,
                                verbose  = TRUE) {
  
  # check if each component of the list is valid, not NULL, not NA etc.
  valid <- !is.null(l_list) & class(l_list) != "NULL" & class(l_list) != "logical" & class(l_list) == "list" & length(l_list) != 0 & sum(!is.na(l_list)) == length(l_list) & sum(!sapply(l_list, is.null)) == length(l_list)
  
  
  if (valid == TRUE & verbose == TRUE) {
    print("This is a valid list")
  } else if (valid == FALSE & err_stop == TRUE) {
    stop("This is not a valid list. At least one element in the list does not contain information.")
  } else if (valid == FALSE & verbose == TRUE){
    warning("This is not a valid list. At least one element in the list does not contain information.")
  }
  
  
} # close the function

