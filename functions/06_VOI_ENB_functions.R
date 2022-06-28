#' Functions to calculate expected net benefit for the different quadrants 
#' The reference strategy is the reject strategy 
#'
#' \code{calc_enb_approve} implements the decision model used.
#'
#' @param iNB Net benefit (NB) of treatment minus net benefit (NB) of control. 
#' @param n_patients_current All patients until new research results lead to approval.
#' @param n_patients_future All patients as of decision based on new research results.
#' @param c_implementation Irrecoverable costs of implementation following approval. Default = 0; Assumed to be negligible for COVID therapies
#' @return Expected net benefit (ENB) of the overall strategy 


calc_enb_approve <- function(iNB, n_patients_current, n_patients_future, c_implementation = 0){
  
 enb  <- iNB * n_patients_current +
         iNB * n_patients_future -
          c_implementation
 enb <- round(enb)
 return(enb)
  
}

#'
#' \code{calc_enb_AWR} implements the decision model used.
#'
#' @param iNB Net benefit (NB) of treatment minus net benefit (NB) of control.
#' @param rct_ratio Proportion of patients in RCT randomly assigned to new therapy
#' @param sample_size Sample size of the RCT. Also possible to give a vector of sample size 
#' @param EVPSI Expected value of sample information of the sample size. Also possible to give a vector to calculate for multiple sample size
#' @param n_patients_current All patients until new research results lead to approval
#' @param n_patients_future All patients as of decision based on new research results
#' @param c_RCT_fixed Fixed costs of an RCT
#' @param c_RCT_pp Per person costs of an RCT
#' @param c_implementation Irrecoverable costs of implementation following approval. Default = 0; Assumed to be negligible for COVID therapies
#' @param c_reversal Costs required to undo the implementation. Default = 0; Assumed to be negligible for COVID therapies
#' @return Expected net benefit (ENB) of the overall strategy 



calc_enb_AWR <- function(iNB, rct_ratio, sample_size, EVPSI, n_patients_current, n_patients_future, c_RCT_fixed, c_RCT_pp, c_implementation = 0, c_reversal = 0){
 
   enb <- iNB * (n_patients_current - (1 - rct_ratio) * sample_size) +
       (iNB + EVPSI) * n_patients_future -
       (c_RCT_fixed + sample_size * c_RCT_ppo) -
       c_implementation - c_reversal
   enb <- round(enb)
  return(enb)
}



#' \code{calc_enb_OIR} implements the decision model used.
#'
#' @param iNB Net benefit (NB) of treatment minus net benefit (NB) of control.
#' @param rct_ratio Proportion of patients in RCT randomly assigned to new therapy
#' @param sample_size Sample size of the RCT. Also possible to give a vector of sample size 
#' @param EVPSI Expected value of sample information of the sample size. Also possible to give a vector to calculate for multiple sample size
#' @param n_patients_future All patients as of decision based on new research results
#' @param c_RCT_fixed Fixed costs of an RCT
#' @param c_RCT_pp Per person costs of an RCT
#' @return Expected net benefit (ENB) of the overall strategy 

calc_enb_OIR <- function(iNB, rct_ratio, sample_size, EVPSI, n_patients_future, c_RCT_fixed, c_RCT_pp){
  
  enb <-  iNB * rct_ratio * sample_size  +
    ((if_else(iNB <= 0, 0, iNB) + EVPSI) * n_patients_future) -
    (c_RCT_fixed + sample_size * c_RCT_ppo) 
    
  enb <- round(enb)
  
  return(enb)
}
