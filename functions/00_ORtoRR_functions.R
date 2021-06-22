# Odds Ratio to Relative Risk @EK newly added by stijntje

# Example https://clincalc.com/Stats/ConvertOR.aspx?example
# RR = risk ratio, OR = Odds ratio, P_ref = prevalence outcome in reference group
# RR = prob exposed/prob non-exposed, OR = odds exposed/odds non-exposed

# OR <- 0.51
# P_ref <- 0.93

# RR <- OR / ( (1-P_ref) + (P_ref* OR) ) 
# RR # 0.94



#-------------------------------------------------------------------#
####   Calculates a risk ratio from an odds ratio  ####
#-------------------------------------------------------------------#
#'
#' \code{OR_to_RR} Calculates a risk ratio from an odds ratio
#' 
#' @param OR    The odds ratio
#' @param p_prev_ref  Prevalence of the outcome in the reference group
#' @return 
#' The RR, the risk ratio
#' @import 
#' @export


OR_to_RR <- function(OR, p_prev_ref){
  # Argument 
  RR <- OR / ((1 - p_prev_ref) + (p_prev_ref * OR))
  return(RR)
}

# Zhang J, Yu KF. What's the relative risk? A method of correcting the odds ratio in cohort studies of common outcomes. JAMA. 1998 Nov 18;280(19):1690-1. doi: 10.1001/jama.280.19.1690. PMID: 9832001.