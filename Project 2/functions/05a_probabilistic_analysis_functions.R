# The code used where the PSA function for lognormal was calculated using the mean and high



#' Generate data for the iterations of  probabilistic sensitivity analysis 
#' 
#' \code{make_psa_df} is used to sample the values for each iteration of a probabilistic sensitivity analysis (PSA) and stores it in a dataframe
#' @param param
#' @param n_iter
#' @param seed
#' @keywords Probabilistic sensitivity analysis, PSA, dataframe
#' @section Details:
#' \code{make_psa_df} adds the values for each PSA iteration to the dataframe
#' @return  param_psa A dataframe with the input values of all parameters of each PSA run. 
#' 
# Make a function to make the PSA data set 
make_psa_df <- function(param, n_iter, seed = 123){
  # Arguments:
  ## param: a dataframe with the paramters
  ## n_iter: the number of PSA iterations
  ## seed : seed to be able to reproduce the results
  # Return:
  ## param_psa: dataframe with estimated PSA parameters
  
  set.seed(seed) # set the seed
  
  # Draw samples for PSA
  param_psa         <- as.data.frame(lapply(param, rep, n_iter))
  param_psa$psa_est <- NA
  param_psa$iter    <- rep(1:n_iter, each = nrow(param)) 
  
  for(i in 1:nrow(param)){    # loop over all parameters, in all diseases
    if(param[i, "Distribution"] == "Triangular"){
      param_psa$psa_est[param_psa$Treatment == param$Treatment[i] & 
                          param_psa$Param   == param$Param[i]] <- with(param[i, ], 
                                                                       rtriangle(n = n_iter,
                                                                                 a = Lo_alpha, 
                                                                                 b = Hi_beta, 
                                                                                 c = Med))
      
    }else{ #if distribution is not triangle, check normal
      if(param[i, "Distribution"] == "Normal"){
        param_psa$psa_est[param_psa$Treatment == param$Treatment[i] & 
                            param_psa$Param      == param$Param[i]] <- with(param[i, ], 
                                                                            rnorm(n = n_iter ,
                                                                                  mean = Med,
                                                                                  sd = (Hi_beta - Med)/1.96))
      }else{ #if distribution is not triangle, nor normal, check LogNormal
        if(param[i, "Distribution"] == "Lognormal"){
          param_psa$psa_est[param_psa$Treatment == param$Treatment[i] & 
                              param_psa$Param      == param$Param[i]] <- exp(with(param[i, ], 
                                                                                  rnorm(n = n_iter ,
                                                                                        mean = log(Med),
                                                                                        sd = (log(Hi_beta) - log(Med)) / 1.96)))
        }else{ #if distribution is not triangle/normal/logNormal, check Beta
          if(param[i, "Distribution"] == "Beta"){
            param_psa$psa_est[param_psa$Treatment == param$Treatment[i] & 
                                param_psa$Param      == param$Param[i]] <- with(param[i, ],
                                                                                rbeta(n = n_iter,
                                                                                      shape1 = Lo_alpha,
                                                                                      shape2 = Hi_beta)) 
          }else{ #if distribution is not triangle/normal/lognormal/beta, check uniform
            if(param[i, "Distribution"] == "Uniform"){
              param_psa$psa_est[param_psa$Treatment == param$Treatment[i] & 
                                  param_psa$Param      == param$Param[i]] <- with(param[i, ],
                                                                                  runif(n = n_iter,
                                                                                        min = Lo_alpha,
                                                                                        max = Hi_beta)) 
            }else{ #if distribution is not triangle/normal/lognormal/beta/unif, check NA
              if(param[i, "Distribution"] == "NA"){
                param_psa$psa_est[param_psa$Treatment == param$Treatment[i] & 
                                    param_psa$Param      == param$Param[i]] <- with(param[i, ],
                                                                                    rep(x = Med,
                                                                                          n_iter)) 
              }
        }
      }
    }
      }
    }
  }
  
  return(param_psa)  # Return the parameter values for the PSA runs
}  




#' Generate probabilistic sensitivity analysis dataframe in DARTH style
#' 
#' \code{gen_psa} is used to compute the probabilistic sensitivity analysis (PSA) dataframe in DARTH style.
#' @param df_param_psa
#' @keywords Probabilistic sensitivity analysis, PSA, DARTH
#' @section Details:
#' \code{gen_psa} reorganises the values from the dataframe that is entered in the function to a data frame with values of each paramters of each PSA run.
#' @return  df_psa_input A data frame with the input values of all parameters of each PSA run. Dimension of the data frame are number of PSA simulations * parameters 
#' 

gen_psa <- function(df_param_psa){
  # Argument
  ## df_param_psa: A dataframe with the parameter values for the PSA
  # Return
  ## df_psa_input: The input values of each PSA iteration in a dataframe of size n_sim x n_parameters
  
  param_names <- unique(df_param_psa$Param)
  n_sim <- length(unique(df_param_psa$iter))
  
  # Create a matrix to store the results 
  df_psa_input <- matrix(data = NA, nrow = n_sim, ncol = length(param_names),
                         dimnames = list(c(1:n_sim), param_names))
  
  for(p in 1:length(param_names)){
    df_psa_input[, p] <- df_param_psa$psa_est[df_param_psa$Param == param_names[p]]
  }
  
  df_psa_input <- as.data.frame(df_psa_input) # Make a dataframe from the matrix 
  
  return(df_psa_input) # Return the psa input dataframe
}


### Make a function to do some calculations for the parameters based on the sampled list
est_param_m_Param_to_list <- function(l_list){
  with(as.list(l_list), {
    ## Costs
    
    # Cost of being alive after the hospitalization period
    c_Alive_post <- c_Healthcare * n_YpC 
    # Treatment Cost
    p_Public_insurance                <- 1 - p_Private_insurance 
    c_Trt                             <- ((c_Trt_private * p_Private_insurance) + 
                                            (c_Trt_public  * p_Public_insurance)) * n_Trt
    # Cost of hospital stay noTrt group
    c_H_noTrt <- c_Hospital * LOS_noTrt  
    c_I_noTrt <- (c_I_vent * p_vent + c_I_noVent * (1 - p_vent)) * LOS_noTrt 
    
    # Cost of hospital stay Trt group
    c_H_Trt <- c_Hospital * LOS_Trt   
    c_I_Trt <- (c_I_vent * p_vent + c_I_noVent * (1 - p_vent)) * LOS_Trt
    
    # Cost in hospital cycle
    c_Hospital_mix     <- (p_IC * c_I_noTrt) +  ((1 - p_IC) * c_H_noTrt)
    c_Hospital_mix_trt <- (p_IC * c_I_Trt)   +  ((1 - p_IC) * c_H_Trt)
    
    # Treatment effect
    t1               <- t_Die_mj  # time majority dies
    t2               <- n_days_timespan1 - t1 # time RCT - majority dies
    t3               <- n_DpC - t1 - t2       # remaining time
    
    # Test if adds up to total days in a cycle
    if (sum(t1 + t2 + t3) != n_DpC) {warning("Splitting the first cycle in terms did   not work")}
    
    ###### RETURN OUTPUT  
    l_list$c_Alive_post       <- c_Alive_post
    l_list$p_Public_insurance <- p_Public_insurance
    l_list$c_Trt <- c_Trt
    l_list$c_H_noTrt <- c_H_noTrt
    l_list$c_I_noTrt <- c_I_noTrt
    l_list$c_H_Trt <- c_H_Trt
    l_list$c_I_Trt <- c_I_Trt
    l_list$c_Hospital_mix <- c_Hospital_mix
    l_list$c_Hospital_mix_trt <- c_Hospital_mix_trt
    l_list$t1 <- t1
    l_list$t2 <- t2
    l_list$t3 <- t3
    out <- as.list(l_list)
    return(out) 
  })
}
## end
