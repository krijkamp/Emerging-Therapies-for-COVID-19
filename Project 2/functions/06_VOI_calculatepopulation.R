calculate_population <- function(study = study ){
  
    l_param_trt[[study]]$startcount         <- if(l_dates[study] < min(df_hosp_usa$date)) {(min(df_hosp_usa$date))} else {(l_dates[study])}
    l_param_trt[[study]]$endcount           <- max(df_hosp_usa$date)
    l_param_trt[[study]]$startcount_trial   <- l_param_trt[[study]]$startcount
    l_param_trt[[study]]$endcount_trial     <- as.Date(mondate(l_param_trt[[study]]$startcount_trial) + 3) #
    
    l_hosp       <- calculate_hospitalizations(df_hospitalisation = df_hosp_usa, 
                                               startcount = l_param_trt[[study]]$startcount, 
                                               endcount   = l_param_trt[[study]]$endcount)
    l_hosp_trial <- calculate_hospitalizations(df_hospitalisation = df_hosp_usa, 
                                               startcount = l_param_trt[[study]]$startcount_trial, 
                                               endcount   = l_param_trt[[study]]$endcount_trial) 
    
    l_out <- list(n_H_year = l_hosp$pop_mean,
                  n_H_trial = l_hosp_trial$pop_mean)
    
    return(l_out)
      
}




