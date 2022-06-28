# Functions for the hospitalizaion data

# Dataset for hospitalization: Reference_hospitalization_all_locs.csv [current projection]
# Relevant parameters:
# location_name: name country or subnational location
# admis_mean: Mean hospital admissions
# admis_lower: Lower uncertainty bound of hospital admissions per day
# admis_upper: Upper uncertainty bound of hospital admissions per day
# new_ICU_mean: Mean number of new people going to the ICU by day
# new_ICU_lower: Lower uncertainty bound of nr of people going to ICU per day
# new_ICU_upper: Upper uncertainty bound of nr of people going to ICU per day

modify_df_hospitalization <- function(df_hospitalisation, country = "United States of America"){
  # Overwrite dates of admission "text" as "date" datatype
  df_hospitalisation$date <- as.Date(df_hospitalisation$date)
  # Create dataframe subset for country of interest. Currently only USA.
  df_hosp_country <- param_hosp[df_hospitalisation$location_name == country, ]
  
  return(df_hosp_country) # the country specific hospitalizaton data
}

calculate_hospitalizations <- function(df_hospitalisation, startcount = "2020-12-11", endcount = "2020-04-01" ){
  
  row_startcount_hosp <- which(grepl(startcount, df_hospitalisation$date))
  row_endcount_hosp   <- which(grepl(endcount, df_hospitalisation$date))
  
  # Select timeperiod over which number of hospitalizations should be calculated
  hosp_timeperiod_usa <- df_hospitalisation[c(row_startcount_hosp:row_endcount_hosp), ]
  
  # Calculate expected number of people to still be hospitalized from today   onward
  pop_mean  <- sum(hosp_timeperiod_usa$admis_mean)
  pop_lower <- sum(hosp_timeperiod_usa$admis_lower) 
  pop_upper <- sum(hosp_timeperiod_usa$admis_upper) 
  
  l_out <- list(pop_mean = pop_mean,
                pop_lower = pop_lower,
                pop_upper = pop_upper)
  
  return(l_out)
  
}