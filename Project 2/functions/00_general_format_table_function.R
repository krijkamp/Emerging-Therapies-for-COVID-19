#----------------------------------------------------------------------------#
#' Format the summary table
#'
#' \code{format_table_summary} formats the summary table for the base-case and PSA results 
#'
#' @param table_summary a dataframe object - table with summary results
#' @param v_names_subet a vector - with the full names of the treatments to select
#' @param round_factor  a logical value to round the values to Millions and Thousands. Default = FALSE
#' @return              a dataframe object - formatted summary table
#' 
#' 


format_table_summary  <- function(table_summary, v_names_subset = NULL, EVPPI = FALSE, round_factor = FALSE) {
  
  if(EVPPI == FALSE){
    colnames(table_summary)[colnames(table_summary) 
                            %in% c("Cost-effective", 
                                   "Incr cost Rx",
                                   "Incr effect Rx", 
                                   "ICER",
                                   "Incr NMB",
                                   "Incr NHB")] <- 
      c("Is treatment cost-effective?", 
        "Incremental Costs ($)", 
        "Incremental QALYs",
        "ICER ($/QALY)",
        "Incremental net monetary benefit ($)",
        "Incremental net health benefit (QALY)"
      )
  } else {
  colnames(table_summary)[colnames(table_summary) 
                      %in% c("Cost-effective", 
                             "Incr cost Rx",
                             "Incr effect Rx", 
                             "ICER",
                             "Incr NMB",
                             "Incr NHB",
                             "EVPPI",
                             "Future patients",
                             "Current patients")] <- 
    c("Is treatment cost-effective?", 
      "Incremental Costs ($)", 
      "Incremental QALYs",
      "ICER ($/QALY)",
      "Incremental net monetary benefit ($)",
      "Incremental net health benefit (QALY)",
      "EVPPI",
      "Future patients",
      "Current patients"
      )
  }
  # make all columns that are numbers numeric 
  table_summary[, -1] <- data.frame(lapply(table_summary[, -1], as.numeric))
  
  
  
  table_summary$`Incremental QALYs`       <- round(table_summary$`Incremental QALYs`, 3)
  table_summary$`Incremental Costs ($)`   <- round(table_summary$`Incremental Costs ($)`, 3)
  table_summary$`ICER ($/QALY)`           <- round(table_summary$`ICER ($/QALY)`, 0)
  
  table_summary$`Incremental net monetary benefit ($)`  <- round(table_summary$`Incremental net monetary benefit ($)`, 0)
  table_summary$`Incremental net health benefit (QALY)` <- round(table_summary$`Incremental net health benefit (QALY)`, 3)
  
  #table_summary[is.na(table_summary)] <- "-"
  
  if(EVPPI == TRUE) {
  table_summary$`Future patients`        <- round(table_summary$`Future patients`)
  table_summary$`Current patients`       <- round(table_summary$`Current patients`)
  table_summary$`EVPPI`                  <- round(table_summary$`EVPPI`)
  }
  
  if(length(v_names_subset)>0){
    # select the values
    table_summary <- table_summary[v_names_subset, ]
    rownames(table_summary) <- remove_date_from_name(rownames(table_summary))
  }
  
  # Order based on ICER
  table_summary <- table_summary[order(table_summary$`ICER ($/QALY)`), ]
  
  v_new_order <-    c(which(table_summary[, "Is treatment cost-effective?"] == "Yes*"),
                      which(table_summary[, "Is treatment cost-effective?"] == "Yes"),
                      which(table_summary[, "Is treatment cost-effective?"] == "Yes**"),
                      which(table_summary[, "Is treatment cost-effective?"] == "No"),
                      which(table_summary[, "Is treatment cost-effective?"] == "No*"),
                      which(table_summary[, "Is treatment cost-effective?"] == "No**"))
  
  table_summary <- table_summary[v_new_order, ]
  
  return(table_summary)
}



#----------------------------------------------------------------------------#
#' Format the summary table
#'
#' \code{format_table_summary} formats the summary table for the base-case and PSA results 
#'
#' @param table_summary a dataframe object - table with summary results
#' @param v_names_subet a vector - with the full names of the treatments to select
#' @param round_factor  a logical value to round the values to Millions and Thousands. Default = FALSE
#' @return              a dataframe object - formatted summary table
#' 
#' 


format_table_LY_summary  <- function(table_summary, v_names_subset = NULL, EVPPI = FALSE, round_factor = FALSE) {
  
  if(EVPPI == FALSE){
    colnames(table_summary)[colnames(table_summary) 
                            %in% c("Cost-effective", 
                                   "Incr cost Rx",
                                   "Incr effect Rx", 
                                   "ICER",
                                   "Incr NMB",
                                   "Incr NHB")] <- 
      c("Is treatment cost-effective?", 
        "Incremental Costs ($)", 
        "Incremental LYs",
        "ICER ($/LY)",
        "Incremental net monetary benefit ($)",
        "Incremental net health benefit (LY)"
      )
  } else {
    colnames(table_summary)[colnames(table_summary) 
                            %in% c("Cost-effective", 
                                   "Incr cost Rx",
                                   "Incr effect Rx", 
                                   "ICER",
                                   "Incr NMB",
                                   "Incr NHB",
                                   "EVPPI",
                                   "Future patients",
                                   "Current patients")] <- 
      c("Is treatment cost-effective?", 
        "Incremental Costs ($)", 
        "Incremental LYs",
        "ICER ($/LY)",
        "Incremental net monetary benefit ($)",
        "Incremental net health benefit (LY)",
        "EVPPI",
        "Future patients",
        "Current patients"
      )
  }
  # make all columns that are numbers numeric 
  table_summary[, -1] <- data.frame(lapply(table_summary[, -1], as.numeric))
  
  
  
  table_summary$`Incremental LYs`       <- round(table_summary$`Incremental LYs`, 3)
  table_summary$`Incremental Costs ($)`   <- round(table_summary$`Incremental Costs ($)`, 3)
  table_summary$`ICER ($/LY)`           <- round(table_summary$`ICER ($/LY)`, 0)
  
  table_summary$`Incremental net monetary benefit ($)`  <- round(table_summary$`Incremental net monetary benefit ($)`, 0)
  table_summary$`Incremental net health benefit (LY)`   <- round(table_summary$`Incremental net health benefit (LY)`, 3)
  
  #table_summary[is.na(table_summary)] <- "-"
  
  if(EVPPI == TRUE) {
    table_summary$`Future patients`        <- round(table_summary$`Future patients`)
    table_summary$`Current patients`       <- round(table_summary$`Current patients`)
    table_summary$`EVPPI`                  <- round(table_summary$`EVPPI`)
  }
  
  if(length(v_names_subset)>0){
    # select the values
    table_summary <- table_summary[v_names_subset, ]
    rownames(table_summary) <- remove_date_from_name(rownames(table_summary))
  }
  
  # Order based on ICER
  table_summary <- table_summary[order(table_summary$`ICER ($/LY)`), ]
  
  v_new_order <-    c(which(table_summary[, "Is treatment cost-effective?"] == "Yes*"),
                      which(table_summary[, "Is treatment cost-effective?"] == "Yes"),
                      which(table_summary[, "Is treatment cost-effective?"] == "Yes**"),
                      which(table_summary[, "Is treatment cost-effective?"] == "No"),
                      which(table_summary[, "Is treatment cost-effective?"] == "No*"),
                      which(table_summary[, "Is treatment cost-effective?"] == "No**"))
  
  table_summary <- table_summary[v_new_order, ]
  
  return(table_summary)
}



# Divide by Million $
#for (k in c("EVPPI")){
#  df_summary_cea_QALY_PSA_plot[k, !is.na(df_summary_cea_QALY_PSA_plot[k, ])] <- round(as.numeric#(df_summary_cea_QALY_PSA_plot[k, !is.na(df_summary_cea_QALY_PSA_plot[k, ])]) / 1e6)
#  n_new_name <- paste(k, "(Million $)", sep =" ")
#  rownames(df_summary_cea_QALY_PSA_plot)[rownames(df_summary_cea_QALY_PSA_plot) == k] <- n_new_name
#}


#----------------------------------------------------------------------------#
#' Format the summary table
#'
#' \code{format_table_summary} formats the CEA table.
#'
#' @param table_summary a dataframe object - table with summary results
#' @param v_names_subet a vector - with the full names of the treatments to select
#' @param round_factor  a logical value to round the values to Millions and Thousands. Default = FALSE
#' @return              a dataframe object - formatted summary table
#' 

format_table_final <- function(table2format, v_round = NULL, v_thousand = NULL, v_million = NULL, v_billion = NULL, integer = FALSE){
  
  if(!is.null(v_thousand)){
    for(r in v_thousand){
      table2format[r, ] <- round(as.numeric(table2format[r, ])/1e3, if_else(r %in% v_round, 0, if_else(integer == TRUE, 0, 3)))
      rownames(table2format)[rownames(table2format) 
                              %in% r] <- paste(r, "(Thousand)")
    }
  }

  if(!is.null(v_million)){
    for(r in v_million){
      table2format[r, ] <- round(as.numeric(table2format[r, ])/1e6, if_else(r %in% v_round, 0, if_else(integer == TRUE, 0, 3)))
      rownames(table2format)[rownames(table2format) 
                            %in% r] <- paste(r, "(Million)")
    }
  }
  
  if(!is.null(v_billion)){
    for(r in v_billion){
      table2format[r, ] <- round(as.numeric( table2format[r, ])/1e9, if_else(r %in% v_round, 0, if_else(integer == TRUE, 0, 3)))
      rownames(table2format)[rownames(table2format) 
                            %in% r] <- paste(r, "(Billion)")
    }
  }
  
  if(!is.null(v_round)){
    for(r in v_round[v_round %nin% c(v_thousand, v_million, v_billion)]){
      table2format[r, ] <- round(as.numeric(table2format[r, ]), 0)
    }
  }
  

 
  return(table2format) 
}



number_ticks <- function(n) {
  function(limits) {
    pretty(limits, n + 1)
  }
}

  


