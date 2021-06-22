### Survival ###

# Create a list, matrix and dataframe to store the results
generate_df_survival <- function(l_out_ce, l_names_drugs){
  #l_out_ce:  The list with all the output including the trace list 
  
  l_list <- l_out_ce$l_param_trt_f
  with(as.list(l_list), { 
  
  n_str <- length(v_names_str) # number of strategies
  n_t    <- l_list$n_t
  
  m_os <- matrix(data = NA, 
                 nrow = n_t + 1, ncol = length(names(l_names_drugs)),
                 dimnames = list(paste("cycle", 0:n_t, sep = " "), names(l_names_drugs)))
  
  # Create a data frame 
  df_survival <-  data.frame(treatment = rep(names(l_names_drugs), 
                                             each = (n_t + 1) * length(v_names_str)), 
                              strategy = rep(rep(c("trt", "notrt"), 
                                                each = (n_t + 1)), length(names(l_names_drugs))),
                                 cycle = rep(rep(0:n_t), 
                                         times = length(names(l_names_drugs)) * length(v_names_str)),
                                    OS = NA)  
  
  for (k in v_names_str){ # loop for the number or strategies 
    
    for (n in names(l_list)){ # loop for the number of treatments 
      
      drug     <- n
      strategy <- k
      l_trace <- l_out_ce[[drug]]$l_trace
      l_M     <- l_trace[[strategy]]

      # calculate the overall survival (OS) probability 
      v_os <- 1 - l_M[, "D"]
      
      m_os[, n] <- v_os  # store the vector in the drug specific column 
      
      # create a dataframe with the overall survival
      df_survival$OS[df_survival$treatment == drug & df_survival$strategy == strategy] <- v_os
      
    }
  }
  return(df_survival)
  })
}

plot_survival <- function(df_survival, cycle_zoom = FALS){
  # df_survival: dataframe with the overall survival for one drug
  # cycle_zoom: zoom in to the first x cycles, default is FALSE. If is is a number zoom in up to this cycle
  
    
    if(cycle_zoom == FALSE){
      # create a simple plot showing the OS of both treatments 
      plot <- ggplot(data = df_surv_loop, aes(x = cycle, y = OS, group = strategy)) + 
        geom_line(aes(linetype = strategy, color = strategy)) +
        ylab ("Overall survival") + 
        ggtitle(paste("Overall survival", df_surv_loop$treatment[n], seb = ""))
    }
    
    if (cycle_zoom > 0){
      # create a simple plot showing the OS of both treatments for the first 5 cycle
      plot <- ggplot(data = df_surv_loop, aes(x = cycle, y = OS, group = strategy)) +  
        geom_line(aes(linetype = strategy, color = strategy)) + xlim(0, cycle_zoom) +
        ylab ("Overall survival") + 
        ggtitle(paste("Overall survival - first 5 cycles", df_surv_loop$treatment[n], seb = ""))
    }
  
  return(plot)
}


### Life expectancy ###

generate_df_life_expectancy <- function(df_survival, v_names_str){
  # Create a matrix to store the life expectancy of each treatment
  # NOTE: this survival estimate does not use a weighted LE for the first cycle but assumes that all patients are alive the entire cycle. This explains the slight difference between the LE estimate using this survival code and the LY in the CEA.
  m_LE <- matrix(data = NA, ncol = 2, nrow = length(unique(df_survival$treatment)))
  rownames(m_LE) <- c(unique(df_survival$treatment))
  colnames(m_LE) <- c(v_names_str)
  
  df_le <- data.frame(treatment = c(sort(rep(unique(df_survival$treatment), 2))),
                      strategy = c(rep(v_names_str, length(unique(df_survival$treatment)))),
                      LE = NA)
  
  
  for (n in unique(df_survival$treatment)){
    df_surv_loop <- filter(df_survival, treatment == n)
    
    df_notrt <- df_surv_loop %>% 
      filter(strategy == v_names_str[1]) 
    
    df_trt <- df_surv_loop %>% 
      filter(strategy == v_names_str[2]) 
    
    # summing probability of OS over time  (i.e. average life expectancy in in days in a 1 year horizon)
    v_le_notrt <- sum(df_notrt$OS) * n_YpC   # summing probability of OS over time  (i.e. average life expectancy in days in a 1 year horizon)
    v_le_trt <- sum(df_trt$OS) * n_YpC   # summing probability of OS over time  (i.e. average life expectancy in days in a 1 year horizon)
    
    m_LE[n, v_names_str[2]] <- v_le_trt
    m_LE[n, v_names_str[1]] <- v_le_notrt
    
    df_le$LE[df_le$treatment == n & df_le$strategy == v_names_str[1]] <- v_le_notrt
    df_le$LE[df_le$treatment == n & df_le$strategy == v_names_str[2]] <- v_le_trt
    
  }
  m_LE <- as.data.frame(m_LE)
  
  output <- list(m_LE  = m_LE,
                 df_le = df_le)
  return(output)
}


plot_life_expectancy <- function(df_le_plot){
  # A function to plot the LE
  # df_le_plot: a dataframe with columns treatment, strategy and LE for one drug
  # 
  
  n <- df_le_plot$treatment[1]
    
  # create a simple plot showing the OS of both treatments 
  plot <- ggplot(data = df_le_plot, aes(x = strategy, y = LE, group = strategy, fill = strategy)) + 
    geom_col() +
    ggtitle(paste("Life expectancy", n, sep = " ")) 
    
  
  return(plot) # return the plot
}

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



