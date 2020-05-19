#####################################################################
## Make transition array (without transition from pre- to postop) ##
#####################################################################

make_m_trans <- function(number_states = NULL,
                         state_names = NULL,
                         number_cycles = NULL, 
                         parameters = NULL,
                         cbsdata = NULL){
# # Arguments:  
  ## number_states: number of states
  ## state_names: names of the states
  ## n_cycles: number of cycles 
  ## parameters: parameters in the model
  ## cbsdata: data of the cbs containing survival probabilities
# Returns: 
  # m_trans: gives the transition probability array
  
  # Construct transition probability array
  m_trans     <- array(data = 0, dim = c(number_states, number_states, n_cycles), 
                       dimnames = list(state_names, state_names))
  
  # standard probability to die (from CBS data) 
  # based on the average age of the patient. The parameters is used of all cycles
  # this probability is added to the disease specific mortality
  p_die <- rep(cbsdata[cbsdata$age >= parameters["Age"], ]$prob, each = 52)
  
  # Fill the array, first health state is the health state the individual start
  # the second health state is the health state the individual transitions to
  # during the cycle 
  
  # From Preop
  m_trans["Preop", "Dead", ]              <- parameters["Surv_no_tx"] + p_die
  m_trans["Preop", "Preop", ]             <- 1 - m_trans["Preop", "Dead", ] 
 
  # From Dead 
  m_trans["Dead", "Dead",  ]              <- 1
  m_trans["Dead", c("Preop", "Postop"), ] <- 0

  # From Postop
  m_trans["Postop", "Dead", ]             <- parameters["Surv_tx"] + p_die
  m_trans["Postop", "Preop", ]            <- 0
  m_trans["Postop", "Postop", ]           <- 1 - m_trans["Postop", "Dead", ]
  
  # For oncological problems: set survival with treatment after tumor_dbl_time equal to survival without treatment. To find the first cycle were the survival is the same we round the tumor doubling time to the lowest integer value using ceiling. From that cycle to the final cycle the transition probabilities after surgery are the same as the probabilities before surgery
  
  if("Time_noeff_Surv" %in% names(parameters)){

    m_trans["Postop", "Dead", c(ceiling(parameters["Time_noeff_Surv"]) : number_cycles)] <- 
     m_trans["Preop", "Dead", c(ceiling(parameters["Time_noeff_Surv"]) : number_cycles)]
    
    # calculate the remaining transition probability to stay in the Postop state
    m_trans["Postop", "Postop", c(ceiling(parameters["Time_noeff_Surv"]) : number_cycles)] <- 1 -        m_trans["Postop", "Dead", c(ceiling(parameters["Time_noeff_Surv"]) : number_cycles)]
  }
  
  ##BUILD IN ERROR MESSAGE to check if transition matrix is valid (i_e_, each row should add up to 1)
  valid <- apply(m_trans, 3, function(x) sum(rowSums(x)) == number_states)
  if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(number_cycles)))) {
    invalid_cycles <- which(valid == FALSE)
    if(length (invalid_cycles) > 10){
      stop(paste("This is not a valid transition Matrix. Issued in more that", as.character(length(invalid_cycles)), "cycles", sep = " ")) 
    }
    else {   
      stop(print("This is not a valid transition Matrix. Issues in: "), 
           print(paste("Cycles", as.character(invalid_cycles), sep = ""))) }
}
    
  ##BUILD IN ERROR MESSAGE 
  # This message checks for unlikely high probabilities to die 
  if(TRUE %in% (c(m_trans["Preop", "Dead", ], m_trans["Postop","Dead", ]) > 0.9)){
    print(paste("There transition probabilities to the Dead state higher than 0.9"))
    break
  }
  
  return(m_trans) # return the transition probability array
  
}

######################################################
## Include delay to operation in array              ##
######################################################
trans_operation <- function(trans_matrix = NULL, weeks_until_op = NULL, number_cycles = NULL){
# Arguments
  # trans_matrix:   the transition array created by make_m_trans
  # weeks_until_op: the number of weeks until the surgery (operation)
  # number_cycles:  the total number of cycles 
# Returns:
  # trans_operation: returns the transition probability array 
                   # taking operation delay into account
  
  
  # Before delay
  # For the duration of the delay, no transition from Preop to Postop is possible 
  
  number_states <- dim(trans_matrix)[1]  # count the number of health states
  
  # If delay = 999, make the weeks_until_op as large as n.cycles
  if(weeks_until_op==999){
    weeks_until_op <- number_cycles-1
  }
  
  trans_matrix["Preop",
               "Postop", 
               1:weeks_until_op] <- 0 
 
   # After the duration of delay, those that are still alive in Preop and 
   # don't die during that cycle will get an operation. 
  trans_matrix["Preop",
               "Preop", 
               1:weeks_until_op] <- 1 - trans_matrix["Preop", 
                                                     "Dead", 
                                                     1:weeks_until_op]
  # After delay
  trans_matrix["Preop", 
               "Postop", 
               (weeks_until_op + 1):n_cycles] <- 1 - trans_matrix["Preop",
                                                                  "Dead",
                                                                  (weeks_until_op + 1):number_cycles] 
  # after these weeks no-one is in the Preop state anymore
  # and the Preop - Preop probability is zeor
  trans_matrix["Preop", 
               "Preop", 
               (weeks_until_op + 1):number_cycles] <- 0
  
  # Check if transition matrix is valid (i_e_, each row should add up to 1)
  valid <- apply(trans_matrix, 3, function(x) sum(rowSums(x)) == number_states)
  if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(number_cycles)))) {
    invalid_cycles <- which(valid == FALSE)
    if(length (invalid_cycles) > 10){
      stop(paste("This is not a valid transition Matrix. Issued in more that", as.character(length(invalid_cycles)), "cycles", sep = " ")) 
    }
    else {   
      stop(print("This is not a valid transition Matrix. Issues in: "), 
           print(paste("Cycles", as.character(invalid_cycles), sep = ""))) }
  }
  return(trans_matrix)  # Return the transition array
}


######################################################
## Calculate the derivative / slope for QALY loss   ##
######################################################
# a function to calculate the derivative 
calculateDerivative <- function(df_pooled, plot = FALSE, plot_ind = FALSE, cum = FALSE, folder = "figures"){
# Arguments:
  ## df_pooled:     pooled results data
  ## plot:          do we have to make a plot (default FALSE) 
  ## plot_ind:      do we like to have individual plots (default is FALSE)
  ## cum:           do we like to have the cumulative derivative (default FALSE)
  ## folder:        folder to store the values
# Return:
  ## df_derivative: a data frame with the values 
  #   stores figures in the specified folder
  
  pop_names <- as.character(sort(unique(df_pooled$Label)))
  delay <- unique(df_pooled$delay)
  
  difference <- matrix(NA, nrow = length(pop_names), ncol = length(delay) - 1,
                       dimnames = list(pop_names, c(paste("period", 1:(length(delay)-1), sep = " ")) ))
  
  
  df_derivative <-  data.frame(Label = sort(rep(pop_names, length(delay) - 1)),
                               weeks   = rep(delay[-length(delay)], length(pop_names)))
  
  for (d in 1:length(pop_names)){
    data <- df_pooled[df_pooled$Label == pop_names[d]]
    difference[d, ]  <- diff(data$QALY_med)/diff(data$delay) # store it in a matrix 
    
    df_derivative$derivative[df_derivative$Label == pop_names[d]] <- diff(data$QALY_med)/diff(data$delay) # Store in a data frame 
    
    df_derivative$derivative_lo[df_derivative$Label == pop_names[d]] <- (diff(data$QALY_med) - diff(data$QALY_lo)) / diff(data$delay) # Store in a data frame 
    
    df_derivative$derivative_hi[df_derivative$Label == pop_names[d]] <- (diff(data$QALY_med) + diff(data$QALY_hi)) / diff(data$delay) # Store in a data frame   
    
    ## IF we like to add the code to make the figures for each Population we can add this here 
    ## PLACEHOLDER ### 
    
  }
  

  if (plot == TRUE){
    
    plot_derivative <- ggplot(df_derivative, aes(x = weeks, y = derivative, ymin = derivative_lo, ymax =  derivative_hi)) + 
      geom_ribbon(alpha = 0.4) +
      geom_line(cex = 2) +
      facet_wrap(~ Label) +
      xlab("Delay (weeks)") +
      ylab("QALY loss/week") +
      theme_bw() +
      theme(text = element_text(size = 16), 
            strip.text = element_text(size = 16))+
      scale_y_continuous(limits=c(-0.03,0.001))
  }
  return(list(plot_derivative, df_derivative)) 
}

################################################
## Function to plot all results of the model  ##
################################################

plotPopulationOutcomes <- function(data, folder = "figures/", size_cm=10){
  # Arguments:
  ## data: pooled results
  ## folder: folder where figures are stored
  
  pop_names <- sort(unique(data$Label))

  # Create a figure of the QALYs with confidence interval of the PSA
  for (d in 1:length(pop_names)) {
    avg_y <- median(as.numeric(data$QALY_med[data$Label == pop_names[d]]))
    plot_QALY <- ggplot(data[data$Label == pop_names[d],], aes(x = delay, y = QALY_med, ymin = QALY_lo, ymax = QALY_hi)) + 
      geom_ribbon(alpha = 0.4) +
      geom_line(cex = 2) +
      facet_wrap(~ Label) +
      xlab("Delay (weeks)") +
      ylab("QALY") +
      scale_y_continuous(limits = c(avg_y - 10,
                                    avg_y + 10))+
      theme_bw() +
      theme(text = element_text(size = 16), 
            strip.text = element_text(size = 16))
    p <- gsub(' ', "_", pop_names[d]) # Replace the spaces with an underscore to save the file
    ggsave(filename = paste(folder,p,"_QALY.png", sep = ""),
           plot = plot_QALY, device = "png", 
           height=size_cm, width=size_cm, units = "cm")
  }
}



