# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)
# Petros Pechlivanoglou, PhD (6) 

# In collaboration of: 		
# 1 Center for Research and Teaching in Economics (CIDE), Drug Policy Program, Mexico
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada


# used for sensitivity analysis

calculate_ce_out <- function (l_params_all, n_wtp = 100000) {
  with(as.list(l_params_all), {
    
    # Static characteristics
    v_x      <- runif(n_i, min = 0.95, max = 1.05) # treatment effect modifier at baseline                                         
    v_age0   <- sample(x = dist_Age$age, prob = dist_Age$prop, size = n_i, replace = TRUE) # sample from age distribution an initial age for every individual
    df_X     <- data.frame(ID = 1:n_i, x = v_x, Age = v_age0)
    
    
    #### 05.1 Probability function ####
    # The Probs function that updates the transition probabilities of every cycle is shown below.
    
    Probs <- function(M_t, df_X, v_Ts, t) { 
      # Arguments:
      # M_t: health state occupied  at cycle t (character vector)
      # df_X: dataframe with individual caracteristics
      # v_Ts: time an individual is sick
      # t:     current cycle 
      # Returns: 
      #   transition probabilities for that cycle
      
      m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  # create matrix of state transition probabilities
      rownames(m_p_t) <-  v_n                               # give the state names to the rows
      
      # lookup baseline probability and rate of dying based on individual characteristics
      p_HD_all <- inner_join(df_X, p_mort, by = c("Age"))
      p_HD     <- p_HD_all[M_t == "H","p_HD"]
      
      # update the v_p with the appropriate probabilities   
      m_p_t[, M_t == "H"]  <- rbind(1 - p_HS1 - p_HD, p_HS1, 0, p_HD)                             # transition probabilities when healthy
      m_p_t[, M_t == "S1"] <- rbind(p_S1H, 1 - p_S1H - p_S1S2 - p_S1D[v_Ts], p_S1S2, p_S1D[v_Ts]) # transition probabilities when sick
      m_p_t[, M_t == "S2"] <- rbind(0, 0, 1 - p_S2D, p_S2D)                                       # transition probabilities when sicker
      m_p_t[, M_t == "D"]  <- rbind(0, 0, 0, 1)                                                   # transition probabilities when dead   
      return(t(m_p_t))
    }       
    
    #### 05.2 Cost function ####
    # The Costs function estimates the costs at every cycle.
    
    Costs <- function (M_t, Trt = FALSE) {
      # M_t: health state occupied by individual i at cycle t (character variable)
      # Trt:  is the individual being treated? (default is FALSE) 
      
      c_t <- 0                                 # by default the cost for everyone is zero 
      c_t[M_t == "H"]  <- c_H                  # update the cost if healthy
      c_t[M_t == "S1"] <- c_S1 + c_Trt * Trt   # update the cost if sick conditional on treatment
      c_t[M_t == "S2"] <- c_S2 + c_Trt * Trt   # update the cost if sicker conditional on treatment
      c_t[M_t == "D"]  <- c_D                  # update the cost if dead
      
      return(c_t)        		                   # return the costs
    }
    
    #### 05.3 Health outcome function ####
    # The Effs function to update the utilities at every cycle.
    
    Effs <- function (M_t, df_X, Trt = FALSE, cl = 1) {
      # M_t: health state occupied by individual i at cycle t (character variable)
      # df_Pop: inidividual characteristics inclusing Age, Sex and the effect mofifier of the treatment effect
      # Trt:  is the individual treated? (default is FALSE) 
      # cl:   cycle length (default is 1)
      
      u_t <- 0                                       # by default the utility for everyone is zero
      u_t[M_t == "H"]  <- u_H                        # update the utility if healthy
      u_t[M_t == "S1" & Trt == FALSE] <- u_S1        # update the utility if sick
      u_t[M_t == "S1" & Trt == TRUE]  <- u_Trt * df_X$x[M_t == "S1"]  # update the utility if sick but on treatment (adjust for individual effect modifier) 
      u_t[M_t == "S2"] <- u_S2                       # update the utility if sicker
      u_t[M_t == "D"]  <- u_D                        # update the utility if dead
      
      QALYs <-  u_t * cl            # calculate the QALYs during cycle t
      return(QALYs)                 # return the QALYs
    }
    
    #### 06 Run Microsimulation ####
    MicroSim <- function(n_i, df_X , Trt = FALSE, seed = 1) {
      # Arguments:  
      # n_i:     number of individuals
      # df_X     data frame with individual data 
      ## Age      age of the individuals
      ## Sex      sex of the indivuduals 
      ## x        effect modifier  
      # Trt:     is this the individual receiving treatment? (default is FALSE)
      # seed:    defauls is 1
      
      set.seed(seed) # set the seed
      
      # create three matrices called m_M, m_C and m_E
      # number of rows is equal to the n_i, the number of columns is equal to n_t  (the initial state and all the n_t cycles)
      # m_M is used to store the health state information over time for every individual
      # m_C is used to store the costs information over time for evey individual
      # m_E is used to store the effects information over time for every individual
      
      m_M <- m_C <- m_E <- m_Ts <-  matrix(nrow = n_i, ncol = n_t + 1, 
                                           dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                           paste("cycle", 0:n_t, sep = " ")))  
      
      m_M [, 1] <- v_M_init    # initial health state at cycle 0 for individual i
      v_Ts      <- v_Ts_init   # initialize time since illnes onset for individual i
      
      m_C[, 1]  <- Costs(m_M[, 1], Trt)         # calculate costs per individual during cycle 0
      m_E[, 1]  <- Effs (m_M[, 1], df_X, Trt)   # calculate QALYs per individual during cycle 0
      
      # open a loop for time running cycles 1 to n_t 
      for (t in 1:n_t) {
        v_p <- Probs(m_M[, t], df_X, v_Ts, t)             # calculate the transition probabilities for the cycle based on  health state t
        m_M[, t + 1]  <- samplev(v_p, 1)                  # sample the current health state and store that state in matrix m_M 
        m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt)         # calculate costs per individual during cycle t + 1
        m_E[, t + 1]  <- Effs(m_M[, t + 1], df_X, Trt)    # calculate QALYs per individual during cycle t + 1
        
        v_Ts <- if_else(m_M[, t + 1] == "S1", v_Ts + 1, 0) # update time since illness onset for t + 1 
        df_X$Age[m_M[,t + 1] != "D"]  <- df_X$Age[m_M[, t + 1] != "D"] + 1
        
        
      } # close the loop for the time points 
      
      # calculate  
      tc <- m_C %*% v_dwc    # total (discounted) cost per individual
      te <- m_E %*% v_dwe    # total (discounted) QALYs per individual 
      tc_hat <- mean(tc)     # average (discounted) cost 
      te_hat <- mean(te)     # average (discounted) QALYs
      
      # store the results from the simulation in a list
      results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, tc_hat = tc_hat, te_hat = te_hat)   
      return(results)  # return the results
    } # end of the MicroSim function  
    
    ### Run the simulation for both no treatment and treatment options
    outcomes_no_trt  <- MicroSim(n_i, df_X, Trt = FALSE, seed = 1)
    outcomes_trt     <- MicroSim(n_i, df_X, Trt = TRUE, seed = 1)
    
    ## Vector with total discounted mean Costs and QALYs
    v_tc_d      <- c(outcomes_no_trt$tc_hat, outcomes_trt$tc_hat)
    v_tu_d      <- c(outcomes_no_trt$te_hat, outcomes_trt$te_hat)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d     <- v_tu_d * n_wtp - v_tc_d
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_tu_d,
                        NMB      = v_nmb_d)
    return(df_ce)
  }
  )
}


samplev <- function(m_Probs, m) {
  # Arguments
  # m_Probs: matrix with probabilities (n.i * n.s)
  # m:       number of states than need to be sampled per individual  
  # Return
  # ran:    n.i x m matrix filled with sampled health state(s) per individual
  
  d <- dim(m_Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
  n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
  k <- d[2]          # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m_Probs)[[2]]  # extract the names of the health states considered for sampling
  if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k       # create a sequence from 1:k (number of health states considered)
  # create a matrix 
  ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
  U <- t(m_Probs)    # transposed m_Probs matrix n.i x n.s --> n.s x n.i 
  
  for(i in 2:k) {    # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
  }
  if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
    un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
    ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  }
  ran # return the new health state per individual n.i x m
} # close the function 


# plot health state trace
plot_m_TR <- function(m_M, title = "Health state trace" ) {
  # plot the distribution of the Pop across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  matplot(m_TR, type = "l", main = title, col= 1:n_states,
          ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  legend("topright", v_n, col = 1:n_states,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)
  
}

#-----------------------------------------------------------------------------------------------#
#### R function to extract the parameters of a beta distribution from mean and st. deviation ####
#-----------------------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
betaPar <- function(m, s) {
  a <- m * ((m * (1 - m) / s ^ 2) - 1)
  b <- (1 - m) * ((m * (1 - m) / s ^ 2) - 1)
  list(a = a, b = b)
}

beta_params <- function(mean, sigma) {
  alpha <- ((1 - mean) / sigma ^ 2 - 1 / mean) * mean ^ 2
  beta  <- alpha * (1 / mean - 1)
  params <- list(alpha = alpha, beta = beta)
  return(params)
}

#-------------------------------------------------------------------------------------------------#
#### R function to extract the parameters of a gamma distribution from mean and st. deviation  ####
#-------------------------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
gammaPar <- function(m, s) {   
  # m: mean  
  # s: standard deviation 
  shape <- m ^ 2 / s ^ 2
  scale <- s ^ 2 / m
  list(shape = shape, scale = scale)
}

#-----------------------------------------------------------------------------------------------#
#### R function to create Markov Trace ####
#-----------------------------------------------------------------------------------------------#
CalculateMarkovTrace <- function(m_Trans, v_MT1, n_s, v_n, n_t){
  # Arguments
  # m_Trans Transition probability matrix (control or intervention)
  # v_MT1   Starting cohort allocation
  # n_s    No. of health states
  # v_n    Vector with health state names
  # n_t:    No. of cycles
  # Return
  # m_Trace The cohort trace matrix 
  
  m_Trace <- matrix(NA, nrow = n_t + 1, ncol = n_s,
                    dimnames = list(paste("Cycle", 0:n_t, sep = " "), v_n))
  # Creates a trace matrix for the allocation of the cohort in each cycle for each health state
  m_Trace[1, ] <- v_MT1
  # First cycle is the starting allocation of the cohort
  for(i in 2:(n_t + 1)) {
    m_Trace[i, ] <- t(m_Trace[i - 1, ]) %*% m_Trans
  }
  # Fills the rows of the trace matrix by multiplying the first row 
  # of the trace matrix with the transition probability matrix
  m_Trace
  # Function returns the trace matrix
}
#-----------------------------------------------------------------------------------------------#
#### R function to plot Markov Trace ####
#-----------------------------------------------------------------------------------------------#
PlotTrace <- function(trace, xlab, title, txtsize = 12) {
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
  
  trace <- data.frame(time = seq(1, (nrow(trace))), trace)
  trace <- melt(trace, id.vars = "time")
  plot_trace <- ggplot(trace, aes(x = time, y = value, colour = variable)) +
    geom_line() +
    scale_colour_hue("States", l = 50) +
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

#-----------------------------------------------------------------------------------------------#
#### R function to plot Markov Trace with set limits x/y axis  ####
#-----------------------------------------------------------------------------------------------#
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
  
  trace <- data.frame(time = seq(1, (nrow(trace))), trace)
  trace <- melt(trace, id.vars = "time")
  plot_trace <- ggplot(trace, aes(x = time, y = value, group = variable)) +
    geom_line(aes(linetype = variable, color = variable)) +
    scale_colour_hue("States", l = 50) +
    scale_linetype_discrete("States")+
    scale_y_continuous(limits=c(0,1)) +
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


#-----------------------------------------------------------------------------------------------#
#### R function to create distribution  ####
#-----------------------------------------------------------------------------------------------#

# Creating beta distribution function
beta_mom <- function(mean, var){
  # Arguments
  # mean: mean of the probability
  # var: variance of the probability
  # Returns
  # beta distribution for the probability
  term <- mean * (1 - mean) / var - 1
  alpha <- mean * term
  beta <- (1 - mean) * term
  if (var >= mean * (1 - mean)) stop("var must be less than mean * (1 - mean)")
  return(list(alpha = alpha, beta = beta))
}

beta_mom(0.8, 0.1)

# Creating a lognormal distributions
lnorm_mom <- function(mean, sd){
  # Arguments
  # mean : mean of the relative risk (in our case)
  # sd : standard deviation of the relative risk
  # Returns
  # lognormal distribution for the probability
  if (mean > 0){
    sigma2 <- log((sd ^ 2 + mean ^ 2) /mean ^ 2)
    mu <- log(mean) - 1/2 * sigma2
  } else{
    stop("Mean must be positive")
  }
  return(list(mu = mu, sigma2 = sigma2))
}

#-----------------------------------------------------------------------------------------------#
#### R function to change probabilities to rates and rates to probabilities  ####
#-----------------------------------------------------------------------------------------------#
ProbToRate <- function(p){
  # argument
  # p : the probability
  # Retunrs:
  # r : rate
  r <- -log(1 - p)
  r
}

RateToProb <- function(r, t){
  p <- 1 - exp(-r * t)
  return(p)
}
#-----------------------------------------------------------------------------------------------#
#### R function to generate net monetary and net heatlh benefit  ####
#-----------------------------------------------------------------------------------------------#

# Function for net health benefit (NBM)
calculateNHB <- function(effectiveness, costs, WTP){
  NHB <- effectiveness - costs / WTP
  return(NHB)
}

# Function for net monetary benefit
calculateNMB <- function(effectiveness, costs, WTP) {
  NMB <- effectiveness * WTP - costs
  return(NMB)
}

#-----------------------------------------------------------------------------------------------#
#### R function that converts a VAS score into a standard gamble equivalent  ####
#-----------------------------------------------------------------------------------------------#
convertVAStoUtility <- function(vas_score, r){
  # Arguments
  ## vas_score: as reported by individual 
  ## r: conversion factor ranging between 1.6 and 2.3
  # Returns 
  ## The utility values 
  utility <- 1 - (1 - vas_score/100) ^ r
  return(utility)
}

#-----------------------------------------------------------------------------------------------#
#### R function for Value of Information Analysis  ####
#-----------------------------------------------------------------------------------------------#
#### Formatting functions ####
# Run them all before continuing!
# Function for number of axis ticks in ggplot2 graphs
number_ticks <- function(n) {function(limits) pretty(limits, n)} 
# Total population affected by the decision
TotPop <- function(time, prev, incid, disc = 0){
  # Computes total population afected by technology
  #
  # Args:
  #   time:  vector with time points defining technology lifetime
  #   prev:  present prevalence
  #   incid: incidence
  #   disc:  discount factor; deafult = 0.
  #
  # Returns:
  #   tot.pop: total population afected by technology over technology lifetime
  #  
  # Technology Life Time, the last entry of vector `time`
  LT            <- time[length(time)]
  # Vector with population afected by the technolgy at each time point
  pop_time      <- c(prev, rep(incid, (length(time)-1))) 
  # Vector with present value of population afected by the technolgy at each time point
  disc_pop_time <- pop_time/(1+disc)^time
  # Total population afected by the technology
  tot_pop <-sum(disc_pop_time)
}
# Cost of Research
CostRes <- function(fixed.cost = 0, 
                    samp_size, 
                    cost_per_patient, 
                    INMB, 
                    clin_trial = TRUE, n_arms = 2){
  # Computes the cost of collecting information (i.e., through a research study)
  #
  # Args:
  #   fixed_cost:       fixed cost of collecting information
  #                     (e.g., fixed cost of a clinical trial); default = 0
  #   samp_size:               vector with sample sizes
  #   cost_per_patient: cost per patient in research study
  #   INMB:             Incremental Net Monetary Benefit
  #   clin_trial:       indicator whether calculation is for a clinical trial;
  #                     default = TRUE
  #   n_arms:           Number of arms in research study design; default = 2
  #
  # Returns:
  #   cost_res: vector with the total cost of collecting information for each simple size
  #
  if (clin_trial){
    Cost_Res <- fixed_cost + n_arms*samp_size*cost_per_patient + samp_size*INMB
  } else { # E.g., cohort study
    Cost_Res <- fixed_cost + samp_size*cost_per_patient
  }
  return(Cost_Res)
}


#-----------------------------------------------------------------------------------------------#
#### R function that converts a VAS score into a standard gamble equivalent  ####
#-----------------------------------------------------------------------------------------------#

# Function that converts a VAS score into a standard gamble equivalent
convertVAStoUtility <- function(vas_score, r){
  # Arguments
  ## vas_score: as reported by individual 
  ## r: conversion factor ranging between 1.6 and 2.3
  # Returns 
  ## The utility values 
  utility <- 1 - (1 - vas_score/100) ^ r
  return(utility)
}

beta_params <- function(mean, sigma) {
  alpha <- ((1 - mean) / sigma ^ 2 - 1 / mean) * mean ^ 2
  beta  <- alpha * (1 / mean - 1)
  params <- list(alpha = alpha, beta = beta)
  return(params)
}


#-----------------------------------------------------------------------------------------------#
#### R function to change probabilities to rates and rates to probabilities  ####
#-----------------------------------------------------------------------------------------------#
ProbToRate <- function(p){
  # argument
  # p : the probability
  # Retunrs:
  # r : rate
  r <- -log(1 - p)
  r
}

RateToProb <- function(r, t){
  p <- 1 - exp(-r * t)
  return(p)
}

