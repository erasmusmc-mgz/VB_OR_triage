
# Load functions 
source("R/app0_functions.R")  # Load a function to install packages
source("R/functions.R")       # Load general functions useful for Markov models 
source("R/functions_PSA.R")   # Load the PSA function  
source("R/model.R")           # Code of the main model

# Only required once, 
# uncommend when running the code for the first time
# v_packages_to_install <- c("ggplot2", "data.table", "triangle")
# install_and_load(v_packages_to_install)
# install_github("DARTH-git/dampack", force = TRUE)

# Load packages
library(ggplot2)
library(data.table)
library(triangle)
library(dampack)


# Load data
param <- data.frame(readxl::read_xlsx("Data/Model parameters.xlsx"))
cbs   <- read.csv("Data/CBS lifetable.csv", sep = ";")

# document which paramaters do not have a unit and will be deleted
param_NA <- param[which(is.na(param$Unit)), ]

# Remove parameters without unit
param <- param[!is.na(param$Unit), ]

# Make variables numeric
numvar <- c("Med", "Lo", "Hi")
for (i in numvar) {param[, i] <- as.numeric(param[, i])}

#############################
## Survival to right units ##
#############################

# ## median survival times --> hazard to die 
 param[param$Unit == "Median survival time (weeks)", c("Med", "Lo", "Hi")]  <- log(2)/param[param$Unit == "Median survival time (weeks)", c("Med", "Hi", "Lo")]
# explanation convertion: https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Survival_Parameter_Conversion_Tool.pdf

## 10-year survival rate --> probability to die per week
param[param$Unit == "Probability 10-year survival", c("Med", "Lo", "Hi")]  <- RateToProb(ProbToRate(1 - param[param$Unit == "Probability 10-year survival", c("Med", "Hi", "Lo")]), 1/(10 * 52))

## 5-year survival rate --> probability to die per week
param[param$Unit == "Probability 5-year survival", c("Med", "Lo", "Hi")]  <- RateToProb(ProbToRate(1 - param[param$Unit == "Probability 5-year survival", c("Med", "Hi", "Lo")]), 1/(5 * 52))

## 3-year survival rate --> probability to die per week
param[param$Unit == "Probability 3-year survival", c("Med", "Lo", "Hi")]  <- RateToProb(ProbToRate(1 - param[param$Unit == "Probability 3-year survival", c("Med", "Hi", "Lo")]), 1/(3 * 52))


## 1-year survival rate --> probability to die per week
param[param$Unit == "Probability 1-year survival", c("Med", "Lo", "Hi")]  <- RateToProb(ProbToRate(1 - param[param$Unit == "Probability 1-year survival", c("Med", "Hi", "Lo")]), 1/(1 * 52))

## 30-day survival rate --> probability to die per week
param[param$Unit == "Probability 30 day survival", c("Med", "Lo", "Hi")]  <- RateToProb(ProbToRate(1 - param[param$Unit == "Probability 30 day survival", c("Med", "Hi", "Lo")]), 1/30)


## Mortality rate per person year --> probability to die per week
param[param$Unit == "Mortality rate per person-year", c("Med","Lo","Hi")] <- param[param$Unit == "Mortality rate per person-year", c("Med","Lo","Hi")]/52  

#### CBS data ####
#CBS data: prob per year --> prob per week
cbs$Man_kans   <- RateToProb(ProbToRate(p = cbs$Man_kans),   1/52)
cbs$Vrouw_kans <- RateToProb(ProbToRate(p = cbs$Vrouw_kans), 1/52)

#CBS data: prob per sex --> average prob
cbs$prob <- (cbs$Man_kans + cbs$Vrouw_kans) / 2

#CBS data: Leeftijd --> age
cbs$age <- as.numeric(substr(cbs$Leeftijd, start = 1, stop = 2))


######################
## Other units      ##
######################
## Days --> weeks
param[param$Unit == "Days", c("Med", "Lo", "Hi")]  <- param[param$Unit == "Days", c("Med", "Lo", "Hi")]/7

## Months --> weeks
param[param$Unit == "Months", c("Med", "Lo", "Hi")]  <- param[param$Unit == "Months", c("Med", "Lo", "Hi")]/12*52


### Model structure 
state_names <- c("Preop", "Postop", "Dead")  # names of health states
n_s         <- length(state_names)           # lenght of the health states

# Delay in weeks
delay <- seq(from = 2, to = 52, by = 10)

#Also add a no-operation state (delay is infinite)
delay <- c(delay, 999)

# Create a dataframe to store results in
n_iter <- 100 # Number of iterations in the PSA


# Create a data frame to store the results (df_res)
# expand grid makes a date frame for all combinations of the supplied vectors
df_res <- expand.grid(iter = 1:n_iter,                  # PSA iteration
                      Pop = unique(param$Population),   # Disease name 
                      delay = delay,                    # time of delay
                      QALY = NA,                        # QALYs
                      LY = NA,                          # LY: qol=1
                      AAC = NA,
                      AAC_ly = NA,# Area above the curve (the loss in QALYs)
                      AAC_delay = NA,
                      AAC_delay_ly = NA)                # AAC/delay       

# Add the surgeries/labels per population in the intervention column
df_res$Intervention <- param$Intervention[match(df_res$Pop, param$Population)]
df_res$Label        <- param$Label[match(df_res$Pop, param$Population)]
df_res$Source       <- param$Source[match(df_res$Pop, param$Population)]

pop_names <- sort(unique(param$Population))  # substract the disease names and order them in alphabetic order 

# Make a dataframe with PSA parameters 
param_psa <- make_psa_df(param = param, n_iter = n_iter, seed = 19)

# Remove the unit column, since these units are now all adjusted to weekly probabilities or age in years
param_psa <- subset(param_psa, select = -c(Unit)) # Store the weekly probabilities in the dataframe param


# loop over every disease 
for (d in pop_names){
  for(it in 1:n_iter){   # for every PSA iteration
    
    # PSA estimates
    # Extract vector of parameters
    p_vector        <- param_psa$psa_est[param_psa$iter == it & param_psa$Population == d]
    # give the corresponding parameter names
    names(p_vector) <- param_psa$Param  [param_psa$iter == it & param_psa$Population == d]
    
    ##BUILD IN ERROR MESSAGE 
    # The we make use of 9 parameters. In case we have more then 9 parameter we have two populations with the same name.
    if(length(p_vector) > 9){
      print(paste("There are populations with the same name (",d,")"))
      break
      }
     
    # Round age to the nearest number
    p_vector["Age"] <- round(p_vector["Age"], 0)
    
    # Time horizon (in weeks)
    n_years     <- as.vector(100 - p_vector["Age"])  # time horizon in years
    n_cycles    <- 52 * n_years           # number of model cycles
    
    
    # create the transition matrix
    m_trans <- make_m_trans(number_states = n_s,
                            state_names = state_names,
                            number_cycles = n_cycles, 
                            parameters = p_vector,
                            cbsdata = cbs)
    
    # rewards vectors (useful for state rewards)
    v_utility <- numeric(n_s)
    v_utility <- c(p_vector["QoL_no_tx"], p_vector["QoL_tx"], 0)
    names(v_utility) <- state_names
    
    # The utility vector, where there is no gain in QoL of the surgery
    v_utility_noeff           <- v_utility
    v_utility_noeff["Postop"] <- v_utility_noeff["Preop"]
    
    
    # Start the senario analysis for different delays. 
    for(i in 1:length(delay)){
      
      ## Include delay in transition matrix
      m_trans <- trans_operation(trans_matrix   = m_trans,
                                 weeks_until_op = delay[i],
                                 number_cycles   = n_cycles)
      
      
      # create the Markov trace for the individuals 
      m_trace <- matrix(NA, nrow = n_cycles + 1, ncol = n_s,
                        dimnames = list(paste("cycle", 0:n_cycles, sep = " "), state_names)) 
      
      # initialize the start health state
      m_trace[1, ] <- c(1, 0, 0)  # all individuals start out in the Preop health states
      
      
      ## Calculate Markov trace for this senario of delay
      for (t in 1:n_cycles){     # loop through the number of cycles
        m_trace[t + 1, ] <- t(m_trace[t, ]) %*% m_trans[, , t]    # estimate the Markov trace 
        # for the next cycle (t + 1)
      } # close the loop for the markov model 
      
      # Calculate the QALYs
      # Vector with the total effects (QALYs per cycle)
      if("Time_noeff_QoL" %in% names(p_vector) &   # If the effect on QoL is lost after a particular time
         delay[i] > p_vector["Time_noeff_QoL"]){ # And if the delay is larger than the time untill the effect is lost
        v_tu  <- (m_trace   %*%  v_utility_noeff)  # use the no-effect vector
      }else{
        v_tu  <- (m_trace   %*%  v_utility)        # use the utility vector wíth effect
      }
      
      # Make discount vector 
      d_e_effect <- 0.015 # discount rate of 0.015 
      times <- seq(0, n_cycles/52, by = 1/52)  # factor to adjust discount weights for cycle length. Every cycle is 3 months
      v_dwe      <- 1 / ((1 + d_e_effect) ^ (times)) # make a vector with discount weight
      
      # Calculate the total QALYs
      QALY <- (t(v_tu) %*%  v_dwe) / n_years
      
      # Store discounted total QALY result 
      df_res[df_res$iter == it & df_res$Pop == d & df_res$delay == delay[i], "QALY"] <- QALY[1]
      
      # Calculate the LYs
      v_tu_ly <- m_trace %*% c(1,1,0)
      
      # Apply discount
      LY      <- (t(v_tu_ly) %*%  v_dwe) / n_years
      
      #Store in results
      df_res[df_res$iter == it & df_res$Pop == d & df_res$delay == delay[i], "LY"] <- LY
      
      # Calculate the area above the curve for QALYs 
      area.total <- df_res[df_res$Pop == d & df_res$iter == it, ]$QALY[1] * df_res[df_res$Pop == d & df_res$iter == it, ]$delay[i]
      auc <- DescTools::AUC(x    = df_res[df_res$Pop == d & df_res$iter==it, ]$delay, 
                            y    = df_res[df_res$Pop == d & df_res$iter==it, ]$QALY, 
                            from = df_res[df_res$Pop == d & df_res$iter==it, ]$delay[1], 
                            to   = df_res[df_res$Pop == d & df_res$iter==it, ]$delay[i], 
                            method = "spline") # Calculate the area under the curve 
      aac <- area.total - auc # area above the curve 
      
      # Store the results
      df_res[df_res$iter==it & df_res$Pop == d & df_res$delay == delay[i], "AAC"] <- aac
      df_res[df_res$iter==it & df_res$Pop == d & df_res$delay == delay[i], "AAC_delay"] <- aac/delay[i]
      
      # Calculate the area above the curve for LYs 
      area.totally <- df_res[df_res$Pop == d & df_res$iter == it, ]$LY[1] * df_res[df_res$Pop == d & df_res$iter == it, ]$delay[i]
      aucly <- DescTools::AUC(x    = df_res[df_res$Pop == d & df_res$iter==it, ]$delay, 
                            y    = df_res[df_res$Pop == d & df_res$iter==it, ]$LY, 
                            from = df_res[df_res$Pop == d & df_res$iter==it, ]$delay[1], 
                            to   = df_res[df_res$Pop == d & df_res$iter==it, ]$delay[i], 
                            method = "spline") # Calculate the area under the curve 
      aacly <- area.totally - aucly # area above the curve 
      
      # Store the results
      df_res[df_res$iter==it & df_res$Pop == d & df_res$delay == delay[i], "AAC_ly"] <- aacly
      df_res[df_res$iter==it & df_res$Pop == d & df_res$delay == delay[i], "AAC_delay_ly"] <- aacly/delay[i]
    }
    
    # Display simulation progress
    if(it/(n_iter/10) == round(it/(n_iter/10), 0)) { # display progress every 10%
      cat('\r', paste(it/n_iter * 100, "% PSA of disease", d, "done & ", round(which(pop_names == d)/length(pop_names) * 100, 0), "% total PSA", sep = " "))
    }
  } # close the number if PSA iterations
} # close the loop for diseases


df_res$AAC_delay[is.nan(df_res$AAC_delay)] <- 0 # Divided by 0 is not relevant --> put to 0
df_res$AAC_delay_ly[is.nan(df_res$AAC_delay_ly)] <- 0 # Divided by 0 is not relevant --> put to 0


# Pool results PSA
df_res <- data.table(df_res)
results_pooled <- df_res[,.(QALY_med      = median(QALY),
                            QALY_lo       = quantile(QALY, probs = 0.025),
                            QALY_hi       = quantile(QALY, probs = 0.975),
                            LY_med        = median(LY),
                            LY_lo         = quantile(LY, probs = 0.025),
                            LY_hi         = quantile(LY, probs = 0.975),
                            AAC_med       = median(AAC),
                            AAC_lo        = quantile(AAC, probs = 0.025),
                            AAC_hi        = quantile(AAC, probs = 0.975),
                            AAC_delay_med = median(AAC_delay),
                            AAC_delay_lo  = quantile(AAC_delay, probs = 0.025),
                            AAC_delay_hi  = quantile(AAC_delay, probs = 0.975),
                            AAC_delay_ly_med = median(AAC_delay_ly),
                            AAC_delay_ly_lo  = quantile(AAC_delay_ly, probs = 0.025),
                            AAC_delay_ly_hi  = quantile(AAC_delay_ly, probs = 0.975)),
                         by = .(Label, delay)]


# save files in the output folder 
save(df_res,         file = "output/res_psa.Rdata")
save(results_pooled, file = "output/psa_pooled.Rdata")
save(param,          file = "output/input_param.Rdata")
save(param_psa,      file = "output/psa_parameters.RData")

# plot results of all diseases seperately
plotPopulationOutcomes(data = results_pooled[results_pooled$delay!=999,],
                       folder = "figures/QALY_per_pop/",
                       size_cm = 15) 

# Calculate the QALY loss per week (the derivative)

# for(p in unique(results_pooled$Label)){
#   res_derivative <- calculateDerivative(df_pooled = results_pooled[results_pooled$Label==p&
#                                                                      results_pooled$delay!=999,,], 
#                                         plot = TRUE, plot_ind = FALSE, cum = FALSE, folder = "figures")
#   # Make a figure of that derivative
#   p <- gsub(' ', "_", p) # Replace the space with an unscore
#   ggsave(filename = paste("figures/QALY_per_pop/",p,"_derivatives.png", sep = ""), plot = res_derivative[[1]], device = "png")
# }





  