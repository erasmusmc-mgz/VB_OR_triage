

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
    if(param[i, "Distribution"] == "Triangle"){
      param_psa$psa_est[param_psa$Population == param$Population[i] & 
                        param_psa$Param   == param$Param[i]] <- with(param[i, ], 
                                   rtriangle(n = n_iter,
                                             a = Lo, 
                                             b = Hi, 
                                             c = Med))
      
    }else{ #if distribution is not triangle
      if(param[i, "Distribution"] == "Normal"){
        param_psa$psa_est[param_psa$Population == param$Population[i] & 
                          param_psa$Param      == param$Param[i]] <- with(param[i, ], 
                                                                         rnorm(n = n_iter ,
                                                                               mean = Med,
                                                                               sd = (Hi - Med)/1.96))
      }else{ #if distribution is not triangle, nor normal
        if(param[i, "Distribution"] == "Lognormal"){
          param_psa$psa_est[param_psa$Population == param$Population[i] & 
                            param_psa$Param      == param$Param[i]] <- exp(with(param[i, ], 
                                                                           rnorm(n = n_iter ,
                                                                                 mean = log(Med),
                                                                                 sd = (log(Hi) - log(Med)) / 1.96)))
        }else{ #if distribution is not triangle/normal/lognormal
          if(param[i,"Distribution"] == "Beta"){
            param_psa$psa_est[param_psa$Population == param$Population[i] & 
                                param_psa$Param      == param$Param[i]] <- with(param[i, ],
                                                                                rbeta(n = n_iter,
                                                                                 shape1 = Lo,
                                                                                 shape2 = Hi)) 
          }
        }
      }
    }
  }
  ## For diseases with no Surv_notx, calculate this from the Surv_tx and treatment effect
  for(p in param$Population[param$Distribution == "Calc_from_tx_eff"]){
    #In case we have a hazard ratio
    if(param$Unit[param$Population == p&param$Param == "Tx_eff"] == "HR"){
      
      # Since HR = hazard_tx/hazard_notx, harzard_notx = hazard_tx/HR
      h_tx <- ProbToRate(param_psa$psa_est[param_psa$Population == p &
                                  param_psa$Param == "Surv_tx"])
      hrs  <- param_psa$psa_est[param_psa$Population == p &
                                  param_psa$Param == "Tx_eff"]
      param_psa$psa_est[param_psa$Population == p &
                          param_psa$Param == "Surv_no_tx"] <- RateToProb(r = h_tx/hrs, t = 1) 
    }
    #In case we have a relative risk
    if(param$Unit[param$Population == p & param$Param == "Tx_eff"] =="RR"){
      
      # Since RR = prob_tx/prob_notx, prob_notx = prob_tx/RR
      s_tx <- param_psa$psa_est[param_psa$Population == p&
                                  param_psa$Param == "Surv_tx"]
      rrs  <- param_psa$psa_est[param_psa$Population == p &
                                  param_psa$Param == "Tx_eff"]
      param_psa$psa_est[param_psa$Population == p &
                          param_psa$Param == "Surv_no_tx"] <- s_tx/rrs 
    }
    #In case we have a odds ratio
    if(param$Unit[param$Population == p & param$Param == "Tx_eff"] =="OR"){
      
      # Since OR = odds_tx/odds_notx, odds_notx = odds_tx/HR
      s_tx <- param_psa$psa_est[param_psa$Population == p &
                                  param_psa$Param == "Surv_tx"]
      o_tx <- s_tx/(1-s_tx)
      ors  <- param_psa$psa_est[param_psa$Population == p &
                                  param_psa$Param == "Tx_eff"]
      
      o_notx <- o_tx/ors
      param_psa$psa_est[param_psa$Population == p &
                          param_psa$Param == "Surv_no_tx"] <- o_notx/(1 + o_notx)
    } 
  }
  
  # If Utilities are equal pre/post intervention (QoL_tx vs QoL_no_tx), make them equal
  ## First: identify which populations have equal QoL pre/post
  which_pop_df <- param[param$Param%in%c("QoL_tx", "QoL_no_tx"),c("Population","Param","Med")]
  which_pop    <- rep(0, nrow(which_pop_df))
  for(i in seq(2, length(which_pop))){if(which_pop_df$Med[i] == which_pop_df$Med[i-1]){which_pop[i] <- 1}}
  
  ## Second: set the estimates of the PSA for QoL_no_tx to QoL_tx
  for(pop in which_pop_df$Population[which_pop == 1]){
    param_psa$psa_est[
      param_psa$Population == pop &
        param_psa$Param   == "QoL_no_tx"] <- # Replace by the values of the PSA estimates of QoL_tx  
      param_psa$psa_est[param_psa$Population == pop &
                          param_psa$Param == "QoL_tx"]
   } 

  return(param_psa)  # Return the parameter values for the PSA runs
}  
