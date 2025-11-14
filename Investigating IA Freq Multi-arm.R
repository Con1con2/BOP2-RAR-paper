## IA Freq multi-arm (and RAR comparisons)
# updated

library(Rcpp) # Load package 'Rcpp'

library(MyPackageV2) # contains RCPP code for parallilisation 

library(devtools)

library(microbenchmark)

library(doParallel)  # Parallelisation
library(foreach)  # Parallelisation

library(roxygen2)

library(ggplot2)


library(viridis) # colours

library(ggpubr)


Multi_arm_BOP2_V2 = function(treats_eff, treats_tox, n = 60*length(treats_eff), IAn = c(n), control = FALSE, phi_e = 0.45, phi_t = 0.3, threshold_type = 1, gamma = 0.98, lamda = 0.81, exact = TRUE){
  ## Multi Arm BOP2 aims to extend BOP2 to the multi arm setting, like in Mulier et al 2024
  
  # treats_ ss a vector of treatments. If control = TRUE, then treats[1] will be the control 
  # treats_eff is the vector of efficacy (ex 0.7 is the probability of a desirable outcome)
  # treats_tox  is the vector of toxicity (ex 0.6 is the probabilty of toxic response)
  # So treats_eff = c(0.1,0.2,0.3) and treats_tox = c(0.5,0.4,0.3) correspond to 3 treatments
  # One with 0.1 success 0.5 tox,then 0.2 success 0.4 tox, then 0.3 success and tox
  # IAn is a vector of possible IA
  
  # n will mean how many patients in total
  # IAn is a vector of target IA
  # control = TRUE/FALSE determines whether there is a control or not. If TRUE, then treats[1] is a control
  # phi_e, phi_f are the efficacy/futility thresholds for the uncontrolled setting
  # lamda are gamma are parameters that are used to control the type one error whilst maximising power 
  
  # For ease of calculation, we will take a Dirichlet(1,1,1,1) prior 
  # (this allows us to use Beta_comp_exact)
  
  
  num_treat = length(treats_eff) # number of treatments
  continue = rep(1,num_treat)  # Will be 1 if we continue with a treatment, 0 futile (no test for efficacy)
  
  early_stop = numeric(num_treat) # Recording whether stopping for futility
  
  final_sample_size = numeric(num_treat) # where the final sample sizes will go 
  
  # Where the responses for efficacy and toxicity will go
  Responses_eff = vector("list",num_treat)
  Responses_tox = vector("list",num_treat)
  
  per_arm = IAn[1] %/% num_treat
  
  remainder = IAn[1] %% num_treat
  
  sampling_vec = rep(1:num_treat,per_arm)
  
  #browser()
  
  if (remainder > 0){
    # Randomly assigning the leftovers
    remainder_allo = rmultinom(1,remainder, rep(1/num_treat, num_treat))
    
    for (i in 1:num_treat){
      sampling_vec = append(sampling_vec, rep(i,remainder_allo[i]))  # allocating the remainders
    }
    
    
  }
  
  
  sampling_vec = sample(sampling_vec) # Randomises the order, now this is the sampling order
  
  for (i in sampling_vec){  # Initial bit of samples for every treatment before IA
    Responses_eff[[i]] = append(Responses_eff[[i]],rbinom(1, 1, treats_eff[i]))
    Responses_tox[[i]] = append(Responses_tox[[i]],rbinom(1, 1, treats_tox[i]))
  }
  
  ## First IA assignment done
  
  # Without control 
  if (control == FALSE){
    
    if (threshold_type == 1){  # There are varying threshold types the paper gives
      
      #Threshold type 1 is most similar to the original BOP2 paper. 
      
      threshold = 1 - (lamda * ((IAn[1] / n)^gamma))  #
    }
    
    
    for (i in 1:num_treat){
      
      ## Temp values to do the beta calcs with 
      temp_n = length(Responses_eff[[i]])
      temp_eff = sum(Responses_eff[[i]])
      temp_tox = sum(Responses_tox[[i]])
      
      #print(threshold)
      #print(temp_n)
      #print(temp_eff)
      #print(temp_tox)
      
       # Standard beta comp against a default
      
      if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
        
        
        continue[i] = 0  # This has either exceeded futility or toxicity
        
        if (length(IAn) != 1){
          early_stop[i] = 1 # Early stopped
        }
        
        if (identical(continue,numeric(num_treat)) ){  # If every arm has stopped
          
          # Quickly generating the OCs as the trial is over
          
          final_sample_size = numeric(num_treat)
          
          for (k in 1:num_treat){
            final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
          }
          
          return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
        }
        
      }
      
    }
    
    if (length(IAn) == 1){ # NO IA, SO TRIAL ENDS
      final_sample_size = numeric(num_treat)
      
      for (k in 1:num_treat){
        final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
      }
      
      return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
    }
    
    
    for (i in (2:length(IAn))){
      
      
      to_sample = numeric(0)
      for (j in 1:num_treat){
        if (continue[j] == 1){
          to_sample = append(to_sample, j) # which arms we should continue sampling from
        } 
      }
      
      per_arm = (IAn[i] - IAn[i-1]) %/% length(to_sample)
      
      remainder = (IAn[i] - IAn[i-1]) %% length(to_sample)
      
      sampling_vec = rep(to_sample,per_arm)
      
      
      if (remainder > 0){
        # Randomly assigning the leftovers
        remainder_allo = rmultinom(1,remainder, rep(1/length(to_sample), length(to_sample)))
        
        for (j in 1:num_treat){
          sampling_vec = append(sampling_vec, rep(i,to_sample[remainder_allo[i]]))  # allocating the remainders
        }
        
        
      }
      
      sampling_vec = sample(sampling_vec) # Randomises the order, now this is the sampling order
      
      for (j in sampling_vec){  # Initial bit of samples for every treatment before IA
        Responses_eff[[i]] = append(Responses_eff[[i]],rbinom(1, 1, treats_eff[i]))
        Responses_tox[[i]] = append(Responses_tox[[i]],rbinom(1, 1, treats_tox[i]))
      }
      # Sampling done
      
      ## Onto the IA
      
      
      for (i in 1:num_treat){  # Looking at every treatment
        
        if (continue[i] == 1){  # We only want to analyse the arms being updated
          
          ## Temp values to do the beta calcs with 
          temp_n = length(Responses_eff[[i]])
          temp_eff = sum(Responses_eff[[i]])
          temp_tox = sum(Responses_tox[[i]])
          
          
          if (threshold_type == 1){
            
            
            threshold = 1 - (lamda * ((IAn[i] / n)^gamma))   #
          }
          
          #print(threshold)
          #print(temp_n)
          #print(temp_eff)
          #print(temp_tox)
          
          
          
          
          if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
            
            continue[i] = 0  # This has either exceeded futility or toxicity
            
            if (j != num_treat){
              early_stop[i] = 1 # Early stopped
            }
            
            if (identical(continue,numeric(num_treat)) ){  # If every arm has stopped
              
              # Quickly generating the OCs as the trial is over
              
              final_sample_size = numeric(num_treat)
              
              for (k in 1:num_treat){
                final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
              }
              
              return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
            }
            
          }
        }  
      }
      
      
    }
    
        
    
        
    
    for (k in 1:num_treat){
      final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
    }
      
      
      
  }
    
  
  ## With control 
  else{
    
    if (threshold_type == 1){  # There are varying threshold types the paper gives
      
      #Threshold type 1 is most similar to the original BOP2 paper. 
      
      threshold = 1 - (lamda * ((IAn[1] / n)^gamma))
    }
    
    con_n = length(Responses_eff[[1]])  ## Control results
    con_eff = sum(Responses_eff[[1]])
    con_tox = sum(Responses_tox[[1]])
    
    
    for (i in 2:num_treat){
      
      ## Temp values to do the beta calcs with 
      temp_n = length(Responses_eff[[i]])
      temp_eff = sum(Responses_eff[[i]])
      temp_tox = sum(Responses_tox[[i]])
      
      #print(threshold)
      #print(temp_n)
      #print(temp_eff)
      #print(temp_tox)
      
      if (exact){  # If exact, the computation is a bit quicker, but the priors have to be changed
        if ( rcpp_exact_beta( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  |  rcpp_exact_beta( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
          
          continue[i] = 0  # This has either exceeded futility or toxicity
          
          if (length(IAn) != 1){
            early_stop[i] = 1 # Early stopped
          }
          
          if (identical(continue,numeric(num_treat)) ){  # If every arm has stopped
            
            # Quickly generating the OCs as the trial is over
            
            final_sample_size = numeric(num_treat)
            
            for (k in 1:num_treat){
              final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
            }
            
            return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
          }
          
        }
        
      }
      else{
        if ( Beta_comp( 0.45 + temp_eff, 0.55 + temp_n - temp_eff, 0.45 + con_eff, 0.55 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 0.3 + temp_tox, 0.7 + temp_n - temp_tox, 0.3 + con_tox, 0.7 + con_n - con_tox) > threshold ){
          
          continue[i] = 0  # This has either exceeded futility or toxicity
          
          if (j != num_treat){
            early_stop[i] = 1 # Early stopped
          }
          
          if (identical(continue,numeric(num_treat)) ){  # If every arm has stopped
            
            # Quickly generating the OCs as the trial is over
            
            final_sample_size = numeric(num_treat)
            
            for (k in 1:num_treat){
              final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
            }
            
            return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
          }
          
        }
        
      }
    }
    
    if (length(IAn) == 1){ # NO IA, SO TRIAL ENDS
      final_sample_size = numeric(num_treat)
      
      for (k in 1:num_treat){
        final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
      }
      
      return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
    }
      
      for (i in (2:length(IAn))){
        to_sample = numeric(0)
        for (j in 1:num_treat){
          if (continue[j] == 1){
            to_sample = append(to_sample, j) # which arms we should continue sampling from
          } 
        }
        
        per_arm = (IAn[i] - IAn[i-1]) %/% length(to_sample)
        
        remainder = (IAn[i] - IAn[i-1]) %% length(to_sample)
        
        sampling_vec = rep(to_sample,per_arm)
        
        if (remainder > 0){
          # Randomly assigning the leftovers
          remainder_allo = rmultinom(1,remainder, rep(1/length(to_sample), length(to_sample)))
          
          for (j in 1:num_treat){
            sampling_vec = append(sampling_vec, rep(i,to_sample[remainder_allo[i]]))  # allocating the remainders
          }
          
          
        }
        
        sampling_vec = sample(sampling_vec) # Randomises the order, now this is the sampling order
        
        for (j in sampling_vec){  # Initial bit of samples for every treatment before IA
          Responses_eff[[i]] = append(Responses_eff[[i]],rbinom(1, 1, treats_eff[i]))
          Responses_tox[[i]] = append(Responses_tox[[i]],rbinom(1, 1, treats_tox[i]))
        }
        # Sampling done
        
        ## Onto the IA
        
        for (j in 2:num_treat){  # Looking at every treatment, not the control
          
          if (continue[j] == 1){  # We only want to analyse the arms being updated
            
            ## Temp values to do the beta calcs with 
            temp_n = length(Responses_eff[[j]])
            temp_eff = sum(Responses_eff[[j]])
            temp_tox = sum(Responses_tox[[j]])
            
            
            if (threshold_type == 1){
              
              
              threshold = 1 - (lamda * ((IAn[i] / n)^gamma))
            }
            
            #print(threshold)
            #print(temp_n)
            #print(temp_eff)
            #print(temp_tox)
            
            
            
            
            if (exact){
              if ( rcpp_exact_beta( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  |  rcpp_exact_beta( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
                
                continue[i] = 0  # This has either exceeded futility or toxicity
                
                if (j != num_treat){
                  early_stop[i] = 1 # Early stopped
                }
                
                if (identical(continue,numeric(num_treat)) ){  # If every arm has stopped
                  
                  # Quickly generating the OCs as the trial is over
                  
                  final_sample_size = numeric(num_treat)
                  
                  for (k in 1:num_treat){
                    final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
                  }
                  
                  return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
                }
                
              }
              
            }
            else{
              if ( Beta_comp( 0.45 + temp_eff, 0.55 + temp_n - temp_eff, 0.45 + con_eff, 0.55 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 0.3 + temp_tox, 0.7 + temp_n - temp_tox, 0.3 + con_tox, 0.7 + con_n - con_tox) > threshold ){
                
                continue[i] = 0  # This has either exceeded futility or toxicity
                
                if (j != num_treat){
                  early_stop[i] = 1 # Early stopped
                }
                
                if (identical(continue,numeric(num_treat)) ){  # If every arm has stopped
                  
                  # Quickly generating the OCs as the trial is over
                  
                  final_sample_size = numeric(num_treat)
                  
                  for (k in 1:num_treat){
                    final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE
                  }
                  
                  return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
                }
                
              }
              
            }
          }
          
        }
        
      }
      
      
    for (k in 1:num_treat){
      final_sample_size[k] = length(Responses_eff[[k]]) # SAMPLE SIZE 
    }
    
      
    }
    
    
    
  
  
  
  # Results 
  # print(Responses_eff)
  # print(Responses_tox)  
  
  to_return = list( early_stop, final_sample_size, continue ) # The three variables of interst
  # Need these to replicate their results
  
  return(to_return)
  
}

Multi_arm_BOP2_Max_PA = function(treats_eff, treats_tox, n_pa = 60, IAn_pa = c(n), control = FALSE, phi_e = 0.45, phi_t = 0.3, gamma = 0.71, lamda = 0.83, PairedSeed = FALSE){
  ## Same as above but this one focuses on per-arm. 
  # In particular, n_pa, and IAn_pa refer to the maximum and IA choices per arm. 
  # So n_pa = 60, IAn = c(15, 30, 45, 60) would match the Mulier paper. 
  
  
  if (PairedSeed != FALSE){  
    set.seed(PairedSeed) # sets seed for paired comparisons
  }
  
  
  num_treat = length(treats_eff)
  
  num_IA = length(IAn_pa) # Num of IAs
  continue = rep(1,num_treat) # Indicates whether we should continue 
  
  Responses_eff = vector("list",num_treat)
  Responses_tox = vector("list",num_treat)
  
  # OCs that will be filled in later
  early_stop = numeric(num_treat)
  final_sample_size = numeric(num_treat)
  
  
  for (i in 1:num_IA){
    # Going per IA
    if (i == 1){  # WOrks out how many to sample per IA
      to_sample = IAn_pa[1]
    }
    else{
      to_sample = IAn_pa[i] - IAn_pa[i-1]
    }
    
    
    for (j in which(continue == 1)){
      Responses_eff[[j]] = append(Responses_eff[[j]],rbinom(to_sample, 1, treats_eff[j]))
      Responses_tox[[j]] = append(Responses_tox[[j]],rbinom(to_sample, 1, treats_tox[j]))
      
    }
    # Sampling done
    
    ## IA
    if (control){
      # Control beta comparisons
      
      con_n = length(Responses_eff[[1]])  ## Control results
      con_eff = sum(Responses_eff[[1]])
      con_tox = sum(Responses_tox[[1]])
      
      for (j in 2:num_treat){
        # Skipping the control
        
        if (continue[j] == 1){
          # Performs the IA on arms still recreuiting
          
          
          
            
          ## Temp values to do the beta calcs with 
          temp_n = length(Responses_eff[[j]])
          temp_eff = sum(Responses_eff[[j]])
          temp_tox = sum(Responses_tox[[j]])
            
          #print(threshold)
          #print(temp_n)
          #print(temp_eff)
          #print(temp_tox)
          
          threshold = 1 - (lamda * ((IAn_pa[i] / n_pa)^gamma))   # threshold based on how far through 
          
            
          
          if ( 1 - rcpp_exact_beta( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  | 1 - rcpp_exact_beta( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
          
            continue[j] = 0
            
            if (i != num_IA){
              early_stop[j] = 1  # Early stops if rejection happens before end 
              
            }
            final_sample_size[j] = temp_n # Records final sample size
                  
          }
          
          
        }
        
      }
      # Finished this round of IA
      if (identical(c(1,rep(0,num_treat - 1)), continue)){
        # If all of the non control arms have finished their runs
        
        if (i != num_IA){
          early_stop[1] = 1 # Early stops if rejection happens before end 
          
        }
        
        
        final_sample_size[1] = temp_n
        
        return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
        
      }
      
      # Loops back
    }
    else{
      # Uncontrolled analysis
      for (j in 1:num_treat){
        # Not skipping a control
        
        if (continue[j] == 1){
          # Performs the IA on arms still recreuiting
          
          temp_n = length(Responses_eff[[j]])
          temp_eff = sum(Responses_eff[[j]])
          temp_tox = sum(Responses_tox[[j]])
          
            
            
          threshold = 1 - (lamda * ((IAn_pa[i] / n_pa)^gamma))   # threshold based on how far through 
          
          #print(threshold)
          #print(temp_n)
          #print(temp_eff)
          #print(temp_tox)
          
          
          
          
          if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
            # If we get a rejection
            
            continue[j] = 0
            
            if (i != num_IA){
              early_stop[j] = 1  # Early stops if rejection happens before end 
              
            }
            final_sample_size[j] = temp_n # Records final sample size
          
          }
        }
        
      }
      # Finished this round of IA
      
      if (identical(c(rep(0,num_treat)), continue)){
        # If all of the arms have finished their runs
        
        
        return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
        
      }
      
      # Loops back
      
      
    }
    
  }
  
  # If we get here, that means that we've reached the end and there is at least one
  # significant non-control arm
  
  for (i in 1:num_treat){
    # all that needs updating are the final sample sizes 
    if (final_sample_size[i] == 0){
      final_sample_size[i] = n_pa
    }
  }
  
  
  
  return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
  
  
}

Multi_arm_BOP2_V3 = function(treats_eff, treats_tox, n = 60*length(treats_eff), IAn = c(n), control = FALSE, phi_e = 0.45, phi_t = 0.3, gamma = 0.71, lamda = 0.83, PairedSeed = FALSE){
  ## This uses the skeleton of Multi_arm_BOP2_Max_PA to try and improve Multi_arm_BOP2_V2
  # As Multi_arm_BOP2_Max_PA has a better skeleton, all we need to do is change the allo section
  
  if (PairedSeed != FALSE){  
    set.seed(PairedSeed) # sets seed for paired comparisons
  }
  
  num_treat = length(treats_eff)
  
  num_IA = length(IAn) # Num of IAs
  continue = rep(1,num_treat) # Indicates whether we should continue 
  
  Responses_eff = vector("list",num_treat)
  Responses_tox = vector("list",num_treat)
  
  # OCs that will be filled in later
  early_stop = numeric(num_treat)
  final_sample_size = numeric(num_treat)
  
  sample_order = numeric(0) # NEW, records the sampling order (interest of RAR)
  
  
  for (i in 1:num_IA){
    # Going per IA
    if (i == 1){  # Works out how many to sample per IA
      
      to_sample = IAn[1] # how many will be sampled in total
      
      # denotes what each arm gets min 
      each_arm = to_sample %/% num_treat
      
      # What will be used to choose the samples
      sample_vec = rep(each_arm, num_treat)
      
      if (each_arm != IAn[1]){
        # There are leftovers
        extra = rmultinom(1,to_sample %% num_treat, rep(1/num_treat, num_treat))
        sample_vec = extra + sample_vec
      }
      
      ## Now just randomisaing an order, (As this is PBD, won't grow too apart here)
      
      for (k in 1:num_treat){
        sample_order = append(sample_order, rep(k, sample_vec[k]))
      }
      sample_order = sample(sample_order)
      
      
    }
    else{
      to_sample = IAn[i] - IAn[i-1] # how many will be sampled in total
      
      # This is where the RAR goes, for the BRAR version. For now, it is just the same as above
      
      active_arms = sum(continue) # We only care about the remaining arms
      
      # denotes what each arm gets min 
      each_arm = to_sample %/% active_arms
      
      # What will be used to choose the samples
      sample_vec = rep(each_arm, active_arms)
      
      if (each_arm != (IAn[i] - IAn[i-1])){
        # There are leftovers
        extra = rmultinom(1,to_sample %% active_arms, rep(1/active_arms, active_arms))
        sample_vec = extra + sample_vec
      }
      
      ## Now just randomisaing an order, (As this is PBD, won't grow too apart here)
      
      new_sample_order = numeric(0)
      
      for (k in 1:active_arms){
        new_sample_order = append(new_sample_order, rep(which(continue == 1)[k], sample_vec[k]))
      }
      new_sample_order = sample(new_sample_order)
      
      sample_order = append(sample_order, new_sample_order)
      
      #browser()
      
    }
    
    counter = 1  # Need a counter to accurately reference sample_Vec
    for (j in which(continue == 1)){
       
      Responses_eff[[j]] = append(Responses_eff[[j]],rbinom(sample_vec[counter], 1, treats_eff[j]))
      Responses_tox[[j]] = append(Responses_tox[[j]],rbinom(sample_vec[counter], 1, treats_tox[j]))
      
      counter = counter + 1
    }
    # Sampling done
    
    ## IA
    if (control){
      # Control beta comparisons
      
      con_n = length(Responses_eff[[1]])  ## Control results
      con_eff = sum(Responses_eff[[1]])
      con_tox = sum(Responses_tox[[1]])
      
      for (j in 2:num_treat){
        # Skipping the control
        
        if (continue[j] == 1){
          # Performs the IA on arms still recreuiting
          
          
          
          
          ## Temp values to do the beta calcs with 
          temp_n = length(Responses_eff[[j]])
          temp_eff = sum(Responses_eff[[j]])
          temp_tox = sum(Responses_tox[[j]])
          
          #print(threshold)
          #print(temp_n)
          #print(temp_eff)
          #print(temp_tox)
          
          threshold = 1 - (lamda * ((IAn[i] / n)^gamma))   # threshold based on how far through 
          
          #browser()
          
          if ( 1 - rcpp_exact_beta( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  | 1 - rcpp_exact_beta( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
            
            #browser()
            
            continue[j] = 0
            
            if (i != num_IA){
              early_stop[j] = 1  # Early stops if rejection happens before end 
              
            }
            final_sample_size[j] = temp_n # Records final sample size
            
          }
          
          
        }
        
      }
      # Finished this round of IA
      if (identical(c(1,rep(0,num_treat - 1)), continue)){
        # If all of the non control arms have finished their runs
        
        if (i != num_IA){
          early_stop[1] = 1 # Early stops if rejection happens before end 
          
        }
        
        
        final_sample_size[1] = temp_n
        
        return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
        
      }
      
      # Loops back
    }
    else{
      # Uncontrolled analysis
      for (j in 1:num_treat){
        # Not skipping a control
        
        if (continue[j] == 1){
          # Performs the IA on arms still recreuiting
          
          temp_n = length(Responses_eff[[j]])
          temp_eff = sum(Responses_eff[[j]])
          temp_tox = sum(Responses_tox[[j]])
          
          
          
          threshold = 1 - (lamda * ((IAn[i] / n)^gamma))   # threshold based on how far through 
          
          #print(threshold)
          #print(temp_n)
          #print(temp_eff)
          #print(temp_tox)
          
          
          
          
          if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
            # If we get a rejection
            
            continue[j] = 0
            
            if (i != num_IA){
              early_stop[j] = 1  # Early stops if rejection happens before end 
              
            }
            final_sample_size[j] = temp_n # Records final sample size
            
          }
        }
        
      }
      # Finished this round of IA
      
      if (identical(c(rep(0,num_treat)), continue)){
        # If all of the arms have finished their runs

        
        return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
        
      }
      
      # Loops back
      
      
    }
    
  }
  
  # If we get here, that means that we've reached the end and there is at least one
  # significant non-control arm
  
  for (i in 1:num_treat){
    # all that needs updating are the final sample sizes 
    if (final_sample_size[i] == 0){
      final_sample_size[i] = length(Responses_eff[[i]])
    }
  }
  
  
  
  return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
  
  
}


# Allo-prob will alternate between the Trippa scheme and the maximum scheme.
allo_prob_brar_V2 = function(successes, toxics, current_ns, N, control = FALSE, trippa = TRUE, phi_e = 0.45, phi_t = 0.3, prior_alpha = rep(phi_e,length(successes)), prior_beta = rep(1 - phi_e,length(successes)),prior_alpha_tox = rep(phi_t,length(successes)), prior_beta_tox = rep(1 - phi_t,length(successes))){
  # Goal: To be given a vector of success, toxic, current sample size per arm, and maximum sample size to get allocation probs
  # Successes is vector of successes
  # toxics is vector of toxics
  # current_ns us a vector of the current sample sizes per arm
  # N is the maxmium possible sample size
  # control determines whether it is a control or not.
  # Trippa determines whether we use Trippa randomisation or not. 
  
  
  
  # Useful for later calculations
  beta_eff = current_ns - successes
  beta_tox = current_ns - toxics
  
  num_arms = length(successes)
  
  if (num_arms == 1){
    return(1)  ## Don't need to do any maths, only one probability makes sense for one arm
  }
  
  if (trippa){
    
    ## Each non control arm works out the probability it is better than
    # predefined values, then gets normalised, then the control gets fixed (but this can be done 
    # outside of the function potentially)
    
    if (control){
      
      trippa_prob = numeric(num_arms)
      
      c = sum(current_ns)/(2*N)  # Tuning factor
      
      for (i in 2:num_arms){
        trippa_prob[i] = (1 - pbeta(phi_e, successes[i] + phi_e, beta_eff[i] + 1 - phi_e ))*( pbeta(phi_t, toxics[i] + phi_t, beta_tox[i] + 1 - phi_t ) )
        
        trippa_prob[i] = trippa_prob[i]^c
        
      }
      
      trippa_prob = trippa_prob/sum(trippa_prob)
      
      trippa_prob[1] = ((exp(max(current_ns[2:num_arms]) - current_ns[1]))/3)^c
      
      trippa_prob = trippa_prob/sum(trippa_prob)
      
    }
    else{
      
      trippa_prob = numeric(num_arms)
      
      c = sum(current_ns)/(2*N)  # Tuning factor
      
      for (i in 1:num_arms){
        trippa_prob[i] = 1 - pbeta(phi_e, successes[i] + phi_e, beta_eff[i] + 1 - phi_e )
      
        trippa_prob[i] = trippa_prob[i]^c
        
      }
      
      trippa_prob = trippa_prob/sum(trippa_prob)
      
    }
    
    return(trippa_prob)
    
  }
  else{
    
    eff_probs = numeric(num_arms)
    tox_probs = numeric(num_arms)
    
    looping_suc = rep(successes, 2) + 1  # so we can loop around and make the following code easier
    looping_tox = rep(toxics,2) + 1 # Beta(1,1) prior
    
    looping_suc_b = rep(beta_eff,2) + 1  # (Beta(1,1) priors to avoid 0 errors)
    looping_tox_b = rep(beta_tox,2) + 1 
    
    if (num_arms == 3){
      # If we're working with 3 arms
      for (i in 1:num_arms){
        eff_probs[i] = rcpp_exact_beta_3_arm(looping_suc[i], looping_suc_b[i],looping_suc[i+1], looping_suc_b[i+1],looping_suc[i+2], looping_suc_b[i+2])
        tox_probs[i] = rcpp_exact_beta_3_arm(looping_tox[i], looping_tox_b[i],looping_tox[i+1], looping_tox_b[i+1],looping_tox[i+2], looping_tox_b[i+2])
      }
      
    
    
    }
    else if (num_arms == 2){
      # If we're working with 2 arms
      for (i in 1:num_arms){
        eff_probs[i] = 1 - rcpp_exact_beta(looping_suc[i], looping_suc_b[i],looping_suc[i+1], looping_suc_b[i+1])
        tox_probs[i] = 1 - rcpp_exact_beta(looping_tox[i], looping_tox_b[i],looping_tox[i+1], looping_tox_b[i+1])
      }
      
    }
    else{
      # If we're working with 4 arms
      
      for (i in 1:num_arms){
        eff_probs[i] = rcpp_exact_beta_4_arm(looping_suc[i], looping_suc_b[i],looping_suc[i+1], looping_suc_b[i+1],looping_suc[i+2], looping_suc_b[i+2],looping_suc[i+3], looping_suc_b[i+3])
        tox_probs[i] = rcpp_exact_beta_4_arm(looping_tox[i], looping_tox_b[i],looping_tox[i+1], looping_tox_b[i+1],looping_tox[i+2], looping_tox_b[i+2],looping_tox[i+3], looping_tox_b[i+3])
      }
      
    }
    
    ## We have the max probabilities, now we combine and return
    
    no_tox_probs = 1 - tox_probs
    
    good_probs = eff_probs * no_tox_probs  # Could weight these differently here
    
    for (i in 1:num_arms){
      
      if (good_probs[i] < 0.01){
        good_probs[i] = 0.01  # Prevents any probability becoming 0, with volatility of max method
      }
    }
    
    good_probs = good_probs / sum(good_probs)
    
    c = sum(current_ns)/(2*N)  # Tuning factor
    
    
    
    final_prob = numeric(num_arms)

    
    ## This can give errors
    for (i in 1:num_arms){
      
      final_prob[i] = good_probs[i] ^ c
    }
    
    final_prob = final_prob / sum(final_prob)
    
    #browser()
    
    return(final_prob)
  }
  
  
  
  
}

Multi_arm_BOP2_brar_V2 = function(treats_eff, treats_tox, n = 180, IAn = c(n), control = FALSE, phi_e = 0.45, phi_t = 0.3, lamda = 0.8, gamma = 0.94, trippa = FALSE, PairedSeed = FALSE){
  ## This uses the skeleton of Multi_arm_BOP2_Max_PA to try and improve Multi_arm_BOP2_V2
  # As Multi_arm_BOP2_Max_PA has a better skeleton, all we need to do is change the allo section
  
  if (PairedSeed != FALSE){  
    set.seed(PairedSeed) # sets seed for paired comparisons
  }
  
  num_treat = length(treats_eff)
  
  num_IA = length(IAn) # Num of IAs
  continue = rep(1,num_treat) # Indicates whether we should continue 
  
  Responses_eff = vector("list",num_treat)
  Responses_tox = vector("list",num_treat)
  
  # OCs that will be filled in later
  early_stop = numeric(num_treat)
  final_sample_size = numeric(num_treat)
  
  sample_order = numeric(0) # NEW, records the sampling order (interest of RAR)
  
  
  for (i in 1:num_IA){
    # Going per IA
    if (i == 1){  # Works out how many to sample per IA
      
      to_sample = IAn[1] # how many will be sampled in total
      
      # denotes what each arm gets min 
      each_arm = to_sample %/% num_treat
      
      # What will be used to choose the samples
      sample_vec = rep(each_arm, num_treat)
      
      if (each_arm != IAn[1]){
        # There are leftovers
        extra = rmultinom(1,to_sample %% num_treat, rep(1/num_treat, num_treat))
        sample_vec = extra + sample_vec
      }
      
      ## Now just randomisaing an order, (As this is PBD, won't grow too apart here)
      
      for (k in 1:num_treat){
        sample_order = append(sample_order, rep(k, sample_vec[k]))
      }
      sample_order = sample(sample_order)
      
      
    }
    else{
      to_sample = IAn[i] - IAn[i-1] # how many will be sampled in total
      
      # This is where the RAR goes, for the BRAR version.
      
      active_arms = sum(continue) # We only care about the remaining arms
      
      # Vectors necessary for allo_prob_brar
      temp_successes = numeric(active_arms)
      temp_tox = numeric(active_arms)
      temp_n = numeric(active_arms)
      
      # Counter is just so we can cycle through the temps to fill
      counter = 0
      
      for (k in which(continue == 1)){
        
        counter = counter + 1
        
        temp_successes[counter] = sum(Responses_eff[[k]])
        temp_tox[counter] = sum(Responses_tox[[k]])
        temp_n[counter] = length(Responses_eff[[k]])
        
        
      }
      
      # ADD CASE FOR WHEN 2 ARM 
      
      allo_prob = allo_prob_brar_V2(temp_successes, temp_tox, temp_n, control = control, trippa = trippa, N = n)
      # denotes what each arm gets min 
      
      
      ## PBD, if want to go fully coin based randomisation this is the line you change
      each_arm = floor(to_sample * allo_prob)  
      
      # What will be used to choose the samples
      sample_vec = each_arm
      
      if (sum(each_arm) != to_sample){
        # There are leftovers
        extra = rmultinom(1,to_sample - sum(each_arm), allo_prob)
        sample_vec = extra + sample_vec
      }
      
      ## Now just randomisaing an order, (As this is PBD, won't grow too apart here)
      
      new_sample_order = numeric(0)
      
      for (k in 1:active_arms){
        new_sample_order = append(new_sample_order, rep(which(continue == 1)[k], sample_vec[k]))
      }
      new_sample_order = sample(new_sample_order)
      
      sample_order = append(sample_order, new_sample_order)
      
      #browser()
      
    }
    
    counter = 1  # Need a counter to accurately reference sample_Vec
    for (j in which(continue == 1)){
      
      Responses_eff[[j]] = append(Responses_eff[[j]],rbinom(sample_vec[counter], 1, treats_eff[j]))
      Responses_tox[[j]] = append(Responses_tox[[j]],rbinom(sample_vec[counter], 1, treats_tox[j]))
      
      counter = counter + 1
    }
    # Sampling done
    
    ## IA
    if (control){
      # Control beta comparisons
      
      con_n = length(Responses_eff[[1]])  ## Control results
      con_eff = sum(Responses_eff[[1]])
      con_tox = sum(Responses_tox[[1]])
      
      for (j in 2:num_treat){
        # Skipping the control
        
        if (continue[j] == 1){
          # Performs the IA on arms still recreuiting
          
          #browser()
          
          
          threshold = 1 - (lamda * ((IAn[i] / n)^gamma))   # threshold based on how far through 
          
          ## Temp values to do the beta calcs with 
          temp_n = length(Responses_eff[[j]])
          temp_eff = sum(Responses_eff[[j]])
          temp_tox = sum(Responses_tox[[j]])
          
          #print(threshold)
          #print(temp_n)
          #print(temp_eff)
          #print(temp_tox)
          
          #browser()
          
          
          if ( 1 - rcpp_exact_beta( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  | 1 - rcpp_exact_beta( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
            
            
            continue[j] = 0
            
            #browser()
            
            if (i != num_IA){
              early_stop[j] = 1  # Early stops if rejection happens before end 
              
            }
            final_sample_size[j] = temp_n # Records final sample size
            
          }
          
          
        }
        
      }
      # Finished this round of IA
      if (identical(c(1,rep(0,num_treat - 1)), continue)){
        # If all of the non control arms have finished their runs
        
        if (i != num_IA){
          early_stop[1] = 1 # Early stops if rejection happens before end 
          
        }
        
        
        final_sample_size[1] = temp_n
        
        return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
        
      }
      
      # Loops back
    }
    else{
      # Uncontrolled analysis
      for (j in 1:num_treat){
        # Not skipping a control
        
        if (continue[j] == 1){
          # Performs the IA on arms still recreuiting
          
          temp_n = length(Responses_eff[[j]])
          temp_eff = sum(Responses_eff[[j]])
          temp_tox = sum(Responses_tox[[j]])
          
          
          
          threshold = 1 - (lamda * ((IAn[i] / n)^gamma))   # threshold based on how far through 
          
          #print(threshold)
          #print(temp_n)
          #print(temp_eff)
          #print(temp_tox)
          
          
          
          
          if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
            # If we get a rejection
            
            continue[j] = 0
            
            if (i != num_IA){
              early_stop[j] = 1  # Early stops if rejection happens before end 
              
            }
            final_sample_size[j] = temp_n # Records final sample size
            
          }
        }
        
      }
      # Finished this round of IA
      
      if (identical(c(rep(0,num_treat)), continue)){
        # If all of the arms have finished their runs
        
        
        return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
        
      }
      
      # Loops back
      
      
    }
    
  }
  
  # If we get here, that means that we've reached the end and there is at least one
  # significant non-control arm
  
  for (i in 1:num_treat){
    # all that needs updating are the final sample sizes 
    if (final_sample_size[i] == 0){
      final_sample_size[i] = length(Responses_eff[[i]])
    }
  }
  
  
  
  return(list( early_stop, final_sample_size, continue ))  # So ends the trial here
  
  
}

OC_gen = function(treats_eff, treats_tox, n = 180, IAn = c(n), lambda = 0.9, gamma = 0.9, RAR = FALSE, trippa = FALSE, n_sim = 10000, control = TRUE, max_pa = FALSE, PairedSeed = FALSE){
  # Generate OCs from above functions
  # Want ESS, power, Early stopping rate
  
  n_arms = length(treats_eff)
  
  Power = numeric(n_arms)
  ESS = numeric(n_arms)
  Early_Stop = numeric(n_arms)
  FWER = 0
  
  
  
  
  
  if (RAR){
    for (i in 1:n_sim){
      
      temp_data = Multi_arm_BOP2_brar_V2(treats_eff = treats_eff, treats_tox = treats_tox, n = n, IAn = IAn, control = control, trippa = trippa, lamda = lambda, gamma = gamma, PairedSeed = PairedSeed)
      
      Power = Power + temp_data[[3]]/n_sim
      ESS = ESS + temp_data[[2]]/n_sim
      Early_Stop = Early_Stop + temp_data[[1]]/n_sim
      
      if (control){
        temp_data[[c(3,1)]] = 0 # Makes the control 0 to not screw over FWER readings 
      }
      
      FWER = FWER + max(temp_data[[3]])/n_sim
      
      if (PairedSeed != FALSE){  
        PairedSeed = PairedSeed + 1 # This makes sure that every sim of the trial uses the same seed, same data generation
          # whilst also making sure that the data points are still i.i.d to one another(as different seed)
      }
      
      
    }
    
    
  }
  else{
    
    for (i in 1:n_sim){
        
      if (max_pa){
        temp_data = Multi_arm_BOP2_Max_PA(treats_eff = treats_eff, treats_tox = treats_tox, n_pa = n, IAn_pa = IAn, control = control, lamda = lambda, gamma = gamma, PairedSeed = PairedSeed)
      }else{
        temp_data = Multi_arm_BOP2_V3(treats_eff = treats_eff, treats_tox = treats_tox, n = n, IAn = IAn, control = control, lamda = lambda, gamma = gamma, PairedSeed = PairedSeed)
      }
      
      
      Power = Power + temp_data[[3]]/n_sim
      ESS = ESS + temp_data[[2]]/n_sim
      Early_Stop = Early_Stop + temp_data[[1]]/n_sim
      
      if (control){
        temp_data[[c(3,1)]] = 0 # Makes the control 0 to not screw over FWER readings 
      }
      
      
      FWER = FWER + max(temp_data[[3]])/n_sim
      
      #browser()
      
      if (PairedSeed != FALSE){  
        PairedSeed = PairedSeed + 1 # This makes sure that every sim of the trial uses the same seed, same data generation
        # whilst also making sure that the data points are still i.i.d to one another(as different seed)
      }
      
    }
    
    
  }
  
  return(list(Early_Stop,ESS,Power, FWER))
  
}

grid_search_flex_multi = function(control_eff = c(0.45,0.45,0.45), control_tox = c(0.3,0.3,0.3), treats_eff = c(0.45,0.45,0.6), treats_tox = c(0.3,0.3,0.2), n = 180, IAn = c(n), control = FALSE, RAR = FALSE, trippa = FALSE, n_sim = 10000, target_FWER = 0.1, max_pa = FALSE){
  # Grid searches uses OC_gen. Will use parallelisation. 
  # Control_eff and control_tox will be used to control the FWER
  
  pot_lamda = seq(0.7 + 0.1*RAR,1,0.01 + 0.01*RAR) # Can vary the parameters to try and get the right kinda area, going off table S9 from supp material
  pot_gamma = seq(0.7 + 0.1*RAR,1,0.01 + 0.01*RAR)
  
  temp_len = length(pot_lamda)  # Before the length gets updated
  
  # For easy looping in Foreach 
  pot_lamda = rep(pot_lamda, length(pot_gamma))
  pot_gamma = rep(pot_gamma, each = temp_len)
  
  par_power = numeric(length(pot_gamma))
  
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl) 
  
  par_power = foreach( i = 1:length(pot_gamma), .combine = c, .export = c("OC_gen","Multi_arm_BOP2_brar_V2","allo_prob_brar_V2","Multi_arm_BOP2_V3", "Multi_arm_BOP2_Max_PA", "rcpp_exact_beta_4_arm","rcpp_exact_beta_3_arm","rcpp_exact_beta")) %dopar%{
    
    # Grid searching power
    if (OC_gen(control_eff, control_tox, n = n, IAn = IAn, lambda = pot_lamda[i], gamma = pot_gamma[i], RAR = RAR, trippa = trippa, control = control, n_sim = n_sim, max_pa = max_pa)[[4]] < target_FWER){
      return(OC_gen(treats_eff, treats_tox, n = n, IAn = IAn, lambda = pot_lamda[i], gamma = pot_gamma[i], RAR = RAR, trippa = trippa, control = control, n_sim = n_sim, max_pa = max_pa)[[c(3, length(treats_eff))]])
    }
    else{
      return(0)
    }
    
    
  }
  
  stopCluster(cl)
  
  max_index = which.max(par_power)
  
  results = c("Lambda" = pot_lamda[max_index], "Gamma" = pot_gamma[max_index], "Power" = par_power[max_index])
  
  print(results)
  
  return(results)
  
  
}

## RESULTS
# RAR FALSE, no IA, control = FALSE, 0.83 lambda, 0.71 gamma, 0.74 power
# RAR TRUE, control = FALSE, Trippa = FALSE c(60,120,180), 0.8 lambda, 1 gamma, 0.84 power
# RAR TRUE, control = FALSE, Trippa = TRUE, c(60,120,180), 0.8 lambda, 0.94 gamma, 0.83 power
# RAR FALSE, no IA, control = TRUE, 0.77 lambda, 0.86 gamma, 0.57 power
# RAR TRUE, IA c(60,120,180), Trippa = FALSE, 0.8 lambda, 1 gamma, 0.53 power
# RAR TRUE, IA c(60,120,180), Trippa = TRUE, 0.8 lambda, 1 gamma, 0.526 power

allo_prob_plot = function(control = FALSE){
  # Will plot how the allocation probability changes as efficacy increases versus the two different schemes
  
  trippa_data = numeric(31)
  trippa_data_control= numeric(31)
  
  max_data = numeric(31)
  max_data_control = numeric(31)
  
  for (i in 0:30){
    trippa_temp = allo_prob_brar_V2(c(13,13,i),c(9,9,9), current_ns = c(30,30,30), control = control, trippa = TRUE, N = 180)
    
    trippa_data[i+1] = trippa_temp[3]
    trippa_data_control[i+1] = trippa_temp[1]
    
    max_temp = allo_prob_brar_V2(c(13,13,i),c(9,9,9), current_ns = c(30,30,30), control = control, trippa = FALSE, N = 180)
    
    max_data[i+1] = max_temp[3]
    max_data_control[i+1] = max_temp[1]
    
  }
  
  par(mfrow = c(2,1))
  
  plot(x = 0:30, y = trippa_data, main = "Trippa scheme allocation probability plot.", col = "blue", ylim = c(0,1), ylab = "Probability", xlab = "Number of successes in experimental treatment")
  points(x = 0:30, y = trippa_data_control, col = "red", pch = 2)
  lines(x = rep(13,11), y = 0:10/10, col = "green")
  legend("topleft",legend = c("Control","Treatment"),pch = c(2,1), col = c("red","blue"))
  
  plot(x = 0:30, y = max_data, col = "blue", main = "Max scheme allocation probability plot.", ylim = c(0,1), ylab = "Probability", xlab = "Number of successes in experimental treatment")
  points(x = 0:30, y = max_data_control, col = "red", pch = 2)
  lines(x = rep(13,11), y = 0:10/10, col = "green")
  legend("topleft",legend = c("Control","Treatment"),pch = c(2,1), col = c("red","blue"))
  
  #browser()
  
  par(mfrow = c(1,1))
  
}

One_ESS_plot_multi = function(treats_eff, treats_tox, n = 180, IAn = c(n), lambda = 0.83, gamma = 0.71, RAR = FALSE, trippa = FALSE, n_sim = 10000, control = FALSE, max_pa = FALSE, PairedSeed = FALSE, FWER = FALSE){
  # Will plot ESS for multi-arm cases.
  # Will just see how the ESS changes for each arm as you change the IA placement
  
  num_arms = length(treats_eff)
  
  iter = n %/% num_arms
  remain = n %% num_arms
  
  IA_list = vector("list", iter ) 
  Power_vec = numeric(iter)
  FWER_vec = numeric(iter)
  
  # To get the IAs
  for (i in 1:iter ){
    temp_data = OC_gen(treats_eff, treats_tox, n, IAn = c(remain + i*num_arms, n), lambda, gamma, RAR, trippa, n_sim, control, max_pa, PairedSeed = PairedSeed)

    IA_list[[i]] = temp_data[[2]]
    Power_vec[i] = temp_data[[c(3,3)]]  # only gives power when last element
    FWER_vec[i] = temp_data[[4]]
  }
  
  ESS_plot = numeric(iter)
  
  # First IA
  for (i in 1:iter){
    ESS_plot[i] = IA_list[[c(i,1)]]
  }
  
  #browser()
  
  if (TRUE){
    
    
    for (j in 2:num_arms){
      ## Plotting the other ESS, have to extract them from the list first
      
      for (i in 1:iter){
        
        ESS_plot = append(ESS_plot,IA_list[[c(i,j)]])
      }
    }
    
    #browser()
  
    plot_df = data.frame(ESS_plot,Power_plot = rep(Power_vec,num_arms), IA_placement = rep(seq(remain + num_arms,n,num_arms),num_arms), treat = factor(rep(1:num_arms,each = length(Power_vec))))
    
    
    ggplot_ess = ggplot(plot_df, aes(y = ESS_plot, x = IA_placement)) +
      geom_line(aes(y = ESS_plot, colour = treat, linetype = treat), linewidth = 1.5) +
      labs( title = "ESS per arm versus Interim placement - Multi arm.", x = "IA placement", y = "ESS per arm", colour = "Arms", linetype = "Arms") +
      scale_color_viridis(discrete=TRUE, option="viridis")
    
    
    print(ggplot_ess)
    
    ggplot_power = ggplot(plot_df, aes(y = Power_plot, x = IA_placement)) +
      geom_line(aes(y = Power_plot)) +
      labs( title = "Least power versus Interim placement - Multi arm", x = "IA placement", y = "Least power")
    print(ggplot_power)
    
    to_print_plot = ggarrange(ggplot_ess, ggplot_power,
              ncol = 1, nrow = 2, common.legend = TRUE)
    
    print(to_print_plot)
    
    if (FWER){
      plot_df = data.frame(ESS_plot,FWER_plot = rep(FWER_vec,num_arms), IA_placement = rep(seq(remain + num_arms,n,num_arms),num_arms), treat = factor(rep(1:num_arms,each = length(Power_vec))))
      
      ggplot_power = ggplot(plot_df, aes(y = FWER_plot, x = IA_placement)) +
        geom_line(aes(y = FWER_plot)) +
        labs( title = "FWER versus Interim placement - Multi arm", x = "IA placement", y = "FWER")
      print(ggplot_power)
      
    }
  
  }
  
  if (FALSE){ # regular R plot
    
    
    plot(y = Power_vec, x = seq(remain + num_arms,n,num_arms), main = "Least Power versus Interim placement", xlab = "IA placement", ylab = "Least Power", type = "l")
    
    plot(y = ESS_plot, x = seq(remain + num_arms,n,num_arms), main = "ESS versus Interim placement", xlab = "IA placement", ylab = "ESS", type = "l", ylim = c(30,100), col = 1)
    for (j in 2:num_arms){
      ## Plotting the other ESS, have to extract them from the list first
      
      for (i in 1:iter){
        
        ESS_plot[i] = IA_list[[c(i,j)]]
      }
      
      lines(x = seq(remain + num_arms,n,num_arms), y = ESS_plot, col = j)
      
    }
  
    
    #Legend
    legend("topright",
           legend = treats_eff,
           col = 1:num_arms,
           bty = "n",
           lty = c(1,1,1),
           pt.cex = 2,
           cex = 1.2,
           text.col = "black",
           horiz = F ,
           inset = c(0.1, 0.1))
  
  }
}



## Reproducible plots

## Figure 7

Figure_7 = FALSE

if (Figure_7){
  set.seed(2025)
  
  One_ESS_plot_multi(c(0.45,0.45,0.6),c(0.3,0.3,0.2), PairedSeed = 2025)
  
}

Appendix_for_Figure_7 = FALSE
if (Appendix_for_Figure_7){
  set.seed(2025)
  
  One_ESS_plot_multi(c(0.45,0.45,0.45),c(0.3,0.3,0.3), PairedSeed = 2025, FWER = TRUE)
  
}


## Table 3, Table 4, Table 5, Table 6 are each in the multi arm setting

Table_3 = FALSE  # DONE
Table_4 = FALSE # DONE
Table_5 = FALSE # DONE
Table_6 = FALSE # DONE


## All lambda and gamma below obtained via grid search at type one error 0.1
if (Table_3){
  set.seed(2025)
  
  # Mulier
  print(OC_gen(c(0.45,0.45,0.45,0.6),c(0.3,0.3,0.3,0.2), n = 60, IAn = c(15,30,45,60), control = TRUE, RAR = FALSE, trippa = FALSE, max_pa = TRUE, lambda = 0.78, gamma = 0.92  ))
  
  # Reallo 
  print(OC_gen(c(0.45,0.45,0.45,0.6),c(0.3,0.3,0.3,0.2), n = 240, IAn = c(60,120,180,240), control = TRUE, RAR = FALSE, trippa = FALSE, max_pa = FALSE, lambda = 0.76, gamma = 0.82 ))
  
  # Trippa
  print(OC_gen(c(0.45,0.45,0.45,0.6),c(0.3,0.3,0.3,0.2), n = 240, IAn = c(60,120,180,240), control = TRUE, RAR = TRUE, trippa = TRUE, max_pa = FALSE, lambda = 0.8, gamma = 1  ))
  
  # Max
  print(OC_gen(c(0.45,0.45,0.45,0.6),c(0.3,0.3,0.3,0.2), n = 240, IAn = c(60,120,180,240), control = TRUE, RAR = TRUE, trippa = FALSE, max_pa = FALSE, lambda = 0.8, gamma = 1  ))
  
}
if (Table_4){
  set.seed(2025)
  
  # Mulier
  print(OC_gen(c(0.45,0.45,0.6),c(0.3,0.3,0.2), n = 60, IAn = c(20,40,60), control = FALSE, RAR = FALSE, trippa = FALSE, max_pa = TRUE, lambda = 0.78, gamma = 0.98  ))
  
  # Reallo 
  print(OC_gen(c(0.45,0.45,0.6),c(0.3,0.3,0.2), n = 180, IAn = c(60,120,180), control = FALSE, RAR = FALSE, trippa = FALSE, max_pa = FALSE, lambda = 0.79, gamma = 0.9))
  
  # Trippa
  print(OC_gen(c(0.45,0.45,0.6),c(0.3,0.3,0.2), n = 180, IAn = c(60,120,180), control = FALSE, RAR = TRUE, trippa = TRUE, max_pa = FALSE, lambda = 0.8, gamma =1))
  
  # Max
  print(OC_gen(c(0.45,0.45,0.6),c(0.3,0.3,0.2), n = 180, IAn = c(60,120,180), control = FALSE, RAR = TRUE, trippa = FALSE, max_pa = FALSE, lambda = 0.8, gamma = 0.96))
}
if (Table_5){
  set.seed(2025)
  
  # Mulier
  print(OC_gen(c(0.45,0.55,0.6,0.65),c(0.3,0.25,0.2,0.15), n = 60, IAn = c(15,30,45,60), control = TRUE, RAR = FALSE, trippa = FALSE, max_pa = TRUE, lambda = 0.78, gamma = 0.92 ))
  
  # Reallo 
  print(OC_gen(c(0.45,0.55,0.6,0.65),c(0.3,0.25,0.2,0.15), n = 240, IAn = c(60,120,180,240), control = TRUE, RAR = FALSE, trippa = FALSE, max_pa = FALSE, lambda = 0.76, gamma = 0.82))
  
  # Trippa
  print(OC_gen(c(0.45,0.55,0.6,0.65),c(0.3,0.25,0.2,0.15), n = 240, IAn = c(60,120,180,240), control = TRUE, RAR = TRUE, trippa = TRUE, max_pa = FALSE, lambda = 0.8, gamma = 1 ))
  
  # Max
  print(OC_gen(c(0.45,0.55,0.6,0.65),c(0.3,0.25,0.2,0.15), n = 240, IAn = c(60,120,180,240), control = TRUE, RAR = TRUE, trippa = FALSE, max_pa = FALSE, lambda = 0.8, gamma = 1 ))
}
if (Table_6){
  set.seed(2025)
  
  # Mulier
  print(OC_gen(c(0.55,0.6,0.65),c(0.25,0.2,0.15), n = 60, IAn = c(20,40,60), control = FALSE, RAR = FALSE, trippa = FALSE, max_pa = TRUE, lambda = 0.78, gamma = 0.98 ))
  
  # Reallo 
  print(OC_gen(c(0.55,0.6,0.65),c(0.25,0.2,0.15), n = 180, IAn = c(60,120,180), control = FALSE, RAR = FALSE, trippa = FALSE, max_pa = FALSE, lambda = 0.79, gamma = 0.9))
  
  # Trippa
  print(OC_gen(c(0.55,0.6,0.65),c(0.25,0.2,0.15), n = 180, IAn = c(60,120,180), control = FALSE, RAR = TRUE, trippa = TRUE, max_pa = FALSE, lambda = 0.8, gamma =1))
  
  # Max
  print(OC_gen(c(0.55,0.6,0.65),c(0.25,0.2,0.15), n = 180, IAn = c(60,120,180), control = FALSE, RAR = TRUE, trippa = FALSE, max_pa = FALSE, lambda = 0.8, gamma = 0.96))
}


## Figure 12 and 13 aren't stochastic, just by calling allo_prob_plot
Figure_12_13 = TRUE
if(Figure_12_13){
  
  
  allo_prob_plot(control = TRUE)
  
}

## Figure 15

Figure_15 = FALSE

if (Figure_15){
  set.seed(2025)
  
  One_ESS_plot_multi(c(0.6,0.45,0.45),c(0.2,0.3,0.3), RAR = TRUE, trippa = FALSE)
  
}


## Figure 16

Figure_16 = FALSE

if (Figure_16){
  set.seed(2025)
  
  One_ESS_plot_multi(c(0.6,0.45,0.45),c(0.2,0.3,0.3), RAR = TRUE, trippa = TRUE)
  
}

Appendix_2 = FALSE
if(Appendix_2){
  set.seed(2025)
  
  One_ESS_plot_multi(c(0.45,0.45,0.6),c(0.3,0.3,0.2), lambda = 0.83, gamma = 0.71, RAR = FALSE)
  
}
