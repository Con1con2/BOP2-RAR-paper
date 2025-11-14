## BOP2 2 report reproducible
# Leave the booleans as they are and run the script for reproducible results.
# Edit the booleans and arguments to tailor them as you need
# updated

library(tidyverse) # For plotting
library(ggplot2) # For RAR (And other) plots
library(ggthemes) # For colours
library(knitr) # For tables
library(kableExtra) # For tables
library(doParallel)  # Parallelisation
library(doRNG) #  Reproducible in parallel
library(foreach)  # Parallelisation
library(rBeta2009) # For quicker beta 
#library(webshot2) # For saving the kables
# library(microbenchmark) # Not needed to run, but useful to compare speed of different functions


## These booleans determine whether you want to optimise the parameters again
# Both are two arms are optimised for H0 0.2, H1 0.4

# Optimise the Balanced Randomisation (Standard BOP2) procedure two arm
Opt_BR = FALSE

# Optimise the Bayesian Response Adaptive Randomisation procedure two arm
Opt_BRAR = FALSE


# Both Multi arms are optimised for H0: 0.45,0.3 three arm, H1: one arm 0.6,0.2 two arm 0.45,0.3
# Optimsies BOP2 multi arm parameters
Opt_BR_Multi_Arm = FALSE

# Optimised BOP2 Multi Arm BRAR parameters
Opt_BRAR_Multi_Arm = FALSE

## These figure variables represent whether that figure will be produced or not

Fig_2_4 = FALSE

Fig_5 = FALSE

Fig_6_7 = TRUE

Fig_8 = FALSE

Fig_9 = FALSE

Fig_10 = FALSE


## Figures made for the 2nd draft of the report


New_Figure_8 = FALSE

Appendix_1 =  FALSE

Appendix_2 = FALSE





## Base functions for BOP2 2 ARM

Beta_comp = function(A1,B1, A2, B2,n = 1000){
  ## Beta_comp is a function finds the probability than one beta is larger than another beta via simulation (empirical c.d.f). 
  # A1 and B1 are the alpha and beta parameters of first Beta that you want to see is larger
  # A2 and B2 are the same for the second beta
  # n is the number of trials to run
  
  Beta1 = rbeta(n,A1,B1)
  Beta2 = rbeta(n,A2,B2)  # Simulates the different betas
  
  p = mean(Beta1>Beta2) # Empirical c.d.f. calculation
  
  return(p)
}

Beta_comp_exact = function(A1,B1, A2, B2){
  ## Beta_comp_exact does an exact Beta computation comparing two beta variables when the prior is an uniformative 1,1
  ## If the parameters are not integers, this will not work. 
  ## Derived on https://www.evanmiller.org/bayesian-ab-testing.html
  
  # In testing, using this function instead of beta_comp performed analysis ~ 10 times faster.
  
  # A1, B1 are the alpha and beta parameters for the beta we want to see is higher. A2 and B2 are the parameters for the other beta
  
  if( identical(as.integer(c(A1,B1, A2, B2)), c(A1,B1, A2, B2))){  # We need integer priors
    print("Error: not integer priors")
    return()
  }
  
  prob = 0
  
  for (i in 0:(A2-1)){
    prob = prob + beta(A1 + i, B1 + B2)/((B2 + i)*(beta(i+1,B2)*beta(A1,B1)))
  }
  
  return(1- prob)  # This returns the probability that the A1,B1 beta variable is bigger than the A2,B2 beta variable
}

Threshold_prob = function(lamda, gamma, n_total, n_current){
  ## Function responsible for the threshold probabilities for early stopping
  # gamma and lamda are the same variables in the BOP2 paper
  # n_total is the total sample size (of the arm)
  # n_current is the current sample size (of the arm)
  
  
  Futility = lamda * (n_current/n_total)^gamma
  Super = 2 * pnorm ( (qnorm ( (1 + lamda)/2))/sqrt(n_current/n_total)) - 1
  return(c(Futility, Super))
}

## NOTE! n is DIFFERENT for Sim_trial and Sim_trial_RAR. 
# In Sim_trial, n (and IAn) represent the sample size PER ARM
# In Sim_trial_RAR and after, n (and IAn) represents the TOTAL sample size
Sim_trial = function(control, treat, n, IAn, prior_a = 1, prior_b = 1, gamma = 0.92, lamda = 0.9, Norm_approx = FALSE, test_type = "BOP2_FE"){
  ## Sim_trial runs a 2 arm BOP2 design. Starting off with a simple Binomial endpoint, uninformative prior.
  # control is the ORR of the control treatment
  # treat is the ORR of the target treatment =
  # n is the total sample size per arm
  # IAn is the number of samples needed before an Interim Analysis (carried out at equal intervals)
  # RAR determines whether to do Response Adaptive Randomisation or not. RAR = False means that assignment will be random in a 1:1 ratio.
  # prior_a and prior_b are the priors of the beta distribution, starting with an uniformative 0.5,0.5
  # gamma and lamda are variables that appear in the threshold probabilities
  # Norm approx determines whether we use a normal approximation for the Beta calculations
  # test_type determines how we test our interim analysis, BOP2 is the intended BOP2 deisgn, OF, Pocock, and GF are the other boundaries to test again. _FE at the end stands for futility effifacy, as opposed to just futility.
  
  
  # We take a beta(0.5,0.5) as our uniformative prior, and binomial with p for our data | parameter p
  
  data_control = numeric(n)
  data_treat = numeric(n) # Creating our empty data sets
  
  Num_IA = 1 # records how many IAs have been carried out
  
  while(Num_IA <= (n%/%IAn + 1) ){
    
    
    if ((n - Num_IA*IAn) >= 0) {  # If our sample size allows us to simulate n more
      
      data_control[(1+ (Num_IA-1)*IAn):(Num_IA*IAn)] = rbinom(IAn, 1, control)
      data_treat[(1+ (Num_IA-1)*IAn):(Num_IA*IAn)] = rbinom(IAn, 1, treat)  # Indexing looks complicated, but is just making sure the new data goes in the right spot
      
      
    }
    else{
      
      if (n%%IAn != 0){ # if there is anything more to be simulated
        data_control[ (n - n%%IAn):n] = rbinom(n%%IAn, 1, control)  
        data_treat[ (n - n%%IAn):n] = rbinom(n%%IAn, 1, treat) # Simulating the remaining slots
      }
      ## Final decision rule
      # If efficate
      # return true
      # else
      # return false
      
      
      A1 = prior_a + sum(data_treat)  # Updating the posterior
      B1 = prior_b + (n - sum(data_treat))
      
      A2 = prior_a + sum(data_control)  # Updating the posterior
      B2 = prior_b + (n - sum(data_control))
      
      if (Norm_approx == FALSE){
        
        
        #IA_test = Beta_comp(A1,B1,A2,B2,100000)
        
        IA_test = Beta_comp_exact(A1,B1,A2,B2)
        
      }
      else{   # This normal approximation is primarly used for computational efficiency, such as the grid search.
        Mu_mean = A2/(A2+B2) - A1/(A1+B1)  # This way aroud because pnorm works out the lower tail, so the probability that the second variable is bigger than the first
        Sigma_var = A1*B1/((A1+B1)^2*(A1+B1+1)) + A2*B2/((A2+B2)^2*(A2+B2+1))
        IA_test = pnorm(0, Mu_mean, sqrt(Sigma_var))
      }
      
      t = 1 # Useful for calculations below 
      a = 0.1 # BOP2's chosen type 1 error, necessary for calculating thresholds
      
      if (test_type == "BOP2_FE"){
        ## BOP2 double threshold
        Probs = Threshold_prob(lamda,gamma,n,n)  # Vector of threshold probabilities, [1] is futility, [2] is super. Same when trial is finished. 
      }
      else if (test_type == "BOP2"){
        #BOP2 futility threshold
        
        Probs = Threshold_prob(lamda,gamma,n,n)
        # As this is the final test, there is a "success" threshold
      }
      else if (test_type == "OF_FE"){
        # OF double threshold
        
        temp_p = 2 - 2*pnorm( qnorm((1-a/2))/sqrt(t))
        
        Probs = c(temp_p, 1-temp_p)
        
      }
      else if (test_type == "OF"){
        # OF futility threshold
        
        temp_p = 2 - 2*pnorm( qnorm((1-a/2))/sqrt(t))
        
        Probs = c(temp_p, 1-temp_p)  # Both superior and futility because of final test
      }
      else if (test_type == "Pocock_FE"){
        # Pocock double threshold
        
        temp_p = a * log( 1 + (exp(1)-1)*t)
        
        Probs = c(temp_p, 1-temp_p) 
      } 
      else if (test_type == "Pocock"){
        # Pocock futility threshold
        
        temp_p = a * log( 1 + (exp(1)-1)*t)
        
        Probs = c(temp_p, 1-temp_p) # Both superior and futility because of final test
      }
      else if (test_type == "GF_FE"){
        # GF double threshold
        
        Probs = c(a * (1- exp(2*t))/(1-exp(2)), 1 - a * (1- exp(4*t))/(1-exp(4)))
        
        
      }
      else if (test_type == "GF"){
        # GF futility threshold
        
        
        Probs = c(a * (1- exp(2*t))/(1-exp(2)), 1 - a * (1- exp(4*t))/(1-exp(4)))
        
        
      }
      
      #print(Probs)
      #print(IA_test)
      
      if (IA_test > Probs[2]){  # Threshold can be changed to be more adaptive
        return(c("Reject Null" = TRUE, "Sample size" = n*2))  # Final supereority check, n*2 for total sample size (as n per arm)
      } else {
        return(c("Reject Null" = FALSE, "Sample size" = n*2)) # Didn't pass a superoirty check so futile.
      }
    }
    
    # IA decision rule 
    
    # if futility
    # return false
    # else if efficate
    # return true
    # else
    # Num_IA = Num_IA + 1
    
    
    
    A1 = prior_a + sum(data_treat)  # Updating the posterior
    B1 = prior_b + (Num_IA*IAn - sum(data_treat))
    
    A2 = prior_a + sum(data_control)  # Updating the posterior
    B2 = prior_b + (Num_IA*IAn - sum(data_control))
    
    # print(data_treat)
    # print(data_control) # testing to see data set
    
    # IA_test = Beta_comp(A1,B1,A2,B2,100000) # test statistic for IA analysis
    
    IA_test = Beta_comp_exact(A1,B1,A2,B2) # test statistic for IA analysis
    
    t = Num_IA*IAn/n # Useful for calculations below 
    a = 0.1 # BOP2's chosen type 1 error, necessary for calculating thresholds
    
    if (test_type == "BOP2_FE"){
      ## BOP2 double threshold
      Probs = Threshold_prob(lamda,gamma,n,Num_IA*IAn)  # Vector of threshold probabilities, [1] is futility, [2] is super. Same when trial is finished. 
    }
    else if (test_type == "BOP2"){
      #BOP2 futility threshold
      
      Probs = Threshold_prob(lamda,gamma,n,Num_IA*IAn)
      Probs[2] = 1
      # No efficacy threshold
    }
    else if (test_type == "OF_FE"){
      # OF double threshold
      
      temp_p = 2 - 2*pnorm( qnorm((1-a/2))/sqrt(t))
      
      Probs = c(temp_p, 1-temp_p)
      
    }
    else if (test_type == "OF"){
      # OF futility threshold
      
      temp_p = 2 - 2*pnorm( qnorm((1-a/2))/sqrt(t))
      
      Probs = c(temp_p, 1)  # No efficacy threshold
    }
    else if (test_type == "Pocock_FE"){
      # Pocock double threshold
      
      temp_p = a * log( 1 + (exp(1)-1)*t)
      
      Probs = c(temp_p, 1-temp_p) 
    } 
    else if (test_type == "Pocock"){
      # Pocock futility threshold
      
      temp_p = a * log( 1 + (exp(1)-1)*t)
      
      Probs = c(temp_p, 1) # No efficacy threshold
    }
    else if (test_type == "GF_FE"){
      # GF double threshold
      
      Probs = c(a * (1- exp(2*t))/(1-exp(2)), 1 - a * (1- exp(4*t))/(1-exp(4)))
      
      
    }
    else if (test_type == "GF"){
      # GF futility threshold
      
      
      Probs = c(a * (1- exp(2*t))/(1-exp(2)), 1) # No efficacy threshold
      
      
    } 
    
    #print(IA_test) #testing to see proba
    #print(Probs) #testing to see thresholds
    
    
    if (IA_test > Probs[2]){  # Threshold can be changed to be more adaptive
      return(c("Reject Null" = TRUE, "Sample size" = Num_IA*IAn*2))  # STOPPING EARLY FOR SUPERIORITY *2 for total sample size (as n per arm)
    } else if (IA_test < Probs[1]){
      return(c("Reject Null" = FALSE, "Sample size" = Num_IA*IAn*2)) # STOPPING EARLY FOR FUTILITY
    } else{
      Num_IA = Num_IA + 1 # ONTO NEXT IA
    }
    # print(Num_IA)
    # print(n%%IAn + 1) testing to check it's looping
  }
}

Sim_trial_RAR = function(control, treat, n, IAn, burn = IAn, prior_a = 1, prior_b = 1, gamma = 0.88, lamda = 0.96, Norm_approx = FALSE, print_plots = TRUE, return_power = FALSE){
  
  ## Sim_trial_RAR runs a 2 arm BOP2 design with RAR Thompson Sampling. Starting off with a simple Binomial endpoint, uninformative prior.
  # control is the ORR of the control treatment
  # treat is the ORR of the target treatment =
  # n is the total sample size
  # IAn is the number of samples needed before an Interim Analysis (carried out at equal intervals)
  # burn is how many samples we should do 1:1 equal randomisation before implemening RAR. Note that the Interim Analysis will not begin until the burn has finished.Burn needs to be at least 1 to avoid errors
  # prior_a and prior_b are the priors of the beta distribution, starting with an uniformative 0.5,0.5
  # gamma and lamda are variables that appear in the threshold probabilities
  # Norm approx determines whether we use a normal approximation for the Beta calculations
  
  
  # We take a beta(0.5,0.5) as our uniformative prior, and binomial with p for our data | parameter p
  
  # Slightly differently from Sim_trial, n is now the total sample size, rather than the size of each arm. Likewise, after taking IAn samples we will perform an
  # interim analysis, regardless of how the treatments were assigned. Finally, rather than some kinda of fixed block random assignment we had in the Sim_trial ( in every
  # 20 samples there would be 10 to control and 10 to treatment), this will instead use a 1:1 ratio of equal assignment before RAR begins. This is because this
  # segues better into RAR 
  
  plot_prob = rep(0.5, burn) # The assignment probability we will keep track of (for plotting later. It is already known that it is 0.5 for the burn period.
  
  count = 0
  C_count = numeric(0)
  T_count = numeric(0)  # Variables to keep track of the number of patients assigned to the trial, and then the control and treatment respectively 
  
  data_control = numeric(0)
  data_treat = numeric(0) # Creating our empty data sets
  
  for (i in 1:burn){  # Going through our burn period before adaptive designs take place
    
    if (rbinom(1,1,0.5) == 1){   #ER
      # Control
      data_control[length(data_control) + 1] = rbinom(1,1,control)  # Simulating control treatment
      
      if (i == 1){  # Base case is necessary to avoid errors with the recursion
        C_count[1] = 1
        T_count[1] = 0
      }
      else{
        C_count[i] = C_count[i-1] + 1
        T_count[i] = T_count[i-1]
      }
    }
    else{  
      # Treatment
      data_treat[length(data_treat) + 1] = rbinom(1,1,treat)   # Simulating treatment treatment
      
      if (i == 1){   # Base case is necessary to avoid errors with the recursion
        C_count[1] = 0
        T_count[1] = 1
      }
      else{
        C_count[i] = C_count[i-1] 
        T_count[i] = T_count[i-1] + 1
      }
      
    }
    count = count + 1
    
  }
  
  ### IA Anlysis
  # Work out posterior
  # Post good enough, accept
  # Post not good enough, reject
  # O/w continue
  
  while (TRUE){   # Loop until the trial finishes, at n =  count there is a mandatory break so the loop will break eventually
    
    A1 = prior_a + sum(data_treat)  # Updating the posterior
    B1 = prior_b + (tail(T_count, n=1) - sum(data_treat))
    
    A2 = prior_a + sum(data_control)  # Updating the posterior
    B2 = prior_b + (tail(C_count, n=1) - sum(data_control))
    
    if (Norm_approx == FALSE){
      
      ## One of these is exact, one is empirical c.d.f
      
      # IA_test = Beta_comp(A1,B1,A2,B2,100000)
      IA_test = Beta_comp_exact(A1,B1,A2,B2)
      
    }
    else{   # This normal approximation is primarly used for computational efficiency, such as the grid search.
      Mu_mean = A2/(A2+B2) - A1/(A1+B1)  # This way aroud because pnorm works out the lower tail, so the probability that the second variable is bigger than the first
      Sigma_var = A1*B1/((A1+B1)^2*(A1+B1+1)) + A2*B2/((A2+B2)^2*(A2+B2+1))
      IA_test = pnorm(0, Mu_mean, sqrt(Sigma_var))
    }
    
    Probs = Threshold_prob(lamda, gamma, n, count)  # Vector of threshold probabilities, [1] is futility, [2] is super. Same when trial is finished. 
    # Note that these parameters have been tuned using the non RAR version, there is no reason they should hold here
    
    
    
    # print(data_treat)
    # print(data_control)
    
    if (IA_test > Probs[2]){  # Adaptive threshold
      
      if (return_power){
        return(TRUE)
      }
      
      count_df = data.frame(Patient = 1:length(T_count),Treatment = T_count, Control = C_count, plot_prob)
      count_df = pivot_longer(count_df,cols = Treatment:Control, names_to = "Arm", values_to = "Count")  # Tidying the data a bit
      
      # Various data that can be interesting 
      # print(Probs)
      # print(IA_test)
      # print(count_df)
      # print(data_treat)
      # print(data_control)
      
      
      
      
      if (print_plots == TRUE){
        
        probability_plot = ggplot(count_df, aes( y = plot_prob, x = Patient)) +
          geom_line() +
          labs(
            title = "Allocation Probability",
            y = "Probability of treatment arm."
          ) +
          scale_y_continuous(breaks = seq(0,1,0.1), limits=c(0, 1))
        
        
        
        count_plot = ggplot(count_df, aes(fill=Arm, y= Count, x=Patient)) + 
          geom_bar(position="dodge", stat="identity") +
          labs(
            title = "Patient Allocation"  
          ) +
          scale_color_colorblind()
        
        print(probability_plot)
        print(count_plot)
      }
      
      To_return = vector("list",2)
      To_return[[1]] = TRUE
      To_return[[2]] = c("Successes" = sum(data_treat) + sum(data_control), "Number in treatment arm" = length(data_treat), "Number of patients overall" = count )
      # Successes, num in treatment arm, and total number overall
      
      return(To_return)  # Supereority check
      break
    } else if (IA_test < Probs[1] | n == count ) {   # Futility threshold or end of trial  
      
      if (return_power){
        return(FALSE)
      }
      
      count_df = data.frame(Patient = 1:length(T_count),Treatment = T_count, Control = C_count, plot_prob)
      count_df = pivot_longer(count_df,cols = Treatment:Control, names_to = "Arm", values_to = "Count")  # Tidying the data a bit
      
      # Various data that can be interesting 
      # print(Probs)
      # print(IA_test)
      # print(count_df)
      # print(data_treat)
      # print(data_control)
      
      
      
      
      if (print_plots == TRUE){
        
        probability_plot = ggplot(count_df, aes( y = plot_prob, x = Patient)) +
          geom_line() +
          labs(
            title = "Allocation Probability",
            y = "Probability of treatment arm."
          ) +
          scale_y_continuous(breaks = seq(0,1,0.1), limits=c(0, 1))
        
        
        
        count_plot = ggplot(count_df, aes(fill=Arm, y= Count, x=Patient)) + 
          geom_bar(position="dodge", stat="identity") +
          labs(
            title = "Patient Allocation"  
          ) +
          scale_color_colorblind()
        
        print(probability_plot)
        print(count_plot)
      }
      
      
      
      To_return = vector("list",2)
      To_return[[1]] = FALSE
      To_return[[2]] = c("Successes" = sum(data_treat) + sum(data_control), "Number in treatment arm" = length(data_treat), "Number of patients overall" = count )
      # Successes, num in treatment arm, and total number overall
      
      return(To_return) # Was found futile or not superior at the end of trial
      break
    }
    
    
    
    for (i in 1:(min(IAn, n - count))){  # The min means that we won't exceed n total sample size
      
      A1 = prior_a + sum(data_treat)  # Updating the posterior
      B1 = prior_b + (tail(T_count, n=1) - sum(data_treat))
      
      A2 = prior_a + sum(data_control)  # Updating the posterior
      B2 = prior_b + (tail(C_count,n=1) - sum(data_control))
      
      
      
      if (Norm_approx == FALSE){
        
        ## One of these is exact, one is an empirical c.d.f.
        
        # allo_prob = Beta_comp(A1,B1,A2,B2,100000)  # Allocation prob
        allo_prob = Beta_comp_exact(A1,B1,A2,B2)
        
      }
      else{   # This normal approximation is primarly used for computational efficiency, such as the grid search.
        Mu_mean = A2/(A2+B2) - A1/(A1+B1)  # This way around because pnorm works out the lower tail, so the probability that the second variable is bigger than the first
        Sigma_var = A1*B1/((A1+B1)^2*(A1+B1+1)) + A2*B2/((A2+B2)^2*(A2+B2+1))
        allo_prob = pnorm(0, Mu_mean, sqrt(Sigma_var))
      }
      
      
      
      
      # Read papers about TS
      
      # BRAR with tuning
      allo_prob_tune = (1-allo_prob)^(count/(2*n))  # Common way of smoothing the allocation probability. The 1-p is because the previous calculations take the probability that the treatment is better than the control, but the following code wants the opposite
      allo_prob_comp = (allo_prob)^(count/(2*n))   # Compliment
      
      allo_prob = allo_prob_tune/(allo_prob_tune + allo_prob_comp)
      
      plot_prob[length(plot_prob)+1] = (1 - allo_prob)
      
      # Clipping with untuned BRAR
      # allo_prob = max(c(1 - allo_prob, 0.2), na.rm = TRUE)   # na.rm is needed here as since we don't currently control the probabilty, it can quickly converge to 1 and lead to degenerate results
      
      # plot_prob[length(plot_prob)+1] = min(c(1 - allo_prob,0.8), na.rm= TRUE)
      
      
      if (rbinom(1,1,allo_prob) == 1){   # RAR
        # Control
        data_control[length(data_control) + 1] = rbinom(1,1,control)  # Simulating control treatment
        
        
      }
      else{  
        # Treatment
        data_treat[length(data_treat) + 1] = rbinom(1,1,treat)   # Simulating treatment treatment
        
        
      }
      count = count + 1
      C_count = append(C_count,length(data_control))
      T_count = append(T_count,length(data_treat))
    }
  }
  
  
}

Sim_trial_RAR_PBD = function(control, treat, n, IAn, burn = IAn, prior_a = 1, prior_b = 1, gamma = 0.88, lamda = 0.96, Norm_approx = FALSE, print_plots = TRUE, return_power = FALSE){
  
  ## Sim_trial_RAR_PBD runs a 2 arm BOP2 design with RAR Thompson Sampling. Starting off with a simple Binomial endpoint, uninformative prior.
  ## The difference with sim_trial RAR is that the burn in period will now have a Permuted Block Design randomisation, rather than complete random selection
  
  # control is the ORR of the control treatment
  # treat is the ORR of the target treatment =
  # n is the total sample size
  # IAn is the number of samples needed before an Interim Analysis (carried out at equal intervals)
  # burn is how many samples we should do 1:1 equal randomisation before implemening RAR. Note that the Interim Analysis will not begin until the burn has finished.Burn needs to be at least 1 to avoid errors
  ## NEED BURN TO BE EVEN!
  # prior_a and prior_b are the priors of the beta distribution, starting with an uniformative 0.5,0.5
  # gamma and lamda are variables that appear in the threshold probabilities
  # Norm approx determines whether we use a normal approximation for the Beta calculations
  
  
  # We take a beta(1,1) as our uniformative prior, and binomial with p for our data | parameter p
  
  # Slightly differently from Sim_trial, n is now the total sample size, rather than the size of each arm. Likewise, after taking IAn samples we will perform an
  # interim analysis, regardless of how the treatments were assigned. Finally, rather than some kinda of fixed block random assignment we had in the Sim_trial ( in every
  # 20 samples there would be 10 to control and 10 to treatment), this will instead use a 1:1 ratio of equal assignment before RAR begins. This is because this
  # segues better into RAR 
  
  
  count = 0
  C_count = numeric(0)
  T_count = numeric(0)  # Variables to keep track of the number of patients assigned to the trial, and then the control and treatment respectively 
  
  data_control = numeric(0)
  data_treat = numeric(0) # Creating our empty data sets
  
  
  
  if (burn %% 2 == 1){
    burn = burn + 1
    # making burn even
  }
  
  plot_prob = rep(0.5, burn) # The assignment probability we will keep track of (for plotting later. It is already known that it is 0.5 for the burn period.
  
  
  PBD = rbinom(burn,1,0.5)  # Permuted Block design
  
  # need burn even 
  while (mean(PBD) != 0.5){
    PBD = rbinom(burn,1,0.5)   # Still gives a randomised sequence, but forces equal allocation for each one
    # Note that this is technically rejection sampling. For large burn ins, this will likely be very slow
  }
  
  
  
  
  for (i in 1:burn){  # Going through our burn period before adaptive designs take place
    
    
    
    
    if (PBD[i] == 1){   #ER
      # Control
      data_control[length(data_control) + 1] = rbinom(1,1,control)  # Simulating control treatment
      
      if (i == 1){  # Base case is necessary to avoid errors with the recursion
        C_count[1] = 1
        T_count[1] = 0
      }
      else{
        C_count[i] = C_count[i-1] + 1
        T_count[i] = T_count[i-1]
      }
    }
    else{  
      # Treatment
      data_treat[length(data_treat) + 1] = rbinom(1,1,treat)   # Simulating treatment treatment
      
      if (i == 1){   # Base case is necessary to avoid errors with the recursion
        C_count[1] = 0
        T_count[1] = 1
      }
      else{
        C_count[i] = C_count[i-1] 
        T_count[i] = T_count[i-1] + 1
      }
      
    }
    count = count + 1
    
  }
  
  ### IA Anlysis
  # Work out posterior
  # Post good enough, accept
  # Post not good enough, reject
  # O/w continue
  
  while (TRUE){   # Loop until the trial finishes, at n =  count there is a mandatory break so the loop will break eventually
    
    A1 = prior_a + sum(data_treat)  # Updating the posterior
    B1 = prior_b + (tail(T_count, n=1) - sum(data_treat))
    
    A2 = prior_a + sum(data_control)  # Updating the posterior
    B2 = prior_b + (tail(C_count, n=1) - sum(data_control))
    
    if (Norm_approx == FALSE){
      
      ## One of these is exact, one is empirical c.d.f
      
      # IA_test = Beta_comp(A1,B1,A2,B2,100000)
      IA_test = Beta_comp_exact(A1,B1,A2,B2)
      
    }
    else{   # This normal approximation is primarly used for computational efficiency, such as the grid search.
      Mu_mean = A2/(A2+B2) - A1/(A1+B1)  # This way aroud because pnorm works out the lower tail, so the probability that the second variable is bigger than the first
      Sigma_var = A1*B1/((A1+B1)^2*(A1+B1+1)) + A2*B2/((A2+B2)^2*(A2+B2+1))
      IA_test = pnorm(0, Mu_mean, sqrt(Sigma_var))
    }
    
    Probs = Threshold_prob(lamda, gamma, n, count)  # Vector of threshold probabilities, [1] is futility, [2] is super. Same when trial is finished. 
    # Note that these parameters have been tuned using the non RAR version, there is no reason they should hold here
    
    
    
    # print(data_treat)
    # print(data_control)
    
    if (IA_test > Probs[2]){  # Adaptive threshold
      
      if (return_power){
        return(TRUE)
      }
      
      count_df = data.frame(Patient = 1:length(T_count),Treatment = T_count, Control = C_count, plot_prob)
      count_df = pivot_longer(count_df,cols = Treatment:Control, names_to = "Arm", values_to = "Count")  # Tidying the data a bit
      
      # Various data that can be interesting 
      # print(Probs)
      # print(IA_test)
      # print(count_df)
      # print(data_treat)
      # print(data_control)
      
      
      
      
      if (print_plots == TRUE){
        
        probability_plot = ggplot(count_df, aes( y = plot_prob, x = Patient)) +
          geom_line() +
          labs(
            title = "Allocation Probability",
            y = "Probability of treatment arm."
          ) +
          scale_y_continuous(breaks = seq(0,1,0.1), limits=c(0, 1))
        
        
        
        count_plot = ggplot(count_df, aes(fill=Arm, y= Count, x=Patient)) + 
          geom_bar(position="dodge", stat="identity") +
          labs(
            title = "Patient Allocation"  
          ) +
          scale_color_colorblind()
        
        print(probability_plot)
        print(count_plot)
      }
      
      To_return = vector("list",2)
      To_return[[1]] = TRUE
      To_return[[2]] = c("Successes" = sum(data_treat) + sum(data_control), "Number in treatment arm" = length(data_treat), "Number of patients overall" = count )
      # Successes, num in treatment arm, and total number overall
      
      return(To_return)  # Supereority check
      break
    } else if (IA_test < Probs[1] | n == count ) {   # Futility threshold or end of trial  
      
      if (return_power){
        return(FALSE)
      }
      
      count_df = data.frame(Patient = 1:length(T_count),Treatment = T_count, Control = C_count, plot_prob)
      count_df = pivot_longer(count_df,cols = Treatment:Control, names_to = "Arm", values_to = "Count")  # Tidying the data a bit
      
      # Various data that can be interesting 
      # print(Probs)
      # print(IA_test)
      # print(count_df)
      # print(data_treat)
      # print(data_control)
      
      
      
      
      if (print_plots == TRUE){
        
        probability_plot = ggplot(count_df, aes( y = plot_prob, x = Patient)) +
          geom_line() +
          labs(
            title = "Allocation Probability",
            y = "Probability of treatment arm."
          ) +
          scale_y_continuous(breaks = seq(0,1,0.1), limits=c(0, 1))
        
        
        
        count_plot = ggplot(count_df, aes(fill=Arm, y= Count, x=Patient)) + 
          geom_bar(position="dodge", stat="identity") +
          labs(
            title = "Patient Allocation"  
          ) +
          scale_color_colorblind()
        
        print(probability_plot)
        print(count_plot)
      }
      
      
      
      To_return = vector("list",2)
      To_return[[1]] = FALSE
      To_return[[2]] = c("Successes" = sum(data_treat) + sum(data_control), "Number in treatment arm" = length(data_treat), "Number of patients overall" = count )
      # Successes, num in treatment arm, and total number overall
      
      return(To_return) # Was found futile or not superior at the end of trial
      break
    }
    
    
    
    for (i in 1:(min(IAn, n - count))){  # The min means that we won't exceed n total sample size
      
      A1 = prior_a + sum(data_treat)  # Updating the posterior
      B1 = prior_b + (tail(T_count, n=1) - sum(data_treat))
      
      A2 = prior_a + sum(data_control)  # Updating the posterior
      B2 = prior_b + (tail(C_count,n=1) - sum(data_control))
      
      
      
      if (Norm_approx == FALSE){
        
        ## One of these is exact, one is an empirical c.d.f.
        
        # allo_prob = Beta_comp(A1,B1,A2,B2,100000)  # Allocation prob
        allo_prob = Beta_comp_exact(A1,B1,A2,B2)
        
      }
      else{   # This normal approximation is primarly used for computational efficiency, such as the grid search.
        Mu_mean = A2/(A2+B2) - A1/(A1+B1)  # This way around because pnorm works out the lower tail, so the probability that the second variable is bigger than the first
        Sigma_var = A1*B1/((A1+B1)^2*(A1+B1+1)) + A2*B2/((A2+B2)^2*(A2+B2+1))
        allo_prob = pnorm(0, Mu_mean, sqrt(Sigma_var))
      }
      
      
      
      
      # Read papers about TS
      
      # BRAR with tuning
      allo_prob_tune = (1-allo_prob)^(count/(2*n))  # Common way of smoothing the allocation probability. The 1-p is because the previous calculations take the probability that the treatment is better than the control, but the following code wants the opposite
      allo_prob_comp = (allo_prob)^(count/(2*n))   # Compliment
      
      allo_prob = allo_prob_tune/(allo_prob_tune + allo_prob_comp)
      
      plot_prob[length(plot_prob)+1] = (1 - allo_prob)
      
      # Clipping with untuned BRAR
      # allo_prob = max(c(1 - allo_prob, 0.2), na.rm = TRUE)   # na.rm is needed here as since we don't currently control the probabilty, it can quickly converge to 1 and lead to degenerate results
      
      # plot_prob[length(plot_prob)+1] = min(c(1 - allo_prob,0.8), na.rm= TRUE)
      
      
      if (rbinom(1,1,allo_prob) == 1){   # RAR
        # Control
        data_control[length(data_control) + 1] = rbinom(1,1,control)  # Simulating control treatment
        
        
      }
      else{  
        # Treatment
        data_treat[length(data_treat) + 1] = rbinom(1,1,treat)   # Simulating treatment treatment
        
        
      }
      count = count + 1
      C_count = append(C_count,length(data_control))
      T_count = append(T_count,length(data_treat))
    }
  }
  
  
}

OC_sim_trial_RAR = function(control, treat, n, IAn, Sim_n = 1000,gamma = 0.92, lamda = 0.9, PBD = TRUE, return_power = FALSE){
  ## OC_sim_trial_RAR finds Operating Characterestics of interest for BOP2 RAR design
  
  # control, treat, n, IAn are the same as in Sim_trial_RAR
  # Sim_n is how many simulations to carry out 
  
  # We want to find out the number of successes (Var), the proportion in treatment arm (+var), the average number of patients, and the probability that at least 55% of the final number of control or treat
  Successes = numeric(Sim_n)
  Treat_prop = numeric(Sim_n)
  Num_patients = numeric(Sim_n)
  Num_rejects = numeric(Sim_n)
  
  for (i in 1:Sim_n){
    
    if (PBD){  # If we want to use Permuted Block Design for our burn in 
      
      temp_data = Sim_trial_RAR_PBD(control, treat, n, IAn, print_plots = FALSE, gamma = gamma, lamda = lamda, return_power = return_power )
    }
    else{
      
      temp_data = Sim_trial_RAR(control, treat, n, IAn, print_plots = FALSE, gamma = gamma, lamda = lamda, return_power = return_power)
    }
    
    if (return_power){
      Num_rejects[i] = temp_data[[1]]
    } else{
      
      Num_rejects[i] = temp_data[[1]]
      Successes[i] = temp_data[[c(2,1)]]
      Treat_prop[i] = temp_data[[c(2,2)]]/temp_data[[c(2,3)]]
      Num_patients[i] = temp_data[[c(2,3)]]
    }
  }
  
  if (return_power){
    return(mean(Num_rejects))
  }
  
  Treat_55_prob = mean(Treat_prop >= 0.55)  # Probability of treatment arm being above 55% proportion
  Control_55_prob = mean(Treat_prop <= 0.45) # Probability of control arm being above 55% proportion
  
  Power = mean(Num_rejects)
  Mean_suc = mean(Successes)
  sd_suc = sqrt(var(Successes))
  Mean_treat = mean(Treat_prop)
  sd_treat = sqrt(var(Treat_prop))
  Mean_patient = mean(Num_patients)
  sd_patient = sqrt(var(Num_patients))
  
  return( c("Power" = Power, "Mean success" = Mean_suc, "S.d. Successes" = sd_suc, "Mean prop to treatment" = Mean_treat, "S.d. of prop of treatment" = sd_treat, "Mean number of patients " = Mean_patient, "S.d. of patients" = sd_patient, "Probability of 55% treatment" = Treat_55_prob, "Probability of 55% control" = Control_55_prob))
  
}

Complete_summary_sim_trial_RAR = function(gamma = 0.92, lamda = 0.9, PBD = TRUE, cond = FALSE, one_by_one = TRUE){
  ## This function aims to do a complete summary of the operating characteretics of the BOP2 trial with RAR implemented once it's been optimsed for some gamma, lamda
  # It will cover a range of treatment and controls, and will use parallelisation to do so efficiently
  # 0.1:0.9 for treatment and 0.1:0.9 for the control, using 80n and 20 Interim analysis (what lamda, gamma default values are optimised for)
  
  # I want a plot showing proportion of treatment, probability of 55% to treatment, expected successes, and power
  
  # 9 treatments and controls so it fits nicely onto a square grid
  controls = 1:9 / 10
  treats = controls
  
  par_results = vector("list",81)
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  
  par_results = foreach( i = 0:80, .export = c("Sim_trial_RAR","Beta_comp_exact","Threshold_prob","OC_sim_trial_RAR", "pivot_longer", "ggplot","Sim_trial_RAR_PBD"))  %dopar% {
    
    # What control and treatment we're going to be looking at. +1 so we don't need to fix 0 indexing
    
    c_i = 1 + i%/%9  # How many whole 9s go into i
    t_i = 1 + i%%9 # Remainder
    
    # So, the order it goes is (control = 0.1, treatment = 0.1), treatment increases until 0.9, and then it resets to (0.2,0.1), and treatment increases again
    
    if (cond){
      par_results[[i+1]] = cond_post_burn_plots(controls[c_i],treats[t_i],80,20, gamma = gamma, lamda = lamda, PBD = PBD)
    }
    else{
      par_results[[i+1]] = OC_sim_trial_RAR(controls[c_i],treats[t_i],80,20, gamma = gamma, lamda = lamda, PBD = PBD)
    }
  }
  
  stopCluster(cl) # deactivates
  
  if (one_by_one == FALSE){
    par(mfrow = c(3,3))
  }
  else{
    par(mfrow = c(1,1))
  }
  
  for (i in 1:9){
    
    temp_power = numeric(9)
    
    for (j in 1:9){
      temp_power[j] = par_results[[c(9*(i-1) + j,1)]] 
    }
    
    
    plot(x = treats, y= temp_power, main= c("Power against control:",controls[i]), type = "l", xlab = "Treatments", ylab = "Power", ylim = c(0,1), xlim=c(0,1), sub = "Red (0.1), Blue null")
    abline(h = 0.1, col = "red", lwd=1, lty=2)
    abline(v = controls[i], col = "blue")
  }
  
  for (i in 1:9){
    
    temp_prop = numeric(9)
    
    for (j in 1:9){
      temp_prop[j] = par_results[[c(9*(i-1) + j,4)]]
    }
    
    #browser()
    
    barplot(names.arg = treats, width = rep(0.1,9), height= temp_prop, main= c("Proportion to treatment, control:",controls[i]), xlab = "Treatments", ylab = "Proportion in treatment", ylim = c(0,1), xlim=c(0,1), sub = "Red 50%")
    
    lines(x = c(0.07,2),y = c(0.5,0.5), col = "red")
    #abline(v = controls[i], col = "blue")
  }
  
  for (i in 1:9){
    
    temp_prob = numeric(9)
    
    for (j in 1:9){
      temp_prob[j] = par_results[[c(9*(i-1) + j,8)]]
    }
    
    
    plot(x = treats, y= temp_prob, main= c("Prob 55% in treatment arm, control:",controls[i]), type = "l", xlab = "Treatments", ylab = "Proportion in treatment", ylim = c(0,1), xlim=c(0,1), sub = "Blue null, red 50%")
    
    abline(h = 0.5, col = "red")
    abline(v = controls[i], col = "blue")
  }
  
  for (i in 1:9){
    
    temp_patients = numeric(9)
    temp_success = numeric(9)
    
    for (j in 1:9){
      temp_patients[j] = par_results[[c(9*(i-1) + j,6)]]
      temp_success[j] = par_results[[c(9*(i-1) + j,2)]]
    }
    
    
    plot(x = treats, y= temp_patients, main= c("Num Patients + successes, control:",controls[i]), type = "l", xlab = "Treatments", ylab = "Count", ylim = c(0,100), xlim=c(0,1), sub = "Red success, green max")
    
    max_treat = (0.5*(controls[i] + treats ) + 0.5*abs(controls[i] - treats))  # Maximum between control and treatment
    
    avr_treat = 0.5*(controls[i] + treats)
    
    # print(max_treat) #testing
    
    lines(x = treats, y = temp_patients * max_treat , col = "green")  # Theoretical maxmium expected successes
    
    lines(x = treats, y = temp_success, col = "red")  # Actual successes
    
    lines(x = treats, y = temp_patients * avr_treat, col = "purple" )
  }
  
  par(mfrow = c(1,1))
  
  par(mfrow = c(2,1))
  
  for (i in c(1,5)){
    
    temp_prop = numeric(9)
    
    for (j in 1:9){
      temp_prop[j] = par_results[[c(9*(i-1) + j,4)]]
    }
    
    
    barplot(names.arg = treats, width = rep(0.1,9), height= temp_prop, main= c("Proportion of patients allocated to the experimental treatment versus treatment success rate."), xlab = "Treatment success rate", ylab = "Proportion in treatment", ylim = c(0,1), xlim=c(0,1))
    
    lines(x = c(0.07,2),y = c(0.5,0.5), col = "red")
    #abline(v = controls[i], col = "blue")
  }
  
  browser()
  
  par(mfrow = c(1,1))
  
}

Sim_trial_repeats = function(control, treat, n, IAn, n_sim = 1000, test_type = "BOP2_FE", gamma = 0.88, lamda = 0.96 ){
  ## Repeats Sim_trial n_sim times and returns mean rejection and mean sample size
  
  Rejects =  numeric(n_sim)
  Samples = numeric(n_sim)
  
  
  
  for ( i in 1:n_sim){
    temp_data = Sim_trial(control, treat,n,IAn, gamma = gamma, lamda = lamda, test_type = test_type)
    
    Rejects[i] = temp_data[1]
    Samples[i] = temp_data[2]
    
  }
  
  
  return(c(mean(Rejects),mean(Samples)))
  
}

BOP2_2ARM_data_parallel = function(gamma = 0.92, lamda = 0.9){
  ## This function aims to recreate the data from the BOP2 2 arm paper, in particular comparing power and expected sample size with 3 other classical tests
  ## The function has the same functionality as BOP2_2ARM_data but uses parallelisation for efficient running, and can be much quicker
  
  # The gamma and lamda were found using grid_search for 0.1,0.3,20,80 (what this is testing)
  
  ## Initialising everything
  Controls = rep(c(0.2,0.25,0.3), each = 6)
  
  Treats = rep(c(0.1,0.2,0.25,0.3,0.4,0.5),3)
  
  Boundaries = c("OF","Pocock","GF","BOP2","OF_FE","Pocock_FE","GF_FE","BOP2_FE")
  
  n_sim = 10000
  
  n = 40
  
  IAn = 10
  
  #Dataframe for just the futility thresholds
  
  Data_Futility = data.frame("Control" = Controls, "Treatment" =  Treats, "OF Null rejection" = numeric(18), "OF expected sample size" = numeric(18),
                             "Pocock Null rejection" = numeric(18), "Pocock expected sample size" = numeric(18),
                             "GF Null rejection" = numeric(18), "GF expected sample size" = numeric(18),
                             "BOP2 Null rejection" = numeric(18), "BOP2 expected sample size" = numeric(18)
  )
  
  #Dataframe for the futility and success thresholds
  Data_Fut_Treat = data.frame("Control" = Controls, "Treatment" =  Treats, "OF Null rejection" = numeric(18), "OF expected sample size" = numeric(18),
                              "Pocock Null rejection" = numeric(18), "Pocock expected sample size" = numeric(18),
                              "GF Null rejection" = numeric(18), "GF expected sample size" = numeric(18),
                              "BOP2 Null rejection" = numeric(18), "BOP2 expected sample size" = numeric(18)
  )
  
  ## Parallelisation set up 
  
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl) 
  
  
  Results_par = vector("list", length(Controls) )  # The foreach target
  
  Results_par = foreach (i = 1:length(Controls), .export = c("Sim_trial","Beta_comp_exact","Threshold_prob","Sim_trial_repeats" ))%dopar% {
    
    temp_boundary = vector("list",length(Boundaries))
    
    for(k in 1:length(Boundaries)){
      
      
      temp_boundary[[k]] = Sim_trial_repeats(Controls[i],Treats[i],n,IAn,n_sim, gamma = gamma, lamda = lamda, test_type = Boundaries[k])  # simulating data
      #print(c(Controls[i],Treats[i],i,temp_data))  WARNING Print error checking in parallel foreach loops slows the loops down and can mess with results
    }
    
    
    Results_par[[i]] = temp_boundary
  }
  
  stopCluster(cl)
  
  # Transcribing the data from a big list to readable data frames.  
  for (i in 1:length(Controls)){
    for (k in 1:length(Boundaries)){
      if (k<5){
        Data_Futility[[1+2*k]][i] = Results_par[[c(i,k,1)]]
        Data_Futility[[2+2*k]][i] = Results_par[[c(i,k,2)]]    # The first index represents the control/treatment, the second represents the boundary, the third represents null reject/sample size 
      }
      else{
        Data_Fut_Treat[[1+2*(k-4)]][i] = Results_par[[c(i,k,1)]]
        Data_Fut_Treat[[2+2*(k-4)]][i] = Results_par[[c(i,k,2)]]    # The first index represents the control/treatment, the second represents the boundary, the third represents null reject/sample size 
      }
    }
  }
  
  to_return = vector("list",2)
  to_return[[1]] = Data_Futility
  to_return[[2]] = Data_Fut_Treat
  return(to_return)
}

BOP2_rar_size_test =  function(gamma = 0.92){
  ## This function will just output a table of type one errors for the BOP2 test with RAR
  
  controls = 1:9/10
  
  control_par_loop = rep(controls, 10)  # This is a handy vector to cycle through in the foreach loop 
  
  par_results = numeric(81) # What the paralllisation will go into
  
  lamda_seq = seq(0.9,0.99,0.01)  # The lamda we want to look at
  
  lamda_seq = rep(lamda_seq, each = 9) # Making it useful for the foreach loop
  
  #Parallaleisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl) 
  
  
  # We want to do all of the lamda = 0.9 first at the different controls, then 0.91,0.92, ext
  
  par_results = foreach (i = 1:90, .export = c("Sim_trial_RAR","Beta_comp_exact","Threshold_prob","OC_sim_trial_RAR","pivot_longer" ))%dopar% {
    
    par_results[i] = OC_sim_trial_RAR(control_par_loop[i], control_par_loop[i], 80, 20, Sim_n = 10000, gamma = gamma, lamda = lamda_seq[i]  )[[1]]
    
    return(par_results[i])
    
  }
  
  stopCluster(cl) # Stopping parallisation
  
  
  vector_results = numeric(90)  # A vector to transcribe the parallel list
  
  for ( i in 1:90){
    vector_results[i] = par_results[[i]]
  }
  
  
  table_df = data.frame("Control" = controls, "Type 1 error 0.9 lambda" = vector_results[1:9],
                        "Type 1 error 0.91 lambda" = vector_results[10:18],
                        "Type 1 error 0.92 lambda" = vector_results[19:27],
                        "Type 1 error 0.93 lambda" = vector_results[28:36],
                        "Type 1 error 0.94 lambda" = vector_results[37:45],
                        "Type 1 error 0.95 lambda" = vector_results[46:54],
                        "Type 1 error 0.96 lambda" = vector_results[55:63],
                        "Type 1 error 0.97 lambda" = vector_results[64:72],
                        "Type 1 error 0.98 lambda" = vector_results[73:81],
                        "Type 1 error 0.99 lambda" = vector_results[82:90])
  
  kable(table_df, digits = 3, format = "html", row.names = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover"),
                  full_width = T,
                  font_size = 15,
                  position = "left")
  
}

power_plot_RAR_v_ER = function(gamma = 0.92, lamda = 0.9, n = 80, IAn = 20, n_sim = 10000, PBD = TRUE){
  ## power_plot_RAR_v_ER is a function that calculates power across a range of 
  # null and alternate hypothesis for both the RAR trial and the BR trial (at
  # fixed lambda gamma) and then compares them in a plot. 
  
  controls = 1:9 / 10
  treats = controls
  
  par_results = vector("list",81)
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  
  par_results = foreach( i = 0:80, .export = c("Sim_trial_RAR","Beta_comp_exact","Threshold_prob","OC_sim_trial_RAR", "pivot_longer", "ggplot","Sim_trial_RAR_PBD","Sim_trial", "Sim_trial_repeats"))  %dopar% {
    
    # What control and treatment we're going to be looking at. +1 so we don't need to fix 0 indexing
    
    c_i = 1 + i%/%9  # How many whole 9s go into i
    t_i = 1 + i%%9 # Remainder
    
    # So, the order it goes is (control = 0.1, treatment = 0.1), treatment increases until 0.9, and then it resets to (0.2,0.1), and treatment increases again
    
    
    temp_BR_power = Sim_trial_repeats(controls[c_i],treats[t_i], n/2, IAn/2, gamma = gamma, lamda = lamda, n_sim = n_sim)[1] # Balanced/equal randomisation
    
    temp_BRAR_power = OC_sim_trial_RAR(controls[c_i],treats[t_i],n ,IAn , gamma = gamma, lamda = lamda, PBD = PBD, Sim_n = n_sim) # BRAR
    
    par_results[[i+1]] = list(temp_BR_power,temp_BRAR_power)
  }
  
  stopCluster(cl) # deactivates
  
  par(mfrow = c(3,3))
  
  for (i in 1:9){
    
    temp_power_BRAR = numeric(9)
    temp_power_BR = numeric(9)
    
    for (j in 1:9){
      temp_power_BRAR[j] = par_results[[c(9*(i-1) + j,2,1)]] 
      temp_power_BR[j] = par_results[[c(9*(i-1) + j,1)]]
    }
    
    
    plot(x = treats, y= temp_power_BRAR, main= c("Power against control:",controls[i]), type = "l", xlab = "Treatments", ylab = "Power", ylim = c(0,1), xlim=c(0,1), sub = "Red (0.1), Blue null, Green Balanced Randomisation, Black BRAR")
    lines(x = treats, y= temp_power_BR, col = "green")
    abline(h = 0.1, col = "red", lwd=1, lty=2)
    abline(v = controls[i], col = "blue")
  }
  
  
  par(mfrow = c(1,1))
}

cond_post_burn_plots = function(control, treat, n, IAn, Sim_n = 1000,gamma = 0.92, lamda = 0.9, PBD = TRUE){
  
  ## cond_post_burn_plots is mainly a copy of OC_sim_trial_RAR 
  ## The difference is in the objective: cond_post_burn_plots only aims to get conditional plots after the burn in period
  
  
  # control, treat, n, IAn are the same as in Sim_trial_RAR
  # Sim_n is how many simulations to carry out 
  
  # We want to find out the number of successes (Var), the proportion in treatment arm (+var), the average number of patients, and the probability that at least 55% of the final number of control or treat
  Successes = numeric(Sim_n)
  Treat_prop = numeric(Sim_n)
  Num_patients = numeric(Sim_n)
  Num_rejects = numeric(Sim_n)
  Num_post_burn_in = 0 
  
  for (i in 1:Sim_n){
    
    if (PBD){  # If we want to use Permuted Block Design for our burn in 
      
      temp_data = Sim_trial_RAR_PBD(control, treat, n, IAn, print_plots = FALSE, gamma = gamma, lamda = lamda)
    }
    else{
      
      temp_data = Sim_trial_RAR(control, treat, n, IAn, print_plots = FALSE, gamma = gamma, lamda = lamda)
    }
    
    
    if ( temp_data[[c(2,3)]] == IAn){ # If it terminated at the burn in
      
      # Adjusting the data so it doesn't get recorded
      
      Num_rejects[i] = NaN
      Successes[i] = NaN
      Treat_prop[i] = NaN
      Num_patients[i] = NaN
    }
    
    else{
      
      Num_rejects[i] = temp_data[[1]]
      Successes[i] = temp_data[[c(2,1)]]
      Treat_prop[i] = temp_data[[c(2,2)]]/temp_data[[c(2,3)]]
      Num_patients[i] = temp_data[[c(2,3)]]
      
      Num_post_burn_in = Num_post_burn_in  + 1
      
    }
    
    
  }
  
  Treat_55_prob = mean(Treat_prop >= 0.55, na.rm = TRUE)  # Probability of treatment arm being above 55% proportion
  Control_55_prob = mean(Treat_prop <= 0.45, na.rm = TRUE) # Probability of control arm being above 55% proportion
  
  
  
  Power = mean(Num_rejects, na.rm = TRUE)
  Mean_suc = mean(Successes, na.rm = TRUE)
  sd_suc = sqrt(var(Successes, na.rm = TRUE))
  Mean_treat = mean(Treat_prop, na.rm = TRUE)
  sd_treat = sqrt(var(Treat_prop, na.rm = TRUE))
  Mean_patient = mean(Num_patients, na.rm = TRUE)
  sd_patient = sqrt(var(Num_patients, na.rm = TRUE))
  
  return( c("Power" = Power, "Mean success" = Mean_suc, "S.d. Successes" = sd_suc, "Mean prop to treatment" = Mean_treat, "S.d. of prop of treatment" = sd_treat, "Mean number of patients " = Mean_patient, "S.d. of patients" = sd_patient, "Probability of 55% treatment" = Treat_55_prob, "Probability of 55% control" = Control_55_prob, "Number of trials passing burn in" = Num_post_burn_in))
  
}

parameter_heatmaps = function(control = 0.2, RAR = FALSE, PBD = TRUE, n_sim = 10000){
  ## This function investigates the effect of parameters on type one, power
  # and how varying them can have an effect
  
  gamma = seq(0.1,1,0.05)
  lambda = seq(0.1,1,0.05)
  
  temp_len = length(gamma)
  gamma = rep(gamma, each = length(lambda))
  lambda = rep(lambda, temp_len)
  
  # Representing type one errors
  alpha = numeric(length(lambda))
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  
  alpha = foreach(i = 1:length(gamma), .export = c("Sim_trial_RAR","Beta_comp_exact","Threshold_prob","OC_sim_trial_RAR", "pivot_longer", "ggplot","Sim_trial_RAR_PBD","Sim_trial", "Sim_trial_repeats"), .combine = c)  %dopar% {
    
    
    return(Sim_trial_repeats(control, control, 80 , 20, gamma = gamma[i], lamda = lambda[i], n_sim = n_sim)[1])
  }
  
  stopCluster(cl) # deactivates
  
  heat_data = data.frame(lambda,gamma,alpha)
  
  var_subtitle = sprintf("Control is %s RAR is %f",control, RAR)
  
  ggp <- ggplot(heat_data, aes(lambda, gamma)) +                           # Create heatmap with ggplot2
    geom_tile(aes(fill = alpha)) +
    labs(
      title = "Type one error heatmap",
      subtitle = var_subtitle
    )
  print(ggp)
  
}

type_one_error_heatmap = function(opt_gamma = c(), opt_lambda = c(), controls = seq(0.1,0.9,0.1), RAR = FALSE, PBD = TRUE, n_sim = 10000){
  ## This function will investigate how fitting for one control can effect the 
  # type one error if the control is actually something else: testing robustness.
  
  # opt_gamma, opt_lambda are the optimal gamma and lambda parameters
  # controls is the vector of controls they have been optimised at. 
  
  type_one_error = numeric(length(controls)^2)
  
  for (i in 1:length(controls)){
    # The fitted control
    for (j in 1:length(controls)){
      # The actual control
      if (RAR){
        type_one_error[ length(controls)*(i-1) + j] = OC_sim_trial_RAR(controls[i], controls[i], 80 , 20, gamma = opt_gamma[j], lamda = opt_lambda[j], Sim_n = n_sim, PBD = PBD)[1]
      }
      else{
        type_one_error[ length(controls)*(i-1) + j] = Sim_trial_repeats(controls[i], controls[i], 40 , 10, gamma = opt_gamma[j], lamda = opt_lambda[j], n_sim = n_sim)[1]
      }
    }
    
    
    
  }
  
  # Making it into a dataframe
  
  Fitted_Controls = rep(controls, each = length(controls))
  
  Actual_Controls = rep(controls, length(controls))
  
  heat_data = data.frame(Fitted_Controls, Actual_Controls, type_one_error)
  
  
  var_subtitle = sprintf("RAR is %f", RAR)
  
  #browser()
  
  # heatmap
  ggp <- ggplot(heat_data, aes(Fitted_Controls, Actual_Controls)) +                           # Create heatmap with ggplot2
    geom_tile(aes(fill = type_one_error)) +
    labs(
      title = "Fitted parameters type one error with other controls",
      subtitle = var_subtitle
    ) +
    scale_fill_gradient2(low="navy", mid="white", high="red", 
                         midpoint=0.1, limits=c(0.05,0.15))
  print(ggp)
  
  
  
  
  
}

Grid_search_par =  function(null, alt, size, n = 80, IAn = 20, n_sim = 1000, RAR = FALSE, PBD = TRUE, return_power = TRUE){
  ## Grid search function to find a lamda and gamma that controls for size and optimises power
  # null is what null treatment hypothesis should be taken
  # alt is what alterbative hypothesis for the treatment should be taken
  # size is the target type 1 error
  
  # The extra RAR bit makes it so the more intense RAR computations have to be done less,
  # at the cost of some accuracy
  pot_lamda = seq(0.7 + 0.1*RAR,1,0.01 + 0.01*RAR) # Can vary the parameters to try and get the right kinda area, going off table S9 from supp material
  pot_gamma = seq(0.7 + 0.1*RAR,1,0.01 + 0.01*RAR)
  
  temp_len = length(pot_lamda)  # Before the length gets updated
  
  # For easy looping in Foreach 
  pot_lamda = rep(pot_lamda, length(pot_gamma))
  pot_gamma = rep(pot_gamma, each = temp_len)
  
  par_power = numeric(length(pot_gamma))
  
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  
  par_power = foreach( i = 1:length(pot_gamma), .combine = c, .export = c("Sim_trial_RAR","Beta_comp_exact","Threshold_prob","OC_sim_trial_RAR", "pivot_longer", "ggplot","Sim_trial_RAR_PBD","Sim_trial", "Sim_trial_repeats")) %dopar%{
    
    if (RAR){
      # Do RAR power stuff
      
      if ( OC_sim_trial_RAR(null, null, n = n, IAn =  IAn, Sim_n = n_sim, gamma = pot_gamma[i], lamda = pot_lamda[i], PBD = PBD, return_power = return_power ) < 0.1){
        return(OC_sim_trial_RAR(null, alt, n = n, IAn =  IAn, Sim_n = n_sim, gamma = pot_gamma[i], lamda = pot_lamda[i], PBD = PBD, return_power = return_power))
      }
      else{
        return(0)
      }
      
      
    }
    else{
      # Do BR power stuff
      
      if ( Sim_trial_repeats(null, null, n = n/2, IAn =  IAn/2, n_sim = n_sim, gamma = pot_gamma[i], lamda = pot_lamda[i] )[1] < 0.1){
        return(Sim_trial_repeats(null, alt, n = n/2, IAn =  IAn/2, n_sim = n_sim, gamma = pot_gamma[i], lamda = pot_lamda[i])[1])
      }
      else{
        return(0)
      }
      
    }
    
    
  }
  
  stopCluster(cl) # deactivates
  
  max_index = which.max(par_power)
  
  results = c("Lambda" = pot_lamda[max_index], "Gamma" = pot_gamma[max_index], "Power" = par_power[max_index])
  
  print(results)
  
  return(results)
  
  
}

optimal_vec = function(control = seq(0.1,0.9,0.1), treats = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9,0.95), n_sim = 10000, RAR = FALSE, PBD = TRUE, return_power = TRUE){
  # This function takes a vector of controls and treatments and then optimises at each one
  
  if ( length(control) != length(treats) ){
    return("Mismatching control and treat length")
  }
  
  n = length(control)
  
  # Vectors for the optimal parameters to go into 
  opt_lambda = numeric(n)
  opt_gamma = numeric(n)
  
  for (i in 1:n){
    temp_para = Grid_search_par(control[i],treats[i],0.1, n_sim = n_sim, RAR = RAR, PBD = PBD, return_power = return_power)
    
    opt_lambda[i] = temp_para[1]
    opt_gamma[i] = temp_para[2]
  }
  
  return(list(opt_lambda,opt_gamma))
  
}

## Base functions for BOP2 Multi-arm (uncontrolled)
# H0 is normally eff 0.45,tox 0.3 per arm, H1 is normally eff 0.6, tox 0.2 per arm.  


Multi_arm_BOP2 = function(treats_eff, treats_tox, n = 60, IAn = 15, control = FALSE, phi_e = 0.45, phi_t = 0.3, threshold_type = 1, gamma = 0.98, lamda = 0.81, exact = TRUE){
  ## Multi Arm BOP2 aims to extend BOP2 to the multi arm setting, like in Mulier et al 2024
  
  # treats_ ss a vector of treatments. If control = TRUE, then treats[1] will be the control 
  # treats_eff is the vector of efficacy (ex 0.7 is the probability of a desirable outcome)
  # treats_tox  is the vector of toxicity (ex 0.6 is the probabilty of toxic response)
  # So treats_eff = c(0.1,0.2,0.3) and treats_tox = c(0.5,0.4,0.3) correspond to 3 treatments
  # One with 0.1 success 0.5 tox,then 0.2 success 0.4 tox, then 0.3 success and tox
  
  # Like the multi arm paper, n means n patients per arm
  # IAn is after how many patients to perform an interim analysis
  # The default values are the ones they used for their testing
  # control = TRUE/FALSE determines whether there is a control or not. If TRUE, then treats[1] is a control
  # phi_e, phi_f are the efficacy/futility thresholds for the uncontrolled setting
  # lamda are gamma are parameters that are used to control the type one error whilst maximising power 
  
  
  
  # For ease of calculation, we will take a Dirichlet(1,1,1,1) prior 
  # (this allows us to use Beta_comp_exact)
  
  
  num_treat = length(treats_eff)
  continue = rep(1,num_treat)  # Will be 1 if we continue with a treatment, 0 futile (no test for efficacy)
  
  # Where the responses for efficacy and toxicity will go
  Responses_eff = vector("list",num_treat)
  Responses_tox = vector("list",num_treat)
  
  for (i in 1:num_treat){  # Initial bit of samples for every treatment before IA
    Responses_eff[[i]] = rbinom(IAn, 1, treats_eff[i])
    Responses_tox[[i]] = rbinom(IAn, 1, treats_tox[i])
  }
  
  
  
  if (control == FALSE){
    
    if (threshold_type == 1){  # There are varying threshold types the paper gives
      
      #Threshold type 1 is most similar to the original BOP2 paper. 
      
      threshold = 1 - (lamda * ((IAn / n)^gamma))
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
      
      
      if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
        
        continue[i] = 0  # This has either exceeded futility or toxicity
        
      }
      
    }
    
    remain = n - IAn # Number of max samples left
    
    remain_IAn = IAn # Number of samples left until the next IA
    
    while( remain > 0 ){  # Will continue the sampling until we're done
      
      # Only want to simulate those that are still eligible
      for (i in 1:num_treat){
        if (continue[i] == 1){
          
          # Simulating new variables
          Responses_eff[[i]] = append(Responses_eff[[i]],rbinom(1, 1, treats_eff[i]))
          Responses_tox[[i]] = append(Responses_tox[[i]],rbinom(1, 1, treats_tox[i]))
          
          
        }
      }
      
      remain_IAn = remain_IAn - 1
      
      remain =  remain - 1
      
      if (remain_IAn == 0 & remain != 0){   # Interim analysis, if remain = 0 skips to final analysis
        
        remain_IAn = IAn
        
        for (i in 1:num_treat){  # Looking at every treatment
          
          if (continue[i] == 1){  # We only want to analyse the arms being updated
            
            ## Temp values to do the beta calcs with 
            temp_n = length(Responses_eff[[i]])
            temp_eff = sum(Responses_eff[[i]])
            temp_tox = sum(Responses_tox[[i]])
            
            
            if (threshold_type == 1){
              
              
              threshold = 1 - (lamda * ((temp_n / n)^gamma))
            }
            
            #print(threshold)
            #print(temp_n)
            #print(temp_eff)
            #print(temp_tox)
            
            
            
            
            if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
              
              continue[i] = 0  # This has either exceeded futility or toxicity
              
            }
          }  
        }
        
      }
      
      
      
      
    }
    
    early_stop = rep(TRUE, num_treat)  # Recording whether a treatment stopped early
    
    final_sample_size = numeric(num_treat) # Recording the samples size of each treatment
    
    # One final test for efficacy/toxicity
    
    for (i in 1:num_treat){  # Looking at every treatment
      
      final_sample_size[i] = length(Responses_eff[[i]])
      
      if (continue[i] == 1){  # We only want to analyse the arms being updated
        
        early_stop[i] = FALSE # If it has made it to the final analysis, it has not stopped early
        
        ## Temp values to do the beta calcs with 
        temp_n = length(Responses_eff[[i]])
        temp_eff = sum(Responses_eff[[i]])
        temp_tox = sum(Responses_tox[[i]])
        
        if (threshold_type == 1){  # There are varying threshold types the paper gives
          
          #Threshold type 1 is most similar to the original BOP2 paper. 
          
          threshold = 1 - lamda 
        }
        
        
        #print(threshold)
        #print(temp_n)
        #print(temp_eff)
        #print(temp_tox)
        #browser()
        
        if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
          
          continue[i] = 0  # This has either exceeded futility or toxicity
          
        }
      }  
    }
    
    
  }
  
  else{
    
    if (threshold_type == 1){  # There are varying threshold types the paper gives
      
      #Threshold type 1 is most similar to the original BOP2 paper. 
      
      threshold = 1 - (lamda * ((IAn / n)^gamma))
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
        if ( Beta_comp_exact( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
          
          continue[i] = 0  # This has either exceeded futility or toxicity
          
        }
        
      }
      else{
        if ( Beta_comp( 0.45 + temp_eff, 0.55 + temp_n - temp_eff, 0.45 + con_eff, 0.55 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 0.3 + temp_tox, 0.7 + temp_n - temp_tox, 0.3 + con_tox, 0.7 + con_n - con_tox) > threshold ){
          
          continue[i] = 0  # This has either exceeded futility or toxicity
          
        }
        
      }
    }
    
    
    remain = n - IAn # Number of max samples left
    
    remain_IAn = IAn # Number of samples left until the next IA
    
    while( remain > 0 ){  # Will continue the sampling until we're done
      
      # Only want to simulate those that are still eligible
      for (i in 1:num_treat){
        if (continue[i] == 1){
          
          # Simulating new variables
          Responses_eff[[i]] = append(Responses_eff[[i]],rbinom(1, 1, treats_eff[i]))
          Responses_tox[[i]] = append(Responses_tox[[i]],rbinom(1, 1, treats_tox[i]))
          
          
        }
      }
      
      remain_IAn = remain_IAn - 1
      
      remain =  remain - 1
      
      if (remain_IAn == 0 & remain != 0){   # Interim analysis, if remain = 0 skips to final analysis
        
        remain_IAn = IAn
        
        con_n = length(Responses_eff[[1]])  ## Control results
        con_eff = sum(Responses_eff[[1]])
        con_tox = sum(Responses_tox[[1]])
        
        for (i in 2:num_treat){  # Looking at every treatment, not the control
          
          if (continue[i] == 1){  # We only want to analyse the arms being updated
            
            ## Temp values to do the beta calcs with 
            temp_n = length(Responses_eff[[i]])
            temp_eff = sum(Responses_eff[[i]])
            temp_tox = sum(Responses_tox[[i]])
            
            
            if (threshold_type == 1){
              
              
              threshold = 1 - (lamda * ((temp_n / n)^gamma))
            }
            
            #print(threshold)
            #print(temp_n)
            #print(temp_eff)
            #print(temp_tox)
            
            
            
            
            if (exact){
              if ( Beta_comp_exact( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
                
                continue[i] = 0  # This has either exceeded futility or toxicity
                
              }
              
            }
            else{
              if ( Beta_comp( 0.45 + temp_eff, 0.55 + temp_n - temp_eff, 0.45 + con_eff, 0.55 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 0.3 + temp_tox, 0.7 + temp_n - temp_tox, 0.3 + con_tox, 0.7 + con_n - con_tox) > threshold ){
                
                continue[i] = 0  # This has either exceeded futility or toxicity
                
              }
              
            }
          }
          
        }
        
      }
      
      
      
      
    }
    
    early_stop = rep(TRUE, num_treat)  # Recording whether a treatment stopped early
    
    final_sample_size = numeric(num_treat) # Recording the samples size of each treatment
    
    # Manually inputting the control results
    early_stop[1] = FALSE
    final_sample_size[1] = length(Responses_eff[[1]])
    
    
    
    # One final test for efficacy/toxicity
    con_n = length(Responses_eff[[1]])  ## Control results
    con_eff = sum(Responses_eff[[1]])
    con_tox = sum(Responses_tox[[1]])
    
    
    for (i in 2:num_treat){  # Looking at every treatment, not control
      
      final_sample_size[i] = length(Responses_eff[[i]])
      
      if (continue[i] == 1){  # We only want to analyse the arms being updated
        
        early_stop[i] = FALSE # If it has made it to the final analysis, it has not stopped early
        
        ## Temp values to do the beta calcs with 
        temp_n = length(Responses_eff[[i]])
        temp_eff = sum(Responses_eff[[i]])
        temp_tox = sum(Responses_tox[[i]])
        
        if (threshold_type == 1){  # There are varying threshold types the paper gives
          
          #Threshold type 1 is most similar to the original BOP2 paper. 
          
          threshold = 1 - lamda 
        }
        
        
        #print(threshold)
        #print(temp_n)
        #print(temp_eff)
        #print(temp_tox)
        #browser()
        
        if (exact){
          if ( Beta_comp_exact( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
            
            continue[i] = 0  # This has either exceeded futility or toxicity
            
          }
          
        }
        else{
          if ( Beta_comp( 0.45 + temp_eff, 0.55 + temp_n - temp_eff, 0.45 + con_eff, 0.55 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 0.3 + temp_tox, 0.7 + temp_n - temp_tox, 0.3 + con_tox, 0.7 + con_n - con_tox) > threshold ){
            
            continue[i] = 0  # This has either exceeded futility or toxicity
            
          }
          
        }
      }  
    }
    
    
  }
  
  
  # Results 
  # print(Responses_eff)
  # print(Responses_tox)  
  
  to_return = list( early_stop, final_sample_size, continue ) # The three variables of interst
  # Need these to replicate their results
  
  return(to_return)
  
}

type_one_error_exact = function(lambda = 0.81,gamma = 0.98, p_eff = 0.45, p_tox = 0.3, least_power = FALSE, phi_e = 0.45, phi_t = 0.3){
  ## This function returns an exact type one error for the lambda and gamma given
  # It calculates it via working out every feasible run and then working out the prob
  
  # p_eff and p_tox and the probability of an efficacious and toxic result respectively
  # least_power is boolean indicating whether you want the power returned rather than the FWER
  
  # Phi_e and Phi_t are the beta prior parameters as well as the null scenario
  
  # Threshold function
  C_n = 1 - lambda * ((1:4)/4)^gamma 
  
  threshold_eff = numeric(4)
  threshold_tox = numeric(4)
  
  
  # Calculating the efficacy threshold values
  for( j in 1:4){
    i = 0
    
    
    while (pbeta(phi_e,phi_e + i, 1 - phi_e + 15*j - i) > C_n[j]) {
      i = i + 1
      
      if (1 - phi_e + 15*j - i < 0){
        # Error catching, grid search can make the lambda gamma 1
        i = i - 1
        break
      }
      
    }
    
    threshold_eff[j] = i 
    #print(i)
  }
  
  # Caclulating the toxicity threshold values
  for( j in 1:4){
    i = 0
    
    while (pbeta(phi_t,phi_t + 15*j - i, 1 - phi_t  + i) < 1 - C_n[j]) {
      i = i + 1
      
      if ( 1 - phi_t + 15*j - i < 0){
        # Error catching, grid search can make the lambda gamma 1
        i = i - 1
        break
      }
      
    }
    threshold_tox[j] = i 
    #print( 15*j - i)
  }
  
  # This quadruple for loop determines whether the run is feasible for effiacy,
  # and if so adds the probability of it happening 
  eff_prob = 0
  for (i in 0:15){
    for (j in 0:15){
      for (k in 0:15){
        for (r in 0:15){
          if ((i >= threshold_eff[1]) & (i + j >= threshold_eff[2]) & (i + j + k >= threshold_eff[3]) & (i + j + k + r >= threshold_eff[4])){
            eff_prob = eff_prob + (pbinom(i,15,p_eff) - pbinom((i-1),15,p_eff))*(pbinom(j,15,p_eff) - pbinom(j-1,15,p_eff))*(pbinom(k,15,p_eff) - pbinom(k-1,15,p_eff))*(pbinom(r,15,p_eff) - pbinom(r-1,15,p_eff))
          }
        }
      }
    }
  }
  # Same for toxicity
  tox_prob = 0
  for (i in 0:15){
    for (j in 0:15){
      for (k in 0:15){
        for (r in 0:15){
          if ((i <= 15 - threshold_tox[1]) & (i + j <= 30 - threshold_tox[2]) & (i + j + k <= 45 - threshold_tox[3]) & (i + j + k + r <= 60 - threshold_tox[4])){
            tox_prob = tox_prob + (pbinom(i,15,p_tox) - pbinom((i-1),15,p_tox))*(pbinom(j,15,p_tox) - pbinom(j-1,15,p_tox))*(pbinom(k,15,p_tox) - pbinom(k-1,15,p_tox))*(pbinom(r,15,p_tox) - pbinom(r-1,15,p_tox))
          }
        }
      }
    }
  }
  
  # The run is successful if both eff and tox weren't rejected, so both need to be met
  prob = eff_prob*tox_prob
  
  if (least_power){
    # This returns the least power
    return(prob)
  }
  # Otherwise it calculates the FWER 
  
  
  # As each arm is i.i.d , then the FWER is just the product of them all
  exact_type_one_error = 1 - (1 - prob)^3
  
  return(exact_type_one_error)
}

OC_exact = function(lambda = 0.81,gamma = 0.98, p_eff = 0.45, p_tox = 0.3, least_power = FALSE, phi_e = 0.45, phi_t = 0.3){
  ## This function returns an exact type one error for the lambda and gamma given
  # It calculates it via working out every feasible run and then working out the prob
  
  # p_eff and p_tox and the probability of an efficacious and toxic result respectively
  # least_power is boolean indicating whether you want the power returned rather than the FWER
  
  
  # Threshold function
  C_n = 1 - lambda * ((1:4)/4)^gamma 
  
  threshold_eff = numeric(4)
  threshold_tox = numeric(4)
  
  
  # Calculating the efficacy threshold values
  for( j in 1:4){
    i = 0
    
    while (pbeta(phi_e,phi_e + i, 1 - phi_e + 15*j - i) > C_n[j]) {
      i = i + 1
    }
    
    threshold_eff[j] = i 
    #print(i)
  }
  
  # Caclulating the toxicity threshold values
  for( j in 1:4){
    i = 0
    
    while (pbeta(phi_t,phi_t + 15*j - i, 1 - phi_t  + i) < 1 - C_n[j]) {
      i = i + 1
    }
    threshold_tox[j] = i 
    #print( 15*j - i)
  }
  
  # This quadruple for loop determines whether the run is feasible for effiacy,
  # and if so adds the probability of it happening 
  
  # Type one error/power
  eff_prob = 0
  
  # The probability of stopping at 15,30,45, and 60 respectively for eff
  # Can be used to calculate the ESS (dot product with 15,30,45,60)
  # Can be used to calculate early stopping (sum of first 3)
  p_x_eff = c(0,0,0,0)
  
  # Used to avoid multiplicity in probabilities
  no_repeats = numeric(0)
  
  for (i in 0:15){
    for (j in 0:15){
      for (k in 0:15){
        for (r in 0:15){
          if ((i >= threshold_eff[1]) & (i + j >= threshold_eff[2]) & (i + j + k >= threshold_eff[3]) & (i + j + k + r >= threshold_eff[4])){
            # Temp as we're adding it to two things so saves recalculating
            temp_prob = (pbinom(i,15,p_eff) - pbinom((i-1),15,p_eff))*(pbinom(j,15,p_eff) - pbinom(j-1,15,p_eff))*(pbinom(k,15,p_eff) - pbinom(k-1,15,p_eff))*(pbinom(r,15,p_eff) - pbinom(r-1,15,p_eff))
            
            eff_prob = eff_prob + temp_prob
            
            # Adding to probability of ending at 4
            p_x_eff[4] = p_x_eff[4] + temp_prob
            
            
            # The below series of if/elses looks complicated, but ultimately is just
            # Summing up the probabilities of the situations where it ended at 15,
            # 30, 45, and 60 without rejecting the null.
          } else{
            
            if (i < threshold_eff[1]){
              # If ends at 15
              
              # Checks whether we've already calculated this probability
              if(! (i %in% no_repeats)){
                
                p_x_eff[1] = p_x_eff[1] + (pbinom(i,15,p_eff) - pbinom((i-1),15,p_eff))
                
                no_repeats = append(no_repeats,i)
              }
            } else {
              
              if ((i + j < threshold_eff[2])){
                # If ends at 30
                
                # Checking whether i+j is in no_repeats, +100 added for no overlap, same reasosning for 20
                if(! ((i + 20*j + 100) %in% no_repeats)){
                  
                  p_x_eff[2] = p_x_eff[2] + (pbinom(i,15,p_eff) - pbinom((i-1),15,p_eff))*(pbinom(j,15,p_eff) - pbinom(j-1,15,p_eff))
                  
                  no_repeats = append(no_repeats,i + 20*j + 100)
                  
                }
              } else{
                
                if ((i + j + k < threshold_eff[3])){
                  # If ends at 45
                  
                  # Checking where i+j+k in no_repeats, +1000 added for no overlap, same reasosing for 20,500
                  if(! ((i + 20*j + 500*k + 1000) %in% no_repeats)){
                    
                    p_x_eff[3] = p_x_eff[3] + (pbinom(i,15,p_eff) - pbinom((i-1),15,p_eff))*(pbinom(j,15,p_eff) - pbinom(j-1,15,p_eff))*(pbinom(k,15,p_eff) - pbinom(k-1,15,p_eff))
                    
                    no_repeats = append(no_repeats,i + 20*j + 500*k + 1000)
                    
                  }
                } else{
                  # If ends at 60
                  
                  p_x_eff[4] = p_x_eff[4] + (pbinom(i,15,p_eff) - pbinom((i-1),15,p_eff))*(pbinom(j,15,p_eff) - pbinom(j-1,15,p_eff))*(pbinom(k,15,p_eff) - pbinom(k-1,15,p_eff))*(pbinom(r,15,p_eff) - pbinom(r-1,15,p_eff))
                  
                }
              }
            }
          }
        }
      }
    }
  }
  # Same for toxicity
  tox_prob = 0
  
  # Used to avoid multiplicity in probabilities
  no_repeats = numeric(0)
  
  p_x_tox = c(0,0,0,0)
  for (i in 0:15){
    for (j in 0:15){
      for (k in 0:15){
        for (r in 0:15){
          if ((i <= 15 - threshold_tox[1]) & (i + j <= 30 - threshold_tox[2]) & (i + j + k <= 45 - threshold_tox[3]) & (i + j + k + r <= 60 - threshold_tox[4])){
            
            # Temp as we're adding it to two things so saves calculating it twice
            temp_prob = (pbinom(i,15,p_tox) - pbinom((i-1),15,p_tox))*(pbinom(j,15,p_tox) - pbinom(j-1,15,p_tox))*(pbinom(k,15,p_tox) - pbinom(k-1,15,p_tox))*(pbinom(r,15,p_tox) - pbinom(r-1,15,p_tox))
            
            tox_prob = tox_prob + temp_prob
            
            # Adding to probability of ending at 4
            p_x_tox[4] = p_x_tox[4] + temp_prob
            
            
            # The below series of if/elses looks complicated, but ultimately is just
            # Summing up the probabilities of the situations where it ended at 15,
            # 30, 45, and 60 without rejecting the null.
            
          } else{
            
            if (i > 15 - threshold_tox[1]){
              # If ends at 15
              
              # Checks whether we've already calculated this probability
              if(! (i %in% no_repeats)){
                
                p_x_tox[1] = p_x_tox[1] + (pbinom(i,15,p_tox) - pbinom((i-1),15,p_tox))
                
                no_repeats = append(no_repeats,i)
                
              }
            } else {
              
              
              if ((i + j > 30 - threshold_tox[2])){
                # If ends at 30
                
                #+100 added for no overlap, same for 20*, so that 1,2, is different form 2,1
                if(! ((i+20*j+100) %in% no_repeats)){
                  
                  p_x_tox[2] = p_x_tox[2] + (pbinom(i,15,p_tox) - pbinom((i-1),15,p_tox))*(pbinom(j,15,p_tox) - pbinom(j-1,15,p_tox))
                  
                  no_repeats = append(no_repeats,i+20*j+100)
                  
                }
                
              } else{
                
                if ((i + j + k > 45 - threshold_tox[3])){
                  # If ends at 45
                  
                  # +1000 added for no overlap, 20,500 (15^2 + 1) added to ensure no overlap
                  if(! ((i+20*j+500*k+1000) %in% no_repeats)){
                    
                    p_x_tox[3] = p_x_tox[3] + (pbinom(i,15,p_tox) - pbinom((i-1),15,p_tox))*(pbinom(j,15,p_tox) - pbinom(j-1,15,p_tox))*(pbinom(k,15,p_tox) - pbinom(k-1,15,p_tox))
                    
                    no_repeats = append(no_repeats,i+20*j+500*k+1000)
                    
                  }
                } else{
                  # If ends at 60
                  
                  p_x_tox[4] = p_x_tox[4] + (pbinom(i,15,p_tox) - pbinom((i-1),15,p_tox))*(pbinom(j,15,p_tox) - pbinom(j-1,15,p_tox))*(pbinom(k,15,p_tox) - pbinom(k-1,15,p_tox))*(pbinom(r,15,p_tox) - pbinom(r-1,15,p_tox))
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  # The run is successful if both eff and tox weren't rejected, so both need to be met
  prob = eff_prob*tox_prob
  
  # The vector of probabilities for ending at 15, 30, 45, 60
  p_x_joint = numeric(4)
  
  # Ends at 15
  p_x_joint[1] = p_x_eff[1] + p_x_tox[1] - p_x_eff[1]*p_x_tox[1]
  
  # Ends at 30 (doesn't end at 15 and then ends at 30)
  p_x_joint[2] = p_x_tox[2]*(1-p_x_eff[1]) + p_x_eff[2]*(1-p_x_tox[1]) - p_x_tox[2]*(p_x_eff[2])
  
  # Ends at 60 (both make it to 60)
  p_x_joint[4] = p_x_tox[4]*p_x_eff[4]
  
  # Ends at 45 (what is left) (could be made more rigorous by checking it equals the derived value like in 30's case)
  p_x_joint[3] = 1 - sum(p_x_joint)
  
  #browser()  # For testing
  
  # Expected sample size
  ESS = sum(p_x_joint*c(15,30,45,60))
  
  # Early stopping
  ES = sum(p_x_joint*c(1,1,1,0))
  
  if (least_power){
    # This returns the least power
    return(prob)
  }
  # Otherwise it calculates the FWER 
  
  
  # As each arm is i.i.d , then the FWER is just the product of them all
  # But we don't want the FWER for power of each arm
  exact_type_one_error = prob
  
  return(c("Power" = exact_type_one_error,"ESS" = ESS, "ES" = ES))
}

Grid_Search_Multi_Arm = function(target_FWER = 0.1, n_sim = 10000, exact = TRUE){
  # Does a grid search to opimtise parameters 
  
  # exact is a boolean that indicates whether we want to use the exact type one error or not
  
  pot_lamda = seq(0.6,1,0.01)
  pot_gamma = seq(0.6,1,0.01)
  
  # Need a temp_n here as we overwrite pot_lamda next line but still need its length
  temp_n = length(pot_lamda)
  
  pot_lamda = rep(pot_lamda, length(pot_lamda))
  pot_gamma = rep(pot_gamma, each = temp_n )  # Vectors for foreach to cycle through 
  
  # Control and treatment values to test against
  controls_eff = c(0.45,0.45,0.45)
  controls_tox = c(0.3,0.3,0.3)
  
  treats_eff = c(0.45,0.45,0.6)  # Slightly more efficate
  treats_tox = c(0.3,0.3,0.2)  # Slightly less toxic 
  
  # Parallilisation stuff
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  # Note that this will record the least power: the power in the Least Favourable Circumstance where only one arm is efficate
  par_power = numeric(length(pot_lamda))
  
  par_power = foreach( i = 1:length(pot_lamda), .combine = c, .export = c("Multi_arm_BOP2", "Beta_comp_exact", "Beta_comp", "type_one_error_exact")) %dopar%{
    
    #for( i in 1:length(pot_lamda)){   # For testing purposes
    # Note, the FWER and power will be examined using the paper's values
    if (exact){
      # NOTE EXACT CURRENTLY USES MAGIC NUMBERS 0.45,0.3 AND 0.6,0.2
      
      
      if ( type_one_error_exact(lambda = pot_lamda[i], gamma = pot_gamma[i], p_eff = 0.45, p_tox = 0.3) < 0.1){
        par_power[i] = (type_one_error_exact(lambda = pot_lamda[i], gamma = pot_gamma[i], p_eff = 0.6, p_tox = 0.2, least_power = TRUE))
      } else{
        par_power[i] = 0
      }
      
      
    } else{
      # NON exact version, classical
      
      FWER  = numeric(n_sim)
      
      for (j in 1:n_sim){
        FWER[j] = max(Multi_arm_BOP2(controls_eff,controls_tox, lamda = pot_lamda[i], gamma = pot_gamma[i])[[3]]) # Will record 1 if any passes
      }
      
      if (mean(FWER) <= target_FWER){  # If it meets the FWER requirements
        
        temp_power = numeric(n_sim)  # SOmething to store the successes in before passing them to par_power
        
        for (j in 1:n_sim){
          temp = Multi_arm_BOP2(treats_eff,treats_tox,lamda = pot_lamda[i], gamma = pot_gamma[i])[[3]]  # Temp so reading is easier, as we go into an if next line
          
          if ( temp[3] == 1 ){
            temp_power[j] = 1
          }
        }
        
        par_power[i] = mean(temp_power)
        
        
        
      }
      else{
        par_power[i] = 0
      }
      
    }
    
  }
  
  stopCluster(cl) # deactivates
  
  #print(par_power)
  
  max_index = which.max(par_power)
  
  to_print = c("Lamda" = pot_lamda[max_index], "Gamma" = pot_gamma[max_index], "Power" = max(par_power))
  
  print(to_print)
  
  return(to_print)
}

# Could probably also be made exact like the type one error for things like sample size and early ending propotion
OC_Multi_Arm_BOP = function(treat_eff, treat_tox, n = 60, IAn = 15, control = FALSE, gamma = 0.98, lamda = 0.81, n_sim = 10000, threshold_type = 1, RAR = FALSE, Return_FWER = FALSE, Return_Least_Power = FALSE, exact = TRUE){
  ## OC_Multi_Arm_BOP gets the operating characteristics of the Multi Arm BOP, namely FWER, rejection of null, sample size, and whether there was early stopping
  
  # treat_eff is the vector of treatments
  # treat_tox is the vector of toxicitites
  # n is the max number of samples PER ARM
  # IAN is the number of samples before an IA PER ARM
  # control determines whether trea_eff[1] and control_eff[1] are controls or not
  # gamma and lamda are the threshold parameters
  # n_sim is how many simulations we're doing
  # RAR specifies whether we're looking at the OC of the RAR version or not (NOTE, this changes the defintion of n and IAn)
  
  # Initialising vectors to be filled
  FWER_vec = numeric(n_sim) # Where we will record the FWER (sometimes not applicable)
  
  Power_vec = vector("list", length(treat_eff))  # Recording the power/type 1 error of the different arms
  
  Early_stop = vector("list", length(treat_eff)) # Recording the proportion of early stopping in each arm
  
  Sample_size = vector("list", length(treat_eff)) # Recording average sample size
  
  for (i in 1:length(treat_eff)){
    # Initialising the lists so that each item is a vector that can be easily added to 
    
    Power_vec[[i]] = numeric(1)
    Early_stop[[i]] = numeric(1)
    Sample_size[[i]] = numeric(1)
  }
  
  
  if (exact & !(RAR)){
    ## Only works for non RAR case 
    for (i in 1:length(treat_eff)){
      
      temp_OC = OC_exact(p_eff = treat_eff[i], p_tox = treat_tox[i])
      
      Power_vec[[i]] = temp_OC[1]
      Early_stop[[i]] = temp_OC[3]
      Sample_size[[i]] = temp_OC[2]
      
      
      
    }
    
    # Assuming that the 1st treatment is the control case
    
    FWER_vec = 1 - (1 - Power_vec[[i]])^3
    
      # Could potentially be parallised
  }
  else{
  
    for( i in 1:n_sim){
      
      if (RAR){
        temp_data = Multi_arm_BOP2_brar(treats_eff = treat_eff, treats_tox = treat_tox, n = n , IAn = IAn, control = control, lamda = lamda, gamma = gamma, threshold_type = threshold_type) 
      }
      else{
        temp_data = Multi_arm_BOP2(treats_eff = treat_eff, treats_tox = treat_tox, n = n , IAn = IAn, control = control, lamda = lamda, gamma = gamma, threshold_type = threshold_type)
      }
      
      if (control){
        temp_data[[c(3,1)]] = 0 # Ignoring the control in FWER calculations
      }
      
      FWER_vec[i] = max(temp_data[[3]])  # If there is a positive result, then it is a breach of FWER in the null setting
      
      for (j in 1:length(treat_eff)){
        # updating the lists 
        
        Power_vec[[c(j,i) ]] = temp_data[[c(3,j)]]
        Early_stop[[c(j,i)]] = temp_data[[c(1,j)]]
        Sample_size[[c(j,i)]] = temp_data[[c(2,j)]]
      }
      
      
    }
  
  
  }
  
  
  
  
  if (Return_FWER == TRUE){  # If we just want the FWER or least power, we can just return those
    return(mean(FWER_vec))
  }
  else if (Return_Least_Power == TRUE){
    return(mean(Power_vec[[length(treat_eff)]]))  ## NOT GENERALISED, ONLY LOOKS AT DESIGNS WHERE THE TREATMENT IS last
  }
  else{
    
    print("FWER")
    print(mean(FWER_vec))
    for (j in 1:length(treat_eff)){
      # printing results
      
      print(c("Treatment" = treat_eff[j], "Toxic" = treat_tox[j]))
      
      print(c("Power" = mean(Power_vec[[j]]))) 
      print(c("Early Stop" = mean(Early_stop[[j]])))
      print(c("Sample size" = mean(Sample_size[[j]])))
    }
    
    # Needs to be generalised somehow, probably incorportaing it into the above loop, but works for reproduction purposes
    if (control){
      to_return_list = list(c(mean(Power_vec[[1]]),mean(Power_vec[[2]]),mean(Power_vec[[3]]),mean(Power_vec[[4]] )),
                            c(mean(Early_stop[[1]]),mean(Early_stop[[2]]),mean(Early_stop[[3]]),mean(Early_stop[[4]])),
                            c(mean(Sample_size[[1]]),mean(Sample_size[[2]]),mean(Sample_size[[3]]),mean(Sample_size[[4]])))  
    }
    else{
      to_return_list = list(c(mean(Power_vec[[1]]),mean(Power_vec[[2]]),mean(Power_vec[[3]])),
                            c(mean(Early_stop[[1]]),mean(Early_stop[[2]]),mean(Early_stop[[3]])),
                            c(mean(Sample_size[[1]]),mean(Sample_size[[2]]),mean(Sample_size[[3]])))      
    }
    
  }
  
  return(to_return_list)
  
  
}

match_data_table = function(gamma = 0.98, lamda = 0.81){
  ## Aims to match the data from the BOP2 multi arm paper in table 1
  
  # Control and treatment efficacies and toxicities
  C_e = 0.45
  C_t = 0.3
  T_e = 0.6
  T_t = 0.2
  
  # And their extra efficacious scenario 5
  T_plus_e = 0.65
  
  
  # Data stores for the five situations (consdiered a for loop, but really not saving much efficiency with just 5 examples)
  temp_data_1 = OC_Multi_Arm_BOP(c(C_e,C_e,C_e),c(C_t,C_t,C_t),gamma = gamma, lamda = lamda)
  
  temp_data_2 = OC_Multi_Arm_BOP(c(C_e,C_e,T_e),c(C_t,C_t,T_t),gamma = gamma, lamda = lamda)
  
  temp_data_3 = OC_Multi_Arm_BOP(c(C_e,T_e,T_e),c(C_t,T_t,T_t),gamma = gamma, lamda = lamda)
  
  temp_data_4 = OC_Multi_Arm_BOP(c(T_e,T_e,T_e),c(T_t,T_t,T_t),gamma = gamma, lamda = lamda)
  
  temp_data_5 = OC_Multi_Arm_BOP(c(T_plus_e,T_plus_e,T_plus_e),c(T_t,T_t,T_t),gamma = gamma, lamda = lamda)
  
  
  results = data.frame("Efficacy" = c(C_e,C_e,C_e,C_e,C_e,T_e,C_e,T_e,T_e,T_e,T_e,T_e,T_plus_e,T_plus_e,T_plus_e),
                       "Toxicity" = c(C_t,C_t,C_t,C_t,C_t,T_t,C_t,T_t,T_t,T_t,T_t,T_t,T_t,T_t,T_t),
                       "Reject Null rate" = c(temp_data_1[[1]],temp_data_2[[1]],temp_data_3[[1]],temp_data_4[[1]],temp_data_5[[1]]),
                       "Early Stop rare" = c(temp_data_1[[2]],temp_data_2[[2]],temp_data_3[[2]],temp_data_4[[2]],temp_data_5[[2]]),
                       "Average sample size" = c(temp_data_1[[3]],temp_data_2[[3]],temp_data_3[[3]],temp_data_4[[3]],temp_data_5[[3]]))
  
  return(results)
}

repeated_control_prob = function(Alpha, Beta, Alpha_tox, Beta_tox, c = 0.5, control = FALSE, phi_e = 0.45, phi_t = 0.3, exact = TRUE){
  ## Similar to max_beta_calc, this function aims to generate a prob to used for RAR
  # It will calculate the probability that better than the control arm/ control value
  # Then it will calculate allocation probabilites by normalising these against each other
  
  
  # The variables are the same as referenced in other functions
  # c is the tuning parameter for BRAR 
  
  if (control){  # If we have a control arm (taking the first treat as control)
    
    n =  length(Alpha)
    
    prob_vec = numeric(n) #efficacy vector
    tox_vec = numeric(n) # toxicity vector
    
    
    
    
    if (exact){
      
      for( i in 1:n){  # Note that this also does the control arm
        
        
        prob_vec[i] = 1 - Beta_comp_exact(Alpha[1],Beta[1],Alpha[i],Beta[i])
        tox_vec[i] = 1 - Beta_comp_exact(Alpha_tox[1],Beta_tox[1],Alpha_tox[i],Beta_tox[i])
        
      }
      
    }
    else{
      
      
      for( i in 1:n){
        prob_vec[i] = 1 - Beta_comp(Alpha[1],Beta[1],Alpha[i],Beta[i], 10000)
        tox_vec[i] = 1 - Beta_comp(Alpha_tox[1],Beta_tox[1],Alpha_tox[i],Beta_tox[i],10000)
        
      }
      
      
    }
    
  }
  else{  # Otherwise
    
    n =  length(Alpha)
    
    prob_vec = numeric(n)
    tox_vec = numeric(n)
    
    for( i in 1:n){
      
      prob_vec[i] = 1 - pbeta(phi_e,Alpha[i],Beta[i])
      tox_vec[i] = pbeta(phi_t,Alpha_tox[i],Beta_tox[i])
      
    }
  }
  
  return_prob_vec = prob_vec * tox_vec  # Assuming independence, probability of desirable result 
  
  return_prob_vec = return_prob_vec^c  # tuning
  
  return_prob_vec = return_prob_vec / sum(return_prob_vec) # Normalising into probabilities
  
  
  return(return_prob_vec)
}

allo_prob_brar = function(successes, toxics, n_arms, N, current_n, control = FALSE, phi_e = 0.45, phi_t = 0.3, prior_alpha = rep(phi_e,length(successes)), prior_beta = rep(1 - phi_e,length(successes)),prior_alpha_tox = rep(phi_t,length(successes)), prior_beta_tox = rep(1 - phi_t,length(successes))){
  ## allo_prob_brar determines the allocation probabulities of a multi arm BRAR design akin to Trippa et al 2012
  # Success is a vector of successes per arm
  # n_arms is a vector of current intake per arm
  # N is the maximum sample size (overall)
  # control is whether successes[1] and n_arms[1] is a control or not
  # phi_t and phi_e are standard
  # prior_alpha and priob_beta are the priors for alpha and beta
  
  # Right now just handles single endpoint, but can be extended
  
  
  n = length(successes)
  
  return_prob = numeric(n) # the vector of probabilities that will ultimately be returned
  
  #i = sum(n_arms) # The i in i/2N for BRAR tuning replaced by current_n to account for arms dropping out 
  
  if (control){  # So first index of successes and n_arms ism control
    
    if (length(successes) == 1){  # If there is only the control arm, just return 1, no need to do calculations
      return_prob = 1
    }
    else{
      Alpha = successes[2:n]
      
      Beta = n_arms[2:n] - successes[2:n]
      
      Alpha_tox = toxics[2:n]
      
      Beta_tox = n_arms[2:n] - toxics[2:n]
      
      return_prob[2:n] = repeated_control_prob(prior_alpha[2:n] + Alpha, prior_beta[2:n] + Beta, prior_alpha_tox[2:n] + Alpha_tox, prior_beta_tox[2:n] + Beta_tox,   c = current_n/(2*N), phi_e = phi_e, phi_t = phi_t)
      
      eta = 0.25*(current_n/N)
      
      return_prob[1] = 1/(n-1) * (exp( max(n_arms[2:n]) - n_arms[1]))^eta  # As in Trippa et al 2012, making sure the control doesnt fall behind the best treatment
      
      return_prob = return_prob / sum(return_prob) # Normalising the probabilities
      
    }
    
    
    ## Error catching/debugging
    #if (NA %in% return_prob){
    #  browser()
    #  
    #}
    
  }
  
  else{
    Alpha = successes
    
    Beta = n_arms - successes
    
    Alpha_tox = toxics
    
    Beta_tox = n_arms - toxics
    
    return_prob = repeated_control_prob( prior_alpha + Alpha, prior_beta + Beta, prior_alpha_tox + Alpha_tox, prior_beta_tox + Beta_tox, c = current_n/(2*N), phi_e = phi_e, phi_t = phi_t)
  }
  
  return(return_prob)
  
  
}

Multi_arm_BOP2_brar = function(treats_eff, treats_tox, n = 180, IAn = 45, burn = IAn, control = FALSE, phi_e = 0.45, phi_t = 0.3, threshold_type = 1, lamda = 0.77, gamma = 0.94, exact = TRUE){
  ## Multi Arm BOP2 BRAR aims to extend BOP2 to the multi arm setting, like in Mulier et al 2024, but this time with BRAR using allo_prob_brar
  
  # treats_ ss a vector of treatments. If control = TRUE, then treats[1] will be the control 
  # treats_eff is the vector of efficacy (ex 0.7 is the probability of a desirable outcome)
  # treats_tox  is the vector of toxicity (ex 0.6 is the probabilty of toxic response)
  # So treats_eff = c(0.1,0.2,0.3) and treats_tox = c(0.5,0.4,0.3) correspond to 3 treatments
  # One with 0.1 success 0.5 tox,then 0.2 success 0.4 tox, then 0.3 success and tox
  
  # n is now the maxmimum number of patients
  # IAn is after how many patients to perform an interim analysis
  # burn is how mnay patients there are before BRAR begins
  # The default values are the ones they used for their testing
  # control = TRUE/FALSE determines whether there is a control or not. If TRUE, then treats[1] is a control
  # phi_e, phi_f are the efficacy/futility thresholds for the uncontrolled setting
  # lamda are gamma are parameters that are used to control the type one error whilst maximising power
  # return_fwer and return_least_power determine whether you want those exclusively returned or not
  
  
  
  # For ease of calculation, we will take a Dirichlet(1,1,1,1) prior 
  # (this allows us to use Beta_comp_exact)
  
  
  num_treat = length(treats_eff)
  continue = rep(1,num_treat)  # Will be 1 if we continue with a treatment, 0 futile (no test for efficacy)
  
  # Where the responses for efficacy and toxicity will go
  Responses_eff = vector("list",num_treat)
  Responses_tox = vector("list",num_treat)
  
  ## PBD randomisastion 
  
  PBD_rando = rep(1:num_treat, ceiling(burn/num_treat))
  
  PBD_rando = sample(PBD_rando)  # Randomisation of burn in period done via PBD
  
  for( j in 1:length(PBD_rando)){
    
    i = PBD_rando[j] # Which treatment we're allocating to currently
    
    Responses_eff[[i]] = append(Responses_eff[[i]],rbinom(1, 1, treats_eff[i]))
    Responses_tox[[i]] = append(Responses_tox[[i]],rbinom(1, 1, treats_tox[i]))
  }
  
  
  if (control == FALSE){
    
    if (threshold_type == 1){  # There are varying threshold types the paper gives
      
      #Threshold type 1 is most similar to the original BOP2 paper. 
      
      threshold = 1 - (lamda * ((IAn / n)^gamma))
    }
    
    Alpha = rep(0,num_treat)
    
    Beta = rep(0,num_treat)
    
    Alpha_tox = rep(0,num_treat)
    
    Beta_tox = rep(0,num_treat)
    
    
    for (i in 1:num_treat){    # Keep this, standard arm
      
      ## Temp values to do the beta calcs with 
      temp_n = length(Responses_eff[[i]])
      temp_eff = sum(Responses_eff[[i]])
      temp_tox = sum(Responses_tox[[i]])
      
      Alpha[i] = temp_eff 
      Beta[i] = temp_n - temp_eff  # Alpha and Beta useful for later BRAR calculations
      
      Alpha_tox[i] = temp_tox 
      Beta_tox[i] = temp_n - temp_tox  # Alpha and Beta useful for later BRAR calculations
      
      
      #print(threshold)
      #print(temp_n)
      #print(temp_eff)
      #print(temp_tox)
      
      
      if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
        
        continue[i] = 0  # This has either exceeded futility or toxicity
        
      }
      
    }
    
    remain_IAn = IAn
    
    remain = n - length(PBD_rando)
    
    ## If we want, we can either update the allocation probabilities fully sequentially, or whenever there is an interim analysis, after each interim for now
    
    current_n = n - remain
    
    Alpha = Alpha[as.logical(continue)] # Only looking at eligible arms
    
    Beta = Beta[as.logical(continue)]
    
    Alpha_tox = Alpha_tox[as.logical(continue)] # Only looking at eligible arms
    
    Beta_tox = Beta_tox[as.logical(continue)]
    
    vec_prob = allo_prob_brar(successes = Alpha, toxics = Alpha_tox, n_arms = Beta + Alpha, N = n, current_n =  current_n, control = control)
    
    cdf_prob = cumsum(vec_prob)
    
    #print(cdf_prob)
    
    
    while(remain > 0){
      
      if (max(continue) == 0){
        break
      }
      ## Continuing simulating using BRAR now
      
      temp_p = runif(1)
      
      simu_target = which(cdf_prob > temp_p)[1]  # Which should be simulated
      
      sum_continue = cumsum(continue)
      
      simu_target = which(sum_continue == simu_target)[1] ## Making sure the target is the correct active arm
      
      ##Testing
      #browser()
      
      #print(Responses_eff[[simu_target]])
      #print(simu_target)
      
      Responses_eff[[simu_target]] = append(Responses_eff[[simu_target]],rbinom(1, 1, treats_eff[simu_target]))
      Responses_tox[[simu_target]] = append(Responses_tox[[simu_target]],rbinom(1, 1, treats_tox[simu_target]))
      
      remain_IAn = remain_IAn - 1
      
      remain =  remain - 1
      
      if (remain_IAn == 0 & remain != 0){   # Interim analysis, if remain = 0 skips to final analysis
        
        remain_IAn = IAn
        
        for (i in 1:num_treat){  # Looking at every treatment
          
          if (continue[i] == 1){  # We only want to analyse the arms being updated
            
            ## Temp values to do the beta calcs with 
            temp_n = length(Responses_eff[[i]])
            temp_eff = sum(Responses_eff[[i]])
            temp_tox = sum(Responses_tox[[i]])
            
            
            Alpha[i] = temp_eff 
            Beta[i] = temp_n - temp_eff  # Alpha and Beta useful for later BRAR calculations
            
            Alpha_tox[i] = temp_tox 
            Beta_tox[i] = temp_n - temp_tox  # Alpha and Beta useful for later BRAR calculations
            
            
            if (threshold_type == 1){
              
              
              threshold = 1 - (lamda * (((n-remain) / n)^gamma))
            }
            
            #print(threshold)
            #print(temp_n)
            #print(temp_eff)
            #print(temp_tox)
            
            
            
            
            if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
              
              continue[i] = 0  # This has either exceeded futility or toxicity
              
            }
            
          } 
        }
        
        ## Update allocation probabilities 
        Alpha = Alpha[as.logical(continue)] # Only looking at eligible arms
        
        Beta = Beta[as.logical(continue)]
        
        Alpha_tox = Alpha_tox[as.logical(continue)] # Only looking at eligible arms
        
        Beta_tox = Beta_tox[as.logical(continue)]
        
        vec_prob = allo_prob_brar(successes = Alpha, toxics = Alpha_tox, n_arms = Beta + Alpha, N = n, current_n =  current_n, control =  control)
        
        cdf_prob = cumsum(vec_prob)
        
        #print(cdf_prob)
        
      }
      
      
      
      
    }
    
    early_stop = rep(TRUE, num_treat)  # Recording whether a treatment stopped early
    
    final_sample_size = numeric(num_treat) # Recording the samples size of each treatment
    
    # One final test for efficacy/toxicity
    
    for (i in 1:num_treat){  # Looking at every treatment
      
      final_sample_size[i] = length(Responses_eff[[i]])
      
      if (continue[i] == 1){  # We only want to analyse the arms being updated
        
        early_stop[i] = FALSE # If it has made it to the final analysis, it has not stopped early
        
        ## Temp values to do the beta calcs with 
        temp_n = length(Responses_eff[[i]])
        temp_eff = sum(Responses_eff[[i]])
        temp_tox = sum(Responses_tox[[i]])
        
        #print(temp_eff)
        #print(temp_tox)
        
        
        if (threshold_type == 1){  # There are varying threshold types the paper gives
          
          #Threshold type 1 is most similar to the original BOP2 paper. 
          
          threshold = 1 - lamda 
        }
        
        
        #print(threshold)
        #print(temp_n)
        #print(temp_eff)
        #print(temp_tox)
        #browser()
        
        if ((pbeta(phi_e, 0.45 + temp_eff, 0.55 + temp_n - temp_eff) > threshold) | (pbeta(phi_t, 0.3 + temp_tox, 0.7 + temp_n - temp_tox) < 1 - threshold)  ){
          
          continue[i] = 0  # This has either exceeded futility or toxicity
          
        }
      }  
    }
    
    
  }
  
  else{
    
    if (threshold_type == 1){  # There are varying threshold types the paper gives
      
      #Threshold type 1 is most similar to the original BOP2 paper. 
      
      threshold = 1 - (lamda * ((IAn / n)^gamma))
    }
    
    Alpha = rep(0,num_treat)
    
    Beta = rep(0,num_treat)
    
    Alpha_tox = rep(0,num_treat)
    
    Beta_tox = rep(0,num_treat)
    
    con_n = length(Responses_eff[[1]])  ## Control results
    con_eff = sum(Responses_eff[[1]])
    con_tox = sum(Responses_tox[[1]])
    
    ## Manually filling these in as the control doesn't get compared itself
    Alpha[1] = con_eff
    Beta[1] = con_n - con_eff
    
    Alpha_tox[1] = con_tox
    Beta_tox[1] = con_n - con_eff
    
    for (i in 2:num_treat){    # Keep this, standard arm
      
      ## Temp values to do the beta calcs with 
      temp_n = length(Responses_eff[[i]])
      temp_eff = sum(Responses_eff[[i]])
      temp_tox = sum(Responses_tox[[i]])
      
      Alpha[i] = temp_eff 
      Beta[i] = temp_n - temp_eff  # Alpha and Beta useful for later BRAR calculations
      
      Alpha_tox[i] = temp_tox 
      Beta_tox[i] = temp_n - temp_tox  # Alpha and Beta useful for later BRAR calculations
      
      
      #print(threshold)
      #print(temp_n)
      #print(temp_eff)
      #print(temp_tox)
      
      
      if (exact){
        if ( Beta_comp_exact( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
          
          continue[i] = 0  # This has either exceeded futility or toxicity
          
        }
        
      }
      else{
        if ( Beta_comp( 0.45 + temp_eff, 0.55 + temp_n - temp_eff, 0.45 + con_eff, 0.55 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 0.3 + temp_tox, 0.7 + temp_n - temp_tox, 0.3 + con_tox, 0.7 + con_n - con_tox) > threshold ){
          
          continue[i] = 0  # This has either exceeded futility or toxicity
          
        }
        
      }
      
    }
    
    remain_IAn = IAn
    
    remain = n - length(PBD_rando)
    
    ## If we want, we can either update the allocation probabilities fully sequentially, or whenever there is an interim analysis, after each interim for now
    
    current_n = n - remain
    
    Alpha = Alpha[as.logical(continue)] # Only looking at eligible arms
    
    Beta = Beta[as.logical(continue)]
    
    Alpha_tox = Alpha_tox[as.logical(continue)] # Only looking at eligible arms
    
    Beta_tox = Beta_tox[as.logical(continue)]
    
    vec_prob = allo_prob_brar(successes = Alpha, toxics = Alpha_tox, n_arms = Beta + Alpha, N = n, current_n =  current_n, control = control)
    
    #print(vec_prob)
    
    cdf_prob = cumsum(vec_prob)
    
    #print(cdf_prob)
    
    
    while(remain > 0){
      
      if (max(continue[2:num_treat]) == 0){
        break
      }
      ## Continuing simulating using BRAR now
      
      temp_p = runif(1)
      
      simu_target = which(cdf_prob > temp_p)[1]  # Which should be simulated
      
      sum_continue = cumsum(continue)
      
      simu_target = which(sum_continue == simu_target)[1] ## Making sure the target is the correct active arm
      
      ##Testing
      
      #print(Responses_eff[[simu_target]])
      #print(simu_target)
      
      #browser()
      
      Responses_eff[[simu_target]] = append(Responses_eff[[simu_target]],rbinom(1, 1, treats_eff[simu_target]))
      Responses_tox[[simu_target]] = append(Responses_tox[[simu_target]],rbinom(1, 1, treats_tox[simu_target]))
      
      remain_IAn = remain_IAn - 1
      
      remain =  remain - 1
      
      if (remain_IAn == 0 & remain != 0){   # Interim analysis, if remain = 0 skips to final analysis
        
        remain_IAn = IAn
        
        con_n = length(Responses_eff[[1]])  ## Control results
        con_eff = sum(Responses_eff[[1]])
        con_tox = sum(Responses_tox[[1]])
        
        ## Manually filling these in as the control doesn't get compared itself
        Alpha[1] = con_eff
        Beta[1] = con_n - con_eff
        
        Alpha_tox[1] = con_tox
        Beta_tox[1] = con_n - con_eff
        
        for (i in 2:num_treat){  # Looking at every treatment
          
          if (continue[i] == 1){  # We only want to analyse the arms being updated
            
            ## Temp values to do the beta calcs with 
            temp_n = length(Responses_eff[[i]])
            temp_eff = sum(Responses_eff[[i]])
            temp_tox = sum(Responses_tox[[i]])
            
            
            Alpha[i] = temp_eff 
            Beta[i] = temp_n - temp_eff  # Alpha and Beta useful for later BRAR calculations
            
            Alpha_tox[i] = temp_tox 
            Beta_tox[i] = temp_n - temp_tox  # Alpha and Beta useful for later BRAR calculations
            
            
            if (threshold_type == 1){
              
              
              threshold = 1 - (lamda * (((n-remain) / n)^gamma))
            }
            
            #print(threshold)
            #print(temp_n)
            #print(temp_eff)
            #print(temp_tox)
            
            
            
            
            if (exact){
              if ( Beta_comp_exact( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
                
                continue[i] = 0  # This has either exceeded futility or toxicity
                
              }
              
            }
            else{
              if ( Beta_comp( 0.45 + temp_eff, 0.55 + temp_n - temp_eff, 0.45 + con_eff, 0.55 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 0.3 + temp_tox, 0.7 + temp_n - temp_tox, 0.3 + con_tox, 0.7 + con_n - con_tox) > threshold ){
                
                continue[i] = 0  # This has either exceeded futility or toxicity
                
              }
              
            }
            
          } 
        }
        
        ## Update allocation probabilities 
        Alpha = Alpha[as.logical(continue)] # Only looking at eligible arms
        
        Beta = Beta[as.logical(continue)]
        
        Alpha_tox = Alpha_tox[as.logical(continue)] # Only looking at eligible arms
        
        Beta_tox = Beta_tox[as.logical(continue)]
        
        vec_prob = allo_prob_brar(successes = Alpha, toxics = Alpha_tox, n_arms = Beta + Alpha, N = n, current_n =  current_n, control =  control)
        
        #print(vec_prob)
        
        cdf_prob = cumsum(vec_prob)
        
        #print(cdf_prob)
        
      }
      
      
      
      
    }
    
    early_stop = rep(TRUE, num_treat)  # Recording whether a treatment stopped early
    
    final_sample_size = numeric(num_treat) # Recording the samples size of each treatment
    
    # One final test for efficacy/toxicity
    
    con_n = length(Responses_eff[[1]])  ## Control results
    con_eff = sum(Responses_eff[[1]])
    con_tox = sum(Responses_tox[[1]])
    
    ## Manually filling these in as the control doesn't get compared itself
    Alpha[1] = con_eff
    Beta[1] = con_n - con_eff
    
    Alpha_tox[1] = con_tox
    Beta_tox[1] = con_n - con_eff
    
    final_sample_size[1] = length(Responses_eff[[1]])
    
    early_stop[1] = FALSE
    
    for (i in 2:num_treat){  # Looking at every treatment
      
      final_sample_size[i] = length(Responses_eff[[i]])
      
      if (continue[i] == 1){  # We only want to analyse the arms being updated
        
        early_stop[i] = FALSE # If it has made it to the final analysis, it has not stopped early
        
        ## Temp values to do the beta calcs with 
        temp_n = length(Responses_eff[[i]])
        temp_eff = sum(Responses_eff[[i]])
        temp_tox = sum(Responses_tox[[i]])
        
        #print(temp_eff)
        #print(temp_tox)
        
        
        if (threshold_type == 1){  # There are varying threshold types the paper gives
          
          #Threshold type 1 is most similar to the original BOP2 paper. 
          
          threshold = 1 - lamda 
        }
        
        
        #print(threshold)
        #print(temp_n)
        #print(temp_eff)
        #print(temp_tox)
        #browser()
        
        if (exact){
          if ( Beta_comp_exact( 1 + temp_eff, 1 + temp_n - temp_eff, 1 + con_eff, 1 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 1 + temp_tox, 1 + temp_n - temp_tox, 1 + con_tox, 1 + con_n - con_tox) > threshold ){
            
            continue[i] = 0  # This has either exceeded futility or toxicity
            
          }
          
        }
        else{
          if ( Beta_comp( 0.45 + temp_eff, 0.55 + temp_n - temp_eff, 0.45 + con_eff, 0.55 + con_n - con_eff) < 1 - threshold  |  Beta_comp( 0.3 + temp_tox, 0.7 + temp_n - temp_tox, 0.3 + con_tox, 0.7 + con_n - con_tox) > threshold ){
            
            continue[i] = 0  # This has either exceeded futility or toxicity
            
          }
          
        }
      }  
    }
    
    
  }
  
  # Results 
  #print(Responses_eff)
  #print(Responses_tox)  
  
  to_return = list( early_stop, final_sample_size, continue ) # The three variables of interst
  # Need these to replicate their results
  
  return(to_return)
  
}

Grid_Search_Multi_Arm_brar = function(target_FWER = 0.1, n_sim = 1000){
  # Does a grid search to opimtise parameters , this time for the brar version 
  
  pot_lamda = seq(0.6,1,0.01)
  pot_gamma = seq(0.6,1,0.01)
  
  
  # Need a temp_n here as we overwrite pot_lamda next line but still need its length
  temp_n = length(pot_lamda)
  
  pot_lamda = rep(pot_lamda, length(pot_lamda))
  pot_gamma = rep(pot_gamma, each = temp_n )  # Vectors for foreach to cycle through 
  
  # Control and treatment values to test against
  controls_eff = c(0.45,0.45,0.45)
  controls_tox = c(0.3,0.3,0.3)
  
  treats_eff = c(0.45,0.45,0.6)  # Slightly more efficate
  treats_tox = c(0.3,0.3,0.2)  # Slightly less toxic 
  
  # Parallilisation stuff
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  # Note that this will record the least power: the power in the Least Favourable Circumstance where only one arm is efficate
  par_power = numeric(length(pot_lamda))
  
  par_power = foreach( i = 1:length(pot_lamda), .combine = c, .export = c("Multi_arm_BOP2_brar", "Beta_comp_exact","allo_prob_brar","repeated_control_prob")) %dopar%{
    
    #for( i in 1:length(pot_lamda)){   # For testing purposes
    # Note, the FWER and power will be examined using the paper's values
    
    FWER  = numeric(n_sim)
    
    for (j in 1:n_sim){
      FWER[j] = max(Multi_arm_BOP2_brar(controls_eff,controls_tox, lamda = pot_lamda[i], gamma = pot_gamma[i])[[3]]) # Will record 1 if any passes
    }
    
    if (mean(FWER) <= target_FWER){  # If it meets the FWER requirements
      
      temp_power = numeric(n_sim)  # SOmething to store the successes in before passing them to par_power
      
      for (j in 1:n_sim){
        temp = Multi_arm_BOP2_brar(treats_eff,treats_tox,lamda = pot_lamda[i], gamma = pot_gamma[i])[[3]]  # Temp so reading is easier, as we go into an if next line
        
        if ( temp[3] == 1){  # Least power
          temp_power[j] = 1
        }
      }
      
      par_power[i] = mean(temp_power)
      
      
      
    }
    else{
      par_power[i] = 0
    }
    
  }
  
  stopCluster(cl) # deactivates
  
  print(par_power)
  
  max_index = which.max(par_power)
  
  to_print = c("Lamda" = pot_lamda[max_index], "Gamma" = pot_gamma[max_index])
  
  print(to_print)
  
  return(to_print)
}

data_comparisons = function(n_sim = 10000, control = FALSE, lambda_BR = 0.81, gamma_BR = 0.98, lambda_BRAR = 0.77, gamma_BRAR = 0.94){
  ## Data comparisons will look at the first situations in the BOP2 multi arm trial and replicate them, similar to Multi_arm_BOP2_data , but this time
  ## with a better framework for scalability, incoroporating parallelisation and the BRAR method for comparison
  
  # testing_list contains the treat/tox values we want to test, can easily be scaled for different situations
  testing_list = vector("list",5)  
  
  testing_list[[1]] = list(c(0.45,0.3),c(0.45,0.3),c(0.45,0.3))
  
  testing_list[[2]] = list(c(0.45,0.3),c(0.45,0.3),c(0.6,0.2))
  
  testing_list[[3]] = list(c(0.45,0.3),c(0.6,0.2),c(0.6,0.2))
  
  testing_list[[4]] = list(c(0.6,0.2),c(0.6,0.2),c(0.6,0.2))
  
  testing_list[[5]] = list(c(0.65,0.2),c(0.65,0.2),c(0.65,0.2))
  
  
  # Parallilisation stuff
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  n_tests = length(testing_list)  # How many sets of testing parameters we're doing
  
  par_results = vector("list", n_tests) # For each set we do a BRAR and non BRAR
  
  par_results = foreach( i = 1:n_tests, .export = c("Multi_arm_BOP2_brar", "Beta_comp_exact","allo_prob_brar","repeated_control_prob","OC_Multi_Arm_BOP","Multi_arm_BOP2","OC_exact")) %dopar% {
    
    # getting the treatment and toxics
    temp_trt = c(testing_list[[c(i,1,1)]],testing_list[[c(i,2,1)]],testing_list[[c(i,3,1)]])
    temp_tox = c(testing_list[[c(i,1,2)]],testing_list[[c(i,2,2)]],testing_list[[c(i,3,2)]])
    
    # gamma = 0.98, lamda = 0.81 for Equal Randomisation 
    temp_results_ER = OC_Multi_Arm_BOP(temp_trt,temp_tox,gamma = gamma_BR, lamda = lambda_BR)
    
    # lamda = 0.77, gamma = 0.94 for BRAR
    temp_results_BRAR = OC_Multi_Arm_BOP(temp_trt,temp_tox, n = 180, IAn = 45, lamda = lambda_BRAR, gamma = gamma_BRAR, RAR = TRUE)
    
    to_return = list(temp_results_ER,temp_results_BRAR)  # So we can return it outside of the par loop
    
    par_results[[i]] = to_return
    
  }
  
  
  print(par_results)
  
  return(par_results)
}

heat_map_maker_par = function(RAR = FALSE, control = FALSE, n_sim = 1000, lam_low = 0.1, lam_high = 1, lam_res = 0.1, gam_low = 0.1, gam_high = 1, gam_res = 0.1){
  ## Identical to heat_map_maker but implements parallilisation in an attempt to speed things up
  
  ## heat_map_maker aims to make a heat map of the power of the power of the trial, to hopefully make optimisation easier
  
  # lam_low,high, and res are the lower bound, higher bound, and the jumps for lamdba in the sequences
  # It is the same for gamma 
  
  #Values we want to look at 
  lambda = seq(lam_low,lam_high,lam_res)
  gamma = seq(gam_low,gam_high,gam_res)
  
  # So we only need one for loop#
  temp_l = length(lambda)  # Have to get this now before we go to gamma and lambda has already changed
  lambda = rep(lambda,length(gamma))
  gamma = rep(gamma, each = temp_l)
  
  power_hm = numeric(length(lambda)) # Where the powers will be going
  
  
  # Example data for the uncontrolled/controlled situation
  if (control){
    
    controls_eff = c(0.6,0.6,0.6,0.6)
    controls_tox = c(0.4,0.4,0.4,0.4)
    
    treats_eff = c(0.6,0.6,0.6,0.75)  # Slightly more efficate
    treats_tox = c(0.4,0.4,0.4,0.3)  # Slightly less toxic
    
    
  }
  else{
    
    controls_eff = c(0.45,0.45,0.45)
    controls_tox = c(0.3,0.3,0.3)
    
    treats_eff = c(0.45,0.45,0.6)  # Slightly more efficate
    treats_tox = c(0.3,0.3,0.2)  # Slightly less toxic 
    
  }
  
  #Parallel stuff
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  
  power_hm = foreach(i = 1:length(lambda), .combine = c, .export = c("OC_Multi_Arm_BOP","Multi_arm_BOP2", "Beta_comp_exact", "Beta_comp","Multi_arm_BOP2_brar")) %dopar%  {
    
    
    
    temp_FWER = OC_Multi_Arm_BOP(controls_eff,controls_tox, control = control, n_sim = n_sim, RAR = RAR, gamma = gamma[i], lamda = lambda[i], Return_FWER = TRUE)
    
    if (temp_FWER < 0.1){
      return( OC_Multi_Arm_BOP(treats_eff,treats_tox, control = control, n_sim = n_sim, RAR = RAR, gamma = gamma[i], lamda = lambda[i], Return_Least_Power = TRUE))
    }
    else{
      return(0)
    }
    
  }
  
  stopCluster(cl) # deactivates
  
  heat_data = data.frame(lambda,gamma,power_hm)
  
  max_index = which.max(heat_data$power_hm)
  
  var_subtitle = sprintf("Control is %s RAR is %f",control, RAR)
  
  ggp <- ggplot(heat_data, aes(lambda, gamma)) +                           # Create heatmap with ggplot2
    geom_tile(aes(fill = power_hm)) +
    labs(
      title = "Power heatmap",
      subtitle = var_subtitle
    )
  print(ggp)
  
  to_return = c("Lambda" = heat_data$lambda[max_index], "Gamma" = heat_data$gamma[max_index], "Power" = heat_data$power_hm[max_index] )
  
  print(to_return)
  
  return(to_return)
}


### Calling functions to optimise parameters
# Note that as the threshold probabilities are essentially discretised in outcomes
# (there is only as many outcomes as there are binomial realisaitions)
# Then the optimum for some of these may be non-unique as they don't actually
# Change the threshold

if (Opt_BR){
  temp_optimal_BR = Grid_search_par(0.2,0.4,0.1,n_sim = 10000)
  lambda_BR = temp_optimal_BR[1]
  gamma_BR = temp_optimal_BR[2]
} else{
  # What they were on my run through
  lambda_BR = 0.91
  gamma_BR = 0.93
}


if (Opt_BRAR){
  temp_optimal_BRAR = Grid_search_par(0.2,0.4,0.1,n_sim = 10000, RAR = TRUE)
  lambda_BRAR = temp_optimal_BR[1]
  gamma_BRAR = temp_optimal_BR[2]
} else{
  # What they were on my run through
  
  lambda_BRAR = 0.9
  gamma_BRAR = 0.86
}


if (Opt_BR_Multi_Arm){
  
  temp_results = Grid_Search_Multi_Arm()
  
  lambda_Mult = temp_results[1]
  gamma_Mult = temp_results[2]
  
} else{
  # What they were on my run through
  
  lambda_Mult = 0.78
  gamma_Mult = 0.9
  
}


if (Opt_BRAR_Multi_Arm){
  
  # Note this has a lower n_sim than other functions due to the complexity of the Multi Arm BRAR
  temp_results = Grid_Search_Multi_Arm_brar()
  
  lambda_Mult_BRAR = temp_results[1]
  gamma_Mult_BRAR = temp_results[2]
  
} else{
  # What they were on my run through
  
  lambda_Mult_BRAR = 0.75
  gamma_Mult_BRAR = 0.88
  
}


### Calling functions to produce figures

## Figures 2+4
# This returns a table which directly matches the futility, futility+efficacy data tables
# in the two arm paper

if(Fig_2_4){
  set.seed(2025)
  
  Comparison_Results = BOP2_2ARM_data_parallel(gamma = gamma_BR, lamda = lambda_BR)
  
  
  ## Table making
  
  kable(Comparison_Results[[1]], digits = 3, format = "html", row.names = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover"),
                  full_width = T,
                  font_size = 15,
                  position = "left") 
  
  kable(Comparison_Results[[2]], digits = 3, format = "html", row.names = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover"),
                  full_width = T,
                  font_size = 15,
                  position = "left") 
  
  
  
}



## Figure 5
# Bar chart showcasing the effect that RAR can have (NOTE: Illustrative example, seed may be curated)

if(Fig_5){
  set.seed(2025)
  
  
  Sim_trial_RAR_PBD(0.2,0.4,80,20, gamma = gamma_BRAR, lamda = lambda_BRAR)
}

## Figure 6+7
# Returns graphs showcasing proportion to treatment and power against control (Alongside other graphs)
# FIGURE 5 IN NEW PAPER 

if(Fig_6_7){
  # Parallel so needs this seed setting indeed
  registerDoRNG(2025)
  
  Complete_summary_sim_trial_RAR(gamma = gamma_BRAR, lamda = lambda_BRAR, PBD = TRUE, cond = FALSE)
  
  
}

## Figure 8
# Data matching the multi arm Mulier paper (Table)

if(Fig_8){
  set.seed(2025)
  
  match_data_table( gamma = gamma_Mult, lamda = lambda_Mult)
  
}

## Figure 9 
# Considering multi-arm BRAR + BOP2 data table

if(Fig_9){
  set.seed(2025)
  
  results_for_table = data_comparisons(lambda_BR = lambda_Mult, gamma_BR = gamma_Mult, lambda_BRAR = lambda_Mult_BRAR, gamma_BRAR = gamma_Mult_BRAR)
  
  treats = c("A(0.45,0.3)","B(0.45,0.3)","C(0.45,0.3)",
             "A(0.45,0.3)","B(0.45,0.3)","C(0.6,0.2)",
             "A(0.45,0.3)","B(0.6,0.2)","C(0.6,0.2)",
             "A(0.6,0.2)","B(0.6,0.2)","C(0.6,0.2)",
             "A(0.65,0.2)","B(0.65,0.2)","C(0.65,0.2)"
  )
  
  # data that will be filled in from the results_for_table to make a dataframe
  power = numeric(15)
  power_brar = numeric(15)
  
  sample_size = numeric(15)
  sample_size_brar = numeric(15)
  
  early_stop = numeric(15)
  early_stop_brar = numeric(15)
  
  prop_treat = numeric(15)
  prop_treat_brar = numeric(15)
  
  # decoding the results into a dataframe
  for( i in 1:5){
    power[3*i-2] = results_for_table[[c(i,1,1,1)]]
    power[3*i-1] = results_for_table[[c(i,1,1,2)]]
    power[3*i] = results_for_table[[c(i,1,1,3)]]
    
    power_brar[3*i-2] = results_for_table[[c(i,2,1,1)]]
    power_brar[3*i-1] = results_for_table[[c(i,2,1,2)]]
    power_brar[3*i] = results_for_table[[c(i,2,1,3)]]
    
    
    early_stop[3*i-2] = results_for_table[[c(i,1,2,1)]]
    early_stop[3*i-1] = results_for_table[[c(i,1,2,2)]]
    early_stop[3*i] = results_for_table[[c(i,1,2,3)]]
    
    early_stop_brar[3*i-2] = results_for_table[[c(i,2,2,1)]]
    early_stop_brar[3*i-1] = results_for_table[[c(i,2,2,2)]]
    early_stop_brar[3*i] = results_for_table[[c(i,2,2,3)]]
    
    
    sample_size[3*i-2] = results_for_table[[c(i,1,3,1)]]
    sample_size[3*i-1] = results_for_table[[c(i,1,3,2)]]
    sample_size[3*i] = results_for_table[[c(i,1,3,3)]]
    
    prop_treat[3*i - 2] = sample_size[3*i-2]/(sample_size[3*i-2] + sample_size[3*i-1]+ sample_size[3*i])
    prop_treat[3*i - 1] = sample_size[3*i-1]/(sample_size[3*i-2] + sample_size[3*i-1]+ sample_size[3*i])
    prop_treat[3*i] = sample_size[3*i]/(sample_size[3*i-2] + sample_size[3*i-1]+ sample_size[3*i])
    
    
    sample_size_brar[3*i-2] = results_for_table[[c(i,2,3,1)]]
    sample_size_brar[3*i-1] = results_for_table[[c(i,2,3,2)]]
    sample_size_brar[3*i] = results_for_table[[c(i,2,3,3)]]
    
    prop_treat_brar[3*i - 2] = sample_size_brar[3*i-2]/(sample_size_brar[3*i-2] + sample_size_brar[3*i-1]+ sample_size_brar[3*i])
    prop_treat_brar[3*i - 1] = sample_size_brar[3*i-1]/(sample_size_brar[3*i-2] + sample_size_brar[3*i-1]+ sample_size_brar[3*i])
    prop_treat_brar[3*i] = sample_size_brar[3*i]/(sample_size_brar[3*i-2] + sample_size_brar[3*i-1]+ sample_size_brar[3*i])
    
  }
  
  
  
  
  BRAR_comp_df = data.frame("Treatments" = treats,
                            "Power" = power,
                            "Early Stop proportion" = early_stop,
                            "Sample Size" = sample_size,
                            "Prop to treatmemt" = prop_treat,
                            "Power BRAR" = power_brar,
                            "Early Stop proportion BRAR" = early_stop_brar,
                            "Sample Size BRAR" = sample_size_brar,
                            "Prop to treatment BRAR" = prop_treat_brar)
  
  #Making the table
  kable(BRAR_comp_df, digits = 3, format = "html", row.names = TRUE) %>%
    kable_styling(bootstrap_options = c("hover"),
                  full_width = T,
                  font_size = 15,
                  position = "left") %>%
    row_spec(c(1:3,7:9,13:15), background = "grey", color = "white") %>%
    row_spec(6, bold =  T) %>%
    column_spec(c(3,7), background = "red", color = "white") %>%
    column_spec(c(4,8), background = "green", color = "white") %>%
    column_spec(c(5,9), background = "blue", color = "white") %>%
    column_spec(c(6,10), background = "purple", color = "white")
  
}

## Figure 10
# Sensitivity analysis.

if(Fig_10){
  # New figure 10, showing fitting rather than just flat modifications
  # Parallel so needs this instead
  registerDoRNG(2025)
  
  Optimal_list = optimal_vec(n_sim = 10000)
  opt_lambda = Optimal_list[[1]]
  opt_gamma = Optimal_list[[2]]
  
  type_one_error_heatmap(opt_gamma, opt_lambda, n_sim = 10000)
  
  Optimal_list_rar = optimal_vec(n_sim = 10000, RAR = TRUE, return_power = TRUE)
  
  # Will take a few hours to run, will let run overnight
  
  opt_lambda = Optimal_list_rar[[1]]
  opt_gamma = Optimal_list_rar[[2]]
  
  type_one_error_heatmap(opt_gamma, opt_lambda, n_sim = 10000, RAR = TRUE)
  
  
}


# TABLE 2 IN PAPER
if (New_Figure_8){
  ## This will just make a table comparing the BRAR and BR results
  # Control will be 0.2, treatment will range from 0.1 to 0.4
  
  set.seed(2025)
  
  treat = seq(0.1,0.4,0.1) # List of treatments
  
  Power_BR = numeric(4)
  ESS_BR = numeric(4)
  Prop_BR = rep(0.5,4)  # We already know BR is 50/50
  
  Power_BRAR = numeric(4)
  ESS_BRAR = numeric(4)
  Prop_BRAR = rep(0,4)
  
  for (i in 1:4){
    
    # Temp holders of the data so they can be added to the right place
    temp_BR = Sim_trial_repeats(0.2,i*0.1,40,10, lamda = lambda_BR, gamma = gamma_BR, n_sim = 10000)
    temp_BRAR = OC_sim_trial_RAR(0.2,i*0.1,80,20, lamda = lambda_BRAR, gamma = gamma_BRAR, Sim_n = 10000)
    
    Power_BR[i] = temp_BR[1]
    ESS_BR[i] = temp_BR[2]
    
    Power_BRAR[i] = temp_BRAR[1]
    ESS_BRAR[i] = temp_BRAR[6]
    Prop_BRAR[i] = temp_BRAR[4]
    
  
  }
  
  data_Results = data.frame("Treatments" = treat, 
                            "Power BR" = Power_BR, "ESS BR" = ESS_BR, "Prop BR" = Prop_BR,
                            "Power BRAR" = Power_BRAR, "ESS_BRAR" = ESS_BRAR, "Prop BRAR" = Prop_BRAR)
  
  #Making the table
  kable(data_Results, digits = 3, format = "html", row.names = TRUE) %>%
    kable_styling(bootstrap_options = c("hover"),
                  full_width = T,
                  font_size = 15,
                  position = "left") %>%
    row_spec(c(1,3), background = "grey", color = "white")
  
  
}

if (Appendix_2){
  # This will just make a plot of the threshold values. 
  # No randomness so no seed
  
  n = 80 # Total sample size
  fut_val = numeric(n)
  suc_val = numeric(n)
  
  for ( i in 1:n){
    
    temp_threshold = Threshold_prob(lambda_BR, gamma_BR, n, i)
    
    fut_val[i] = temp_threshold[1]
    suc_val[i] = temp_threshold[2]
    
  }
  
  gg_dataframe = data.frame("Patients" = 1:n,fut_val,suc_val)
  
  threshold_plot = ggplot(gg_dataframe, aes(  x = Patients)) +
    geom_line(y = fut_val, colour = "red") +
    geom_line(y = suc_val, colour = "blue") +
    labs(
      title = "Threshold probabilities",
      y = "Probability needed to meet threshold",
      x = "Patients out of 80",
      subtitle = "0.91 Lambda, 0.93 Gamma. Blue is superiority. Red is futility."
    ) +
    scale_y_continuous(breaks = seq(0,1,0.1), limits=c(0, 1)) +
    guides(colour = guide_legend(title = "Title"))
  
  print(threshold_plot)
 
  
   
}



if(Appendix_1){
  # Illustrative plot about what can happen if you don't use PBD
  # As illustrative, seed may be hand selected
  
  set.seed(2025 +3)
  
  
  Sim_trial_RAR(0.2,0.4,80,20, gamma = gamma_BRAR, lamda = lambda_BRAR)
}


