## Flexible IA placement code
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

par(mfrow = c(1,1))

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

Beta_comp_exact_3_arm = function(A1,B1, A2, B2, A3, B3){
  ## Exact probability that a Beta(A1,B1) is the maxmimum of three betas
  # Once again taken from Miller's blog

  # The big sum component in Miller's algorithm
  big_sum = 0

  for (i in 0:(A2-1)){
    for (j in 0:(A3-1)){

      big_sum = big_sum + beta(i + j + A1, B1 + B2 + B3)/((B2 + i)*beta(1+i,B2)*(B3+j)*beta(1+j,B3)*beta(A1,B1))

    }
  }

  return_prob = 1 - Beta_comp_exact(A2,B2,A1,B1) - Beta_comp_exact(A3,B3,A1,B1) + big_sum

  return(return_prob)


}

Beta_comp_exact_4_arm = function(A1,B1, A2, B2, A3, B3, A4, B4){
  # Same as other beta's this one is for four arms
  # Once again using Miller's algorithm
  # Determines the porbability of wheter the Beta(A1, B1) is the maximum

  # Once again we get a big sum
  big_sum = 0

  for (i in 0:(A2-1)){
    for (j in 0:(A3-1)){
      for (k in 0:(A4-1)){

        big_sum = big_sum + exp(lbeta(A1+i+j+k, B1+B2+B3+B4)
                     - log(B2+i) - log(B3+j) - log(B4+k)
                     - lbeta(1+i, B2) - lbeta(1+j, B3) - lbeta(1+k, B4)
                     - lbeta(A1, B1))


      }
    }
  }


  return_prob = 1 - Beta_comp_exact(A1,B1,A2,B2) - Beta_comp_exact(A1,B1,A3,B3) - Beta_comp_exact(A1,B1,A4,B4) + Beta_comp_exact_3_arm(A1,B1,A2,B2,A3,B3) + Beta_comp_exact_3_arm(A1,B1,A2,B2,A4,B4) + Beta_comp_exact_3_arm(A1,B1,A4,B4,A3,B3) - big_sum

  return(return_prob)

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

# c_exact_beta, just making it collapsible
if (FALSE){
  cppFunction('


   double c_exact_beta(double a1, double b1, double a2, double b2) {

     double prob = 0;

     // Saves doing a beta calc each loop
     double x = R::beta(a1,b1);

     for (int i = 0; i<=(a2-1); i++){

      //double y = R::beta(i+1,b2)

      prob = prob + R::beta(a1 + i, b1 + b2)/((b2 + i)*(R::beta(i+1,b2)*x));

      //Rcout << prob;

     }


     return prob;
  }
  ')
}

#
## NOTE! n is DIFFERENT for Sim_trial and Sim_trial_RAR.
# In Sim_trial, n (and IAn) represent the sample size PER ARM
# In Sim_trial_RAR and after, n (and IAn) represents the TOTAL sample size

## Will be rewriting this so it is consistent with the RAR, n representing the total sample size
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



Sim_trial_V2 = function(control, treat, n, IAn = c(n), prior_a = 1, prior_b = 1, gamma = 0.92, lambda = 0.9, beta_cpp = TRUE, PairedSeed = FALSE){
  ## This function aims to update Sim_trial to give more flexibility with IA placement
  # The main difference is that IAn will be able to take a vector, and that vector will represent the
  # IA placements

  # control is the control effect
  # treat is the treat effect
  # n is the TOTAL sample size
    # so n/2 is the sample size per arm
  # IAn is the vector that represents the placement of the IA.
  # Prior_a,_b represents the priors for theta
  # gamma lambda are the relevant BOP2 parameters
  # beta_cpp is a boolean which represents whether we want to use the c function for beta comp
  
  if (PairedSeed != FALSE){  
    set.seed(PairedSeed) # sets seed for paired comparisons
  }

  data_control = numeric(0)
  data_treat = numeric(0) # Creating our empty data sets


  IAn = sort(IAn) # So we have it in the right order

  num_IA = length(IAn)

  spare_allo = 2 # This will be a number representing where the leftover patient
   # goes if we have an odd analysis point. 1 means last was control, 0 means treatment

  for (i in 1:num_IA){
    if ( i == 1){
      cur_sample = IAn[1]
    }
    else{
      cur_sample = IAn[i] - IAn[i-1]
    }
    if (cur_sample %% 2 == 0){
      # Simulating IAn[i]/2 each

      data_control = append(data_control,rbinom(cur_sample/2, 1, control))
      data_treat = append(data_treat,rbinom(cur_sample/2, 1, treat))


    }
    else {
      # One simulates IAn/2 floor, one up, will be opposite of previous choice if relevant
      if (spare_allo == 2){
        spare_allo = rbinom(1,1,0.5)
      }
      else{
        spare_allo = 1 - spare_allo
      }
      data_control = append(data_control,rbinom(floor(cur_sample/2) + spare_allo, 1, control))
      data_treat = append(data_treat,rbinom(floor(cur_sample/2) + 1 - spare_allo, 1, treat))

    }
    ## Simulation complete

    ## Tests

    A1 = prior_a + sum(data_treat)  # Updating the posterior
    B1 = prior_b + (length(data_treat) - sum(data_treat))

    A2 = prior_a + sum(data_control)  # Updating the posterior
    B2 = prior_b + (length(data_control) - sum(data_control))


    if (beta_cpp){
      # Outsourcing to c, should be quicker
      IA_test = 1 - rcpp_exact_beta(A1,B1,A2,B2)
    }
    else{
      # R calculations
      IA_test = Beta_comp_exact(A1,B1,A2,B2) # test statistic for IA analysis
    }

    # Threshold probabilities
    Probs = Threshold_prob(lambda,gamma,n,IAn[i])

    # browser() # for testing

    # Efficacy and futility boundaries
    if (IA_test > Probs[2]){  # Threshold can be changed to be more adaptive
      return(c("Reject Null" = TRUE, "Sample size" = length(data_treat) + length(data_control)))  # STOPPING EARLY FOR SUPERIORITY *2 for total sample size (as n per arm)
    } else if (IA_test < Probs[1]){
      return(c("Reject Null" = FALSE, "Sample size" = length(data_treat) + length(data_control))) # STOPPING EARLY FOR FUTILITY
    }

    # At the end of the trials (if IAn[-1] = n) then Probs[1] = Probs[2], so no break condition should be needed.

  }

}

Sim_trial_V2_RAR = function(control, treat, n, IAn = c(n), prior_a = 1, prior_b = 1, gamma = 0.92, lambda = 0.9, beta_cpp = TRUE, full_PBD = TRUE, PairedSeed = FALSE){
  ## This function aims to update Sim_trial to give more flexibility with IA placement
  # The main difference is that IAn will be able to take a vector, and that vector will represent the
  # IA placements

  # control is the control effect
  # treat is the treat effect
  # n is the TOTAL sample size
  # so n/2 is the sample size per arm
  # IAn is the vector that represents the placement of the IA.
  # Prior_a,_b represents the priors for theta
  # gamma lambda are the relevant BOP2 parameters
  # beta_cpp is a boolean which represents whether we want to use the c function for beta comp
  # RAR will begin after the first interim analysis
  # full_PBD  determines whether the whole trial is done with PBD, or just the burn in
  
  if (PairedSeed != FALSE){  
    set.seed(PairedSeed) # sets seed for paired comparisons
  }

  data_control = numeric(0)
  data_treat = numeric(0) # Creating our empty data sets


  IAn = sort(IAn) # So we have it in the right order

  num_IA = length(IAn)

  full_seq = numeric(0) # Will record the full sequence of treatments


  for (i in 1:num_IA){
    if ( i == 1){
      cur_sample = IAn[1]
      seq_treat = rep(c(0,1),floor(IAn[1]/2))  # Needs the first IA to be even in this RAR case
    }
    else{
      cur_sample = IAn[i] - IAn[i-1]
      if (full_PBD){
        # PBD allocation

        # Will prioritise the control if it is not a whole number
        #browser() # testing
        seq_treat = c(rep(1,floor(allo_prob*cur_sample)),rep(0, cur_sample - floor(allo_prob*cur_sample)))
        seq_treat = sample(seq_treat)

      }
      else{
        # Coin based allocation

        seq_treat = rbinom(cur_sample,1,allo_prob) # allo_prob being the probability of being allocated to the treatment

      }
    }

    data_control = append(data_control,rbinom(length(seq_treat) - sum(seq_treat), 1, control))
    data_treat = append(data_treat,rbinom(sum(seq_treat), 1, treat))

    full_seq = append(full_seq,seq_treat)

    ## Simulation complete

    ## Tests

    A1 = prior_a + sum(data_treat)  # Updating the posterior
    B1 = prior_b + (length(data_treat) - sum(data_treat))

    A2 = prior_a + sum(data_control)  # Updating the posterior
    B2 = prior_b + (length(data_control) - sum(data_control))


    if (beta_cpp){
      # Outsourcing to c, should be quicker
      IA_test = 1 - rcpp_exact_beta(A1,B1,A2,B2)
    }
    else{
      # R calculations
      IA_test = Beta_comp_exact(A1,B1,A2,B2) # test statistic for IA analysis
    }



    # Threshold probabilities
    Probs = Threshold_prob(lambda,gamma,n,IAn[i])

    # browser() # for testing

    # Efficacy and futility boundaries
    if (IA_test > Probs[2]){  # Threshold can be changed to be more adaptive

      # Interesting OCs
      null_reject = TRUE
      sample_size = length(data_treat) + length(data_control)
      treat_size = length(data_treat)
      treat_response = sum(data_treat)
      control_size = length(data_control)
      control_response = sum(data_control)
      order_treat = full_seq

      to_return = c("Null reject" = null_reject,"Sample size" = sample_size,"Treatment size" = treat_size,"Treatment responses" = treat_response,"Control size" = control_size,"Control responses" = control_response,"Order of treatment" = order_treat)


      return(to_return)  # STOPPING EARLY FOR SUPERIORITY *2 for total sample size (as n per arm)
    } else if (IA_test < Probs[1]){
      # Interesting OCs
      null_reject = FALSE
      sample_size = length(data_treat) + length(data_control)
      treat_size = length(data_treat)
      treat_response = sum(data_treat)
      control_size = length(data_control)
      control_response = sum(data_control)
      order_treat = full_seq

      to_return = c("Null reject" = null_reject,"Sample size" = sample_size,"Treatment size" = treat_size,"Treatment responses" = treat_response,"Control size" = control_size,"Control responses" = control_response,"Order of treatment" = order_treat)


      return(to_return)
    }

    # At the end of the trials (if IAn[-1] = n) then Probs[1] = Probs[2], so no break condition should be needed.
    # THIS MAY NOT BE TRUE IF THE FIRST IA IS ODD

    # If we haven't ended, we can generate the probabilities for the next IA
    allo_prob = IA_test

    allo_prob_tune = (allo_prob)^(IAn[i]/(2*n))  # Common way of smoothing the allocation probability. The 1-p is because the previous calculations take the probability that the treatment is better than the control, but the following code wants the opposite
    allo_prob_comp = (1 - allo_prob)^(IAn[i]/(2*n))   # Compliment

    allo_prob = allo_prob_tune/(allo_prob_tune + allo_prob_comp)


  }


}

plots_IA_results = function(IAn, n_sim = 10000, RAR = TRUE, control = 0.2, treat = 0.4, n = 80, gamma = 0.97, lambda = 0.91, beta_cpp = TRUE, full_PBD = TRUE, return_var = FALSE, return_prop = FALSE, PairedSeed = FALSE, return_95 = FALSE){

  ESS = 0 # Will add up sample sizes
  Power = 0 # Will add up power
  
  if (return_var){  # Returns variance
    ESS_vec = numeric(0)
  }
  if (return_prop){  # Returns prop
    prop_vec = numeric(0)
  }
  if (return_95){ #  Returns the 2.5% and 97.5% quantiles for ESS
    ESS_vec = numeric(0)
    Power_vec = numeric(0)
  }



  for (i in 1:n_sim){
    if (RAR){
      temp_results = Sim_trial_V2_RAR(control, treat, n, IAn, gamma = gamma, lambda = lambda, beta_cpp = beta_cpp, full_PBD = full_PBD, PairedSeed = PairedSeed)
    }
    else{
      temp_results = Sim_trial_V2(control, treat, n, IAn, gamma = gamma, lambda = lambda, beta_cpp = beta_cpp, PairedSeed = PairedSeed)
    }
    
    if (PairedSeed != FALSE){  
        PairedSeed = PairedSeed + 1 # This makes sure that every sim of the trial uses the same seed, same data generation
          # whilst also making sure that the data points are still i.i.d to one another(as different seed)
      }

    #browser()

    ESS = temp_results[2] + ESS
    
    if (return_var){
      ESS_vec = append(ESS_vec, temp_results[2])
    }
    if (return_prop){
      prop_vec = append(prop_vec, temp_results[3]/temp_results[2])
    }
    if (return_95){
      ESS_vec = append(ESS_vec, temp_results[2])
      Power_vec = append(Power_vec, temp_results[1])
    }


    Power = Power + temp_results[1]




  }

  #browser()

  ESS = ESS/n_sim
  Power = Power/n_sim
  
  
  if (return_95){
    
    ESS_mean_vec = numeric(0)
    Power_mean_vec = numeric(0)
    
    for (i in 1:(n_sim/100)){  # averaging every 100 ESS so the quantiles are more meaningful
      low_index = 1 + (i - 1)*100
      high_index = 100 + (i - 1)*100
      
      #browser()
      
      ESS_mean_vec = append(ESS_mean_vec, mean(ESS_vec[low_index:high_index]))
      
      Power_mean_vec = append(Power_mean_vec, mean(Power_vec[low_index:high_index]))
    }
    
    #browser()
    
    sorted_ESS = sort(ESS_mean_vec) # Smallest first
    ESS_low  = sorted_ESS[ceiling(0.025*(n_sim/100))]
    ESS_high = sorted_ESS[floor(0.975*(n_sim/100))]
    
    sorted_power = sort(Power_mean_vec ) # Smallest first
    Power_low  = sorted_power[ceiling(0.025*(n_sim/100))]
    Power_high = sorted_power[floor(0.975*(n_sim/100))]
    

    return(c(Power, ESS, ESS_low, ESS_high, Power_low, Power_high))
  }
  
  if (return_var){
  
    ESS_Var = var(ESS_vec)
    ESS_sd = sqrt(ESS_Var)
    return(ESS_sd)
  
  }
  if (return_prop){
    return(c("Mean prop" = mean(prop_vec),"Var prop" = var(prop_vec)))
  }
  

  return(c(Power, ESS))
}

prop_plot = function(PairedSeed = FALSE){
  ## Does the sd plot of IA that is needed in the paper
  
  IA_index = (1:39)*2
  
  prop = numeric(length(IA_index))
  
  sd_prop = numeric(length(IA_index))
  
  

  
  for (i in 1:length(IA_index)){
    temp_result = plots_IA_results(c(IA_index[i],80), return_prop = TRUE, RAR = TRUE, PairedSeed = PairedSeed)
    prop[i] = temp_result[1]
    sd_prop[i] = sqrt(temp_result[2])
  }
  
  plot_data = data.frame(IA_index = (1:39)*2, Proportion = prop, Sd = sd_prop)
  
  
  #plot(x = (1:39)*2, y = prop, main = "Proportion of patients allocated to the superior treatment against IA placement", xlab = "IA placement", ylab = "Proportion")
  #lines(x = c(0,80), y = c(0.5,0.5), col = "red")
  
  ggplot(plot_data) +
    geom_bar( aes(x=IA_index, y=Proportion), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=IA_index, ymin=Proportion-1.96*Sd, ymax=Proportion+1.96*Sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    labs( title = "Proportion allocated to best treatment against IA placement.", x = "IA placement") +
    coord_cartesian(ylim = c(0.45,0.65)) +
    scale_y_continuous(breaks = seq(0.45, 0.65, by = 0.05))
  
  ggplot_ess = ggplot(plot_data, aes(x=IA_index, y=Proportion)) +
    geom_ribbon(aes(ymin=Proportion-1.96*Sd, ymax=Proportion+1.96*Sd, alpha = 0.5), fill = "grey70") +
    geom_line(aes(y = Proportion)) +
    labs( title = "Proportion allocated to best treatment against IA placement.", x = "IA placement") +
    coord_cartesian(ylim = c(0.45,0.65)) +
    scale_y_continuous(breaks = seq(0.45, 0.65, by = 0.05)) +
    theme(legend.position="none")
  
  print(ggplot_ess)
  
  browser()
  
}

one_IA_ESS_plot = function(RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4, PairedSeed = FALSE, n_sim = 100000){
  ## This function will just look at the effect of a singular IA placement on power and ESS

  ESS_plot = numeric(40)
  Power_plot = numeric(40)
  low_ESS = numeric(40)
  high_ESS = numeric(40)
  low_power = numeric(40)
  high_power = numeric(40)


  # These magic numbers are from optimizing at IAn = c(60,80)
  if (RAR){
    lambda = 0.91
    gamma = 0.97
  }

  else{
    lambda = 0.91
    gamma = 0.94
  }
  
  # This is what is returned
  # return(c(Power, ESS, ESS_low, ESS_high, Power_low, Power_high))

  for (i in seq(2,78,2)){
    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, control = control, treat = treat, PairedSeed = PairedSeed, n_sim = n_sim, return_95 = TRUE)
    Power_plot[i/2] = temp_results[1]
    ESS_plot[i/2] = temp_results[2]
    low_ESS[i/2] = temp_results[3]
    high_ESS[i/2] = temp_results[4]
    low_power[i/2] = temp_results[5]
    high_power[i/2] = temp_results[6]

  }
  
  temp_results = plots_IA_results(c(80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, control = control, treat = treat, PairedSeed = PairedSeed, n_sim = n_sim, return_95 = TRUE)
  Power_plot[40] = temp_results[1]
  ESS_plot[40] = temp_results[2]
  low_ESS[40] = temp_results[3]
  high_ESS[40] = temp_results[4]
  low_power[40] = temp_results[5]
  high_power[40] = temp_results[6]

  if (treat == control){
    ## Nul situation, different bounds
    ylim_ess = c(50,80)
    ylim_pow = c(0.05,0.2)
    red_line = c(0.1,0.1) # The red line is the default power/type one error
  }
  else{
    #alternate situation
    ylim_ess = c(60,80)
    ylim_pow = c(0.71,0.76)
    red_line = c(max(Power_plot),max(Power_plot))
  }
  
  print(min(ESS_plot))

  plot_df = data.frame(ESS_plot, Power_plot, IA_placement = seq(2,80,2), low_ESS, high_ESS, low_power, high_power)
  
  #plot(y = ESS_plot, x = seq(2,80,2), main = "ESS versus Interim placement - Illustrative example", xlab = "IA placement", ylab = "ESS", type = "l", ylim = ylim_ess)
  #lines(x = c(0,80), y = c(80,80), col = "red")
  #Legend
  #legend("bottomright",
  #       legend = c("No interim."),
  #       col = c("red"),
  #       text.col = "black",
  #       lty = 1,
  #       horiz = F ,
  #       inset = c(0.1, 0.1))
  
  #plot(y = Power_plot, x = seq(2,80,2), main = "Power versus Interim placement - Illustrative example", xlab = "IA placement", ylab = "Power", type = "l", ylim = ylim_pow)
  #arrows(seq(2,80,2), Power_plot + 1.96*0.5/sqrt(n_sim),seq(2,80,2), Power_plot - 1.96*0.5/sqrt(n_sim), angle = 90, code = 3, length = 0.05, col = "orange" )  # Error bars
  #lines(x = c(0,80), y = red_line, col = "red")
  
  #browser()
  

  
  ggplot_ess = ggplot(plot_df, aes(y = ESS_plot, x = IA_placement)) +
    geom_ribbon(aes(ymin = low_ESS, ymax = high_ESS), fill = "grey70") +
    geom_line(aes(y = ESS_plot)) +
    labs( title = "ESS versus Interim placement - Illustrative example.", x = "IA placement", y = "ESS")
  
  
  print(ggplot_ess)

  ggplot_power = ggplot(plot_df, aes(y = Power_plot, x = IA_placement)) +
    geom_ribbon(aes(ymin = Power_plot -  1.96*0.5/sqrt(n_sim), ymax = Power_plot + 1.96*0.5/sqrt(n_sim)), fill = "grey70") +
    geom_line(aes(y = Power_plot)) +
    labs( title = "Power versus Interim placement - Illustrative example.", x = "IA placement", y = "Power")
  print(ggplot_power)
  
  ggarrange(ggplot_ess, ggplot_power,
            ncol = 1, nrow = 2)
  
}

flex_grid_search = function(control = 0.2, treat = 0.4, n = 80, IAn = c(n), prior_a = 1, prior_b = 1, beta_cpp = FALSE, n_sim = 10000, RAR = TRUE){
  ## This function just updates our grid search code, will implement Parallilisation
  ## NOTE! PARALLISATION CURRENTLY DOESNT WORK WITH RCPP
  ## Update: RCPP should work with parallilisation now, will test to see if ti gives improvements

  pot_lambda = seq(0.7 ,1,0.01) # Can vary the parameters to try and get the right kinda area, going off table S9 from supp material
  pot_gamma = seq(0.7 ,1,0.01)

  temp_len = length(pot_lambda)  # Before the length gets updated

  # For easy looping in Foreach
  pot_lambda = rep(pot_lambda, length(pot_gamma))
  pot_gamma = rep(pot_gamma, each = temp_len)

  par_power = numeric(length(pot_gamma))

  cores = detectCores()

  cl = makeCluster(cores/2)  ## To not get overloaded

  registerDoParallel(cl)  # Activates

  # Comments represent regular looping, not parallilsed.


  #for (i in 1:length(pot_gamma)){
  par_power = foreach( i = 1:length(pot_gamma), .combine = c, .export = c("plots_IA_results","Sim_trial_V2_RAR","Sim_trial_V2","Threshold_prob","Beta_comp_exact","sourceCpp","rcpp_exact_beta")) %dopar%{




    if ( plots_IA_results(IAn = IAn, n_sim = n_sim, RAR = RAR, control = control, treat = control, n = n, gamma = pot_gamma[i], lambda = pot_lambda[i], beta_cpp = beta_cpp)[1] < 0.1){
      return(plots_IA_results(IAn = IAn, n_sim = n_sim, RAR = RAR, control = control, treat = treat, n = n, gamma = pot_gamma[i], lambda = pot_lambda[i], beta_cpp = beta_cpp)[1])
      #par_power[i] = plots_IA_results(IAn = IAn, n_sim = n_sim, RAR = RAR, control = control, treat = treat, n = n, gamma = pot_gamma[i], lambda = pot_lambda[i])[1]
    }
    else{
      return(0)
      #par_power[i] = 0
    }


  }




  stopCluster(cl) # deactivates

  max_index = which.max(par_power)

  #browser()

  results = c("Lambda" = pot_lambda[max_index], "Gamma" = pot_gamma[max_index], "Power" = par_power[max_index])

  print(results)

  return(results)


}

# Make true if more grid-searching is required.
if (FALSE){
  t1 = Sys.time()
  flex_grid_search(beta_cpp = TRUE, RAR = FALSE, IAn = c(60,80), n_sim = 10000)
  t1 = Sys.time() - t1
}

# Optimise for no stops, then optimise for IAn = c(60,80)
# For no stops, control = 0.2, treat = 0.4, RAR =  TRUE we get lambda = 0.89, gamma = 0.89, power = 0.77
# For IAn = c(60,80), control = 0.2, treat = 0.4, RAR =  TRUE we get lambda = 0.91, gamma = 0.97, power = 0.75

# For no stops, control = 0.2, treat = 0.4, RAR  = FALSE, we get lambda = 0.89, gamma = 0.92, power = 0.77
# (Note that this is the same as BRAR in the non stopping case)
# For IAn = c(60,80), control = 0.2, treat = 0.4, RAR = FALSE we get lambda = 0.91, gamma = 0.94, power = 0.75

# For 100,000 n_sim, we get that it has a power of 0.7646. This will be a bar that will be compared against.

two_IA_plot = function(second_IA, RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4){
  ## This function will look at the effect of the number of interim analyses,
  # looking at whether adding another interim analysis can reduce sample size
  # second_IA should be a vector of the IA analysis you wish to investigate
    # The function will place a IA analysis before each IA in second_IA and look at its effects
    # (hence why it is called second_IA, each IA in there is the second of two IA)

  num_plots = length(second_IA)

  # These magic numbers are from optimizing at IAn = c(60,80), and IAn = c(40,60,80)
  if (RAR){
    lambda_og = 0.91
    gamma_og = 0.97
    lambda = 0.91
    gamma = 0.98
  }

  else{
    lambda_og = 0.91
    gamma_og = 0.94
    lambda = 0.91
    gamma = 0.94
  }



  for (i in second_IA){

    ESS_plot = numeric(length(seq(2,i,2)))
    Power_plot = numeric(length(seq(2,i,2)))


    for (j in seq(2,i,2)){
      temp_results = plots_IA_results(c(j,i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, control = control, treat = treat)
      Power_plot[j/2] = temp_results[1]
      ESS_plot[j/2] = temp_results[2]

    }

    base_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma_og, lambda = lambda_og, control = control, treat = treat)

    plot(y = ESS_plot, x = seq(2,i,2), main = "ESS versus Interim placement - 2 IA", xlab = "IA placement", ylab = "ESS", type = "l", ylim = c(50,80))
    lines(x = c(0,80), y = c(base_results[2],base_results[2]), col = "red", lty = 2)
    lines(x = c(0,80), y = c(80,80), col = "red", lty = 1)
    plot(y = Power_plot, x = seq(2,i,2), main = "Power versus Interim placement - 2 IA", xlab = "IA placement", ylab = "Power", type = "l", ylim = c(0.6,0.8))
    lines(x = c(0,80), y = c(base_results[1],base_results[1]), col = "red")
    

  }
}

three_IA_min_ESS = function(lambda = 0.91, gamma = 0.98, RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4){
  # This function aims to find the best possible IA arrange with three IA, maximising ESS


  opt_IA = 0
  best_ess = 80

  for (i in seq(6,78,2)){
    # third IA
    #print(i)
    for (j in seq(4,i,2)){
      # second IA
      for (k in seq(2,j,2)){
        # third IA
        temp_ess = plots_IA_results(c(k,j,i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, n_sim = 1000, control = control, treat = treat)[2]
        if (temp_ess < best_ess){
          best_ess = temp_ess
          opt_IA = c(k,j,i,80)
        }
      }
    }
  }

  type_one_error = plots_IA_results(opt_IA, RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, n_sim = 1000, control = control, treat = control)[1]
  return(c(opt_IA, best_ess, type_one_error))

}


## Best ESS is 56.8 with  IAn = c(32.000, 44.000, 62.000, 80.000)

two_IA_min_ESS = function(lambda = 0.91, gamma = 0.98, RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4, PairedSeed = FALSE){
  # This function aims to find the best possible IA arrange with three IA, maximising ESS


  opt_IA = 0
  best_ess = 80

  for (i in seq(4,78,2)){
    # Second IA
    #print(i)
    for (j in seq(2,i,2)){
      # First IA
        temp_ess = plots_IA_results(c(j,i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, n_sim = 1000, control = control, treat = treat)[2]
        if (temp_ess < best_ess){
          best_ess = temp_ess
          opt_IA = c(j,i,80)
      }
    }
  }

  type_one_error = plots_IA_results(opt_IA, RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, n_sim = 1000, control = control, treat = control)[1]
  return(c(opt_IA, best_ess, type_one_error))

}

one_IA_min_ESS = function(lambda = 0.91, gamma = 0.98, RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4, PairedSeed = FALSE){
  
  min_ESS = 80
  best_IA = c(40,80)
  
  for (i in 1:40){
    temp_ess = plots_IA_results(IAn = c(2*i,80), RAR = RAR, lambda = lambda, gamma = gamma, full_PBD = full_PBD, control = control, treat = treat)[2]
    if (temp_ess < min_ESS){
      min_ESS = temp_ess
      best_IA = c(2*i,80)
    }
  }
  
  type_one_error = plots_IA_results(IAn = best_IA, RAR = RAR, lambda = lambda, gamma = gamma, full_PBD = full_PBD, control = control, treat = control)[1] # make sure type one error is controlled
  return(c(min_ESS, type_one_error))
}

one_IA_ESS_plot_range_null = function(RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4, PairedSeed = FALSE){
  ## This function will just look at the effect of a singular IA placement on power and ESS
  ## This updated function looks over a range of null hypothesis and reports the reuslts in one plot

  # Will go up and down 0.1 for the null. The relative effect will stay the same.
  # Minus variables refer to going down 0.1, plus variables refer to going up by one

  ESS_plot_minus = numeric(39)
  Power_plot_minus = numeric(39)

  ESS_plot = numeric(39)
  Power_plot = numeric(39)

  ESS_plot_plus = numeric(39)
  Power_plot_plus = numeric(39)


  # These magic numbers are from optimizing at IAn = c(60,80)
  if (RAR){
    lambda_minus = .88
    gamma_minus = .78

    lambda = 0.91
    gamma = 0.97

    lambda_plus = .9
    gamma_plus = 0.75
  }

  else{
    lambda_minus = .9
    gamma_minus = .96

    lambda = 0.91
    gamma = 0.94

    lambda_plus = 0.91
    gamma_plus = 0.97
  }

  for (i in seq(2,78,2)){
    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma_minus, lambda = lambda_minus, control = control - 0.1, treat = treat - 0.1)
    Power_plot_minus[i/2] = temp_results[1]
    ESS_plot_minus[i/2] = temp_results[2]

    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, control = control, treat = treat)
    Power_plot[i/2] = temp_results[1]
    ESS_plot[i/2] = temp_results[2]

    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma_plus, lambda = lambda_plus, control = control + 0.1, treat = treat + 0.1)
    Power_plot_plus[i/2] = temp_results[1]
    ESS_plot_plus[i/2] = temp_results[2]

  }

  if (treat == control){
    ## Nul situation, different bounds
    ylim_ess = c(50,80)
    ylim_pow = c(0.05,0.2)
    red_line = c(0.1,0.1) # The red line is the default power/type one error
  }
  else{
    #alternate situation
    ylim_ess = c(50,80)
    ylim_pow = c(0.3,1)
    red_line = c(0.77,0.77)
  }


  plot(y = ESS_plot, x = seq(2,78,2), main = "ESS versus Interim placement", xlab = "IA placement", ylab = "ESS", type = "l", ylim = ylim_ess)
  lines(x = c(0,80), y = c(80,80), col = "red")
  lines(x = seq(2,78,2), y = ESS_plot_minus, col = "blue")
  lines(x = seq(2,78,2), y = ESS_plot_plus, col = "green")

  #Legend
  legend("bottomleft",
         legend = c("0.1", "0.2","0.3"),
         col = c("blue",
                 "black",
                 "green"),
         pch = c(17,19,21),
         bty = "n",
         pt.cex = 2,
         cex = 1.2,
         text.col = "black",
         horiz = F ,
         inset = c(0.1, 0.1))



  plot(y = Power_plot, x = seq(2,78,2), main = "Power versus Interim placement", xlab = "IA placement", ylab = "Power", type = "l", ylim = ylim_pow)
  lines(x = c(0,80), y = red_line, col = "red")
  lines(x = seq(2,78,2), y = Power_plot_minus, col = "blue")
  lines(x = seq(2,78,2), y = Power_plot_plus, col = "green")

  #Legend
  legend("bottomleft",
         legend = c("0.1", "0.2","0.3"),
         col = c("blue",
                 "black",
                 "green"),
         pch = c(17,19,21),
         bty = "n",
         pt.cex = 2,
         cex = 1.2,
         text.col = "black",
         horiz = F ,
         inset = c(0.1, 0.1))

}

one_IA_ESS_plot_range_alt = function(RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4, PairedSeed = FALSE){
  ## This function will just look at the effect of a singular IA placement on power and ESS
  ## This updated function looks over a range of alt hypothesis and repororts the results in one plot

  ESS_plot_minus = numeric(39)
  Power_plot_minus = numeric(39)

  ESS_plot = numeric(39)
  Power_plot = numeric(39)

  ESS_plot_plus = numeric(39)
  Power_plot_plus = numeric(39)


  # These magic numbers are from optimizing at IAn = c(60,80)
  if (RAR){
    lambda_minus = .9
    gamma_minus = .9

    lambda = 0.91
    gamma = 0.97

    lambda_plus = .89
    gamma_plus = 0.96
  }

  else{
    lambda_minus = .89
    gamma_minus = .1

    lambda = 0.91
    gamma = 0.94

    lambda_plus = 0.9
    gamma_plus = 0.86
  }

  for (i in seq(2,78,2)){
    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma_minus, lambda = lambda_minus, control = control, treat = treat - 0.1)
    Power_plot_minus[i/2] = temp_results[1]
    ESS_plot_minus[i/2] = temp_results[2]

    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, control = control, treat = treat)
    Power_plot[i/2] = temp_results[1]
    ESS_plot[i/2] = temp_results[2]

    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma_plus, lambda = lambda_plus, control = control, treat = treat + 0.1)
    Power_plot_plus[i/2] = temp_results[1]
    ESS_plot_plus[i/2] = temp_results[2]

  }

  if (treat == control){
    ## Nul situation, different bounds
    ylim_ess = c(0,80)
    ylim_pow = c(0.05,0.2)
    red_line = c(0.1,0.1) # The red line is the default power/type one error
  }
  else{
    #alternate situation
    ylim_ess = c(0,80)
    ylim_pow = c(0.3,1)
    red_line = c(0.77,0.77)
  }



  plot(y = ESS_plot, x = seq(2,78,2), main = "ESS versus Interim placement - Varying H1", xlab = "IA placement", ylab = "ESS", type = "l", ylim = ylim_ess)
  lines(x = c(0,80), y = c(80,80), col = "red")
  lines(x = seq(2,78,2), y = ESS_plot_minus, col = "blue")
  lines(x = seq(2,78,2), y = ESS_plot_plus, col = "green")

  #Legend
  legend("bottomleft",
         legend = c("0.3", "0.4","0.5"),
         col = c("blue",
                 "black",
                 "green"),
         lty = c(1,1,1),
         bty = "n",
         pt.cex = 2,
         cex = 1.2,
         text.col = "black",
         horiz = F ,
         inset = c(0.1, 0.1))



  plot(y = Power_plot, x = seq(2,78,2), main = "Power versus Interim placement - Varying H1", xlab = "IA placement", ylab = "Power", type = "l", ylim = ylim_pow)
  lines(x = c(0,80), y = red_line, col = "red")
  lines(x = seq(2,78,2), y = Power_plot_minus, col = "blue")
  lines(x = seq(2,78,2), y = Power_plot_plus, col = "green")

  #Legend
  legend("bottomleft",
         legend = c("0.3", "0.4","0.5"),
         col = c("blue",
                 "black",
                 "green"),
         pch = c(17,19,21),
         bty = "n",
         pt.cex = 2,
         cex = 1.2,
         text.col = "black",
         horiz = F ,
         inset = c(0.1, 0.1))

}

one_IA_ESS_plot_range_null_567 = function(RAR = TRUE, full_PBD = TRUE, control = 0.6, treat = 0.8, PairedSeed = FALSE){
  ## This function will just look at the effect of a singular IA placement on power and ESS
  ## This updated function looks over a range of null hypothesis and reports the reuslts in one plot
  # Same thing as before, but now fitted for higher nulls

  # Will go up and down 0.1 for the null. The relative effect will stay the same.
  # Minus variables refer to going down 0.1, plus variables refer to going up by one

  ESS_plot_minus = numeric(39)
  Power_plot_minus = numeric(39)

  ESS_plot = numeric(39)
  Power_plot = numeric(39)

  ESS_plot_plus = numeric(39)
  Power_plot_plus = numeric(39)


  # These magic numbers are from optimizing at IAn = c(60,80)
  if (RAR){
    lambda_minus = .9
    gamma_minus = .78

    lambda = 0.9
    gamma = 0.82

    lambda_plus = .91
    gamma_plus = 0.82
  }

  else{
    lambda_minus = .91
    gamma_minus = .98

    lambda = 0.91
    gamma = 0.94

    lambda_plus = 0.9
    gamma_plus = 0.77
  }

  for (i in seq(2,78,2)){
    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma_minus, lambda = lambda_minus, control = control - 0.1, treat = treat - 0.1)
    Power_plot_minus[i/2] = temp_results[1]
    ESS_plot_minus[i/2] = temp_results[2]

    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, control = control, treat = treat)
    Power_plot[i/2] = temp_results[1]
    ESS_plot[i/2] = temp_results[2]

    temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma_plus, lambda = lambda_plus, control = control + 0.1, treat = treat + 0.1)
    Power_plot_plus[i/2] = temp_results[1]
    ESS_plot_plus[i/2] = temp_results[2]

  }

  if (treat == control){
    ## Nul situation, different bounds
    ylim_ess = c(50,80)
    ylim_pow = c(0.05,0.2)
    red_line = c(0.1,0.1) # The red line is the default power/type one error
  }
  else{
    #alternate situation
    ylim_ess = c(50,80)
    ylim_pow = c(0.3,1)
    red_line = c(0.77,0.77)
  }


  plot(y = ESS_plot, x = seq(2,78,2), main = "ESS versus Interim placement", xlab = "IA placement", ylab = "ESS", type = "l", ylim = ylim_ess)
  lines(x = c(0,80), y = c(80,80), col = "red")
  lines(x = seq(2,78,2), y = ESS_plot_minus, col = "blue")
  lines(x = seq(2,78,2), y = ESS_plot_plus, col = "green")

  #Legend
  legend("bottomleft",
         legend = c("0.5", "0.6","0.7"),
         col = c("blue",
                 "black",
                 "green"),
         pch = c(17,19,21),
         bty = "n",
         pt.cex = 2,
         cex = 1.2,
         text.col = "black",
         horiz = F ,
         inset = c(0.1, 0.1))



  plot(y = Power_plot, x = seq(2,78,2), main = "Power versus Interim placement", xlab = "IA placement", ylab = "Power", type = "l", ylim = ylim_pow)
  lines(x = c(0,80), y = red_line, col = "red")
  lines(x = seq(2,78,2), y = Power_plot_minus, col = "blue")
  lines(x = seq(2,78,2), y = Power_plot_plus, col = "green")

  #Legend
  legend("bottomleft",
         legend = c("0.5", "0.6","0.7"),
         col = c("blue",
                 "black",
                 "green"),
         pch = c(17,19,21),
         bty = "n",
         pt.cex = 2,
         cex = 1.2,
         text.col = "black",
         horiz = F ,
         inset = c(0.1, 0.1))

}

one_IA_ESS_plot_seven_null = function(RAR = TRUE, full_PBD = TRUE, alt = TRUE, PairedSeed = FALSE){
  ## This function will just look at the effect of a singular IA placement on power and ESS
  # alt determines whether we are in the null or alt hypothesis
  
  # Each vector in the list will contain the data for different nulls 
  ESS_plot = list(numeric(39),numeric(39),numeric(39),numeric(39),numeric(39),numeric(39),numeric(39))
  Power_plot = list(numeric(39),numeric(39),numeric(39),numeric(39),numeric(39),numeric(39),numeric(39))
  
  controls = 1:7/10
  
  if (alt){
    treats = controls + 0.2
  }
  else{
    treats = controls
  }
  
  # These magic numbers are from optimizing at IAn = c(60,80). Will keep these like this for now. The type one error graph will indicate whether things need refitting
  if (RAR){
    ## Needs to be filled in with every optimal
    lambda = rep(0.91,7)
    gamma = rep(0.97,7)
  }
  
  else{
    
    lambda = numeric(7)
    gamma = numeric(7)
    
    lambda[1]= .9
    gamma[1] = .96
    
    lambda[2] = 0.91
    gamma[2] = 0.94
    
    lambda[3] = 0.91
    gamma[3] = 0.97
    
    lambda[4] = 0.9
    gamma[4] = 0.76
    
    lambda[5] = .91
    gamma[5] = .98
    
    lambda[6] = 0.91
    gamma[6] = 0.94
    
    lambda[7] = 0.9
    gamma[7] = 0.77
    
  }
  
  for (j in 1:7){
    for (i in seq(2,78,2)){
      temp_results = plots_IA_results(c(i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma[j], lambda = lambda[j], control = controls[j], treat = treats[j], PairedSeed = PairedSeed)
      Power_plot[[c(j,i/2)]] = temp_results[1]
      ESS_plot[[c(j,i/2)]] = temp_results[2]
      
    }
  }
  
  if (alt){
    ylim_ess = c(50,80)
    ylim_pow = c(0.3,1)
    red_line = c(max(Power_plot[[2]]),max(Power_plot[[2]]))
    
  }
  else{
    ylim_ess = c(50,80)
    ylim_pow = c(0.05,0.2)
    red_line = c(0.1,0.1) # The red line is the default power/type one error
  }

  #browser()
  
  ## ggplot plots
  
  if (TRUE){ #  make true for ggplots
    
    #browser()
    
    ESS_df = c(ESS_plot[[1]],ESS_plot[[2]],ESS_plot[[3]],ESS_plot[[4]],ESS_plot[[5]],ESS_plot[[6]],ESS_plot[[7]])
    Power_df = c(Power_plot[[1]],Power_plot[[2]],Power_plot[[3]],Power_plot[[4]],Power_plot[[5]],Power_plot[[6]],Power_plot[[7]])
    
    plot_df = data.frame(IA_placement = rep(seq(2,78,2),7), ESS = ESS_df, Power = Power_df, Null = factor(rep((1:7)/10, each = 39)))
    
    ggplot_ess = ggplot(plot_df, aes(y = ESS, x = IA_placement)) +
      geom_line(aes(colour = Null, linetype = Null), linewidth = 1.5) +
      labs( title = "ESS versus Interim placement - Varying control effect", x = "IA placement", y = "ESS", colour = "Control success rate", linetype = "Control success rate") +
      scale_color_viridis(discrete=TRUE, option="viridis")
      
      #scale_color_manual(values=c("#000000", "#444444", "#888888","#aaaaaa","#cccccc","#eeeeee","#ffffff")) #greyscale
    
    #print(ggplot_ess)
    
    ggplot_power = ggplot(plot_df, aes(y = Power, x = IA_placement)) +
      geom_line(aes(colour = Null, linetype = Null), linewidth = 1.5) +
      labs( title = "Power versus Interim placement - Varying control effect", x = "IA placement", y = "Power", colour = "Control success rate", linetype = "Control success rate") +
      scale_color_viridis(discrete=TRUE, option="viridis")
      
    
    #scale_color_manual(values=c("#000000", "#444444", "#888888","#aaaaaa","#cccccc","#eeeeee","#ffffff"))
    
    #print(ggplot_power)
    
    return_plot = ggarrange(ggplot_ess, ggplot_power,
              ncol = 1, nrow = 2)
    print(return_plot)
  }
  # regular R plots
  if (FALSE){ # Make true for base R plots
  
    plot(y = ESS_plot[[1]], x = seq(2,78,2), main = "ESS versus Interim placement", xlab = "IA placement", ylab = "ESS", type = "l", ylim = ylim_ess, col = 1)
    for (j in 2:7){
      lines(y = ESS_plot[[j]], x = seq(2,78,2), col = j, type = "l")
    }
    legend(x = "bottomright",          # Position
           legend = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7"),  # Legend texts
           lty = rep(1,7),           # Line types
           col = 1:7,           # Line colors
           lwd = 2)                 # Line width
    
    plot(y = Power_plot[[1]], x = seq(2,78,2), main = "Power versus Interim placement", xlab = "IA placement", ylab = "Power", type = "l", ylim = ylim_pow)
    for (j in 2:7){
      lines(y = Power_plot[[j]], x = seq(2,78,2), col = j, type = "l")
    }
    legend(x = "bottomright",          # Position
           legend = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7"),  # Legend texts
           lty = rep(1,7),           # Line types
           col = 1:7,           # Line colors
           lwd = 2)                 # Line width
  
  }
  
}

max_sample_size_plot = function(lambda = 0.91, gamma = 0.98, RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4, PairedSeed = FALSE, type_one = FALSE){
  ## This does a plot of how the ESS can change based on the maximum sample size, and what effect
  # the IA has as max sample size increases
  
  
  # Sample sizes we're looking at
  Max_sample_sizes = c(40,60,80,100,120,140,160)
  
  ESS_plot = list(numeric(19),numeric(29),numeric(39),numeric(49),numeric(59),numeric(69),numeric(79))

  Power_plot = list(numeric(19),numeric(29),numeric(39),numeric(49),numeric(59),numeric(69),numeric(79))

  

  for (j in 1:7){
    for (i in seq(2,20*j  + 18 ,2)){
      temp_results = plots_IA_results(c(i,Max_sample_sizes[j]), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, control = control, treat = treat, n = Max_sample_sizes[j], PairedSeed = PairedSeed)
      
      ESS_plot[[c(j,i/2)]] = temp_results[2]/Max_sample_sizes[j]
      # We're dividing by Max_sample_size here as we just want to see percentage done
      Power_plot[[c(j,i/2)]] = temp_results[1]
      
    }
  }
  
  
  ylim_ess = c(0.4,1)

  
  #browser()
  
  if (TRUE){ #  make true for ggplots
    
    #browser()
    
    IA_placement = numeric(0)
    Null = numeric(0) # representing max sample size, just called null for ease of copy/paste
    
    for (j in 1:7){
      IA_placement = append(IA_placement,seq(2,20*j  + 18 ,2))
      Null = append(Null,rep(20*(j+1),length(seq(2,20*j  + 18 ,2))))
    }
    
    ESS_df = c(ESS_plot[[1]],ESS_plot[[2]],ESS_plot[[3]],ESS_plot[[4]],ESS_plot[[5]],ESS_plot[[6]],ESS_plot[[7]])

    Power_df = c(Power_plot[[1]],Power_plot[[2]],Power_plot[[3]],Power_plot[[4]],Power_plot[[5]],Power_plot[[6]],Power_plot[[7]])
    
    
    #browser()
    
    plot_df = data.frame(IA_placement = IA_placement/Null, ESS = ESS_df,Power = Power_df, Null = factor(Null))
    
    #browser()
    
    ggplot_ess = ggplot(plot_df, aes(y = ESS, x = IA_placement)) +
      geom_line(aes(colour = Null, linetype = Null), linewidth = 1.5) +
      labs( title = "ESS versus Interim placement - Varying maximum sample size", x = "IA placement through proportion of maximum sample size", y = "ESS (proportion of max sample size)", colour = "Maximum sample size", linetype = "Maximum sample size") +
      scale_color_viridis(discrete=TRUE, option="viridis")
    
    #scale_color_manual(values=c("#000000", "#444444", "#888888","#aaaaaa","#cccccc","#eeeeee","#ffffff")) #greyscale
    
    print(ggplot_ess)
    
    if (type_one == TRUE){
      
      ggplot_power = ggplot(plot_df, aes(y = Power, x = IA_placement)) +
        geom_line(aes(colour = Null, linetype = Null), linewidth = 1.5) +
        labs( title = "Power versus Interim placement - Varying maximum sample size", x = "IA placement through proportion of maximum sample size", y = "Type one error", colour = "Maximum sample size", linetype = "Maximum sample size") +
        scale_color_viridis(discrete=TRUE, option="viridis")
      
      print(ggplot_power)
      
    }
    
  }
  
  if (FALSE){
    plot(y = ESS_plot[[1]], x = seq(2,38,2)/Max_sample_sizes[1], main = "ESS versus Interim placement", xlab = "IA placement", ylab = "ESS", type = "l", ylim = ylim_ess, col = 1)
    for (j in 2:7){
      lines(y = ESS_plot[[j]], x = seq(2,20*j  + 18,2)/Max_sample_sizes[j], col = j, type = "l")
    }
    #browser()
    
    legend(x = "bottomright",          # Position
           legend = Max_sample_sizes,  # Legend texts
           lty = rep(1,7),           # Line types
           col = 1:7,           # Line colors
           lwd = 2)                 # Line width
    
    
  }

  
  
}




## RAR x IA prop stuff

three_IA_min_prop = function(lambda = 0.91, gamma = 0.98, RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4){
  # This function aims to find the best possible IA arrange with three IA, maximising ESS
  
  
  opt_IA = 0
  best_prop = 0
  
  for (i in seq(6,78,2)){
    # third IA
    #print(i)
    for (j in seq(4,i,2)){
      # second IA
      for (k in seq(2,j,2)){
        # third IA
        temp_prop = plots_IA_results(c(k,j,i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, n_sim = 1000, control = control, treat = treat, return_prop = TRUE)[1]
        if (temp_prop > best_prop){
          best_prop = temp_prop
          opt_IA = c(k,j,i,80)
        }
      }
    }
  }
  
  return(c(opt_IA, best_prop))
  
}

two_IA_min_prop = function(lambda = 0.91, gamma = 0.98, RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4, PairedSeed = FALSE){
  # This function aims to find the best possible IA arrange with three IA, maximising ESS
  
  
  opt_IA = 0
  best_prop = 0
  
  for (i in seq(4,78,2)){
    # Second IA
    #print(i)
    for (j in seq(2,i,2)){
      # First IA
      temp_prop = plots_IA_results(c(j,i,80), RAR = RAR, full_PBD = full_PBD, gamma = gamma, lambda = lambda, n_sim = 1000, control = control, treat = treat, return_prop = TRUE)[1]
      if (temp_prop > best_prop){
        best_prop = temp_prop
        opt_IA = c(j,i,80)
      }
    }
  }
  
  return(c(opt_IA, best_prop))
  
}

one_IA_min_prop = function(lambda = 0.91, gamma = 0.98, RAR = TRUE, full_PBD = TRUE, control = 0.2, treat = 0.4, PairedSeed = FALSE){
  
  best_prop = 0
  best_IA = 0
  
  for (i in 1:40){
    temp_prop = plots_IA_results(IAn = c(2*i,80), RAR = RAR, lambda = lambda, gamma = gamma, full_PBD = full_PBD, control = control, treat = treat, return_prop = 1)[1]
    if (temp_prop > best_prop){
      best_prop = temp_prop
      best_IA = c(2*i,80)
    }
  }
  
  return(c(best_prop, best_IA))
}


### REPRODUCIBLE FIGURES

## Figure 2+3

#Prints both figure 2+3

Figure_2 = FALSE
if (Figure_2){
  set.seed(2025)
  one_IA_ESS_plot(RAR = FALSE, PairedSeed = 2025)
}

Appendix_for_Figure_2 = FALSE
if (Appendix_for_Figure_2){
  set.seed(2025)
  one_IA_ESS_plot(RAR = FALSE, PairedSeed = 2025, treat = 0.2, n_sim = 10000)
  # gives type one error
}

## Figure 4+5

Figure_4 = FALSE
if (Figure_4){
  set.seed(2025)
  one_IA_ESS_plot_seven_null(RAR = FALSE, PairedSeed = 2025)
}

Figure_6 = FALSE
if (Figure_6){
  set.seed(2025)
  max_sample_size_plot(RAR = FALSE, PairedSeed = 2025)

  
}


Appendix_for_Figure_6 = FALSE
if (Appendix_for_Figure_6){
  set.seed(2025)
  max_sample_size_plot(RAR = FALSE, PairedSeed = 2025, type_one = TRUE, treat = 0.2)
  # type_one = TRUE (type one error control)
  
}


# Figure 7 is on multi-arm code

## Table 1

Table_1 = FALSE
if (Table_1){
  set.seed(2025)
  
  # One IA
  
  print(one_IA_min_ESS(RAR = FALSE))
  print(plots_IA_results(RAR = FALSE, IAn = c(40,80)))
  print(plots_IA_results(RAR = FALSE, IAn = c(52,80)))
  
  
  # Two IA
  
  print(two_IA_min_ESS(RAR = FALSE))
  print(plots_IA_results(RAR = FALSE, IAn = c(26,54,80)))
  print(plots_IA_results(RAR = FALSE, IAn = c(42,62,80)))
  
  # Three IA
  
  print(three_IA_min_ESS(RAR = FALSE))
  print(plots_IA_results(RAR = FALSE, IAn = c(20,40,60,80)))
  print(plots_IA_results(RAR = FALSE, IAn = c(38,52,66,80)))
  
}

Table_RAR_X_IA = FALSE
if (Table_RAR_X_IA){
  set.seed(2025)
  
  # One IA
  
  print(one_IA_min_prop())
  print(plots_IA_results(IAn = c(40,80), return_prop = TRUE))
  print(plots_IA_results(IAn = c(52,80), return_prop = TRUE))
  
  
  # Two IA
  
  print(two_IA_min_prop())
  print(plots_IA_results(IAn = c(26,54,80), return_prop = TRUE))
  print(plots_IA_results(IAn = c(42,62,80), return_prop = TRUE))
  
  # Three IA
  
  print(three_IA_min_prop())
  print(plots_IA_results(IAn = c(20,40,60,80), return_prop = TRUE))
  print(plots_IA_results(IAn = c(38,52,66,80), return_prop = TRUE))
  
}

## Figure 8, ditch

Figure_8 = FALSE
if (Figure_8){
  set.seed(2025)
  two_IA_plot(60, RAR = FALSE)
}


## Figures 9 through 11 are in BOP2 repoducible, as well as Table 2

## Table 3, Table 4, Table 5, Table 6 are each in the multi arm setting



# Figure 12 and 13 aren't stochastic, in multi arm code


## Figuire 14

Figure_14 = TRUE
if(Figure_14){
  set.seed(2025)
  
  prop_plot(PairedSeed = 2025)
  
}
# Maybe add error bars to proportion plot? Or wait for review

## Figure 15 and 16 are in multi arm.

Appendix_1 = FALSE
if(Appendix_1){
  set.seed(2025)
  
  one_IA_ESS_plot_range_alt(RAR = FALSE)
}


