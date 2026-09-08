library(MyPackageV2) # contains RCPP code for parallilisation #

library(doParallel)  # Parallelisation
library(doRNG) #  Reproducible in parallel
library(foreach)  # Parallelisation

library(ggplot2)
library(ggpubr)

library(microbenchmark)

prob_bin_supe = function(n, y, p_0, p_1){
  # The probability that X_1, a binomila r.v., is greater than X_0, a binomial r.v. with the same n by y 
  # P(X_1 > X_0 + y)
  
  prob = 0 
  
  for (i in y:n){
    
    # Probability that X_1 above i, with the exact probability that x_0 = i - y  
    prob = prob + (1 - pbinom(i,n,p_1))*(pbinom(i - y,n,p_0) - pbinom(i - y - 1,n,p_0))
    
  }
  
  return(prob)
  
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

prob_end_bop2 = function(n_pa, p_0, p_1, lambda, gamma, n_total = 80){
  # Calculates the probability of ending the trial early for supe or futility
  # assuming er, n_pa is patients per arm
  # p_0 null prob
  # p_1 is alt prob
  # lambda and gamma optimisation parameters 
  
  end_early_prob = 0
  
  thresh_prob = Threshold_prob(lambda, gamma, n_total, n_pa*2 )
  
  prob_supe = thresh_prob[2]
  
  prob_fute = thresh_prob[1]
  
  for (i in 0:n_pa){
    for (j in 0:n_pa){
      
      prob = 1 - rcpp_exact_beta(i + 1,n_pa - i + 1,j + 1,n_pa - j + 1)
      
      if (prob < prob_fute | prob > prob_supe ){
        end_early_prob = end_early_prob + (pbinom(i,n_pa,p_1 ) - pbinom(i - 1,n_pa,p_1 ))*(pbinom(j,n_pa,p_0 ) - pbinom(j - 1,n_pa,p_0 ))
      }
      
      #browser()
    }
  }
  
  return(end_early_prob)
  
  
}

both_prob_end_bop2 = function(n_pa, p_0, p_1, lambda, gamma, n_total = 80){
  # Calculates the probability of ending the trial early for supe or futility
  ## specifically, returns both probabilities in a vector
  # assuming er, n_pa is patients per arm
  # p_0 null prob
  # p_1 is alt prob
  # lambda and gamma optimisation parameters 
  # Only works for 0 IA, one final 
  
  fute_early_prob = 0
  supe_early_prob = 0
  
  thresh_prob = Threshold_prob(lambda, gamma, n_total, n_pa*2 )
  
  prob_supe = thresh_prob[2]
  
  prob_fute = thresh_prob[1]
  
  for (i in 0:n_pa){
    for (j in 0:n_pa){
      
      prob = 1 - rcpp_exact_beta(i + 1,n_pa - i + 1,j + 1,n_pa - j + 1)
      
      if (prob < prob_fute){
        fute_early_prob = fute_early_prob + (pbinom(i,n_pa,p_1 ) - pbinom(i - 1,n_pa,p_1 ))*(pbinom(j,n_pa,p_0 ) - pbinom(j - 1,n_pa,p_0 ))
      }
      else if (prob > prob_supe){
        supe_early_prob = supe_early_prob + (pbinom(i,n_pa,p_1 ) - pbinom(i - 1,n_pa,p_1 ))*(pbinom(j,n_pa,p_0 ) - pbinom(j - 1,n_pa,p_0 ))
      }
      
      #browser()
    }
  }
  
  return(c("Futility" = fute_early_prob, "Supe" = supe_early_prob))
  
  
}

both_prob_end_bop2_1IA = function(p_0, p_1, lambda, gamma, n_total = 80, IAn){
  # Calculates the probability of ending the trial early for supe or futility
  ## specifically, returns both probabilities in a vector
  # assuming er, n_pa is patients per arm
  # p_0 null prob
  # p_1 is alt prob
  # lambda and gamma optimisation parameters 
  # ONLY WORKS IF IAn[final] =  n_total!
  # Works for up to 4 IA 
  
  fute_early_prob = numeric(length(IAn))
  supe_early_prob = numeric(length(IAn))
  
  n_1_pa = IAn[1]/2
  
  thresh_prob = Threshold_prob(lambda, gamma, n_total, n_1_pa*2 )
  
  prob_supe = thresh_prob[2]
  
  prob_fute = thresh_prob[1]
  
  for (i in 0:n_1_pa){
    for (j in 0:n_1_pa){
      #print(i)
      #1st IA
      
      thresh_prob = Threshold_prob(lambda, gamma, n_total, n_1_pa*2 )
      
      prob_supe = thresh_prob[2]
      
      prob_fute = thresh_prob[1]
      
      prob = 1 - rcpp_exact_beta(i + 1,n_1_pa - i + 1,j + 1,n_1_pa - j + 1)
      
      if (prob < prob_fute){ # Futility
        fute_early_prob[1] = fute_early_prob[1] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))
      }
      else if (prob > prob_supe){ # Supe
        supe_early_prob[1] = supe_early_prob[1] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))
      }
      else{ # Continue 
        n_2_pa = (IAn[2] - IAn[1])/2
        
        
        
        for (a in 0:n_2_pa){
          for (b in 0:n_2_pa){
            # 2nd IA
            #print(a)
            
            thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[2] )
            
            prob_supe = thresh_prob[2]
            
            prob_fute = thresh_prob[1]
            
            prob = 1 - rcpp_exact_beta(i + a + 1,IAn[2]/2 - i - a+ 1,j + b + 1,IAn[2]/2 - j - b + 1)
            
            if (prob < prob_fute){
              fute_early_prob[2] = fute_early_prob[2] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))
            }
            else if (prob > prob_supe){
              supe_early_prob[2] = supe_early_prob[2] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))
            }
            
            
          }
        }
        
      }
      
      #browser()
    }
  }
  
  return(c("Futility" = fute_early_prob, "Supe" = supe_early_prob))
  
  
}


both_prob_end_bop2_4IA = function(p_0, p_1, lambda, gamma, n_total = 80, IAn){
  # Calculates the probability of ending the trial early for supe or futility
  ## specifically, returns both probabilities in a vector
  # assuming er, n_pa is patients per arm
  # p_0 null prob
  # p_1 is alt prob
  # lambda and gamma optimisation parameters 
  # ONLY WORKS IF IAn[final] =  n_total!
  # Works for up to 4 IA 
  
  fute_early_prob = numeric(length(IAn))
  supe_early_prob = numeric(length(IAn))
  
  n_1_pa = IAn[1]/2
  
  thresh_prob = Threshold_prob(lambda, gamma, n_total, n_1_pa*2 )
  
  prob_supe = thresh_prob[2]
  
  prob_fute = thresh_prob[1]
  
  for (i in 0:n_1_pa){
    for (j in 0:n_1_pa){
      #print(i)
      #1st IA
      
      thresh_prob = Threshold_prob(lambda, gamma, n_total, n_1_pa*2 )
      
      prob_supe = thresh_prob[2]
      
      prob_fute = thresh_prob[1]
      
      prob = 1 - rcpp_exact_beta(i + 1,n_1_pa - i + 1,j + 1,n_1_pa - j + 1)
      
      if (prob < prob_fute){ # Futility
        fute_early_prob[1] = fute_early_prob[1] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))
      }
      else if (prob > prob_supe){ # Supe
        supe_early_prob[1] = supe_early_prob[1] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))
      }
      else{ # Continue 
        n_2_pa = (IAn[2] - IAn[1])/2
        
        
        
        for (a in 0:n_2_pa){
          for (b in 0:n_2_pa){
            # 2nd IA
            #print(a)
            
            thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[2] )
            
            prob_supe = thresh_prob[2]
            
            prob_fute = thresh_prob[1]
            
            prob = 1 - rcpp_exact_beta(i + a + 1,IAn[2]/2 - i - a+ 1,j + b + 1,IAn[2]/2 - j - b + 1)
            
            if (prob < prob_fute){
              fute_early_prob[2] = fute_early_prob[2] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))
            }
            else if (prob > prob_supe){
              supe_early_prob[2] = supe_early_prob[2] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))
            }
            else{
              
              n_3_pa = (IAn[3] - IAn[2])/2
              
              for (c in 0:n_3_pa){
                for (d in 0:n_3_pa){
                  
                  #print(3)
                  
                  thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[3] )
                  
                  prob_supe = thresh_prob[2]
                  
                  prob_fute = thresh_prob[1]
                  
                  prob = 1 - rcpp_exact_beta(i + a + c + 1,IAn[3]/2 - i - a - c + 1,j + b + d + 1,IAn[3]/2 - j - b - d+ 1)
                  
                  if (prob < prob_fute){
                    fute_early_prob[3] = fute_early_prob[3] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))*(pbinom(c,n_3_pa,p_1 ) - pbinom(c - 1,n_3_pa,p_1 ))*(pbinom(d,n_3_pa,p_0 ) - pbinom(d - 1,n_3_pa,p_0 ))
                  }
                  else if (prob > prob_supe){
                    supe_early_prob[3] = supe_early_prob[3] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))*(pbinom(c,n_3_pa,p_1 ) - pbinom(c - 1,n_3_pa,p_1 ))*(pbinom(d,n_3_pa,p_0 ) - pbinom(d - 1,n_3_pa,p_0 ))
                  }
                  else{
                    
                    n_4_pa = (IAn[4] - IAn[3])/2
                    
                    thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[4] )
                    
                    prob_supe = thresh_prob[2]
                    
                    prob_fute = thresh_prob[1]
                    
                    for(e in 0:n_4_pa){
                      for (f in 0:n_4_pa){
                        #print(4)
                        
                        prob = 1 - rcpp_exact_beta(i + a + c + e + 1,IAn[4]/2 - i - a - c - e + 1,j + b + d + f + 1,IAn[4]/2 - j - b - d - f+ 1)
                        
                        if (prob < prob_fute){
                          fute_early_prob[4] = fute_early_prob[4] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))*(pbinom(c,n_3_pa,p_1 ) - pbinom(c - 1,n_3_pa,p_1 ))*(pbinom(d,n_3_pa,p_0 ) - pbinom(d - 1,n_3_pa,p_0 ))*(pbinom(e,n_4_pa,p_1 ) - pbinom(e - 1,n_4_pa,p_1 ))*(pbinom(f,n_4_pa,p_0 ) - pbinom(f - 1,n_4_pa,p_0 ))
                        }
                        else if (prob > prob_supe){
                          supe_early_prob[4] = supe_early_prob[4] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))*(pbinom(c,n_3_pa,p_1 ) - pbinom(c - 1,n_3_pa,p_1 ))*(pbinom(d,n_3_pa,p_0 ) - pbinom(d - 1,n_3_pa,p_0 ))*(pbinom(e,n_4_pa,p_1 ) - pbinom(e - 1,n_4_pa,p_1 ))*(pbinom(f,n_4_pa,p_0 ) - pbinom(f - 1,n_4_pa,p_0 ))
                        }
                        
                        
                        
                      }
                    }
                    
                  }
                  
                  
                }
              }
              
            }
            
            
            
          }
        }
        
      }
      
      #browser()
    }
  }
  
  return(c("Futility" = fute_early_prob, "Supe" = supe_early_prob))
  
  
}

both_prob_end_bop2_4IA_par = function(p_0, p_1, lambda, gamma, n_total = 80, IAn){
  # Calculates the probability of ending the trial early for supe or futility
  ## specifically, returns both probabilities in a vector
  # assuming er, n_pa is patients per arm
  # p_0 null prob
  # p_1 is alt prob
  # lambda and gamma optimisation parameters 
  # ONLY WORKS IF IAn[final] =  n_total!
  
  
  
  n_1_pa = IAn[1]/2
  
  thresh_prob = Threshold_prob(lambda, gamma, n_total, n_1_pa*2 )
  
  prob_supe = thresh_prob[2]
  
  prob_fute = thresh_prob[1]
  
  
  par_results = vector("list",n_1_pa + 1) #  Will be a list, where each element in the list is a vector of 8 values, indicating the prob of stopping at each step
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  
  
  #for (i in 0:n_1_pa){
  par_results = foreach( i = 0:n_1_pa, .export = c("rcpp_exact_beta","Threshold_prob"))  %dopar% {
    
    fute_early_prob = numeric(length(IAn))
    supe_early_prob = numeric(length(IAn))
    
    
    for (j in 0:n_1_pa){
      #print(i)
      #1st IA
      
      thresh_prob = Threshold_prob(lambda, gamma, n_total, n_1_pa*2 )
      
      prob_supe = thresh_prob[2]
      
      prob_fute = thresh_prob[1]
      
      prob = 1 - rcpp_exact_beta(i + 1,n_1_pa - i + 1,j + 1,n_1_pa - j + 1)
      
      if (prob < prob_fute){ # Futility
        fute_early_prob[1] = fute_early_prob[1] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))
      }
      else if (prob > prob_supe){ # Supe
        supe_early_prob[1] = supe_early_prob[1] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))
      }
      else{ # Continue 
        n_2_pa = (IAn[2] - IAn[1])/2
        
        
        
        for (a in 0:n_2_pa){
          for (b in 0:n_2_pa){
            # 2nd IA
            #print(a)
            
            thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[2] )
            
            prob_supe = thresh_prob[2]
            
            prob_fute = thresh_prob[1]
            
            prob = 1 - rcpp_exact_beta(i + a + 1,IAn[2]/2 - i - a+ 1,j + b + 1,IAn[2]/2 - j - b + 1)
            
            if (prob < prob_fute){
              fute_early_prob[2] = fute_early_prob[2] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))
            }
            else if (prob > prob_supe){
              supe_early_prob[2] = supe_early_prob[2] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))
            }
            else{
              
              n_3_pa = (IAn[3] - IAn[2])/2
              
              for (c in 0:n_3_pa){
                for (d in 0:n_3_pa){
                  
                  #print(3)
                  
                  thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[3] )
                  
                  prob_supe = thresh_prob[2]
                  
                  prob_fute = thresh_prob[1]
                  
                  prob = 1 - rcpp_exact_beta(i + a + c + 1,IAn[3]/2 - i - a - c + 1,j + b + d + 1,IAn[3]/2 - j - b - d+ 1)
                  
                  if (prob < prob_fute){
                    fute_early_prob[3] = fute_early_prob[3] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))*(pbinom(c,n_3_pa,p_1 ) - pbinom(c - 1,n_3_pa,p_1 ))*(pbinom(d,n_3_pa,p_0 ) - pbinom(d - 1,n_3_pa,p_0 ))
                  }
                  else if (prob > prob_supe){
                    supe_early_prob[3] = supe_early_prob[3] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))*(pbinom(c,n_3_pa,p_1 ) - pbinom(c - 1,n_3_pa,p_1 ))*(pbinom(d,n_3_pa,p_0 ) - pbinom(d - 1,n_3_pa,p_0 ))
                  }
                  else{
                    
                    n_4_pa = (IAn[4] - IAn[3])/2
                    
                    thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[4] )
                    
                    prob_supe = thresh_prob[2]
                    
                    prob_fute = thresh_prob[1]
                    
                    for(e in 0:n_4_pa){
                      for (f in 0:n_4_pa){
                        #print(4)
                        
                        prob = 1 - rcpp_exact_beta(i + a + c + e + 1,IAn[4]/2 - i - a - c - e + 1,j + b + d + f + 1,IAn[4]/2 - j - b - d - f+ 1)
                        
                        if (prob < prob_fute){
                          fute_early_prob[4] = fute_early_prob[4] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))*(pbinom(c,n_3_pa,p_1 ) - pbinom(c - 1,n_3_pa,p_1 ))*(pbinom(d,n_3_pa,p_0 ) - pbinom(d - 1,n_3_pa,p_0 ))*(pbinom(e,n_4_pa,p_1 ) - pbinom(e - 1,n_4_pa,p_1 ))*(pbinom(f,n_4_pa,p_0 ) - pbinom(f - 1,n_4_pa,p_0 ))
                        }
                        else if (prob > prob_supe){
                          supe_early_prob[4] = supe_early_prob[4] + (pbinom(i,n_1_pa,p_1 ) - pbinom(i - 1,n_1_pa,p_1 ))*(pbinom(j,n_1_pa,p_0 ) - pbinom(j - 1,n_1_pa,p_0 ))*(pbinom(a,n_2_pa,p_1 ) - pbinom(a - 1,n_2_pa,p_1 ))*(pbinom(b,n_2_pa,p_0 ) - pbinom(b - 1,n_2_pa,p_0 ))*(pbinom(c,n_3_pa,p_1 ) - pbinom(c - 1,n_3_pa,p_1 ))*(pbinom(d,n_3_pa,p_0 ) - pbinom(d - 1,n_3_pa,p_0 ))*(pbinom(e,n_4_pa,p_1 ) - pbinom(e - 1,n_4_pa,p_1 ))*(pbinom(f,n_4_pa,p_0 ) - pbinom(f - 1,n_4_pa,p_0 ))
                        }
                        
                        
                        
                      }
                    }
                    
                  }
                  
                  
                }
              }
              
            }
            
            
            
          }
        }
        
      }
      
      #browser()
    }
    
    return(c(fute_early_prob,supe_early_prob))
  }
  
  stopCluster(cl) # deactivates
  
  fute_early_prob_return = numeric(length(IAn))
  supe_early_prob_return = numeric(length(IAn))
  
  for (i in 1:(n_1_pa+1)){
    fute_early_prob_return = fute_early_prob_return + par_results[[i]][1:length(IAn)]
    supe_early_prob_return = supe_early_prob_return + par_results[[i]][(length(IAn)+1):(2*length(IAn))]
    
  }
  
  return(c("Futility" = fute_early_prob_return, "Supe" = supe_early_prob_return))
  
  
}

pre_def_prob = function(p_0 = 0.2, p_1 = 0.4, n_total = 80, IAn = c(20,40,60,80), par = TRUE){
  # This function is independent of BOP2, and just calculates the probability of observing every single value
  # Ex, observing 5 vs 6 successes at the first IA, followed by 10 vs 13 successes at the second IA would have a probability p
  # Question about how to index
  # Defunct, the vector would be too big 
  
  num_IA = length(IAn)
  
  browser()
  
  num_perm = prod(((IAn/2)^2))  # IAn/2 per arm per IA, then multiply everything together
  
  prob_df = data.frame("Control_Success_1" = numeric(num_perm),"Control_Success_2" = numeric(num_perm),"Control_Success_3" = numeric(num_perm),"Control_Success_4" = numeric(num_perm),
                       "Treat_Success_1" = numeric(num_perm),"Treat_Success_2" = numeric(num_perm),"Treat_Success_3" = numeric(num_perm),"Treat_Success_4" = numeric(num_perm),
                       "Probability" = numeric(num_perm))
  count = 0
  
  for(a in IAn[1]/2){
    print(a)
    for (b in IAn[1]/2){
      for(c in IAn[2]/2){
        for (d in IAn[2]/2){
          for(e in IAn[3]/2){
            for (f in IAn[3]/2){
              for(g in IAn[4]/2){
                for (h in IAn[4]/2){
                  
                  count = count + 1
                  
                  temp_prob  = (pbinom(a,IAn[1]/2, p_0) - pbinom(a - 1,IAn[1]/2, p_0))*(pbinom(c,(IAn[2]- IAn[1])/2, p_0) - pbinom(c - 1,(IAn[2]- IAn[1])/2, p_0))*(pbinom(e,(IAn[3]- IAn[2])/2, p_0) - pbinom(e - 1,(IAn[3]- IAn[2])/2, p_0))*(pbinom(g,(IAn[4]- IAn[3])/2, p_0) - pbinom(g - 1,(IAn[4]- IAn[3])/2 , p_0))*
                    (pbinom(b,IAn[1]/2, p_1) - pbinom(b - 1,IAn[1]/2, p_1))*(pbinom(d,(IAn[2]- IAn[1])/2, p_1) - pbinom(d - 1,(IAn[2]- IAn[1])/2, p_1))*(pbinom(f,(IAn[3]- IAn[2])/2, p_1) - pbinom(f - 1,(IAn[3]- IAn[2])/2, p_1))*(pbinom(h,(IAn[4]- IAn[3])/2, p_1) - pbinom(h - 1,(IAn[4]- IAn[3])/2 , p_1))
                  
                  prob_df$Control_Success_1[count] = a
                  prob_df$Control_Success_2[count] = c
                  prob_df$Control_Success_3[count] = e
                  prob_df$Control_Success_4[count] = g
                  
                  prob_df$Treat_Success_1[count] = b
                  prob_df$Treat_Success_1[count] = d
                  prob_df$Treat_Success_1[count] = f
                  prob_df$Treat_Success_1[count] = h
                  
                  prob_df$Probability[count] = temp_prob
                  
                }
              }
              
              
            }
          }
          
        }
      }
      
      
    }
  }
  
  return(prob_df)
  
}

#T = Sys.time()
#answer_vec = both_prob_end_bop2_4IA_par(0.2,0.2,0.91,0.93, IAn = c(20,40,60,80))
#T = Sys.time() - T


BOP2_ess_plot = function(n_min = 20, n_total = 80){
  
  ESS_vec = numeric(0)
  
  for (i in (n_min/2):(n_total/2)){
    ESS_vec = append(ESS_vec, i*2 + (n_total - i*2)*(1 - prob_end_bop2(i, 0.4,0.5,0.91,0.94, n_total = n_total)))
  }
  
  #browser()
  
  plot(x = ((n_min/2):(n_total/2))*2, y = ESS_vec, xlab = "IA placement", ylab = "ESS", main = "H1: 0.4 vs 0.5 treatment")
}

heatmap_exact_plot = function(min_lambda = 0.85, max_lambda = 0.95, lambda_res = 0.005, min_gamma = 0.1, max_gamma = 1, gamma_res = 0.05, p_0 =0.2, p_1 = 0.2, n_total = 80, IAn = c(60,80) ){
  # Plots type one error against lambda and gamma variation
  
  lambda_seq = seq(min_lambda, max_lambda, lambda_res)
  gamma_seq = seq(min_gamma, max_gamma, gamma_res)
  
  length_df = length(lambda_seq)*length(gamma_seq)
  
  temp_n = length(lambda_seq)
  
  lambda_seq = rep(lambda_seq, length(gamma_seq))
  gamma_seq = rep(gamma_seq, each = temp_n)
  
  NR_df = data.frame("NR" = numeric(length_df),"ESS" = numeric(length_df), "Lambda" = lambda_seq, "Gamma" = gamma_seq)
  
  for (i in 1:length_df){
    print(i)
    
    temp_data = both_prob_end_bop2_4IA(p_0, p_1, lambda = lambda_seq[i], gamma = gamma_seq[i], n_total, IAn)
    
    temp_NR = sum(temp_data[(length(IAn)+1):(2*length(IAn))])
    
    temp_ESS = sum(IAn*(temp_data[1:length(IAn)] + temp_data[(length(IAn)+1):(2*length(IAn))]))
    
    NR_df$NR[i] = temp_NR
    
    NR_df$ESS[i] = temp_ESS
    
    #browser()
    
  }
  
  if (p_0 == p_1){
  
    ggp <- ggplot(NR_df, aes(Lambda, Gamma)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = NR)) +
      labs(
        title = "Type one error variation in BOP2 - exact"
      ) +
      scale_fill_gradient2(low="navy", mid="white", high="red", 
                           midpoint=0.1, limits=c(0.05,0.15))
    print(ggp)
    
    high_lim_ESS = max(NR_df$ESS)
    low_lim_ESS = min(NR_df$ESS)
    
    ggp <- ggplot(NR_df, aes(Lambda, Gamma)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = ESS)) +
      labs(
        title = "ESS in H0 variation in BOP2 - exact"
      ) +
      scale_fill_gradient2(low="green",high="purple", limits = c(low_lim_ESS, high_lim_ESS), midpoint = (low_lim_ESS + high_lim_ESS)/2)
    print(ggp)
  
  }
  else{
    
    high_lim_power = max(NR_df$NR)
    low_lim_power = min(NR_df$NR)
    
    ggp <- ggplot(NR_df, aes(Lambda, Gamma)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = NR)) +
      labs(
        title = "Power variation in BOP2 - exact"
      ) +
      scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = (low_lim_power + high_lim_power)/2, limits = c(low_lim_power,high_lim_power))
    print(ggp)
    
    high_lim_ESS = max(NR_df$ESS)
    low_lim_ESS = min(NR_df$ESS)
    
    ggp <- ggplot(NR_df, aes(Lambda, Gamma)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = ESS)) +
      labs(
        title = "ESS in H1 variation in BOP2 - exact"
      ) +
      scale_fill_gradient2(low="green",high="purple", limits = c(low_lim_ESS, high_lim_ESS), midpoint = (low_lim_ESS + high_lim_ESS)/2)
    print(ggp)
    
    #browser()
    
  }
  
  
}

#heatmap_exact_plot()
#heatmap_exact_plot(p_1 = 0.4)

end_calc = function(fute_prob, supe_prob){
  # Gives the probabilities of ending at certain IA
  # currently only takes 4 or less
  
  end_list = vector("list",3)
  
  num_i = length(fute_prob)
  
  if (num_i == 1){
    
    end_list[[1]] = fute_prob
    end_list[[2]] = supe_prob
    end_list[[3]] = fute_prob + supe_prob
    
    
  }
  else if (num_i == 2){
    end_list[[1]] = c(fute_prob[1], (1 - fute_prob[1] - supe_prob[1])*fute_prob[2])
    end_list[[2]] = c(supe_prob[1], (1 - fute_prob[1] - supe_prob[1])*supe_prob[2])
    end_list[[3]] = end_list[[1]] + end_list[[2]]
    
    
    
  }
  else if (num_i == 3){
    end_list[[1]] = c(fute_prob[1], (1 - fute_prob[1] - supe_prob[1])*fute_prob[2], (1 - fute_prob[1] - supe_prob[1])*(1 - fute_prob[2] - supe_prob[2])*fute_prob[3])
    end_list[[2]] = c(supe_prob[1], (1 - fute_prob[1] - supe_prob[1])*supe_prob[2], (1 - fute_prob[1] - supe_prob[1])*(1 - fute_prob[2] - supe_prob[2])*supe_prob[3])
    end_list[[3]] = end_list[[1]] + end_list[[2]]
    
    
    
  }
  else if (num_i == 4){
    end_list[[1]] = c(fute_prob[1], (1 - fute_prob[1] - supe_prob[1])*fute_prob[2], (1 - fute_prob[1] - supe_prob[1])*(1 - fute_prob[2] - supe_prob[2])*fute_prob[3], (1 - fute_prob[1] - supe_prob[1])*(1 - fute_prob[2] - supe_prob[2])*(1 - fute_prob[3] - supe_prob[3])*fute_prob[4])
    end_list[[2]] = c(supe_prob[1], (1 - fute_prob[1] - supe_prob[1])*supe_prob[2], (1 - fute_prob[1] - supe_prob[1])*(1 - fute_prob[2] - supe_prob[2])*supe_prob[3], (1 - fute_prob[1] - supe_prob[1])*(1 - fute_prob[2] - supe_prob[2])*(1 - fute_prob[3] - supe_prob[3])*supe_prob[4])
    end_list[[3]] = end_list[[1]] + end_list[[2]]
    
    
    
  }
  
  
  return(end_list)
  
}

OC_exact_approx_calc = function(IAn = c(20,40,60,80), n_total = 80, p_0 = 0.2, p_1 = 0.4, lambda = 0.91, gamma = 0.93){
  
  supe_prob = numeric(length(IAn))
  fute_prob = numeric(length(IAn))

  
  for(i in 1:length(IAn)){
    temp_prob = both_prob_end_bop2(IAn[i]/2,p_0, p_1, lambda, gamma, n_total)
    
    fute_prob[i] = temp_prob[1]
    supe_prob[i] = temp_prob[2]
    
  }
  
  ESS = 0
  NR = 0
  
  end_probs = end_calc(fute_prob, supe_prob)
  
  NR = sum(end_probs[[2]])
  
  ESS = sum(IAn * end_probs[[3]])
  
  #browser()
  
  return(c("ESS" = ESS, "NR" = NR))
  
}

OC_exact_calc_1IA = function(IAn = c(60,80), n_total = 80, p_0 = 0.2, p_1 = 0.4, lambda = 0.91, gamma = 0.93){
  
  prob = both_prob_end_bop2_1IA(p_0 = p_0, p_1 = p_1, lambda = lambda, gamma = gamma, n_total = n_total, IAn = IAn)
  
  NR = prob[3] + prob[4]
  ESS = (prob[1] + prob[3])*IAn[1] + (prob[2] + prob[4])*IAn[2]
  
  return(c("ESS" = ESS, "NR" = NR))
  
}

exact_optim = function(IAn = c(60,80), n_total = 80, p_0 = 0.2, p_1 = 0.4, max_power = TRUE, max_ESS_H0 = TRUE, max_ESS_H1 = TRUE){
  lambda_vec = seq(0.85,0.925, 0.005)
  gamma_vec = seq(0.05,1, 0.05)
  
  lambda_df_vec = rep(lambda_vec, length(gamma_vec))
  gamma_df_vec = rep(gamma_vec, length(lambda_vec))
  
  n_sim = length(lambda_df_vec)
  
  optim_df = data.frame("Lambda" = lambda_df_vec, "Gamma" = gamma_df_vec, "Power" = numeric(n_sim), "ESS_H1" = numeric(n_sim), "Alpha" = numeric(n_sim), "ESS_H0" = numeric(n_sim))
  
  for(i in 1:n_sim){
    
    temp_H0 = OC_exact_calc_1IA(IAn = IAn, n_total = n_total, p_0 = p_0, p_1 = p_0, lambda = lambda_df_vec[i], gamma = gamma_df_vec[i])
    temp_H1 = OC_exact_calc_1IA(IAn = IAn, n_total = n_total, p_0 = p_0, p_1 = p_1, lambda = lambda_df_vec[i], gamma = gamma_df_vec[i])
    
    optim_df$Power[i] = temp_H1[2]
    optim_df$ESS_H1[i] = temp_H1[1]
    
    optim_df$Alpha[i] = temp_H0[2]
    optim_df$ESS_H0[i] = temp_H0[1]
    
  }
  
  if (max_power){
    #browser()
    
    # We want to optimise power
    index = which.max(optim_df$Power[which(optim_df$Alpha < 0.1)])
    lambda = optim_df$Lambda[which(optim_df$Alpha < 0.1)][index]
    gamma = optim_df$Gamma[which(optim_df$Alpha < 0.1)][index]
    
    print(optim_df[which(optim_df$Alpha < 0.1),][index,])
    
  }
  
  if (max_ESS_H0){
    #browser()
    
    # We want to optimise power
    index = which.min(optim_df$ESS_H0[which(abs(optim_df$Alpha - 0.095) < 0.005)]) # abs here to ensure powe type one error is close to target and thus power reasonable
    lambda = optim_df$Lambda[which(abs(optim_df$Alpha - 0.095) < 0.005)][index]
    gamma = optim_df$Gamma[which(abs(optim_df$Alpha - 0.095) < 0.005)][index]
    
    print(optim_df[which(abs(optim_df$Alpha - 0.095) < 0.005),][index,])
    
  }
  
  if (max_ESS_H1){
    #browser()
    
    # We want to optimise power
    index = which.min(optim_df$ESS_H1[which(abs(optim_df$Alpha - 0.095) < 0.005)]) # abs here to ensure powe type one error is close to target and thus power reasonable
    lambda = optim_df$Lambda[which(abs(optim_df$Alpha - 0.095) < 0.005)][index]
    gamma = optim_df$Gamma[which(abs(optim_df$Alpha - 0.095) < 0.005)][index]
    
    print(optim_df[which(abs(optim_df$Alpha - 0.095) < 0.005),][index,])
    
  }
  
  return(c("Lambda" = lambda, "Gamma" = gamma))
  
  
}

exact_optim_par = function(IAn = c(60,80), n_total = 80, p_0 = 0.2, p_1 = 0.4, max_power = TRUE, max_ESS_H0 = FALSE, max_ESS_H1 = FALSE){
  lambda_vec = seq(0.85,0.925, 0.0025)
  gamma_vec = seq(0.05,1, 0.025)
  
  lambda_df_vec = rep(lambda_vec, length(gamma_vec))
  gamma_df_vec = rep(gamma_vec, length(lambda_vec))
  
  n_sim = length(lambda_df_vec)
  
  optim_df = data.frame("Lambda" = lambda_df_vec, "Gamma" = gamma_df_vec, "Power" = numeric(n_sim), "ESS_H1" = numeric(n_sim), "Alpha" = numeric(n_sim), "ESS_H0" = numeric(n_sim))
  
  par_results = vector("list",n_sim)
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  par_results = foreach( i = 1:n_sim, .export = c("rcpp_exact_beta","Threshold_prob","OC_exact_calc_1IA","both_prob_end_bop2_1IA"))  %dopar% {
    
    #for(i in 1:n_sim){
    
    temp_H0 = OC_exact_calc_1IA(IAn = IAn, n_total = n_total, p_0 = p_0, p_1 = p_0, lambda = lambda_df_vec[i], gamma = gamma_df_vec[i])
    temp_H1 = OC_exact_calc_1IA(IAn = IAn, n_total = n_total, p_0 = p_0, p_1 = p_1, lambda = lambda_df_vec[i], gamma = gamma_df_vec[i])
    
    #optim_df$Power[i] = temp_H1[2]
    #optim_df$ESS_H1[i] = temp_H1[1]
    
    #optim_df$Alpha[i] = temp_H0[2]
    #optim_df$ESS_H0[i] = temp_H0[1]
    
    return(c(temp_H1,temp_H0))
    
  }
  
  stopCluster(cl) # deactivates
  
  for (i in 1:n_sim){
    optim_df$Power[i] = par_results[[i]][2]
    optim_df$ESS_H1[i] = par_results[[i]][1]
    
    optim_df$Alpha[i] = par_results[[i]][4]
    optim_df$ESS_H0[i] = par_results[[i]][3]
    
  }
  
  if (max_power){
    #browser()
    
    # We want to optimise power
    index = which.max(optim_df$Power[which(optim_df$Alpha < 0.1)])
    lambda = optim_df$Lambda[which(optim_df$Alpha < 0.1)][index]
    gamma = optim_df$Gamma[which(optim_df$Alpha < 0.1)][index]
    
    #print(optim_df[which(optim_df$Alpha < 0.1),][index,])
    
  }
  
  if (max_ESS_H0){
    #browser()
    
    # We want to optimise power
    index = which.min(optim_df$ESS_H0[which(abs(optim_df$Alpha - 0.095) < 0.005)]) # abs here to ensure powe type one error is close to target and thus power reasonable
    lambda = optim_df$Lambda[which(abs(optim_df$Alpha - 0.095) < 0.005)][index]
    gamma = optim_df$Gamma[which(abs(optim_df$Alpha - 0.095) < 0.005)][index]
    
    #print(optim_df[which(abs(optim_df$Alpha - 0.095) < 0.005),][index,])
    
  }
  
  if (max_ESS_H1){
    #browser()
    
    # We want to optimise power
    index = which.min(optim_df$ESS_H1[which(abs(optim_df$Alpha - 0.095) < 0.005)]) # abs here to ensure powe type one error is close to target and thus power reasonable
    lambda = optim_df$Lambda[which(abs(optim_df$Alpha - 0.095) < 0.005)][index]
    gamma = optim_df$Gamma[which(abs(optim_df$Alpha - 0.095) < 0.005)][index]
    
    #print(optim_df[which(abs(optim_df$Alpha - 0.095) < 0.005),][index,])
    
  }
  
  return(c("Lambda" = lambda, "Gamma" = gamma))
  
  
}

#exact_optim_par(c(38,80))
#exact_optim_par(c(40,80))

#browser()

vary_plot = function(k){
  p = (1:((1-k)*100))/100
  
  y = sqrt(2*p*(1-p) + k*(1-2*p)+k^2)
  
  plot(x=p,y=y)
  
}

#T = Sys.time()
#exact_optim_par()
#T = Sys.time() - T

#heatmap_exact_plot(IAn = c(40,80), p_1 = 0.4)
#heatmap_exact_plot(IAn = c(38,80), p_1 = 0.4)

exact_one_arm_case = function(p_eff = 0.66, p_tox = 0.29, IAn = c(10,20,35,50), n_total = 50, lambda = 0.675, gamma = 1){
  # Case study for the FRAIL M study, will calculate the probabilities here
  
  num_IA = length(IAn)
  fute_early_prob = numeric(num_IA)
  supe_early_prob = numeric(num_IA)
  

  
  thresh_prob = numeric(2*num_IA)
  
  for(i in 1:num_IA){
    temp_thresh = Threshold_prob(lambda, gamma, n_total, IAn[i])
    
    thresh_prob[i] = temp_thresh[1]
    thresh_prob[i + num_IA] = temp_thresh[2]
    
  }
  
  
  for (a in 0:IAn[1]){
    
    #print(a)
    
    #browser()
    
    #if (supe_early_prob + fute_early_prob != pbinom(a-1, IAn[1],p_eff)){
    #  browser()
    #}
    
    for (b in 0:IAn[1]){
      
      p_eff_comp = 1 - pbeta(0.5,1 + a, 1 + IAn[1] - a) # Probability beta r.v. > 0.5
      p_tox_comp = pbeta(0.35, 1 + b, 1 + IAn[1] - b) # Probability beta r.v. < 0.35
      
      # If p_eff*p_tox is sufficiently high, we end early
      # Else, if it's sufficidely low, we end early
      # Else, we continue to the next IA
      
      p_comp = p_eff_comp*p_tox_comp
      #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[1] )
      
      if (p_comp > thresh_prob[1 + num_IA]){
        
        # End early supe, add probability
        supe_early_prob[1] = supe_early_prob[1] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))
        
      }
      else if (p_comp < thresh_prob[1]){
        
        # End early futility, add probability
        
        fute_early_prob[1] = fute_early_prob[1] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))
        
      }
      else{
        # On to next IA, another layer deep
        pts_2 = IAn[2] - IAn[1] # how many patients are recruited for IA2
        
        for (c in 0:pts_2){
          for (d in 0:pts_2){
            
            p_eff_comp = 1 - pbeta(0.5,1 + a + c, 1 + IAn[2] - a - c) # Probability beta r.v. > 0.5
            p_tox_comp = pbeta(0.35, 1 + b + d, 1 + IAn[2] - b - d) # Probability beta r.v. < 0.35
            
            # If p_eff*p_tox is sufficiently high, we end early
            # Else, if it's sufficidely low, we end early
            # Else, we continue to the next IA
            
            p_comp = p_eff_comp*p_tox_comp
            #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[2] )
            
            if (p_comp > thresh_prob[2 + num_IA]){
              
              # End early supe, add probability
              supe_early_prob[2] = supe_early_prob[2] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))
              
            }
            else if (p_comp < thresh_prob[2]){
              
              # End early futility, add probability
              
              fute_early_prob[2] = fute_early_prob[2] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))
              
            }
            else{
              
              pts_3 = IAn[3] - IAn[2]
              
              for (e in 0:pts_3){
                for (f in 0:pts_3){
                  
                  p_eff_comp = 1 - pbeta(0.5,1 + a + c + e, 1 + IAn[3] - a - c - e) # Probability beta r.v. > 0.5
                  p_tox_comp = pbeta(0.35, 1 + b + d + f, 1 + IAn[3] - b - d - f) # Probability beta r.v. < 0.35
                  
                  # If p_eff*p_tox is sufficiently high, we end early
                  # Else, if it's sufficidely low, we end early
                  # Else, we continue to the next IA
                  
                  p_comp = p_eff_comp*p_tox_comp
                  #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[3] )
                  
                  if (p_comp > thresh_prob[3 + num_IA]){
                    
                    # End early supe, add probability
                    supe_early_prob[3] = supe_early_prob[3] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))*(pbinom(e,pts_3,p_eff ) - pbinom(e - 1,pts_3,p_eff ))*(pbinom(f,pts_3,p_tox ) - pbinom(f - 1,pts_3,p_tox ))
                    
                  }
                  else if (p_comp < thresh_prob[3]){
                    
                    # End early futility, add probability
                    
                    fute_early_prob[3] = fute_early_prob[3] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))*(pbinom(e,pts_3,p_eff ) - pbinom(e - 1,pts_3,p_eff ))*(pbinom(f,pts_3,p_tox ) - pbinom(f - 1,pts_3,p_tox ))
                    
                  }
                  else{
                    
                    pts_4 = IAn[4] - IAn[3]
                    
                    for(g in 0:pts_4){
                      for(h in 0:pts_4){
                        
                        p_eff_comp = 1 - pbeta(0.5,1 + a + c + e + g, 1 + IAn[4] - a - c - e - g) # Probability beta r.v. > 0.5
                        p_tox_comp = pbeta(0.35, 1 + b + d + f + h, 1 + IAn[4] - b - d - f - h) # Probability beta r.v. < 0.35
                        
                        # If p_eff*p_tox is sufficiently high, we end early
                        # Else, if it's sufficidely low, we end early
                        # Else, we continue to the next IA
                        
                        p_comp = p_eff_comp*p_tox_comp
                        #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[4] )
                        
                        if (p_comp > thresh_prob[4 + num_IA]){
                          
                          # End early supe, add probability
                          supe_early_prob[4] = supe_early_prob[4] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))*(pbinom(e,pts_3,p_eff ) - pbinom(e - 1,pts_3,p_eff ))*(pbinom(f,pts_3,p_tox ) - pbinom(f - 1,pts_3,p_tox ))*(pbinom(g,pts_4,p_eff ) - pbinom(g - 1,pts_4,p_eff ))*(pbinom(h,pts_4,p_tox ) - pbinom(h - 1,pts_4,p_tox ))
                          
                        }
                        else if (p_comp < thresh_prob[4]){
                          
                          # End early futility, add probability
                          
                          fute_early_prob[4] = fute_early_prob[4] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))*(pbinom(e,pts_3,p_eff ) - pbinom(e - 1,pts_3,p_eff ))*(pbinom(f,pts_3,p_tox ) - pbinom(f - 1,pts_3,p_tox ))*(pbinom(g,pts_4,p_eff ) - pbinom(g - 1,pts_4,p_eff ))*(pbinom(h,pts_4,p_tox ) - pbinom(h - 1,pts_4,p_tox ))
                          
                        }
                        else{
                          browser()
                        }
                        
                      }
                    }
                  }
                }
                
                
              }
            
            }
            
          }
        }
        
        
        
        
      }
    }
  }
  
  return(c("Futility" = fute_early_prob, "Supe" = supe_early_prob))
}

exact_one_arm_case_par = function(p_eff = 0.66, p_tox = 0.29, IAn = c(10,20,35,50), n_total = 50, lambda = 0.675, gamma = 1){
  # Case study for the FRAIL M study, will calculate the probabilities here
  
  num_IA = length(IAn)
  fute_early_prob = numeric(num_IA)
  supe_early_prob = numeric(num_IA)
  
  thresh_prob = numeric(2*num_IA)
  
  for(i in 1:num_IA){
    temp_thresh = Threshold_prob(lambda, gamma, n_total, IAn[i])
    
    thresh_prob[i] = temp_thresh[1]
    thresh_prob[i + num_IA] = temp_thresh[2]
    
  }
  
  
  par_results = vector("list",(IAn[1]+1))
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  par_results = foreach( a = 0:IAn[1], .export = c("rcpp_exact_beta","Threshold_prob","exact_OC_case","exact_one_arm_case"))  %dopar% {
  
  #for (a in 0:IAn[1]){
    par_return = numeric(2*num_IA)
    #print(a)
    
    #browser()
    
    #if (supe_early_prob + fute_early_prob != pbinom(a-1, IAn[1],p_eff)){
    #  browser()
    #}
    
    for (b in 0:IAn[1]){
      
      p_eff_comp = 1 - pbeta(0.5,1 + a, 1 + IAn[1] - a) # Probability beta r.v. > 0.5
      p_tox_comp = pbeta(0.35, 1 + b, 1 + IAn[1] - b) # Probability beta r.v. < 0.35
      
      # If p_eff*p_tox is sufficiently high, we end early
      # Else, if it's sufficidely low, we end early
      # Else, we continue to the next IA
      
      p_comp = p_eff_comp*p_tox_comp
      #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[1] )
      
      if (p_comp > thresh_prob[1 + num_IA]){
        
        # End early supe, add probability
        par_return[num_IA + 1] = par_return[num_IA + 1] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))
        
        
      }
      else if (p_comp < thresh_prob[1]){
        
        # End early futility, add probability
        
        par_return[1] = par_return[1] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))
        
      }
      else{
        # On to next IA, another layer deep
        pts_2 = IAn[2] - IAn[1] # how many patients are recruited for IA2
        
        for (c in 0:pts_2){
          for (d in 0:pts_2){
            
            p_eff_comp = 1 - pbeta(0.5,1 + a + c, 1 + IAn[2] - a - c) # Probability beta r.v. > 0.5
            p_tox_comp = pbeta(0.35, 1 + b + d, 1 + IAn[2] - b - d) # Probability beta r.v. < 0.35
            
            # If p_eff*p_tox is sufficiently high, we end early
            # Else, if it's sufficidely low, we end early
            # Else, we continue to the next IA
            
            p_comp = p_eff_comp*p_tox_comp
            #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[2] )
            
            if (p_comp > thresh_prob[2 + num_IA]){
              
              # End early supe, add probability
              par_return[num_IA + 2] = par_return[num_IA + 2] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))
              
            }
            else if (p_comp < thresh_prob[2]){
              
              # End early futility, add probability
              
              par_return[2] = par_return[2] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))
              
            }
            else{
              
              pts_3 = IAn[3] - IAn[2]
              
              for (e in 0:pts_3){
                for (f in 0:pts_3){
                  
                  p_eff_comp = 1 - pbeta(0.5,1 + a + c + e, 1 + IAn[3] - a - c - e) # Probability beta r.v. > 0.5
                  p_tox_comp = pbeta(0.35, 1 + b + d + f, 1 + IAn[3] - b - d - f) # Probability beta r.v. < 0.35
                  
                  # If p_eff*p_tox is sufficiently high, we end early
                  # Else, if it's sufficidely low, we end early
                  # Else, we continue to the next IA
                  
                  p_comp = p_eff_comp*p_tox_comp
                  #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[3] )
                  
                  if (p_comp > thresh_prob[3 + num_IA]){
                    
                    # End early supe, add probability
                    par_return[num_IA + 3] = par_return[num_IA + 3] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))*(pbinom(e,pts_3,p_eff ) - pbinom(e - 1,pts_3,p_eff ))*(pbinom(f,pts_3,p_tox ) - pbinom(f - 1,pts_3,p_tox ))
                    
                  }
                  else if (p_comp < thresh_prob[3]){
                    
                    # End early futility, add probability
                    
                    par_return[3] = par_return[3] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))*(pbinom(e,pts_3,p_eff ) - pbinom(e - 1,pts_3,p_eff ))*(pbinom(f,pts_3,p_tox ) - pbinom(f - 1,pts_3,p_tox ))
                    
                  }
                  else{
                    
                    pts_4 = IAn[4] - IAn[3]
                    
                    for(g in 0:pts_4){
                      for(h in 0:pts_4){
                        
                        p_eff_comp = 1 - pbeta(0.5,1 + a + c + e + g, 1 + IAn[4] - a - c - e - g) # Probability beta r.v. > 0.5
                        p_tox_comp = pbeta(0.35, 1 + b + d + f + h, 1 + IAn[4] - b - d - f - h) # Probability beta r.v. < 0.35
                        
                        # If p_eff*p_tox is sufficiently high, we end early
                        # Else, if it's sufficidely low, we end early
                        # Else, we continue to the next IA
                        
                        p_comp = p_eff_comp*p_tox_comp
                        #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[4] )
                        
                        if (p_comp > thresh_prob[4 + num_IA]){
                          
                          # End early supe, add probability
                          par_return[num_IA + 4] = par_return[num_IA + 4] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))*(pbinom(e,pts_3,p_eff ) - pbinom(e - 1,pts_3,p_eff ))*(pbinom(f,pts_3,p_tox ) - pbinom(f - 1,pts_3,p_tox ))*(pbinom(g,pts_4,p_eff ) - pbinom(g - 1,pts_4,p_eff ))*(pbinom(h,pts_4,p_tox ) - pbinom(h - 1,pts_4,p_tox ))
                          
                        }
                        else if (p_comp < thresh_prob[4]){
                          
                          # End early futility, add probability
                          
                          par_return[4] = par_return[4] + (pbinom(a,IAn[1],p_eff ) - pbinom(a - 1,IAn[1],p_eff ))*(pbinom(b,IAn[1],p_tox ) - pbinom(b - 1,IAn[1],p_tox ))*(pbinom(c,pts_2,p_eff ) - pbinom(c - 1,pts_2,p_eff ))*(pbinom(d,pts_2,p_tox ) - pbinom(d - 1,pts_2,p_tox ))*(pbinom(e,pts_3,p_eff ) - pbinom(e - 1,pts_3,p_eff ))*(pbinom(f,pts_3,p_tox ) - pbinom(f - 1,pts_3,p_tox ))*(pbinom(g,pts_4,p_eff ) - pbinom(g - 1,pts_4,p_eff ))*(pbinom(h,pts_4,p_tox ) - pbinom(h - 1,pts_4,p_tox ))
                          
                        }
                        else{
                          #browser()
                        }
                        
                      }
                    }
                  }
                }
                
                
              }
              
            }
            
          }
        }
        
        
        
        
      }
    }
    
    return(par_return)
  }
  
  for (i in 0:IAn[1]){
    fute_early_prob = fute_early_prob + par_results[[i+1]][1:num_IA]
    supe_early_prob = supe_early_prob + par_results[[i+1]][(1+num_IA):(2*num_IA)]
  }
  
  return(c("Futility" = fute_early_prob, "Supe" = supe_early_prob))
}

exact_OC_case = function(p_eff = 0.66, p_tox = 0.29, IAn = c(10,20,35,50), n_total = 50, lambda = 0.675, gamma = 1, par = FALSE){
  
  if (par){
    probs = exact_one_arm_case_par(p_eff = p_eff, p_tox = p_tox, IAn = IAn, n_total = n_total, lambda = lambda, gamma = gamma)
  }
  else{
    probs = exact_one_arm_case(p_eff = p_eff, p_tox = p_tox, IAn = IAn, n_total = n_total, lambda = lambda, gamma = gamma)
  }
  
  #browser()
  
  num_IA = length(IAn)
  
  NR = sum(probs[(num_IA + 1):(2*num_IA)])
  
  End_IA_prob = probs[(num_IA + 1):(2*num_IA)] + probs[1:num_IA]
  
  ESS = sum(IAn*End_IA_prob)
  
  return(c("NR" = NR, "ESS" = ESS))
}

sim_case = function(p_eff = 0.96, p_tox = 0.41, IAn = c(10,20,35,50), n_total = 50, lambda = 0.675, gamma = 1, thresh_prob = c(0.6,0.7,0.8,0.9,0.96,0.94,0.92,0.9)){
  
  eff_vec = rbinom(n_total, 1, p_eff)
  tox_vec = rbinom(n_total, 1, p_tox)
  
  num_IA = length(IAn)
  
  for(i in 1:num_IA){
  
    #thresh_prob = Threshold_prob(lambda, gamma, n_total, IAn[i])
    
    num_succ = sum(eff_vec[1:IAn[i]])
    num_tox = sum(tox_vec[1:IAn[i]])
    
    p_comp = (1 - pbeta(0.5, 1 + num_succ, 1 + IAn[i] - num_succ))*pbeta(0.35, 1 + num_tox, 1 + IAn[i] - num_tox)
    
    if (p_comp > thresh_prob[i + num_IA]){
      
      # End early supe, add probability
      return(c("NR" = 1, "SS" = IAn[i]))
      
    }
    else if (p_comp < thresh_prob[i]){
      
      # End early futility, add probability
      
      return(c("NR" = 0, "SS" = IAn[i]))
      
    }
  
  }
  
}

sim_OC_case = function(p_eff = 0.66, p_tox = 0.29, IAn = c(10,20,35,50), n_total = 50, lambda = 0.675, gamma = 1, n_sim = 10000, paired_seed = FALSE, set_seed = 2025){
  
  NR = 0
  ESS = 0
  
  num_IA = length(IAn)
  
  thresh_prob = numeric(2*num_IA)
  
  for(i in 1:num_IA){
    temp_thresh = Threshold_prob(lambda, gamma, n_total, IAn[i])
    
    thresh_prob[i] = temp_thresh[1]
    thresh_prob[i + num_IA] = temp_thresh[2]
    
  }
  
  for(i in 1:n_sim){
    
    if(paired_seed){
      set.seed(set_seed+1)
    }
    
    temp_result = sim_case(p_eff = p_eff, p_tox = p_tox, IAn = IAn, n_total = n_total, lambda = lambda, gamma = gamma, thresh_prob = thresh_prob)
    
    NR = temp_result[1] + NR
    ESS = temp_result[2] + ESS
    
  }
  
  return_vec = c(NR, ESS)/n_sim
  
  return(return_vec)
  
}

optim_case_par = function(p_eff_0 = 0.5, p_eff_1 = 0.66, p_tox_0 = 0.35, p_tox_1 = 0.29, n_total = 50, IAn = c(10,20,35,50), alpha = 0.05){
  # Optimises the case study problem with exact, parallel calculations
  
  gamma_vec = seq(0.1,1,0.1)
  lambda_vec = seq(0.45,0.7,0.025) # done using heuristics , seeing what gives close to 0.05
  
  lambda_df_vec = rep(lambda_vec, length(gamma_vec))
  gamma_df_vec = rep(gamma_vec, length(lambda_vec))
  
  n_sim = length(lambda_df_vec)
  
  optim_df = data.frame("Lambda" = lambda_df_vec, "Gamma" = gamma_df_vec, "Power" = numeric(n_sim), "ESS_H1" = numeric(n_sim), "Alpha" = numeric(n_sim), "ESS_H0" = numeric(n_sim))
  
  par_results = vector("list",n_sim)
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  par_results = foreach( i = 1:n_sim, .export = c("rcpp_exact_beta","Threshold_prob","exact_OC_case","exact_one_arm_case"))  %dopar% {
    
    temp_H0 = exact_OC_case(IAn = IAn, n_total = n_total, p_eff = p_eff_0, p_tox = p_tox_0 , lambda = lambda_df_vec[i], gamma = gamma_df_vec[i])
    temp_H1 = c(0,0)
    
    if (temp_H0[1] < alpha){
      temp_H1 = exact_OC_case(IAn = IAn, n_total = n_total, p_eff = p_eff_1, p_tox = p_tox_1, lambda = lambda_df_vec[i], gamma = gamma_df_vec[i])
    }
    
    
    
    #optim_df$Power[i] = temp_H1[2]
    #optim_df$ESS_H1[i] = temp_H1[1]
    
    #optim_df$Alpha[i] = temp_H0[2]
    #optim_df$ESS_H0[i] = temp_H0[1]
    
    return(c(temp_H1,temp_H0))
    
  }
  
  stopCluster(cl) # deactivates
  
  for (i in 1:n_sim){
    optim_df$Power[i] = par_results[[i]][1]
    optim_df$ESS_H1[i] = par_results[[i]][2]
    
    optim_df$Alpha[i] = par_results[[i]][3]
    optim_df$ESS_H0[i] = par_results[[i]][4]
    
  }
  
  if (TRUE){
    #browser()
    
    # We want to optimise power
    index = which.max(optim_df$Power)
    lambda = optim_df$Lambda[index]
    gamma = optim_df$Gamma[index]
    
    #print(optim_df[which(optim_df$Alpha < 0.1),][index,])
    
  }
  
  return_vec = c("Lambda" = lambda, "Gamma" = gamma)
  
  print(return_vec)
  return(return_vec)
  
  
}

IA_optimiser_case = function(p_eff = 0.66, p_tox = 0.29, n_total = 50, lambda = 0.675, gamma = 1, n_sim = 10000){
  # Only works with 4 IA
  
  #browser()
  
  num_iter = (n_total-1)*(n_total-2)*(n_total - 3)/6
  
  ESS_vec = numeric(num_iter)
  
  IA_list = vector("list", num_iter)
  
  count = 1
  
  for (a in 3:(n_total-1)){
    #print(a)
    for (b in 2:(a-1)){
      for (c in 1:(b-1)){
        
        # Iterates over all the possible placements for the three IA
        ESS_vec[count] = sim_OC_case(p_eff = p_eff, p_tox = p_tox, n_total = n_total, lambda = lambda, gamma= gamma, n_sim = n_sim, IAn = c(c,b,a,n_total))[2]
        
        IA_list[[count]] = c(c,b,a,n_total)
        
        count = count+1
        
      }
    
    }
  }
  
  index = which.min(ESS_vec)
  
  return(c("IA" = IA_list[[index]],"ESS" = ESS_vec[index]))
  
}

IA_optimiser_case_par = function(p_eff = 0.66, p_tox = 0.29, n_total = 50, lambda = 0.675, gamma = 1, n_sim = 10000){
  # Only works with 4 IA, parallel version
  
  #browser()
  
  num_iter = (n_total-1)*(n_total-2)*(n_total - 3)/6
  
  ESS_vec = numeric(num_iter)
  
  IA_list = vector("list", num_iter)
  
  #count = 1
  
  par_results = vector("list",num_iter)
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  par_results = foreach( a = 3:(n_total-1), .export = c("rcpp_exact_beta","Threshold_prob","sim_OC_case","sim_case"))  %dopar% {
  
    iter_a = (a-2)*(a-1)/2 # how many are in the sub iteration
    to_return_par = vector("list", iter_a)
    count = 1
    
  #for (a in 3:(n_total-1)){
    #print(a)
    for (b in 2:(a-1)){
      for (c in 1:(b-1)){
        
        # Iterates over all the possible placements for the three IA
        temp_ESS = sim_OC_case(p_eff = p_eff, p_tox = p_tox, n_total = n_total, lambda = lambda, gamma= gamma, n_sim = n_sim, IAn = c(c,b,a,n_total))[2]
        
        temp_IA = c(c,b,a,n_total)
        
        to_return_par[[count]] = c(temp_IA,temp_ESS)
        
        count = count+1
        
      }
      
    }
    
    return(to_return_par)
  }
  
  stopCluster(cl) # deactivates
  
  count = 1
  
  #browser()
  for (a in 3:(n_total-1)){
    iter_a = (a-2)*(a-1)/2 # how many are in the sub iteration
  
    
    for (i in 1:iter_a){
      #print(c(a,i))
      ESS_vec[count] = par_results[[a-2]][[i]][5]
      IA_list[[count]] = par_results[[a-2]][[i]][1:4]
      
      count = count + 1
    }
  }
  
  browser()
  
  index = which.min(ESS_vec)
  
  return(c("IA" = IA_list[[index]],"ESS" = ESS_vec[index]))
  
}

IA_optimiser_case_par_3IA = function(p_eff = 0.66, p_tox = 0.29, n_total = 50, lambda = 0.675, gamma = 1, n_sim = 10000){
  # Only works with 3 IA, parallel version
  
  #browser()
  
  num_iter = (n_total-1)*(n_total-2)/2
  
  ESS_vec = numeric(num_iter)
  
  IA_list = vector("list", num_iter)
  
  #count = 1
  
  par_results = vector("list",num_iter)
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  par_results = foreach( a = 2:(n_total-1), .export = c("rcpp_exact_beta","Threshold_prob","sim_OC_case","sim_case"))  %dopar% {
    
    iter_a = a-1 # how many are in the sub iteration
    to_return_par = vector("list", iter_a)
    count = 1
    
    #for (a in 3:(n_total-1)){
    #print(a)
    for (b in 1:(a-1)){
      #for (c in 1:(b-1)){
        
        # Iterates over all the possible placements for the three IA
        temp_ESS = sim_OC_case(p_eff = p_eff, p_tox = p_tox, n_total = n_total, lambda = lambda, gamma= gamma, n_sim = n_sim, IAn = c(b,a,n_total))[2]
        
        temp_IA = c(b,a,n_total)
        
        to_return_par[[count]] = c(temp_IA,temp_ESS)
        
        count = count+1
        
      
      
    }
    
    return(to_return_par)
  }
  
  stopCluster(cl) # deactivates
  
  count = 1
  
  #browser()
  for (a in 2:(n_total-1)){
    iter_a = (a-1) # how many are in the sub iteration
    
    
    for (i in 1:iter_a){
      #print(c(a,i))
      ESS_vec[count] = par_results[[a-1]][[i]][4]
      IA_list[[count]] = par_results[[a-1]][[i]][1:3]
      
      count = count + 1
    }
  }
  
  #browser()
  
  index = which.min(ESS_vec)
  
  return(c("IA" = IA_list[[index]],"ESS" = ESS_vec[index]))
  
}

IA_optimiser_case_par_2IA = function(p_eff = 0.66, p_tox = 0.29, n_total = 50, lambda = 0.675, gamma = 1, n_sim = 10000){
  # Only works with 3 IA, parallel version
  
  #browser()
  
  num_iter = (n_total-1)
  
  ESS_vec = numeric(num_iter)
  
  IA_list = vector("list", num_iter)
  
  #count = 1
  
  par_results = vector("list",num_iter)
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  par_results = foreach( a = 1:(n_total-1), .export = c("rcpp_exact_beta","Threshold_prob","sim_OC_case","sim_case"))  %dopar% {
    
    iter_a = 1 # how many are in the sub iteration
    to_return_par = vector("list", iter_a)
    count = 1
    
    #for (a in 3:(n_total-1)){
    #print(a)
    #for (b in 1:(a-1)){
      #for (c in 1:(b-1)){
      
      # Iterates over all the possible placements for the three IA
      temp_ESS = sim_OC_case(p_eff = p_eff, p_tox = p_tox, n_total = n_total, lambda = lambda, gamma= gamma, n_sim = n_sim, IAn = c(a,n_total))[2]
      
      temp_IA = c(a,n_total)
      
      to_return_par[[count]] = c(temp_IA,temp_ESS)
      
      count = count+1
      
      
      
    
    
    return(to_return_par)
  }
  
  stopCluster(cl) # deactivates
  
  count = 1
  
  #browser()
  for (a in 1:(n_total-1)){
    iter_a = 1 # how many are in the sub iteration
    
    
    for (i in 1:iter_a){
      #print(c(a,i))
      ESS_vec[count] = par_results[[a]][[i]][3]
      IA_list[[count]] = par_results[[a]][[i]][1:2]
      
      count = count + 1
    }
  }
  
  #browser()
  
  index = which.min(ESS_vec)
  
  return(c("IA" = IA_list[[index]],"ESS" = ESS_vec[index]))
  
}

BOP2_normal_sim = function(mu_0 = 0, sigma_0 = 2, mu_1 = 1, sigma_1 = 2, prior_mu = 0, prior_sigma = 1, n_total = 80, IAn = c(20,40,60,80), lambda = 0.9, gamma = 0.9, RAR = FALSE){
  
  p_a_1 = IAn[1]/2 # Patients per arm, first IA
  
  control = rnorm(p_a_1, mu_0, sigma_0)
  treat = rnorm(p_a_1, mu_1, sigma_1)
  
  mu_control = (((sigma_0^2)/p_a_1)*prior_mu + mean(control)*(prior_sigma^2))/((sigma_0^2)/p_a_1 + (prior_sigma^2))
  sigma_control = (((sigma_0^2)/p_a_1) * (prior_sigma^2))/((sigma_0^2)/p_a_1 + (prior_sigma^2))
  
  mu_treat = (((sigma_1^2)/p_a_1)*prior_mu + mean(treat)*(prior_sigma^2))/((sigma_1^2)/p_a_1 + (prior_sigma^2))
  sigma_treat = (((sigma_1^2)/p_a_1) * (prior_sigma^2))/((sigma_1^2)/p_a_1 + (prior_sigma^2))
  
  p_z = 1 - pnorm(0,mean = mu_treat - mu_control, sd = sqrt(sigma_treat + sigma_control)) # Probablility treatment is better than the control, which under H0 has mean 0
  
  temp_thresh = Threshold_prob(lambda, gamma, n_total, IAn[1])
  
  #browser()
  
  if (p_z > temp_thresh[2]){
    # Supe
    return(c("NR" = 1, "ESS" = IAn[1], "Prop" = length(treat)/IAn[1]))
  }
  else if (p_z < temp_thresh[1]) {
    # Fute
    return(c("NR" = 0, "ESS" = IAn[1], "Prop" = length(treat)/IAn[1]))
  }
  else{
    for (i in 2:length(IAn)){
      
      if (RAR){
        allo_prob = p_z^(IAn[i-1]/n_total)/(p_z^(IAn[i-1]/n_total) + (1-p_z)^(IAn[i-1]/n_total))
        
        #browser()
      }
      else{
        allo_prob = 0.5
      }
      
      to_sim = IAn[i] - IAn[i-1]
      
      treat_allo = round(allo_prob*to_sim)
      control_allo = to_sim - treat_allo
      
      control = append(control, rnorm(control_allo, mu_0, sigma_0))
      treat = append(treat,rnorm(treat_allo, mu_1, sigma_1))
      
      n_con = length(control)
      n_treat = length(treat)
      
      #browser()
      
      mu_control = (((sigma_0^2)/n_con)*prior_mu + mean(control)*(prior_sigma^2))/((sigma_0^2)/n_con + (prior_sigma^2))
      sigma_control = (((sigma_0^2)/n_con) * (prior_sigma^2))/((sigma_0^2)/n_con + (prior_sigma^2))
      
      mu_treat = (((sigma_1^2)/n_treat)*prior_mu + mean(treat)*(prior_sigma^2))/((sigma_1^2)/n_treat + (prior_sigma^2))
      sigma_treat = (((sigma_1^2)/n_treat) * (prior_sigma^2))/((sigma_1^2)/n_treat + (prior_sigma^2))
      
      p_z = 1 - pnorm(0,mean = mu_treat - mu_control, sd = sqrt(sigma_treat + sigma_control)) # Probablility treatment is better than the control, which under H0 has mean 0
      
      temp_thresh = Threshold_prob(lambda, gamma, n_total, IAn[i])
      
      if (p_z > temp_thresh[2]){
        # Supe
        return(c("NR" = 1, "ESS" = IAn[i], "Prop" = n_treat/IAn[i]))
      }
      else if (p_z < temp_thresh[1]) {
        # Fute
        
        #browser()
        return(c("NR" = 0, "ESS" = IAn[i], "Prop" = n_treat/IAn[i]))
      }
      
    }
  }
  
  
  
}

OC_gen_norm = function(mu_0 = 0, sigma_0 = 2, mu_1 = 1, sigma_1 = 2, prior_mu = 0, prior_sigma = 1, n_total = 80, IAn = c(20,40,60,80), lambda = 0.9, gamma = 0.9, RAR = FALSE, n_sim = 10000){
  
  NR = numeric(n_sim)
  ESS = numeric(n_sim)
  Prop = numeric(n_sim)
  
  for (i in 1:n_sim){
    
    temp_results = BOP2_normal_sim(mu_0 = mu_0, sigma_0 = sigma_0, mu_1 = mu_1, sigma_1 = sigma_1, prior_mu = prior_mu, prior_sigma = prior_sigma, n_total = n_total, IAn = IAn, lambda = lambda, gamma = gamma, RAR = RAR)
    
    NR[i] = temp_results[1]
    ESS[i] = temp_results[2]
    Prop[i] = temp_results[3]
    
        
  }
  
  sort_prop = sort(Prop)
  LQ = sort_prop[round(0.025*n_sim)]
  UQ = sort_prop[round(0.975*n_sim)]
  
  return(c("NR" = mean(NR), "ESS" = mean(ESS), "Prop" = mean(Prop), "LQ prop" = LQ, "UQ prop" = UQ))
  
}

Norm_optim_par = function(mu_0 = 0, sigma_0 = 2, mu_1 = 1, sigma_1 = 2, prior_mu = 0, prior_sigma = 1, n_total = 80, IAn = c(20,40,60,80), RAR = FALSE, n_sim = 10000, alpha = 0.1){
  
  lambda_vec = seq(0.85,0.95,0.01)
  gamma_vec = seq(0.1,1,0.1)
  
  lambda_vec_long = rep(lambda_vec, each = length(gamma_vec))
  gamma_vec_long = rep(gamma_vec, length(lambda_vec))
  
  n_iter = length(lambda_vec_long)
  
  
  par_results = numeric(n_iter) #  Will be a list, where each element in the list is a vector of 8 values, indicating the prob of stopping at each step
  
  ## Parallelisation
  cores = detectCores()
  
  cl = makeCluster(cores/2)  ## To not get overloaded
  
  registerDoParallel(cl)  # Activates 
  
  
  
  #for (i in 0:n_1_pa){
  par_results = foreach( i = 1:n_iter, .combine = c, .export = c("OC_gen_norm","BOP2_normal_sim","Threshold_prob"))  %dopar% {
   
    temp_NR_H0 = OC_gen_norm(mu_0 = mu_0, sigma_0 = sigma_0, mu_1 = mu_0, sigma_1 = sigma_0, prior_mu = prior_mu, prior_sigma = prior_sigma, n_total = n_total, IAn = IAn, lambda = lambda_vec_long[i], gamma = gamma_vec_long[i], RAR = RAR, n_sim = n_sim)[1]
    
    if (temp_NR_H0 < alpha){
      return( OC_gen_norm(mu_0 = mu_0, sigma_0 = sigma_0, mu_1 = mu_1, sigma_1 = sigma_1, prior_mu = prior_mu, prior_sigma = prior_sigma, n_total = n_total, IAn = IAn, lambda = lambda_vec_long[i], gamma = gamma_vec_long[i], RAR = RAR, n_sim = n_sim)[1])
    }
    else{
      return(0)
    }
    
  }
  
  stopCluster(cl) # deactivates
  
  index = which.max(par_results)
  
  return(c("Lambda" = lambda_vec_long[index], "Gamma" = gamma_vec_long[index], "Power" = par_results[index]))
  
}

plot_1IA_norm = function(){
  # Replicates Figure 1, but with the normal endpoint
  
  Power_plot = numeric(39)
  ESS_plot = numeric(39)
  
  IA_vec = (1:39)*2
  
  for(i in 1:39){
    
    print(i)
    
    temp_para = Norm_optim_par(IAn = c(IA_vec[i],80))
    
    temp_result = OC_gen_norm(IAn =  c(IA_vec[i],80), lambda = temp_para[1], gamma = temp_para[2])
    
    Power_plot[i] = temp_result[1]
    ESS_plot[i] = temp_result[2]
    
  }
  
  plot_df = data.frame(ESS_plot, Power_plot, IA_placement = IA_vec)
  
  
  
  
  ggplot_ess = ggplot(plot_df, aes(y = ESS_plot, x = IA_placement)) +
    #geom_ribbon(aes(ymin = low_ESS, ymax = high_ESS), fill = "grey70") +
    geom_line(aes(y = ESS_plot)) +
    labs( title = "ESS versus Interim placement - Normal endpoint.", x = "IA placement", y = "ESS")
  
  
  print(ggplot_ess)
  
  ggplot_power = ggplot(plot_df, aes(y = Power_plot, x = IA_placement)) +
    #geom_ribbon(aes(ymin = Power_plot -  1.96*0.5/sqrt(n_sim), ymax = Power_plot + 1.96*0.5/sqrt(n_sim)), fill = "grey70") +
    geom_line(aes(y = Power_plot)) +
    labs( title = "Power versus Interim placement - Normal endpoint.", x = "IA placement", y = "Power")
  print(ggplot_power)
  
  ggarrange(ggplot_ess, ggplot_power,
            ncol = 1, nrow = 2)
  
  browser()
  
}

plot_1IA_norm_prop = function(){
  # Replicates Figure 1, but with the normal endpoint
  
  Prop_plot = numeric(39)
  LQ_plot = numeric(39)
  UQ_plot = numeric(39)

  IA_vec = (1:39)*2
  
  for(i in 1:39){
    
    print(i)
    
    temp_para = Norm_optim_par(IAn = c(IA_vec[i],80), RAR = TRUE)
    
    temp_result = OC_gen_norm(IAn =  c(IA_vec[i],80), lambda = temp_para[1], gamma = temp_para[2], RAR = TRUE)
    
    Prop_plot[i] = temp_result[3]
    LQ_plot[i] = temp_result[4]
    UQ_plot[i] = temp_result[5]
    
  }
  
  plot_df = data.frame(Prop_plot, LQ_plot, UQ_plot, IA_vec = IA_vec)
  
  
  ggplot_ess = ggplot(plot_df, aes(x=IA_vec, y=Prop_plot)) +
    geom_ribbon(aes(ymin=LQ_plot, ymax=UQ_plot, alpha = 0.5), fill = "grey70") +
    geom_line(aes(y = Prop_plot)) +
    labs( title = "Proportion to best treatment against IA placement. - Normal", x = "IA placement", y = "Proportion") +
    coord_cartesian(ylim = c(0.45,0.75)) +
    scale_y_continuous(breaks = seq(0.45, 0.75, by = 0.05)) +
    theme(legend.position="none")
  
  print(ggplot_ess)
  
  
  #browser()
  
}


case_study_1IA_plot = function(){
  #Figure 1, but for the case study
  
  Power_plot = numeric(49)
  ESS_plot = numeric(49)
  
  for(i in 1:49){
    
    print(i)
    
    temp_para = optim_case_par(IAn = c(i,50))
    
    temp_result = exact_OC_case(IAn = c(i,50), lambda = temp_para[1], gamma = temp_para[2], par = TRUE)
    
    Power_plot[i] = temp_result[1]
    ESS_plot[i] = temp_result[2]
    
  }
  
  plot_df = data.frame(ESS_plot, Power_plot, IA_placement = 1:49)
  
  
  
  
  ggplot_ess = ggplot(plot_df, aes(y = ESS_plot, x = IA_placement)) +
    #geom_ribbon(aes(ymin = low_ESS, ymax = high_ESS), fill = "grey70") +
    geom_line(aes(y = ESS_plot)) +
    labs( title = "ESS versus Interim placement - Illustrative example.", x = "IA placement", y = "ESS")
  
  
  print(ggplot_ess)
  
  ggplot_power = ggplot(plot_df, aes(y = Power_plot, x = IA_placement)) +
    #geom_ribbon(aes(ymin = Power_plot -  1.96*0.5/sqrt(n_sim), ymax = Power_plot + 1.96*0.5/sqrt(n_sim)), fill = "grey70") +
    geom_line(aes(y = Power_plot)) +
    labs( title = "Power versus Interim placement - Illustrative example.", x = "IA placement", y = "Power")
  print(ggplot_power)
  
  ggarrange(ggplot_ess, ggplot_power,
            ncol = 1, nrow = 2)
  
  #browser()
   
}

case_study_optim_IA = FALSE
if (case_study_optim_IA){
  # Gets the IA which minimses the ESS for H0 and H1 (under the same lambda and gamma that optimises the orignial)
  
  registerDoRNG(2025)
  
  print("H1")
  print(IA_optimiser_case_par(n_sim = 10000))
  
  print("H0")
  print(IA_optimiser_case_par(p_eff = 0.5, p_tox = 0.35, n_sim = 10000))
}

case_study_optim_IA_30 = FALSE
if (case_study_optim_IA_30){
  # Gets the IA which minimses the ESS for H0 and H1 (under the same lambda and gamma that optimises the original)
  # Same as above but using the 30% parameters for lambda gamma
  # This is used for Table 2
  
  registerDoRNG(2025)
  
  print("H1")
  print(IA_optimiser_case_par(n_sim = 10000, lambda = 0.675, gamma = 0.8))
  
  print("H0")
  print(IA_optimiser_case_par(p_eff = 0.5, p_tox = 0.35, n_sim = 10000, lambda = 0.675, gamma = 0.8))
}

table_2_gen = FALSE
if (table_2_gen){
  
  ## One IA 
  # Equal spaced
  
  temp_para = optim_case_par(IAn = c(25,50))
  print("Equal - 1 IA")
  print(exact_OC_case(IAn = c(25,50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  registerDoRNG(2025)
  
  optim_vec = IA_optimiser_case_par_2IA(lambda = temp_para[1], gamma = temp_para[2])
  print("Minimal - 1 IA - equal para")
  print(optim_vec)
  print(exact_OC_case(IAn = c(optim_vec[1],50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  # Waiting period (wait 10)
  
  temp_para = optim_case_par(IAn = c(30,50))
  print("Waiting - 1 IA")
  print(exact_OC_case(IAn = c(30,50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  # Minimise ESS 
  
  registerDoRNG(2025)
  
  optim_vec = IA_optimiser_case_par_2IA(lambda = temp_para[1], gamma = temp_para[2])
  print("Minimal - 1 IA - waiting para")
  print(optim_vec)
  print(exact_OC_case(IAn = c(optim_vec[1],50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  
  
  
  ## Two IA
  # Equal spaced
  
  temp_para = optim_case_par(IAn = c(16,33,50))
  print("Equal - 2 IA")
  print(exact_OC_case(IAn = c(16,33,50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  registerDoRNG(2025)
  
  optim_vec = IA_optimiser_case_par_3IA(lambda = temp_para[1], gamma = temp_para[2])
  print("Minimal - 2 IA - equal")
  print(optim_vec)
  print(exact_OC_case(IAn = c(optim_vec[1],optim_vec[2],50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  # Waiting period (wait 10)
  
  temp_para = optim_case_par(IAn = c(23,36,50))
  print("Waiting - 2 IA")
  print(exact_OC_case(IAn = c(27,39,50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  # Minimise ESS 
  
  registerDoRNG(2025)
  
  optim_vec = IA_optimiser_case_par_3IA(lambda = temp_para[1], gamma = temp_para[2])
  print("Minimal - 2 IA - waiting")
  print(optim_vec)
  print(exact_OC_case(IAn = c(optim_vec[1],optim_vec[2],50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  ## Three IA
  # FRAIL-M spacing
  
  #temp_para = optim_case_par(IAn = c(10,20,35,50))
  #print("Equal - 3 IA")
  #print(exact_OC_case(IAn = c(10,20,35,50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  
  # Waiting period (wait 15)
  
  #temp_para = optim_case_par(IAn = c(15,26,38,50))
  #print("Equal - 3 IA")
  #print(exact_OC_case(IAn = c(15,26,38,50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  # Minimise ESS 
  
  #registerDoRNG(2025)
  
  #print("Minimal ESS - 3 IA")
  #optim_vec = IA_optimiser_case_par(n_sim = 10000, lambda = temp_para[1], gamma = temp_para[2]))
  #print("Minimal - 3 IA")
  #print(optim_vec)
  #print(exact_OC_case(IAn = c(optim_vec[1],optim_vec[2],optim_vec[3],50), par = TRUE, lambda = temp_para[1], gamma = temp_para[2]))
  
  
}

case_study = FALSE
if (case_study){
  set.seed(2025)
  T1 = Sys.time()
  print(exact_OC_case())
  T1 = Sys.time() - T1
  
  T2 = Sys.time()
  print(exact_OC_case(par = TRUE))
  T2 = Sys.time() - T2
  #print(sim_OC_case())
  
  #print(exact_OC_case(0.5,0.35))
  #print(sim_OC_case(0.5,0.35))
  
}

appendix_norm = FALSE
if (appendix_norm){
  set.seed(2025)
  registerDoRNG(2025)
  
  plot_1IA_norm_prop()
  plot_1IA_norm()
  
  
}