## ---------------------------------------------------------------------------
## Exact (non-Monte-Carlo) two-arm BOP2 operating characteristics via a
## forward dynamic program over the joint distribution of cumulative
## per-arm successes.
##----------------------------------------------------------------------------

# Code written by Ayon Mukherjee

## Futility / superiority probability thresholds for the BOP2 boundary.
threshold_prob <- function(lam, gam, n_total, n_current) {
  futility <- lam * (n_current / n_total)^gam
  z <- qnorm((1 + lam) / 2)
  superiority <- 2 * pnorm(z / sqrt(n_current / n_total)) - 1
  c(futility = futility, superiority = superiority)
}

## Exact P(Beta(a1,b1) > Beta(a2,b2)) for integer a2,b2
## a1, b1 may be vectors; a2, b2 must be scalar 
## (comparing a vector of possible treatment-arm posteriors against one realised control-arm posterior).
beta_comp_exact <- function(a1, b1, a2, b2) 
{
  a2 <- as.integer(a2)
  b2 <- as.integer(b2)
  i <- 0:(a2 - 1)
  vapply(seq_along(a1), function(k) {
    log_terms <- lbeta(a1[k] + i, b1[k] + b2) - log(b2 + i) -
      lbeta(i + 1, b2) - lbeta(a1[k], b1[k])
    1 - sum(exp(log_terms))
  }, numeric(1))
}

## Exact two-arm operating characteristics.
##   p0, p1   : true control / treatment response rates
##   lam, gam : BOP2 boundary parameters
##   n_total  : maximum (final-analysis) total sample size
##   IAn      : increasing vector of CUMULATIVE TOTAL patients at each
##              analysis (last element must equal n_total); every entry
##              must be even (equal per-arm allocation)
## Returns reject_prob (power/type-I error), ess (expected total sample size), and per-stage futility/superiority stopping mass.

exact_two_arm <- function(p0, p1, lam, gam, n_total, IAn,
                          prior_a = 1, prior_b = 1) 
{
  num_stages <- length(IAn)
  if (any(IAn %% 2 != 0)) stop("exact_two_arm requires even cumulative IA sizes")
  
  ## active[[paste(i, j)]] = probability mass of trials still running with
  ## i cumulative control successes and j cumulative treatment successes.
  active <- new.env(hash = TRUE)
  assign("0,0", 1.0, envir = active)
  
  prev_n_pa <- 0L
  fute_prob <- numeric(num_stages)
  supe_prob <- numeric(num_stages)
  ess <- 0.0
  
  for (stage in seq_len(num_stages)) {
    cum_n <- IAn[stage]
    n_pa <- cum_n %/% 2L
    inc_pa <- n_pa - prev_n_pa
    
    bounds <- threshold_prob(lam, gam, n_total, cum_n)
    fut_bound <- bounds["futility"]
    sup_bound <- bounds["superiority"]
    
    if (inc_pa > 0) {
      k <- 0:inc_pa
      pmf_ctrl <- dbinom(k, inc_pa, p0)
      pmf_treat <- dbinom(k, inc_pa, p1)
    } else {
      pmf_ctrl <- 1.0
      pmf_treat <- 1.0
    }
    
    new_active <- new.env(hash = TRUE)
    add_mass <- function(env, key, mass) {
      cur <- if (exists(key, envir = env, inherits = FALSE)) get(key, envir = env) else 0.0
      assign(key, cur + mass, envir = env)
    }
    
    states <- ls(active)
    for (key in states) {
      w <- get(key, envir = active)
      if (w <= 0) next
      ij <- as.integer(strsplit(key, ",")[[1]])
      i0 <- ij[1]; j0 <- ij[2]
      
      for (dctrl in 0:inc_pa) {
        pc <- pmf_ctrl[dctrl + 1]
        if (pc <= 0) next
        for (dtreat in 0:inc_pa) {
          pt <- pmf_treat[dtreat + 1]
          if (pt <= 0) next
          mass <- w * pc * pt
          if (mass <= 0) next
          
          i <- i0 + dctrl  # cumulative control successes
          j <- j0 + dtreat # cumulative treatment successes
          
          A1 <- prior_a + j
          B1 <- prior_b + (n_pa - j)
          A2 <- prior_a + i
          B2 <- prior_b + (n_pa - i)
          prob <- beta_comp_exact(A1, B1, A2, B2)
          
          if (prob > sup_bound) {
            supe_prob[stage] <- supe_prob[stage] + mass
            ess <- ess + mass * cum_n
          } else if (prob < fut_bound) {
            fute_prob[stage] <- fute_prob[stage] + mass
            ess <- ess + mass * cum_n
          } else {
            add_mass(new_active, paste(i, j, sep = ","), mass)
          }
        }
      }
    }
    
    active <- new_active
    prev_n_pa <- n_pa
  }
  
  list(
    reject_prob = sum(supe_prob),
    ess = ess,
    fute_prob = fute_prob,
    supe_prob = supe_prob
  )
}

## ---------------------------------------------------------------------------
## grid searches (lambda, gamma)
## candidates: represents the joint distribution as a matrix rather than a list of states, and builds the per-stage binomial increment as a
## convolution (banded) matrix. Numerically identical to exact_two_arm().
## ---------------------------------------------------------------------------

.prob_grid_cache <- new.env()

prob_grid <- function(n_pa, prior_a = 1, prior_b = 1) {
  key <- paste(n_pa, prior_a, prior_b, sep = "_")
  if (!exists(key, envir = .prob_grid_cache)) {
    idx <- 0:n_pa
    A <- prior_a + idx
    B <- prior_b + (n_pa - idx)
    k <- 0:n_pa
    out <- matrix(0.0, n_pa + 1, n_pa + 1)  # rows = A2/B2 (control), cols = A1/B1 (treatment)
    for (row in seq_along(idx)) {
      A2 <- A[row]; B2 <- B[row]
      valid_k <- k[k < A2]
      if (length(valid_k) == 0) {
        out[row, ] <- 1.0
        next
      }
      for (col in seq_along(idx)) {
        A1 <- A[col]; B1 <- B[col]
        log_terms <- lbeta(A1 + valid_k, B1 + B2) - log(B2 + valid_k) -
          lbeta(valid_k + 1, B2) - lbeta(A1, B1)
        out[row, col] <- 1 - sum(exp(log_terms))
      }
    }
    assign(key, out, envir = .prob_grid_cache)
  }
  get(key, envir = .prob_grid_cache)
}

## Banded convolution matrix mapping an old-size distribution through a
## length-L increment pmf into a new-size distribution.
conv_matrix <- function(old_size, pmf, new_size) 
{
  Cm <- matrix(0.0, new_size, old_size)
  L <- length(pmf)
  for (kcol in 0:(old_size - 1)) {
    lo <- kcol
    hi <- min(kcol + L - 1, new_size - 1)
    if (hi >= lo) {
      Cm[(lo:hi) + 1, kcol + 1] <- pmf[1:(hi - lo + 1)]
    }
  }
  Cm
}

run_two_arm_fast <- function(p0, p1, lam, gam, n_total, IAn,
                             prior_a = 1, prior_b = 1) {
  prev_n_pa <- 0L
  max_n_pa <- IAn[length(IAn)] %/% 2L
  active <- matrix(0.0, max_n_pa + 1, max_n_pa + 1)
  active[1, 1] <- 1.0
  ess <- 0.0
  reject_prob <- 0.0
  
  for (cum_n in IAn) {
    n_pa <- cum_n %/% 2L
    inc_pa <- n_pa - prev_n_pa
    bounds <- threshold_prob(lam, gam, n_total, cum_n)
    fut_bound <- bounds["futility"]; sup_bound <- bounds["superiority"]
    
    cur <- active[1:(prev_n_pa + 1), 1:(prev_n_pa + 1), drop = FALSE]
    if (inc_pa > 0) {
      k <- 0:inc_pa
      pmf_c <- dbinom(k, inc_pa, p0)
      pmf_t <- dbinom(k, inc_pa, p1)
      Crow <- conv_matrix(prev_n_pa + 1, pmf_c, n_pa + 1)
      Ccol <- conv_matrix(prev_n_pa + 1, pmf_t, n_pa + 1)
      new_active <- Crow %*% cur %*% t(Ccol)
    } else {
      new_active <- cur
    }
    
    pg <- prob_grid(n_pa, prior_a, prior_b)
    ## pg is indexed [control_state, treatment_state] 
    #new_active is indexed
    
    stop_sup <- pg > sup_bound
    stop_fut <- pg < fut_bound
    stop <- stop_sup | stop_fut
    
    reject_prob <- reject_prob + sum(new_active[stop_sup])
    ess <- ess + cum_n * sum(new_active[stop])
    
    cont <- new_active
    cont[stop] <- 0.0
    active <- matrix(0.0, max_n_pa + 1, max_n_pa + 1)
    active[1:(n_pa + 1), 1:(n_pa + 1)] <- cont
    prev_n_pa <- n_pa
  }
  
  list(reject_prob = reject_prob, ess = ess)
}

## ---------------------------------------------------------------------------
## Part A cross-checks the two implementations above (exact_two_arm, the
## direct hashed-state DP, versus run_two_arm_fast) against each other across several configurations. 
## They are structured completely differently -- one walks
## a list of (control, treatment) success-count states one at a time, the
## other propagates a whole probability matrix via convolution.
## ---------------------------------------------------------------------------

cat("=== Part A: exact_two_arm() vs run_two_arm_fast() ===\n")
check_cases <- list(
  list(p0 = 0.2, p1 = 0.4, lam = 0.91, gam = 0.93, n_total = 80,  IAn = c(60, 80)),
  list(p0 = 0.2, p1 = 0.2, lam = 0.91, gam = 0.93, n_total = 80,  IAn = c(60, 80)),
  list(p0 = 0.2, p1 = 0.4, lam = 0.91, gam = 0.98, n_total = 80,  IAn = c(46, 64, 80)),
  list(p0 = 0.2, p1 = 0.4, lam = 0.91, gam = 0.98, n_total = 80,  IAn = c(20, 40, 60, 80)),
  list(p0 = 0.2, p1 = 0.2, lam = 0.91, gam = 0.98, n_total = 80,  IAn = c(20, 40, 60, 80)),
  list(p0 = 0.2, p1 = 0.4, lam = 0.83, gam = 0.71, n_total = 160, IAn = c(80, 160)),
  list(p0 = 0.1, p1 = 0.3, lam = 0.9,  gam = 0.86, n_total = 80,  IAn = c(20, 40, 60, 80)),
  list(p0 = 0.5, p1 = 0.5, lam = 0.85, gam = 0.9,  n_total = 100, IAn = c(50, 100))
)

max_diff_reject <- 0
max_diff_ess <- 0
for (cs in check_cases) {
  r  <- exact_two_arm(cs$p0, cs$p1, cs$lam, cs$gam, cs$n_total, cs$IAn)
  rf <- run_two_arm_fast(cs$p0, cs$p1, cs$lam, cs$gam, cs$n_total, cs$IAn)
  d_reject <- abs(unname(r$reject_prob) - unname(rf$reject_prob))
  d_ess    <- abs(unname(r$ess) - unname(rf$ess))
  max_diff_reject <- max(max_diff_reject, d_reject)
  max_diff_ess    <- max(max_diff_ess, d_ess)
  cat(sprintf("p0=%.2g p1=%.2g lam=%.2g gam=%.2g IA=%-14s  DP: reject=%.6f ess=%.4f  |  fast: reject=%.6f ess=%.4f  |  diff=%.2e/%.2e\n",
              cs$p0, cs$p1, cs$lam, cs$gam, paste(cs$IAn, collapse = ","),
              r$reject_prob, r$ess, rf$reject_prob, rf$ess, d_reject, d_ess))
}
cat(sprintf("\nMax difference between the two implementations across all cases: reject_prob=%.2e, ess=%.2e\n",
            max_diff_reject, max_diff_ess))

##------------------------------------------------------------------------------------------------
## Part B reproduces (Table 6 / Figure 8's two-IA minimal-ESS configuration, and its type-I-error)
##------------------------------------------------------------------------------------------------
cat("=== Part B: reproducing manuscript numbers directly ===\n")
r <- exact_two_arm(0.2, 0.4, 0.91, 0.98, 80, c(46, 64, 80))
cat(sprintf("Two-IA minimal-ESS config (Table 6 / Figure 8, Beta(1,1) prior): power = %.4f, ess = %.3f  [paper reports power=0.738, ess=60.87]\n",
            r$reject_prob, r$ess))

r0 <- exact_two_arm(0.2, 0.2, 0.91, 0.98, 80, c(46, 64, 80))
cat(sprintf("Same config, type-I error: %.4f  [paper reports 0.097]\n", r0$reject_prob))

r3 <- exact_two_arm(0.2, 0.4, 0.91, 0.98, 80, c(34, 46, 64, 80))
cat(sprintf("Three-IA minimal-ESS config: power = %.4f, ess = %.3f  [paper reports power=0.730, ess=57.61]\n",
            r3$reject_prob, r3$ess))