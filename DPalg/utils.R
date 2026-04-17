################################################################################
## this code contain helper functions for the MCMCs
#  for DP mixture of univariate Normals
#  Kernel:   y | mu ~ N(mu, sigma2)
#  Base:     mu ~ N(mu0, tau20)
#  DP:       G ~ DP(alpha, N(mu0, tau20)) with alpha random
################################################################################

# ---- helper: compute log lik at each iteration ----
# used by _main_to_run1.R and _main_to_run2.R
dp_loglik_trace_mixmodel_normal_normal <- function(Y, fit, hyper = c(1, 0, 1)) {
  sd_y   = sqrt(hyper[1])
  
  # Make sure Y is a plain numeric vector
  if (is.list(Y)) Y = unlist(Y)
  Y = as.numeric(Y)
  
  S = length(fit$c_samples)
  if (is.null(fit$phis[[S]])){return(NULL)}
  
  ll = numeric(S)
  
  for (t in seq_len(S)) {
    c_t    = fit$c_samples[[t]]
    phis_t = fit$phis[[t]]
    
    c_t = as.integer(c_t)
    mu  = phis_t[c_t]  # mean for each observation under current allocation
    
    ll[t] = sum(dnorm(Y, mean = mu, sd = sd_y, log = TRUE))
  }
  
  ll
}

# ---- helper: relabel cluster ids to 1..H and return mapping ----
# used by dp_slice_mixmodel_normal_normal, 
# dp_CRPwithAtoms_mixmodel_normal_normal, and 
# dp_CRPnoAtoms_mixmodel_normal_normal

relabel_to_consecutive <- function(c) {
  labs   = sort(unique(c))
  new_id = match(c, labs) 
  list(c = new_id, labs = labs)
}


estimate_dpmm_density_ss <- function(fit, grid, sigma2, burn = 0, thin = 1,
                                     probs = c(0.025, 0.5, 0.975)) {
  # fit    : output of dp_slice_mixmodel_normal_normal()
  # grid   : numeric vector where the density is evaluated
  # sigma2 : kernel variance used in the model
  # burn   : number of initial MCMC iterations to discard
  # thin   : thinning interval
  
  stopifnot(is.list(fit))
  stopifnot(!is.null(fit$c_samples), !is.null(fit$phis))
  stopifnot(is.numeric(grid), length(grid) >= 1)
  stopifnot(is.numeric(sigma2), length(sigma2) == 1, sigma2 > 0)
  stopifnot(burn >= 0, thin >= 1)
  
  Ttot = length(fit$c_samples)
  keep = seq.int(from = burn + 1, to = Ttot, by = thin)
  if (length(keep) == 0) stop("No MCMC iterations left after burn/thin.")
  
  n = length(fit$c_samples[[keep[1]]])
  #dens = numeric(length(grid))
  tt_count = 0
  dens = matrix(NA, nrow = length(keep), ncol = length(grid))
  
  for (tt in keep) {
    tt_count = tt_count + 1
    H_t    = fit$H[[tt]]
    phis_t = fit$phis[[tt]]
    pi_t   = fit$pi[[tt]]
    
    # occupied components only
    mu_h = c(phis_t[seq_len(H_t)], rnorm(1, 0, 1))
    temp = pi_t[seq_len(H_t)]
    w_h = c(temp, 1-sum(temp)) 
    
    # mixture density at this iteration
    dens_t = numeric(length(grid))
    for (h in seq_len(H_t+1)) {
      dens_t = dens_t + w_h[h] * dnorm(grid, mean = mu_h[h], sd = sqrt(sigma2))
    }
    
    #dens = dens + dens_t
    dens[tt_count,] = dens_t
  }
  
  mean_density   = colMeans(dens)
  median_density = apply(dens, 2, median)
  q_density      = apply(dens, 2, quantile, probs = probs)
  
  
  return(list(
    xgrid = xgrid,
    density_draws  = dens,
    mean_density   = mean_density,
    median_density = median_density,
    quantiles      = q_density,
    probs = probs
  ))
}


estimate_dpmm_density_bgs <- function(fit, grid, sigma2, L, burn = 0, thin = 1, 
                                      probs = c(0.025, 0.5, 0.975)) {
  # fit    : output of dp_slice_mixmodel_normal_normal()
  # grid   : numeric vector where the density is evaluated
  # sigma2 : kernel variance used in the model
  # burn   : number of initial MCMC iterations to discard
  # thin   : thinning interval
  
  stopifnot(is.list(fit))
  stopifnot(!is.null(fit$c_samples), !is.null(fit$phis))
  stopifnot(is.numeric(grid), length(grid) >= 1)
  stopifnot(is.numeric(sigma2), length(sigma2) == 1, sigma2 > 0)
  stopifnot(burn >= 0, thin >= 1)
  
  Ttot = length(fit$c_samples)
  keep = seq.int(from = burn + 1, to = Ttot, by = thin)
  if (length(keep) == 0) stop("No MCMC iterations left after burn/thin.")
  
  n = length(fit$c_samples[[keep[1]]])
  #dens = numeric(length(grid))
  dens = matrix(NA, nrow = length(keep), ncol = length(grid))
  
  tt_count = 0
  for (tt in keep) {
    tt_count = tt_count + 1
    unique_c_labels = unique(fit$c_samples[[tt]])
    phis_t   = c( fit$phis[[tt]][unique_c_labels], rnorm(1, 0, 1))
    temp     = fit$pi[[tt]][unique_c_labels]
    w_h      = c(temp, 1-sum(temp)) 
    
    # mixture density at this iteration
    dens_t = numeric(length(grid))
    for (h in seq_len(length(phis_t))) {
      dens_t = dens_t + w_h[h] * dnorm(grid, mean = phis_t[h], sd = sqrt(sigma2))
    }
    
    #dens = dens + dens_t
    dens[tt_count,] = dens_t
  }
  
  mean_density   = colMeans(dens)
  median_density = apply(dens, 2, median)
  q_density      = apply(dens, 2, quantile, probs = probs)
  
  
  return(list(
    xgrid = xgrid,
    density_draws  = dens,
    mean_density   = mean_density,
    median_density = median_density,
    quantiles      = q_density,
    probs = probs
  ))
}
