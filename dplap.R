# Fits a Laplace approximation for a normal DPM mixture model and draws samples
# from either Laplace or skew-Laplace approximation, returning computation times.
this_dir = dirname(normalizePath(sys.frame(1)$ofile))

#source utilities
source(file.path(this_dir, "utils_dplap.R"))

#main function
dp_lap_mixmodel_normal_normal <- function(y, K = 40, sigma = 1,
                                          m0 = 0, s0 = 1, a = 3, b = 3,
                                          n_samp = 5000,
                                          skew = FALSE,
                                          init_par = rep(0, 2 * K), seed = NULL) {
  nsim = length(y)
  
  if(!is.null(seed)){set.seed(seed)}
  
  if (length(init_par) != 2 * K) {
    stop("Length of 'init_par' must be equal to 2 * K.")
  }
  if (sigma <= 0) stop("'sigma' must be > 0.")
  if (s0 <= 0) stop("'s0' must be > 0.")
  if (a <= 0 || b <= 0) stop("'a' and 'b' must be > 0.")
  if (n_samp <= 0) stop("'n_samp' must be > 0.")
  
  log_post <- function(param) {
    log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  }
  
  nlog_post <- function(param) {
    -log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  }
  
  nlog_grad <- function(param) {
    -grad_log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  }
  
  ## --- optimization ---
  start_total <- Sys.time()
  start_opt <- Sys.time()
  
  optimization <- optim(
    par = init_par,
    fn = nlog_post,
    gr = nlog_grad,
    method = "BFGS",
    hessian = TRUE
  )
  
  end_opt <- Sys.time()
  
  m_lap = optimization$par
  H = optimization$hessian
  p = ncol(H)
  
  cov_lap = tryCatch(
    solve(H),
    error = function(e) solve(H + diag(1e-6, p))
  )
  
  eig = eigen(cov_lap, symmetric = TRUE)
  if (any(eig$values <= 0)) {
    eig$values = pmax(eig$values, 1e-8)
    cov_lap    = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  }
  
  cov_lap_sqrt = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  
  ## --- sampling ---
  start_samp <- Sys.time()
  
  if (skew) {
    samples <- sim_eqc_skew_sym(
      nsim = n_samp,
      mu = m_lap,
      Sigma_sqrt = cov_lap_sqrt,
      log_post = log_post
    )
    approximation <- "skew-Laplace"
  } else {
    samples <- mvtnorm::rmvnorm(
      n = n_samp,
      mean = m_lap,
      sigma = cov_lap
    )
    approximation <- "Laplace"
  }
  
  end_samp <- Sys.time()
  end_total <- Sys.time()
  
  list(
    approximation = approximation,
    optimization = optimization,
    mean = m_lap,
    covariance = cov_lap,
    covariance_sqrt = cov_lap_sqrt,
    samples = samples,
    time = list(
      total = end_total - start_total,
      optimization = end_opt - start_opt,
      sampling = end_samp - start_samp
    )
  )
}