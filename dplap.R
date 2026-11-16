# Fits a Laplace approximation for a normal DPM mixture model and draws samples
# from either Laplace or skew-Laplace approximation, returning computation times.

#main function

dp_lap_mixmodel_normal_normal <- function(y, K = 40, sigma = 1,
                                          m0 = 0, s0 = 1, a = 3, b = 3,
                                          n_samp = 5000,
                                          skew = FALSE,
                                          init_par = rep(0, 2 * K), seed = NULL,
                                          tol_merge = 1e-2, max_passes = 50,
                                          multistart_merge = TRUE, n_starts_merge = 4) {
  nsim = length(y)
  
  if(!is.null(seed)){set.seed(seed)}
  
  if (length(init_par) != 2 * K) {
    stop("Length of 'init_par' must be equal to 2 * K.")
  }
  if (sigma <= 0) stop("'sigma' must be > 0.")
  if (s0 <= 0) stop("'s0' must be > 0.")
  if (a <= 0 || b <= 0) stop("'a' and 'b' must be > 0.")
  if (n_samp <= 0) stop("'n_samp' must be > 0.")
  if (max_passes < 1) stop("'max_passes' must be >= 1.")
  if (n_starts_merge < 1) stop("'n_starts_merge' must be >= 1.")
  
  log_post <- function(param) {
    log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  }
  
  nlog_post <- function(param) {
    -log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  }
  
  nlog_grad <- function(param) {
    -grad_log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  }
  
  #fusion
  merge_and_refit <- function(m_lap, y, K, sigma = 1, m0 = 0, s0 = 1, a = 3, b = 3,
                              tol = 1e-2, maxit = 5000, reltol = 1e-12,
                              multistart_merge = FALSE, n_starts_merge = 5) {
    pi_hat = stick_breaking_weights_from_R(m_lap[1:(K - 1)])
    mu_hat = m_lap[K:(2 * K - 1)]
    tau    = m_lap[2 * K]

    group  = rep(NA_integer_, K)
    centre = numeric(0)
    for (h in order(pi_hat, decreasing = TRUE)) {
      d = abs(centre - mu_hat[h])
      if (length(d) > 0 && min(d) < tol) {
        group[h] = which.min(d)
      } else {
        centre   = c(centre, mu_hat[h])
        group[h] = length(centre)
      }
    }
    K_new = length(centre)
    if (K_new < 2) {
      warning("merge_and_refit: single group.")
      return(list(par = m_lap, K = K, groups = list(seq_len(K)), optim = NULL))
    }
    
    f_grp = factor(group, levels = seq_len(K_new))
    w     = as.vector(tapply(pi_hat, f_grp, sum))
    mu    = as.vector(tapply(pi_hat * mu_hat, f_grp, sum)) / w

    o   = order(w, decreasing = TRUE)
    w   = w[o]; mu = mu[o]
    rem = 1 - c(0, cumsum(w))[1:(K_new - 1)]
    V   = pmin(pmax(w[1:(K_new - 1)] / rem, 1e-10), 1 - 1e-10)
    start = c(qlogis(V), mu, tau)
    
    nlp = function(p) -log_posterior_dpm_random_alpha(p, y, K_new, sigma, m0, s0, a, b)
    ngr = function(p) -grad_log_posterior_dpm_random_alpha(p, y, K_new, sigma, m0, s0, a, b)
    
    if (multistart_merge) {
      if (n_starts_merge > 1) {
        pm2 = rbind(start,
                    matrix(runif((n_starts_merge - 1) * length(start), -3, 3),
                           nrow = n_starts_merge - 1))
      } else {
        pm2 = matrix(start, nrow = 1)
      }
      invisible(capture.output(
        ms2 <- multistart(
          parmat  = pm2,
          fn      = nlp,
          gr      = ngr,
          method  = "BFGS",
          hessian = FALSE,
          control = list(maxit = maxit, reltol = reltol, trace = 0)
        )
      ))
      bst = which.min(ms2$value)
      opt = list(par         = as.numeric(unlist(ms2[bst, seq_along(start)])),
                 value       = ms2$value[bst],
                 convergence = ms2$convergence[bst],
                 best_start  = bst,          # 1 = il punto fuso
                 all_values  = ms2$value)
    } else {
      opt = optim(start, nlp, ngr, method = "BFGS",
                  control = list(maxit = maxit, reltol = reltol))
    }
    
    list(par    = opt$par,
         K      = K_new,
         groups = split(seq_len(K), factor(group, levels = o)),  # componenti fuse (indici della passata precedente)
         optim  = opt)
  }
  
  
  
  ## --- optimization ---
  start_total <- Sys.time()
  start_opt <- Sys.time()
  pm <- matrix(runif(4 * length(init_par), -3, 3), nrow = 4)
  invisible(capture.output(
    ms <- multistart(
      parmat  = pm,
      fn      = nlog_post,
      gr      = nlog_grad,
      method  = "BFGS",
      hessian = FALSE,
      control = list(maxit = 2000, reltol = 1e-12,trace = 0),
      
    )
  ))
  npar <- length(init_par)
  best <- which.min(ms$value)
  # par_hat <- unlist(ms[best, 1:npar])
  # H_hat   <- optimHess(par_hat, fn = nlog_post, gr = nlog_grad)
  #
  # optimization <- optim(
  #   par = rep(0, length(init_par)),
  #   fn = nlog_post,
  #   gr = nlog_grad,
  #   method = "BFGS",
  #   hessian = TRUE,
  #   control = list(
  #     maxit = 2000,
  #     reltol = 1e-12
  #   )
  # )
  
  m_lap =  unlist(ms[best, 1:npar])

  K_cur   = K
  par_cur = as.numeric(m_lap)
  passes  = list()
  for (it in seq_len(max_passes)) {
    ref = merge_and_refit(par_cur, y, K_cur, sigma, m0, s0, a, b,
                          tol = tol_merge,
                          multistart_merge = multistart_merge,
                          n_starts_merge = n_starts_merge)
    passes[[it]] = ref
    par_cur = as.numeric(ref$par)
    if (ref$K == K_cur) break          
    K_cur = ref$K
  }
  K     = K_cur
  m_lap = par_cur
  end_opt <- Sys.time()   
  
  log_post  = function(param)  log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  nlog_post = function(param) -log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  nlog_grad = function(param) -grad_log_posterior_dpm_random_alpha(param, y, K, sigma, m0, s0, a, b)
  
  H = optimHess(m_lap, fn = nlog_post, gr = nlog_grad)
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
    #optimization = optimization,
    mean = m_lap,
    covariance = cov_lap,
    covariance_sqrt = cov_lap_sqrt,
    samples = samples,
    time = list(
      total = end_total - start_total,
      optimization = end_opt - start_opt,
      sampling = end_samp - start_samp
    ),
    K = K,
    merge = passes,            
    n_passes = length(passes)
  )
}