# This file provides utilities for Bayesian inference in Dirichlet Process
# Mixture (DPM) models using Laplace and skew-Laplace approximations of the
# posterior distribution. It includes functions for stick-breaking weight
# construction, skew-symmetric sampling (EQC-based proposals), evaluation of
# the log-posterior and its gradient, computation of mixture densities from
# posterior draws, and posterior density summarization and visualization.

library(matrixStats)#version 1.5.0, package used for the function logSumExp

# Converts a real-valued vector R into stick-breaking weights pi:
# input R = numeric vector of length K-1; 
# output pi = numeric vector of length K with nonnegative entries summing to 1.
stick_breaking_weights_from_R <- function(R) {
  
  V = plogis(R)
  K = length(R) + 1L
  
  logV   = log(V)
  log1mV = log1p(-V)
  prefix = c(0, cumsum(log1mV))
  
  log_pi = c(logV + prefix[-K], prefix[K])
  
  exp(log_pi)
  
}

# Constructs a 2×n matrix of symmetric points around xs:
# input x = numeric vector (point to flip), xs = numeric scalar or vector (same length as x);
# output 2×length(x) matrix with rows (x, 2*xs - x)
eqc_ss <- function(x, xs) {
  rbind(x, 2 * xs - x)
}

# Simulates nsim draws from a skew-symmetric distribution:
# input nsim = number of samples; mu = mean vector (length d);
# Sigma_sqrt = d×d matrix square root of covariance; log_post = function ℝ^d → ℝ;
# output x_sim = nsim × d matrix of samples.
sim_eqc_skew_sym <- function(nsim, mu, Sigma_sqrt, log_post) {
  d = length(mu)
  two_mu = 2 * mu
  
  x_sim = sweep(
    matrix(rnorm(nsim * d), nrow = nsim, ncol = d) %*% Sigma_sqrt,
    2, mu, `+`
  )
  
  x_ref = sweep(-x_sim, 2, two_mu, `+`)
  
  lp1 = vapply(
    seq_len(nsim),
    function(i) log_post(x_sim[i, ]),
    numeric(1)
  )
  
  lp2 = vapply(
    seq_len(nsim),
    function(i) log_post(x_ref[i, ]),
    numeric(1)
  )
  
  choose_ref = log(runif(nsim)) < -log1p(exp(lp1 - lp2))
  x_sim[choose_ref, ] = x_ref[choose_ref, , drop = FALSE]
  
  x_sim
}

# Computes the log-posterior of a Dirichlet process mixture model with
# stick-breaking weights and random concentration parameter:
# input par = numeric vector of length 2K (R, mu, tau);
#       y = data vector; K = number of components; sigma = likelihood sd;
#       m0, s0 = prior mean and sd for mu; a, b = prior parameters for alpha;
# output = scalar log-posterior value.
log_posterior_dpm_random_alpha <- function(par, y, K, sigma, m0, s0, a, b) {
  if (length(par) != 2 * K) {
    stop("Length of 'par' must be equal to 2 * K.")
  }
  
  R     = par[1:(K - 1)]
  mu    = par[K:(2 * K - 1)]
  tau   = par[2 * K]
  alpha = exp(tau)
  
  V      = plogis(R)
  logV   = log(V)
  log1mV = log1p(-V)
  prefix = c(0, cumsum(log1mV))
  log_pi = c(logV + prefix[1:(K - 1)], prefix[K])
  
  log_prior_mu  = sum(dnorm(mu, mean = m0, sd = s0, log = TRUE))
  log_prior_R   = sum(logV + alpha * log1mV)
  log_prior_tau = (K + a - 1) * tau - b * alpha
  
  # log N(y_i | mu_k, sigma)
  # outer(y, mu, "-") gives matrix with entries y_i - mu_k
  z = outer(y, mu, "-") / sigma
  log_dens_mat = -0.5 * z^2 - log(sigma) - 0.5 * log(2 * pi)
  
  # each row i: log sum_k exp(log_pi[k] + log_dens_mat[i,k])
  A = sweep(log_dens_mat, 2, log_pi, "+")
  
  row_max = apply(A, 1, max)
  log_like = sum(row_max + log(rowSums(exp(A - row_max))))
  
  log_like + log_prior_mu + log_prior_R + log_prior_tau
}


# Computes the gradient of the scalar log-posterior for a K-component DPM:
# input par = numeric vector of length 2K (R, mu, tau), y = data vector,
#       K, sigma, m0, s0, a, b = model and prior hyperparameters;
# output = numeric gradient vector of length 2K ordered as (d/dR, d/dmu, d/dtau).
grad_log_posterior_dpm_random_alpha <- function(par, y, K, sigma, m0, s0, a, b) {
  if (length(par) != 2 * K) {
    stop("Length of 'par' must be equal to 2 * K.")
  }
  if (sigma <= 0) stop("'sigma' must be > 0.")
  if (s0 <= 0) stop("'s0' must be > 0.")
  if (a <= 0 || b <= 0) stop("'a' and 'b' must be > 0.")
  
  R     = par[1:(K - 1)]
  mu    = par[K:(2 * K - 1)]
  tau   = par[2 * K]
  alpha = exp(tau)
  
  V  = plogis(R)
  pi = stick_breaking_weights_from_R(R)
  
  dens_mat = outer(y, mu, function(yy, mm) dnorm(yy, mean = mm, sd = sigma))
  mix_dens = as.vector(dens_mat %*% pi)
  
  if (any(!is.finite(mix_dens)) || any(mix_dens <= 0)) {
    return(rep(NaN, length(par)))
  }
  
  resp = sweep(dens_mat, 2, pi, FUN = "*") / mix_dens
  
  # gradient with respect to mu
  grad_mu = -(mu - m0) / (s0^2) + colSums(resp * outer(y, mu, "-")) / (sigma^2)
  
  # tail cumulative sums by row, from right to left
  tail_resp = resp
  if (K >= 2) {
    for (j in (K - 1):1) {
      tail_resp[, j] = rowSums(resp[, j:K, drop = FALSE])
    }
  }
  
  grad_R = numeric(K - 1)
  for (j in 1:(K - 1)) {
    grad_R[j] = 1 - (1 + alpha) * V[j] +
      sum(resp[, j] - V[j] * tail_resp[, j])
  }
  
  grad_tau = alpha * sum(log1p(-V)) + (K + a - 1) - b * alpha
  
  c(grad_R, grad_mu, grad_tau)
}

# Evaluates the DPM mixture density at points x for one parameter draw:
# input x = numeric vector of evaluation points; par = numeric parameter vector
#       of length 2K (or 2K-1 if alpha_random = FALSE); K = number of components;
#       sigma = component standard deviation; alpha_random = whether par includes tau;
# output = numeric vector of length(x) with mixture density values.
dpm_density_from_draw <- function(x, par, K, sigma, alpha_random = TRUE) {
  if (sigma <= 0) stop("'sigma' must be > 0.")
  
  expected_len = if (alpha_random) 2 * K else (2 * K - 1)
  if (length(par) != expected_len) {
    stop("Length of 'par' is not consistent with K.")
  }
  
  R  = par[1:(K - 1)]
  mu = par[K:(2 * K - 1)]
  
  pi = stick_breaking_weights_from_R(R)
  
  dens_mat = outer(x, mu, function(xx, mm) dnorm(xx, mean = mm, sd = sigma))
  as.vector(dens_mat %*% pi)
}

# Computes posterior summaries of a DPM density over a grid:
# input xgrid = numeric vector of evaluation points; draws = M x p matrix of
#       posterior parameter draws; K = number of components; sigma = sd;
#       alpha_random = whether draws include tau; probs = quantile levels;
# output = list with density draws (M x length(xgrid)), mean, median, and
#       pointwise quantiles of the density.
posterior_density_dpm <- function(xgrid, draws, K, sigma, alpha_random = TRUE,
                                  probs = c(0.025, 0.5, 0.975)) {
  if (sigma <= 0) stop("'sigma' must be > 0.")
  
  draws = as.matrix(draws)
  M = nrow(draws)
  if (M == 0) stop("'draws' must not be empty.")
  
  expected_ncol = if (alpha_random) 2 * K else (2 * K - 1)
  if (ncol(draws) != expected_ncol) {
    stop("Number of columns of 'draws' is not consistent with K.")
  }
  
  n_grid = length(xgrid)
  
  dens_draws = t(vapply(
    seq_len(M),
    function(i) {
      dpm_density_from_draw(
        x = xgrid,
        par = draws[i, ],
        K = K,
        sigma = sigma,
        alpha_random = alpha_random
      )
    },
    FUN.VALUE = numeric(n_grid)
  ))
  
  mean_density   = colMeans(dens_draws)
  median_density = apply(dens_draws, 2, median)
  q_density      = apply(dens_draws, 2, quantile, probs = probs)
  
  list(
    xgrid = xgrid,
    density_draws  = dens_draws,
    mean_density   = mean_density,
    median_density = median_density,
    quantiles      = q_density,
    probs = probs
  )
}

# Plots posterior summaries of a DPM density over a grid:
# input xgrid = evaluation points; draws = posterior draws (M x p matrix);
#       K, sigma, alpha_random = model settings; probs = quantile levels;
#       n_show = number of sampled density curves to display;
#       show_* = logical flags for plotting components; labels and title;
# output = (invisible) list returned by posterior_density_dpm().
plot_posterior_density_dpm <- function(xgrid, draws, K, sigma, alpha_random = TRUE,
                                       probs = c(0.025, 0.5, 0.975),
                                       n_show = 30, show_draws = TRUE,
                                       show_mean = TRUE, show_median = FALSE,
                                       show_band = TRUE,
                                       xlab = "y", ylab = "density",
                                       main = NULL) {
  post_dens = posterior_density_dpm(
    xgrid = xgrid,
    draws = draws,
    K = K,
    sigma = sigma,
    alpha_random = alpha_random,
    probs = probs
  )
  
  dens_draws = post_dens$density_draws
  M = nrow(dens_draws)
  
  lower = post_dens$quantiles[1, ]
  upper = post_dens$quantiles[length(probs), ]
  
  ymax = max(dens_draws, post_dens$mean_density)
  ymin = 0
  
  plot(xgrid, post_dens$mean_density, type = "n",
       ylim = c(ymin, ymax),
       xlab = xlab, ylab = ylab, main = main)
  
  if (show_draws) {
    n_show = min(n_show, M)
    idx    = sample(seq_len(M), n_show)
    for (m in idx) {
      lines(xgrid, dens_draws[m, ], lty = 3)
    }
  }
  
  if (show_band) {
    lines(xgrid, lower, lty = 2)
    lines(xgrid, upper, lty = 2)
  }
  
  if (show_median) {
    med_row = which.min(abs(probs - 0.5))
    lines(xgrid, post_dens$quantiles[med_row, ], lty = 3, lwd = 2)
  }
  
  if (show_mean) {
    lines(xgrid, post_dens$mean_density, lwd = 3)
  }
  
  invisible(post_dens)
}


tv_from_draws_1d <- function(x, y, nbreaks = 40) {
  rng = range(c(x, y), finite = TRUE)
  
  # avoid degenerate case
  if (rng[1] == rng[2]) {
    return(0)
  }
  
  breaks = seq(rng[1], rng[2], length.out = nbreaks + 1)
  
  hx = hist(x, breaks = breaks, plot = FALSE) 
  hy = hist(y, breaks = breaks, plot = FALSE)
  
  px = hx$density
  py = hy$density
  0.5 * sum(abs(px - py)) * diff(breaks)[1]
}


compute_pointwise_tv_vs_ss <- function(method_obj, ss_obj, nbreaks = 20) {
  X = method_obj$density$density_draws   # rows = posterior draws, cols = grid points
  Y = ss_obj$density$density_draws
  
  ng = ncol(X)
  
  vapply(seq_len(ng), function(j) {
    tv_from_draws_1d(X[, j], Y[, j], nbreaks = nbreaks)
  }, numeric(1))
}

make_tv_boxplot_df <- function(LAP, Skew_LAP, sample_sizes, xgrid) {
  do.call(rbind, lapply(seq_along(sample_sizes), function(i) {
    rbind(
      data.frame(
        i_count = i,
        n = sample_sizes[i],
        method = "Lap",
        grid_point = seq_along(xgrid),
        x = xgrid,
        TV = LAP[[i]]$TV_grid_vs_SS
      ),
      data.frame(
        i_count = i,
        n = sample_sizes[i],
        method = "Skew-Lap",
        grid_point = seq_along(xgrid),
        x = xgrid,
        TV = Skew_LAP[[i]]$TV_grid_vs_SS
      )
    )
  }))
}

