# Truncated mean-field variational Bayes for the univariate Gaussian DPM
#
# Model (matching dplap.R and DPalg/Slice.R):
#   y_i | z_i = k, mu_k ~ N(mu_k, sigma^2)
#   mu_k ~ N(m0, s0^2)
#   v_k | alpha ~ Beta(1, alpha), k = 1,...,K-1; v_K = 1
#   alpha ~ Gamma(a, b), with b a rate
#
# The coordinate-ascent updates are those of Blei and Jordan (2006), with a
# Gamma variational factor for alpha.  The returned `samples` matrix has the
# same column layout used by posterior_density_dpm():
#   (logit(v_1),...,logit(v_{K-1}), mu_1,...,mu_K, log(alpha)).

.vb_logsumexp_rows <- function(x) {
  xmax <- apply(x, 1L, max)
  xmax + log(rowSums(exp(x - xmax)))
}

.vb_log_stick_weights_from_v <- function(v) {
  # Log weights for v_1,...,v_{K-1}, with v_K fixed at one.
  eps <- .Machine$double.eps
  v <- pmin(pmax(v, eps), 1 - eps)
  log_v <- log(v)
  log_1mv <- log1p(-v)
  prefix <- c(0, cumsum(log_1mv))
  c(log_v + prefix[seq_along(v)], prefix[length(prefix)])
}

.vb_conditional_logpost_mu <- function(mu, y, log_pi, sigma, m0, s0) {
  # The normalising constant of p(mu | v, y) is not required.
  log_kernel <- -0.5 * (outer(y, mu, "-") / sigma)^2 -
    log(sigma) - 0.5 * log(2 * pi)
  log_likelihood <- sum(.vb_logsumexp_rows(
    sweep(log_kernel, 2L, log_pi, "+")
  ))
  log_prior <- sum(dnorm(mu, mean = m0, sd = s0, log = TRUE))
  log_likelihood + log_prior
}

.vb_skew_mu_draws <- function(mu_draws, v_draws, mu_center,
                              y, sigma, m0, s0) {
  # Conditional skew-symmetric correction.  For each sampled v, q(mu) is
  # symmetric around mu_center.  A proposal mu is retained with probability
  #   p(mu | v,y) / {p(mu | v,y) + p(2*mu_center-mu | v,y)},
  # and is reflected otherwise.
  n_samp <- nrow(mu_draws)
  reflected <- logical(n_samp)

  for (s in seq_len(n_samp)) {
    mu_proposal <- mu_draws[s, ]
    mu_reflection <- 2 * mu_center - mu_proposal
    log_pi <- .vb_log_stick_weights_from_v(v_draws[s, ])

    lp_proposal <- .vb_conditional_logpost_mu(
      mu_proposal, y, log_pi, sigma, m0, s0
    )
    lp_reflection <- .vb_conditional_logpost_mu(
      mu_reflection, y, log_pi, sigma, m0, s0
    )

    if (!all(is.finite(c(lp_proposal, lp_reflection)))) {
      stop("Non-finite conditional log-posterior in the mu skew correction.")
    }

    # Stable evaluation of log{p(reflection)/(p(proposal)+p(reflection))}.
    lp_max <- max(lp_proposal, lp_reflection)
    log_denominator <- lp_max + log(
      exp(lp_proposal - lp_max) + exp(lp_reflection - lp_max)
    )
    log_prob_reflection <- lp_reflection - log_denominator

    reflected[s] <- log(runif(1L)) < log_prob_reflection
    if (reflected[s]) mu_draws[s, ] <- mu_reflection
  }

  list(draws = mu_draws, reflected = reflected)
}

.vb_beta_entropy <- function(shape1, shape2) {
  lbeta(shape1, shape2) - (shape1 - 1) * digamma(shape1) -
    (shape2 - 1) * digamma(shape2) +
    (shape1 + shape2 - 2) * digamma(shape1 + shape2)
}

.vb_elbo_normal_dpm <- function(y, sigma, m0, s0, a, b, r,
                                mu_mean, mu_var, v_shape1, v_shape2,
                                alpha_shape, alpha_rate) {
  n <- length(y)
  K <- length(mu_mean)
  sigma2 <- sigma^2
  s02 <- s0^2

  elog_v <- digamma(v_shape1) - digamma(v_shape1 + v_shape2)
  elog_1mv <- digamma(v_shape2) - digamma(v_shape1 + v_shape2)
  elog_pi <- c(elog_v, 0) + c(0, cumsum(elog_1mv))
  elog_alpha <- digamma(alpha_shape) - log(alpha_rate)
  ealpha <- alpha_shape / alpha_rate

  elog_lik <- -0.5 * log(2 * pi * sigma2) -
    0.5 * (outer(y, mu_mean, "-")^2 +
             matrix(mu_var, nrow = n, ncol = K, byrow = TRUE)) / sigma2

  allocation_and_likelihood <- sum(r *
    (sweep(elog_lik, 2L, elog_pi, "+") - log(pmax(r, .Machine$double.xmin))))

  mu_prior <- sum(-0.5 * log(2 * pi * s02) -
                    0.5 * ((mu_mean - m0)^2 + mu_var) / s02)
  mu_entropy <- 0.5 * sum(1 + log(2 * pi * mu_var))

  stick_prior <- (K - 1) * elog_alpha +
    sum((ealpha - 1) * elog_1mv)
  stick_entropy <- sum(.vb_beta_entropy(v_shape1, v_shape2))

  alpha_prior <- a * log(b) - lgamma(a) +
    (a - 1) * elog_alpha - b * ealpha
  alpha_entropy <- alpha_shape - log(alpha_rate) + lgamma(alpha_shape) +
    (1 - alpha_shape) * digamma(alpha_shape)

  allocation_and_likelihood + mu_prior + mu_entropy + stick_prior +
    stick_entropy + alpha_prior + alpha_entropy
}

# Fit the Blei-Jordan truncated variational approximation.
dp_vb_mixmodel_normal_normal <- function(
    y, K = 40, sigma = 1, m0 = 0, s0 = 1, a = 3, b = 3,
    n_samp = 2000, skew_mu = FALSE, max_iter = 10000, tol = 1e-8,
    init = c("kmeans", "random"), seed = NULL, verbose = FALSE) {

  init <- match.arg(init)
  y <- as.numeric(y)
  n <- length(y)

  if (!n || any(!is.finite(y))) stop("'y' must contain finite observations.")
  if (length(K) != 1L || K < 2 || K != as.integer(K))
    stop("'K' must be an integer >= 2.")
  K <- as.integer(K)
  if (sigma <= 0 || s0 <= 0) stop("'sigma' and 's0' must be > 0.")
  if (a <= 0 || b <= 0) stop("'a' and 'b' must be > 0.")
  if (n_samp < 0 || n_samp != as.integer(n_samp))
    stop("'n_samp' must be a nonnegative integer.")
  if (length(skew_mu) != 1L || is.na(skew_mu) || !is.logical(skew_mu))
    stop("'skew_mu' must be TRUE or FALSE.")
  if (skew_mu && n_samp == 0L)
    stop("The mu skew correction requires 'n_samp' > 0.")
  if (max_iter < 1 || max_iter != as.integer(max_iter))
    stop("'max_iter' must be a positive integer.")
  if (tol <= 0) stop("'tol' must be > 0.")
  if (!is.null(seed)) set.seed(seed)

  start_total <- Sys.time()
  start_fit <- Sys.time()

  if (init == "kmeans") {
    K_init <- min(K, n, max(1L, length(unique(y))))
    if (K_init == 1L) {
      cl <- rep.int(1L, n)
    } else {
      cl <- kmeans(y, 
                   centers = min(5,n))$cluster
    }
    r <- matrix(1e-3 / K, nrow = n, ncol = K)
    r[cbind(seq_len(n), cl)] <- 1 - 1e-3 + 1e-3 / K
  } else {
    r <- matrix(rexp(n * K), nrow = n, ncol = K)
    r <- r / rowSums(r)
  }

  sigma2 <- sigma^2
  s02 <- s0^2
  elbo <- numeric(max_iter)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    Nk <- colSums(r)

    mu_var <- 1 / (1 / s02 + Nk / sigma2)
    mu_mean <- mu_var * (m0 / s02 + as.vector(crossprod(r, y)) / sigma2)

    v_shape1 <- 1 + Nk[seq_len(K - 1L)]
    tail_counts <- rev(cumsum(rev(Nk)))[2:K]

    if (iter == 1L) {
      alpha_shape <- a + K - 1
      alpha_rate <- b + sum(digamma(v_shape1 + 1) - digamma(1))
    }
    ealpha <- alpha_shape / alpha_rate
    v_shape2 <- ealpha + tail_counts

    elog_1mv <- digamma(v_shape2) - digamma(v_shape1 + v_shape2)
    alpha_shape <- a + K - 1
    alpha_rate <- b - sum(elog_1mv)

    # Refresh q(v) using the newly updated E_q(alpha).
    v_shape2 <- alpha_shape / alpha_rate + tail_counts
    elog_v <- digamma(v_shape1) - digamma(v_shape1 + v_shape2)
    elog_1mv <- digamma(v_shape2) - digamma(v_shape1 + v_shape2)
    elog_pi <- c(elog_v, 0) + c(0, cumsum(elog_1mv))

    elog_lik <- -0.5 * log(2 * pi * sigma2) -
      0.5 * (outer(y, mu_mean, "-")^2 +
               matrix(mu_var, nrow = n, ncol = K, byrow = TRUE)) / sigma2
    logr <- sweep(elog_lik, 2L, elog_pi, "+")
    logr <- logr - .vb_logsumexp_rows(logr)
    r <- exp(logr)

    elbo[iter] <- .vb_elbo_normal_dpm(
      y, sigma, m0, s0, a, b, r, mu_mean, mu_var,
      v_shape1, v_shape2, alpha_shape, alpha_rate
    )

    if (verbose && (iter == 1L || iter %% 50L == 0L))
      message("VB iteration ", iter, "; ELBO = ", signif(elbo[iter], 10))

    if (iter > 1L) {
      relative_change <- abs(elbo[iter] - elbo[iter - 1L]) /
        (1 + abs(elbo[iter - 1L]))
      if (relative_change < tol) {
        converged <- TRUE
        break
      }
    }
  }

  # One final global update makes the returned factors consistent with the
  # final responsibility matrix.
  Nk <- colSums(r)
  mu_var <- 1 / (1 / s02 + Nk / sigma2)
  mu_mean <- mu_var * (m0 / s02 + as.vector(crossprod(r, y)) / sigma2)
  v_shape1 <- 1 + Nk[seq_len(K - 1L)]
  tail_counts <- rev(cumsum(rev(Nk)))[2:K]
  v_shape2 <- alpha_shape / alpha_rate + tail_counts
  elog_1mv <- digamma(v_shape2) - digamma(v_shape1 + v_shape2)
  alpha_shape <- a + K - 1
  alpha_rate <- b - sum(elog_1mv)
  v_shape2 <- alpha_shape / alpha_rate + tail_counts

  end_fit <- Sys.time()
  start_samp <- Sys.time()

  if (n_samp > 0L) {
    v_draws <- matrix(
      rbeta(n_samp * (K - 1L),
            shape1 = rep(v_shape1, each = n_samp),
            shape2 = rep(v_shape2, each = n_samp)),
      nrow = n_samp, ncol = K - 1L
    )
    eps <- sqrt(.Machine$double.eps)
    R_draws <- qlogis(pmin(pmax(v_draws, eps), 1 - eps))
    mu_draws <- matrix(
      rnorm(n_samp * K, mean = rep(mu_mean, each = n_samp),
            sd = rep(sqrt(mu_var), each = n_samp)),
      nrow = n_samp, ncol = K
    )

    alpha_draws <- rgamma(n_samp, shape = alpha_shape, rate = alpha_rate)

    if (skew_mu) {
      skew_result <- .vb_skew_mu_draws(
        mu_draws = mu_draws,
        v_draws = v_draws,
        mu_center = mu_mean,
        y = y,
        sigma = sigma,
        m0 = m0,
        s0 = s0
      )
      mu_draws <- skew_result$draws
      mu_reflected <- skew_result$reflected
    } else {
      mu_reflected <- rep(FALSE, n_samp)
    }

    samples <- cbind(R_draws, mu_draws, log(alpha_draws))
    colnames(samples) <- c(paste0("R", seq_len(K - 1L)),
                           paste0("mu", seq_len(K)), "tau")
  } else {
    samples <- matrix(numeric(0), nrow = 0L, ncol = 2L * K)
    mu_reflected <- logical(0L)
  }

  end_samp <- Sys.time()

  list(
    approximation = if (skew_mu) {
      "conditionally skewed variational Bayes"
    } else {
      "variational Bayes (Blei-Jordan)"
    },
    K = K,
    responsibilities = r,
    mean = list(mu = mu_mean, alpha = alpha_shape / alpha_rate),
    variational = list(
      mu_mean = mu_mean, mu_var = mu_var,
      v_shape1 = v_shape1, v_shape2 = v_shape2,
      alpha_shape = alpha_shape, alpha_rate = alpha_rate
    ),
    samples = samples,
    skew_mu = list(
      applied = skew_mu,
      reflected = mu_reflected,
      reflection_rate = if (length(mu_reflected)) mean(mu_reflected) else NA_real_
    ),
    elbo = elbo[seq_len(iter)],
    iterations = iter,
    converged = converged,
    time = list(
      total = end_samp - start_total,
      optimization = end_fit - start_fit,
      sampling = end_samp - start_samp
    ),
    call = match.call()
  )
}

# Exact variational posterior mean of the mixture density. 
posterior_density_dpm_vb <- function(
    xgrid, fit, sigma = 1, probs = c(0.025, 0.5, 0.975)) {
  if (sigma <= 0) stop("'sigma' must be > 0.")
  if (is.null(fit$variational)) stop("'fit' is not a dp_vb fit.")

  q <- fit$variational
  K <- length(q$mu_mean)
  ev <- q$v_shape1 / (q$v_shape1 + q$v_shape2)
  e1mv <- q$v_shape2 / (q$v_shape1 + q$v_shape2)
  mean_weights <- c(ev, 1) * c(1, cumprod(e1mv))

  # Integral of N(x | mu, sigma^2) against q(mu).
  component_sd <- sqrt(sigma^2 + q$mu_var)
  dens_mat <- vapply(seq_len(K), function(k)
    dnorm(xgrid, mean = q$mu_mean[k], sd = component_sd[k]),
    numeric(length(xgrid)))
  mean_density <- as.vector(dens_mat %*% mean_weights)

  ans <- list(xgrid = xgrid, mean_density = mean_density,
              mean_weights = mean_weights)

  # If variational draws were requested, return summaries in the same format
  # as posterior_density_dpm(), while retaining the exact mean above.
  if (!is.null(fit$samples) && nrow(fit$samples) > 0L) {
    mc <- posterior_density_dpm(xgrid, fit$samples, sigma,
                                alpha_random = TRUE, probs = probs)
    ans$density_draws <- mc$density_draws
    ans$median_density <- mc$median_density
    ans$quantiles <- mc$quantiles
    ans$probs <- probs
    ans$mc_mean_density <- mc$mean_density

    
    if (isTRUE(fit$skew_mu$applied)) {
      ans$unskewed_vb_mean_density <- mean_density
      ans$mean_density <- mc$mean_density
    }
  }
  ans
}
