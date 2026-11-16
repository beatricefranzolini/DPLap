# Plot the joint posterior of atoms attached to the two largest weights.
# Compatible with the final dp_vb, dp_lap (including skew), and dp_slice outputs.
#
# Rank the weights WITHIN EACH DRAW: x is the atom with the largest weight,
# y is the atom with the second-largest weight. Atoms are not sorted by value.
#
# fits: one fit, or a named list of fits for side-by-side panels.
# burn_slice / thin_slice: used only for slice fits; all approximation draws remain.
# mu0 / tau20: slice base-measure mean and VARIANCE, matching hyper[2:3].
# type: scatter, contour, or both. Contours are KDE density levels, not
#       specified-probability credible regions. adjust controls KDE smoothing.
#
# Slice: include all stored weights, including unoccupied atoms. If the omitted
# tail could change the top two, continue its Beta(1, alpha) stick-breaking law.
# Stop as soon as the remaining mass cannot exceed the current second weight.
# Only atoms that enter the current top two need a draw from the base measure.
# This is post-processing and does not change the slice sampler.
# Tail completion consumes RNG draws; call set.seed() beforehand for repeatability.
#
# Return: a ggplot. Its $data contains the extracted draws, weights, and the
# number of additional tail sticks sampled for each retained iteration.
# R was not available in the authoring environment; this file needs an R run.

plot_top2_atoms <- function(
    fits, burn_slice = 2000L, thin_slice = 1L,
    mu0 = 0, tau20 = 1,
    type = c("both", "contour", "scatter"),
    bins = 8L, adjust = 1, point_alpha = 0.08,
    facet_ncol = 2L, max_tail = 100000L, 
    xlim_manual = c(-3,3), ylim_manual = c(-3,3)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  type <- match.arg(type)
  stopifnot(length(burn_slice) == 1L, is.finite(burn_slice),
            burn_slice >= 0, burn_slice == floor(burn_slice),
            length(thin_slice) == 1L, is.finite(thin_slice),
            thin_slice >= 1, thin_slice == floor(thin_slice),
            length(mu0) == 1L, is.finite(mu0),
            length(tau20) == 1L, is.finite(tau20), tau20 > 0,
            all(is.finite(adjust)), length(adjust) %in% c(1L, 2L),
            all(adjust > 0), length(max_tail) == 1L,
            is.finite(max_tail), max_tail >= 1)

  is_fit <- function(x) {
    is.list(x) && (is.matrix(x$samples) ||
      (is.list(x$pi) && is.list(x$phis)))
  }
  if (is_fit(fits)) {
    label <- if (!is.null(fits$approximation)) fits$approximation else "Slice"
    fits <- setNames(list(fits), label)
  }
  if (!is.list(fits) || !length(fits) ||
      !all(vapply(fits, is_fit, logical(1)))) {
    stop("Supply one fitted object, or a list of fitted objects (e.g. VB[[i]]).")
  }
  labels <- names(fits)
  if (is.null(labels)) labels <- paste("Method", seq_along(fits))
  missing <- is.na(labels) | !nzchar(labels)
  labels[missing] <- paste("Method", which(missing))
  names(fits) <- make.unique(labels)

  extract <- function(fit, method) {
    if (is.matrix(fit$samples)) {
      # Infer the actual K from the stored draws, including Laplace's merging.
      M <- fit$samples
      if (!is.numeric(M) || !nrow(M) || ncol(M) < 4L ||
          ncol(M) %% 2L || any(!is.finite(M))) {
        stop(method, ": expected finite parameter draws with 2*K columns, K >= 2.")
      }
      K <- ncol(M) %/% 2L
      keep <- seq_len(nrow(M))
      ans <- matrix(NA_real_, length(keep), 5L)
      for (s in keep) {
        r <- M[s, seq_len(K - 1L)]
        correction <- log1p(exp(-abs(r)))
        log_v <- -pmax(-r, 0) - correction
        log_1mv <- -pmax(r, 0) - correction
        prefix <- c(0, cumsum(log_1mv))
        log_w <- c(log_v + prefix[seq_len(K - 1L)], prefix[K])
        j <- order(log_w, decreasing = TRUE)[1:2]
        atoms <- M[s, K + seq_len(K) - 1L]
        ans[s, ] <- c(atoms[j], exp(log_w[j]), 0)
      }
    } else {
      total <- length(fit$pi)
      if (length(fit$phis) != total || burn_slice >= total) {
        stop(method, ": incompatible pi/phis lists or no draws left after burn-in.")
      }
      keep <- seq.int(burn_slice + 1L, total, by = thin_slice)
      ans <- matrix(NA_real_, length(keep), 5L)
      for (s in seq_along(keep)) {
        tt <- keep[s]
        w <- as.numeric(fit$pi[[tt]])
        atoms <- as.numeric(fit$phis[[tt]])
        if (!length(w) || length(w) != length(atoms) ||
            any(!is.finite(w)) || any(w < 0) ||
            any(!is.finite(atoms)) || sum(w) <= 0 || sum(w) > 1 + 1e-10) {
          stop(method, ": invalid weights/atoms at iteration ", tt)
        }
        j <- head(order(w, decreasing = TRUE), 2L)
        top_w <- c(w[j], rep(0, 2L - length(j)))
        top_atoms <- c(atoms[j], rep(NA_real_, 2L - length(j)))
        remaining <- max(0, 1 - sum(w))
        steps <- 0L
        # Every still-unseen weight is bounded by the remaining TOTAL mass.
        while (remaining > top_w[2L] || top_w[2L] <= 0) {
          if (remaining <= 0) {
            stop(method, ": fewer than two positive weights and no recoverable ",
                 "tail mass at iteration ", tt)
          }
          alpha <- fit$alpha[tt]
          if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
            stop(method, ": a positive stored alpha is needed at iteration ", tt)
          }
          steps <- steps + 1L
          if (steps > max_tail) {
            stop(method, ": increase max_tail to finish iteration ", tt)
          }
          v <- stats::rbeta(1L, 1, alpha)
          new_w <- remaining * v
          remaining <- remaining * (1 - v)
          if (new_w > top_w[2L]) {
            candidate_w <- c(top_w, new_w)
            candidate_atoms <- c(top_atoms, stats::rnorm(1L, mu0, sqrt(tau20)))
            j <- order(candidate_w, decreasing = TRUE)[1:2]
            top_w <- candidate_w[j]
            top_atoms <- candidate_atoms[j]
          }
        }
        ans[s, ] <- c(top_atoms, top_w, steps)
      }
    }
    data.frame(method = method, draw = keep,
               atom_1 = ans[, 1L], atom_2 = ans[, 2L],
               weight_1 = ans[, 3L], weight_2 = ans[, 4L],
               tail_steps = as.integer(ans[, 5L]))
  }

  d <- do.call(rbind, lapply(seq_along(fits), function(i) {
    extract(fits[[i]], names(fits)[i])
  }))
  rownames(d) <- NULL
  d$method <- factor(d$method, levels = names(fits))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = atom_1, y = atom_2))
  if (type %in% c("scatter", "both")) {
    p <- p + ggplot2::geom_point(alpha = point_alpha, size = 0.6, colour = "grey25")
  }
  if (type %in% c("contour", "both")) {
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Contours require package 'MASS'.")
    valid <- vapply(split(d, d$method), function(z) {
      nrow(z) >= 3L && length(unique(z$atom_1)) > 1L &&
        length(unique(z$atom_2)) > 1L
    }, logical(1))
    good <- d$method %in% names(valid)[valid]
    if (any(good)) {
      p <- p + ggplot2::geom_density_2d(
        data = d[good, , drop = FALSE], bins = bins, adjust = adjust,
        contour_var = "density", colour = "#0072B2", linewidth = 0.6
      )
    }
    if (any(!valid)) {
      warning("Using points only for degenerate/too-small samples: ",
              paste(names(valid)[!valid], collapse = ", "))
      if (type == "contour") {
        p <- p + ggplot2::geom_point(data = d[!good, , drop = FALSE],
                                    alpha = point_alpha, size = 0.6)
      }
    }
  }
  p <- p + ggplot2::labs(
    x = "Atom with largest weight", y = "Atom with second-largest weight",
    subtitle = "Weights ranked separately within each draw"
  ) + ggplot2::coord_equal(xlim = xlim_manual, ylim = ylim_manual, expand = FALSE) + ggplot2::theme_bw(base_size = 12)
  if (length(fits) > 1L) {
    p <- p + ggplot2::facet_wrap(~method, ncol = facet_ncol, scales = "fixed")
  } else {
    p <- p + ggplot2::labs(title = names(fits)[1L])
  }
  p
}

# Example: one sample size, four algorithms.
# source("../plot_top2_atoms.R")
# i <- length(sample_sizes)
# set.seed(123)  # controls any necessary slice-tail completion
# p <- plot_top2_atoms(
#   list(VB = VB[[i]], Laplace = LAP[[i]],
#        `Skew-Laplace` = Skew_LAP[[i]], Slice = SS[[i]]),
#   burn_slice = 2000, type = "both"
# )
# print(p)
# ggplot2::ggsave("top_two_atoms.pdf", p, width = 9, height = 7)
#
# Single fit: plot_top2_atoms(LAP[[i]], type = "contour")
# Extracted values: head(p$data)
# KDE reference: https://ggplot2.tidyverse.org/reference/geom_density_2d.html
