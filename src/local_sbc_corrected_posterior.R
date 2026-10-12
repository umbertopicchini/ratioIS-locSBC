# ============================================================
# Refitted RatioIS localized simulation-based calibration
# ============================================================
#
# PRIMARY IMPLEMENTATION
# ----------------------
# For each fixed lambda >= 1:
#
#   1. Construct the local design g_lambda(theta | y_obs) by inflating the
#      within-component covariance matrices of the already fitted final q2
#      posterior mixture at y_obs:
#
#        g_lambda(theta | y_obs)
#          = sum_k omega_k(y_obs)
#              N(theta; m_k(y_obs), lambda^2 Sigma_k).
#
#   2. Generate an independent LOCAL TRAINING BANK
#
#        theta_n^tr ~ g_lambda,
#        y_n^tr     ~ p(y | theta_n^tr),       n = 1,...,N_local.
#
#   3. Fit ONE new local GLLiM model to the training bank (warm-started from
#      the existing q2 fit by default). Use that local GLLiM to construct the
#      local learned summary T_lambda(y): posterior mean, optional covariance,
#      optional mixture weights, followed by local standardization.
#
#   4. Fit ONE new ratio classifier on
#
#        (theta_n^tr, T_lambda(y_n^tr)),
#
#      using the same feature construction as the main RatioIS classifier by
#      default. Ideally this estimates
#
#        r_lambda(theta,T)
#          = p_Tlambda(T | theta) / m_g_lambda(T).
#
#   5. Freeze the local GLLiM, summary preprocessing, and ratio classifier.
#      Generate an INDEPENDENT LOCAL CHECKING BANK
#
#        theta_i^* ~ g_lambda,
#        y_i^*     ~ p(y | theta_i^*),         i = 1,...,N_local_check.
#
#   6. For each pseudo-observation y_i^*, draw cheap inner particles
#
#        theta_i^(m) ~ g_lambda,               m = 1,...,M_is,
#
#      and evaluate the newly fitted ratio at T_lambda(y_i^*). Since the
#      pseudo-prior and inner IS proposal are both g_lambda,
#
#        w_i^(m) proportional to rhat_lambda(theta_i^(m),T_lambda(y_i^*)).
#
#      Raw weights are normalized, then the SAME deterministic clipping rule
#      used for the reported RatioIS posterior is applied:
#
#        M_c = floor(sqrt(M_is)),
#        c   = M_c-th largest normalized raw weight,
#        w_clip^(m) proportional to min(w_raw^(m), c).
#
#      Hence the primary localized pseudo-posterior is the Eq.-(23)-type
#      clipped weighted empirical measure
#
#        Pi_hat_clip,lambda(dtheta | y_i^*)
#          = sum_m wbar_clip,lambda^(m)(y_i^*) delta_theta_i^(m)(dtheta).
#
#      For PIT computation we use this weighted empirical CDF directly; an
#      additional multinomial resampling step is unnecessary and would only
#      add Monte Carlo noise.
#
#   7. Compute marginal PIT values at theta_i^*. CLIPPED weights are primary;
#      raw-weight PITs are retained as a diagnostic comparison.
#
# INTERPRETATION
# --------------
# This tests whether the GLLiM-summary + ratio-estimation + clipped-IS
# architecture can learn calibrated conditional inference under the chosen
# local experiment g_lambda. It does NOT test whether g_1 itself equals the
# true posterior for y_obs.
#
# Computational cost per lambda:
#   - N_local + N_local_check DGM simulations (default 2*N_local);
#   - one local GLLiM refit;
#   - one local ratio-classifier refit;
#   - N_local_check * M_is cheap g_lambda draws and ratio evaluations.
#
# Expected dependencies sourced before this file:
#   gllim_strict_em.R
#   build_surrogates_v14_two_round_prior_defensive.R
#   ratio_estimation_v10.R
# ============================================================

library(mvtnorm)


# ------------------------------------------------------------
# 0. Small utilities
# ------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.require_local_sbc_function <- function(name) {
  if (!exists(name, mode = "function", inherits = TRUE)) {
    stop(
      "Required function '", name, "' is not available. ",
      "Source build_surrogates_v14_two_round_prior_defensive.R and ",
      "ratio_estimation_v10.R before local_sbc_corrected_posterior_v2.R."
    )
  }
  invisible(TRUE)
}

.normalize_gmm <- function(probs, means, covars, param_names = NULL) {
  probs <- as.numeric(probs)
  if (length(probs) < 1 || any(!is.finite(probs)) || any(probs < 0) || sum(probs) <= 0) {
    stop("Invalid Gaussian-mixture probabilities.")
  }
  probs <- probs / sum(probs)

  means <- as.matrix(means)
  if (ncol(means) != length(probs) && nrow(means) == length(probs)) {
    means <- t(means)
  }
  if (ncol(means) != length(probs)) {
    stop("Mixture means are incompatible with the number of components.")
  }

  p <- nrow(means)
  K <- length(probs)

  if (length(dim(covars)) == 2) {
    if (K != 1) stop("A 2D covariance matrix is only valid for K = 1.")
    covars <- array(as.matrix(covars), dim = c(p, p, 1))
  } else {
    covars <- as.array(covars)
  }

  if (!identical(dim(covars), c(p, p, K))) {
    stop("Mixture covariance array has incompatible dimensions.")
  }

  for (k in seq_len(K)) {
    covars[, , k] <- 0.5 * (covars[, , k] + t(covars[, , k]))
  }

  if (is.null(param_names)) param_names <- rownames(means)
  if (is.null(param_names)) param_names <- paste0("theta", seq_len(p))
  rownames(means) <- param_names

  list(
    K = K,
    dim_theta = p,
    probs = probs,
    means = means,
    covars = covars,
    param_names = param_names
  )
}

.log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# theta_mat is p x N.
gmm_log_density <- function(theta_mat, g_mix) {
  theta_mat <- as.matrix(theta_mat)
  g_mix <- .normalize_gmm(
    g_mix$probs, g_mix$means, g_mix$covars,
    param_names = g_mix$param_names
  )

  if (nrow(theta_mat) != g_mix$dim_theta) {
    if (ncol(theta_mat) == g_mix$dim_theta) theta_mat <- t(theta_mat)
  }
  if (nrow(theta_mat) != g_mix$dim_theta) {
    stop("theta_mat has incompatible dimension in gmm_log_density().")
  }

  N <- ncol(theta_mat)
  K <- g_mix$K
  log_comp <- matrix(-Inf, nrow = K, ncol = N)
  X <- t(theta_mat)

  for (k in seq_len(K)) {
    log_comp[k, ] <- log(g_mix$probs[k]) + mvtnorm::dmvnorm(
      X,
      mean = g_mix$means[, k],
      sigma = g_mix$covars[, , k],
      log = TRUE
    )
  }

  if (K == 1) return(as.numeric(log_comp[1, ]))
  apply(log_comp, 2, .log_sum_exp)
}

# Draw from a GMM, rejecting points outside prior support when prior_pdf is given.
sample_valid_from_gmm <- function(n,
                                  g_mix,
                                  prior_pdf = NULL,
                                  max_attempts = NULL,
                                  verbose = FALSE) {
  g_mix <- .normalize_gmm(
    g_mix$probs, g_mix$means, g_mix$covars,
    param_names = g_mix$param_names
  )

  if (length(n) != 1 || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  n <- as.integer(n)

  if (is.null(max_attempts)) max_attempts <- max(10000L, 100L * n)

  p <- g_mix$dim_theta
  out <- matrix(NA_real_, nrow = p, ncol = n)
  rownames(out) <- g_mix$param_names

  accepted <- 0L
  attempted <- 0L

  while (accepted < n && attempted < max_attempts) {
    remaining <- n - accepted
    batch_n <- min(max_attempts - attempted, max(remaining, ceiling(1.25 * remaining)))

    z <- sample.int(g_mix$K, size = batch_n, replace = TRUE, prob = g_mix$probs)
    cand <- matrix(NA_real_, nrow = p, ncol = batch_n)

    for (k in unique(z)) {
      id <- which(z == k)
      cand[, id] <- t(mvtnorm::rmvnorm(
        length(id),
        mean = g_mix$means[, k],
        sigma = g_mix$covars[, , k]
      ))
    }

    attempted <- attempted + batch_n
    ok <- apply(cand, 2, function(x) all(is.finite(x)))

    if (!is.null(prior_pdf)) {
      ok_idx <- which(ok)
      if (length(ok_idx) > 0) {
        pv <- vapply(ok_idx, function(j) {
          val <- prior_pdf(cand[, j])
          if (length(val) != 1) val <- val[1]
          as.numeric(val)
        }, numeric(1))
        ok[ok_idx] <- is.finite(pv) & pv > 0
      }
    }

    keep <- which(ok)
    if (length(keep) > 0) {
      take <- keep[seq_len(min(length(keep), remaining))]
      jj <- accepted + seq_along(take)
      out[, jj] <- cand[, take, drop = FALSE]
      accepted <- accepted + length(take)
    }

    if (verbose && attempted %% 10000 < batch_n) {
      cat("GMM sampling attempts:", attempted, "accepted:", accepted, "of", n, "\n")
    }
  }

  if (accepted < n) {
    stop(
      "Could not obtain enough prior-valid GMM draws. Accepted ", accepted,
      " of ", n, " after ", attempted, " proposals."
    )
  }

  list(
    theta = out,
    attempts = attempted,
    acceptance_rate = n / attempted
  )
}


# ------------------------------------------------------------
# 1. Recover the already fitted q2(theta | y)
# ------------------------------------------------------------

get_q2_cache <- function(surrogate_out, ratio_out = NULL) {
  if (!is.null(ratio_out) && !is.null(ratio_out$last_gllim_cache)) {
    return(ratio_out$last_gllim_cache)
  }

  .require_local_sbc_function("prepare_gllim_posterior_cache")

  prepare_gllim_posterior_cache(
    fit_obj = surrogate_out$fit1_raw,
    num_param = surrogate_out$num_param,
    dim_data = surrogate_out$dim_data
  )
}


q2_theta_mixture_at_y <- function(y, surrogate_out, ratio_out = NULL, param_names = NULL) {
  .require_local_sbc_function("posterior_mixture_from_input")
  cache <- get_q2_cache(surrogate_out, ratio_out)
  gmm <- posterior_mixture_from_input(cache = cache, x = as.numeric(y))
  .normalize_gmm(gmm$probs, gmm$means, gmm$covars, param_names = param_names)
}


# ------------------------------------------------------------
# 2. Build Eq. (25) directly from q2(theta | y_obs)
# ------------------------------------------------------------

build_g_lambda_from_q2 <- function(surrogate_out,
                                   ratio_out = NULL,
                                   observedData = NULL,
                                   lambda = 1,
                                   param_names = NULL,
                                   ridge = 0) {
  if (length(lambda) != 1 || !is.finite(lambda) || lambda < 1) {
    stop("lambda must be a finite scalar >= 1.")
  }

  if (is.null(observedData)) {
    observedData <- ratio_out$observedData %||% surrogate_out$observedData
  }
  if (is.null(observedData)) stop("observedData is required.")

  q2_obs <- q2_theta_mixture_at_y(
    y = observedData,
    surrogate_out = surrogate_out,
    ratio_out = ratio_out,
    param_names = param_names
  )

  covars <- q2_obs$covars
  for (k in seq_len(q2_obs$K)) {
    S <- q2_obs$covars[, , k]
    S <- 0.5 * (S + t(S))
    if (ridge > 0) S <- S + diag(ridge, q2_obs$dim_theta)
    covars[, , k] <- lambda^2 * S
  }

  out <- .normalize_gmm(
    q2_obs$probs,
    q2_obs$means,
    covars,
    param_names = param_names
  )
  out$lambda <- lambda
  out$base <- "q2(theta | y_obs)"
  out
}


# ------------------------------------------------------------
# 3. Reconstruct exactly the classifier statistic T(y)
# ------------------------------------------------------------

ratio_summary_T <- function(y,
                            surrogate_out,
                            ratio_out = NULL) {
  .require_local_sbc_function("posterior_theta_summary_from_cache")

  cache <- get_q2_cache(surrogate_out, ratio_out)
  cfg <- surrogate_out$summary_config %||% list()

  include_cov <- cfg$include_covariance_in_S1
  if (is.null(include_cov)) include_cov <- TRUE

  include_weights <- cfg$include_mixture_weights_in_S1
  if (is.null(include_weights)) include_weights <- TRUE

  raw <- posterior_theta_summary_from_cache(
    y = as.numeric(y),
    cache = cache,
    num_theta = surrogate_out$num_param,
    include_covariance = include_cov,
    include_mixture_weights = include_weights
  )

  std_info <- surrogate_out$S1_standardization
  if (is.null(std_info) || is.null(std_info$center) || is.null(std_info$scale)) {
    stop("surrogate_out$S1_standardization is required to reconstruct T(y).")
  }

  center <- as.numeric(std_info$center)
  scale <- as.numeric(std_info$scale)
  if (length(raw) != length(center) || length(raw) != length(scale)) {
    stop("Raw T(y) dimension is incompatible with the stored standardization.")
  }

  z <- (as.numeric(raw) - center) / scale


  z
}

validate_ratio_summary_reconstruction <- function(surrogate_out,
                                                  ratio_out = NULL,
                                                  observedData = NULL,
                                                  tol = 1e-6,
                                                  verbose = TRUE) {
  if (is.null(observedData)) {
    observedData <- ratio_out$observedData %||% surrogate_out$observedData
  }
  if (is.null(surrogate_out$S1_observedData)) return(invisible(NA_real_))

  T_rebuilt <- ratio_summary_T(observedData, surrogate_out, ratio_out)
  T_stored <- as.numeric(surrogate_out$S1_observedData)

  if (length(T_rebuilt) != length(T_stored)) {
    stop("Reconstructed T(y_obs) has a different dimension from S1_observedData.")
  }

  err <- max(abs(T_rebuilt - T_stored))
  if (verbose) cat("Max |reconstructed T(y_obs) - stored T(y_obs)| =", err, "\n")
  if (!is.finite(err) || err > tol) {
    warning(
      "The reconstructed classifier summary does not match the stored T(y_obs): max error = ",
      signif(err, 5), ". Check summary preprocessing before interpreting SBC results."
    )
  }
  invisible(err)
}


# ------------------------------------------------------------
# 4. Simulate pseudo-data from the local design
# ------------------------------------------------------------

simulate_local_bank <- function(theta_mat,
                                model,
                                model_param = NULL,
                                verbose = TRUE) {
  theta_mat <- as.matrix(theta_mat)
  n <- ncol(theta_mat)
  y_list <- vector("list", n)

  for (i in seq_len(n)) {
    if (verbose && (i %% 100 == 0 || i == n)) {
      cat("Simulated", i, "of", n, "local SBC datasets\n")
    }

    if (is.null(model_param)) {
      y_i <- model(theta_mat[, i])
    } else {
      y_i <- model(theta_mat[, i], model_param = model_param)
    }
    y_list[[i]] <- as.numeric(y_i)
  }

  y_mat <- do.call(cbind, y_list)
  if (!all(is.finite(y_mat))) stop("Non-finite values found in local SBC simulations.")
  y_mat
}


# ------------------------------------------------------------
# 5. Importance weights for one pseudo-observation
# ------------------------------------------------------------

normalize_log_weights <- function(log_w) {
  log_w <- as.numeric(log_w)
  finite <- is.finite(log_w)
  if (!any(finite)) stop("All importance log-weights are non-finite.")

  m <- max(log_w[finite])
  w <- numeric(length(log_w))
  w[finite] <- exp(log_w[finite] - m)
  w[!is.finite(w)] <- 0

  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) stop("All importance weights are zero after stabilization.")
  w / sw
}

clip_ratiois_weights <- function(weights) {
  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  if (sum(weights) <= 0) stop("Cannot clip zero importance weights.")
  weights <- weights / sum(weights)

  M <- length(weights)
  Nc <- max(1L, floor(sqrt(M)))
  ord <- order(weights, decreasing = TRUE)
  cap <- weights[ord[Nc]]
  w_clip <- pmin(weights, cap)
  w_clip / sum(w_clip)
}

weighted_empirical_pit <- function(theta_star, theta_prop, weights) {
  theta_star <- as.numeric(theta_star)
  theta_prop <- as.matrix(theta_prop)
  weights <- as.numeric(weights)
  weights <- weights / sum(weights)

  if (nrow(theta_prop) != length(theta_star)) {
    stop("theta_star and theta_prop have incompatible dimensions.")
  }
  if (ncol(theta_prop) != length(weights)) {
    stop("theta_prop and weights have incompatible lengths.")
  }

  vapply(seq_along(theta_star), function(j) {
    sum(weights[theta_prop[j, ] <= theta_star[j]])
  }, numeric(1))
}

local_ratio_pseudo_posterior_pit <- function(theta_star,
                                         y_star,
                                         g_mix,
                                         surrogate_out,
                                         ratio_out,
                                         prior_pdf = NULL,
                                         M_is = 5000,
                                         ratio_lambda_choice = NULL,
                                         verbose = FALSE) {
  .require_local_sbc_function("predict_log_ratio_batch")

  if (is.null(ratio_out$ratio_fit)) {
    stop("ratio_out$ratio_fit is required for no-refit localized SBC.")
  }

  # ----------------------------------------------------------
  # Primary construction: the inner proposal is g_lambda itself.
  # ----------------------------------------------------------
  # For fixed y_star,
  #
  #   target(theta) proportional to
  #     g_lambda(theta) * p_T(T_star | theta),
  #
  # while the fitted classifier supplies
  #
  #   r_hat(theta, T_star)
  #     approximately p_T(T_star | theta) / m_train(T_star).
  #
  # Since m_train(T_star) is constant in theta and theta is proposed
  # directly from g_lambda, the self-normalized IS weights reduce to
  #
  #   w(theta) proportional to r_hat(theta, T_star).
  # ----------------------------------------------------------

  prop <- sample_valid_from_gmm(
    n = M_is,
    g_mix = g_mix,
    prior_pdf = prior_pdf,
    verbose = FALSE
  )
  theta_prop <- prop$theta

  # Reconstruct exactly the summary preprocessing used by the fitted classifier.
  T_star <- ratio_summary_T(
    y = y_star,
    surrogate_out = surrogate_out,
    ratio_out = ratio_out
  )

  log_r <- predict_log_ratio_batch(
    ratio_fit = ratio_out$ratio_fit,
    theta_mat = t(theta_prop),
    y_used = T_star,
    s = ratio_lambda_choice
  )

  # Because proposal = g_lambda, the g_lambda density cancels from the
  # importance ratio. The denominator m_train(T_star) also cancels after
  # self-normalization, leaving only the fitted ratio values.
  w_raw <- normalize_log_weights(log_r)
  w_clip <- clip_ratiois_weights(w_raw)

  U_raw <- weighted_empirical_pit(theta_star, theta_prop, w_raw)
  U_clip <- weighted_empirical_pit(theta_star, theta_prop, w_clip)

  finite_lr <- log_r[is.finite(log_r)]

  list(
    U_raw = U_raw,
    U_clipped = U_clip,
    ESS_raw = 1 / sum(w_raw^2),
    ESS_clipped = 1 / sum(w_clip^2),
    max_weight_raw = max(w_raw),
    max_weight_clipped = max(w_clip),
    proposal_acceptance_rate = prop$acceptance_rate,
    proposal_attempts = prop$attempts,
    log_ratio_sd = if (length(finite_lr) > 1) stats::sd(finite_lr) else NA_real_,
    log_ratio_min = if (length(finite_lr) > 0) min(finite_lr) else NA_real_,
    log_ratio_max = if (length(finite_lr) > 0) max(finite_lr) else NA_real_,
    T_star = T_star
  )
}


# ------------------------------------------------------------
# 6. PIT summaries
# ------------------------------------------------------------

summarise_pit <- function(U) {
  U <- as.matrix(U)

  data.frame(
    parameter = colnames(U),
    mean = apply(U, 2, mean, na.rm = TRUE),
    sd = apply(U, 2, sd, na.rm = TRUE),
    q05 = apply(U, 2, stats::quantile, probs = 0.05, na.rm = TRUE),
    q25 = apply(U, 2, stats::quantile, probs = 0.25, na.rm = TRUE),
    q50 = apply(U, 2, stats::quantile, probs = 0.50, na.rm = TRUE),
    q75 = apply(U, 2, stats::quantile, probs = 0.75, na.rm = TRUE),
    q95 = apply(U, 2, stats::quantile, probs = 0.95, na.rm = TRUE),
    ks_pvalue = apply(U, 2, function(u) {
      u <- u[is.finite(u)]
      if (length(u) < 5) return(NA_real_)
      suppressWarnings(stats::ks.test(u, "punif")$p.value)
    }),
    row.names = NULL
  )
}


# ------------------------------------------------------------
# 7. Historical PIT histogram and ordinary ECDF plots
# ------------------------------------------------------------

plot_pit_histograms <- function(U,
                                file = NULL,
                                breaks = 20,
                                main_prefix = "Local SBC PIT",
                                cex.main = 0.8,
                                band_level = 0.95,
                                band_col = grDevices::adjustcolor("grey45", alpha.f = 0.35),
                                hist_col = "grey90",
                                hist_border = "grey30") {
  U <- as.matrix(U)

  if (!is.null(file)) {
    grDevices::pdf(file, width = 8, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  p <- ncol(U)
  nr <- ceiling(sqrt(p))
  nc <- ceiling(p / nr)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nr, nc), mar = c(4, 4, 2.5, 1))

  for (j in seq_len(p)) {
    u <- U[, j]
    u <- u[is.finite(u)]
    N <- length(u)

    h <- hist(u, breaks = breaks, plot = FALSE, xlim = c(0, 1))
    bin_widths <- diff(h$breaks)
    if (max(abs(bin_widths - bin_widths[1])) > 1e-12) {
      stop("Reference-band calculation currently requires equal-width histogram bins.")
    }

    B <- length(h$counts)
    alpha <- 1 - band_level
    p_bin <- 1 / B
    lower_count <- stats::qbinom(alpha / 2, size = N, prob = p_bin)
    upper_count <- stats::qbinom(1 - alpha / 2, size = N, prob = p_bin)
    lower_density <- B * lower_count / N
    upper_density <- B * upper_count / N

    ylim_use <- range(c(0, h$density, upper_density), finite = TRUE)

    plot(
      h,
      freq = FALSE,
      xlim = c(0, 1),
      ylim = ylim_use,
      main = paste(main_prefix, colnames(U)[j]),
      cex.main = cex.main,
      xlab = "PIT value",
      ylab = "Density",
      col = hist_col,
      border = hist_border
    )

    rect(
      xleft = 0, xright = 1,
      ybottom = lower_density, ytop = upper_density,
      border = NA, col = band_col
    )
    plot(h, freq = FALSE, add = TRUE, col = NA, border = hist_border)
    abline(h = 1, lty = 2, lwd = 1.2)
  }

  invisible(NULL)
}

plot_pit_ecdf <- function(U,
                          file = NULL,
                          main_prefix = "Local SBC PIT ECDF",
                          cex.main = 0.8) {
  U <- as.matrix(U)

  if (!is.null(file)) {
    grDevices::pdf(file, width = 8, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  p <- ncol(U)
  nr <- ceiling(sqrt(p))
  nc <- ceiling(p / nr)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nr, nc), mar = c(4, 4, 2.5, 1))

  for (j in seq_len(p)) {
    u <- U[, j]
    u <- u[is.finite(u)]

    plot(
      stats::ecdf(u),
      verticals = TRUE,
      do.points = FALSE,
      xlim = c(0, 1), ylim = c(0, 1),
      main = paste(main_prefix, colnames(U)[j]),
      cex.main = cex.main,
      xlab = "u", ylab = "Empirical CDF"
    )
    abline(0, 1, lty = 2)
  }

  invisible(NULL)
}


# ------------------------------------------------------------
# 8. Simultaneous ECDF bands for continuous PIT values
# ------------------------------------------------------------
# Independent implementation of the simulation-based construction described
# by Sailynoja, Burkner & Vehtari (2022, Statistics and Computing).
# The code does NOT source or copy the authors' implementation.

make_ecdf_simultaneous_band <- function(N,
                                        level = 0.95,
                                        n_sim = 5000,
                                        grid_size = min(N, 400),
                                        seed = 271828) {
  if (length(N) != 1 || !is.finite(N) || N < 5 || N != as.integer(N)) {
    stop("N must be an integer >= 5.")
  }
  if (level <= 0 || level >= 1) stop("level must lie in (0,1).")
  if (n_sim < 100) stop("n_sim should be at least 100.")

  N <- as.integer(N)
  grid_size <- max(20L, min(as.integer(grid_size), N))
  z <- seq(1 / (grid_size + 1), grid_size / (grid_size + 1), length.out = grid_size)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  probs <- diff(c(0, z, 1))
  interval_counts <- stats::rmultinom(n_sim, size = N, prob = probs)
  ecdf_counts <- apply(interval_counts, 2, cumsum)[seq_len(grid_size), , drop = FALSE]
  ecdf_sim <- ecdf_counts / N

  target <- level
  alpha <- 1 - level

  coverage_at_gamma <- function(gamma) {
    lo <- stats::qbinom(gamma / 2, size = N, prob = z) / N
    hi <- stats::qbinom(1 - gamma / 2, size = N, prob = z) / N
    inside <- sweep(ecdf_sim, 1, lo, FUN = ">=") &
      sweep(ecdf_sim, 1, hi, FUN = "<=")
    mean(colSums(!inside) == 0)
  }

  g_lo <- 1e-10
  g_hi <- alpha

  if (coverage_at_gamma(g_hi) >= target) {
    gamma_star <- g_hi
  } else {
    for (iter in seq_len(28)) {
      g_mid <- 0.5 * (g_lo + g_hi)
      cov_mid <- coverage_at_gamma(g_mid)
      if (cov_mid >= target) {
        g_lo <- g_mid
      } else {
        g_hi <- g_mid
      }
    }
    gamma_star <- g_lo
  }

  lower <- stats::qbinom(gamma_star / 2, size = N, prob = z) / N
  upper <- stats::qbinom(1 - gamma_star / 2, size = N, prob = z) / N
  achieved <- coverage_at_gamma(gamma_star)

  list(
    z = c(0, z, 1),
    lower = c(0, lower, 1),
    upper = c(0, upper, 1),
    lower_diff = c(0, lower - z, 0),
    upper_diff = c(0, upper - z, 0),
    level = level,
    gamma = gamma_star,
    achieved_coverage_mc = achieved,
    N = N,
    n_sim = n_sim,
    grid_size = grid_size,
    seed = seed
  )
}

plot_pit_ecdf_simultaneous <- function(U,
                                       file = NULL,
                                       main_prefix = "Local SBC ECDF with simultaneous band",
                                       level = 0.95,
                                       cex.main = 0.7,
                                       n_sim = 5000,
                                       grid_size = NULL,
                                       seed = 271828,
                                       band_col = "grey85") {
  U <- as.matrix(U)

  if (!is.null(file)) {
    grDevices::pdf(file, width = 8, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  p <- ncol(U)
  nr <- ceiling(sqrt(p))
  nc <- ceiling(p / nr)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nr, nc), mar = c(4, 4, 2.5, 1))

  bands <- list()

  for (j in seq_len(p)) {
    u <- U[, j]
    u <- u[is.finite(u)]
    N <- length(u)
    key <- as.character(N)

    if (is.null(bands[[key]])) {
      gs <- if (is.null(grid_size)) min(N, 400) else grid_size
      bands[[key]] <- make_ecdf_simultaneous_band(
        N = N, level = level, n_sim = n_sim,
        grid_size = gs, seed = seed
      )
    }
    band <- bands[[key]]

    plot(
      NA,
      xlim = c(0, 1), ylim = c(0, 1),
      main = paste(main_prefix, colnames(U)[j]),
      cex.main = cex.main,
      xlab = "u", ylab = "Empirical CDF"
    )
    polygon(
      c(band$z, rev(band$z)),
      c(band$lower, rev(band$upper)),
      col = band_col, border = NA
    )
    abline(0, 1, lty = 2)
    Fhat <- stats::ecdf(u)(band$z)
    lines(band$z, Fhat, type = "s", lwd = 1.2)
  }

  invisible(bands)
}

plot_pit_ecdf_difference_simultaneous <- function(U,
                                                  file = NULL,
                                                  main_prefix = "Local SBC ECDF difference",
                                                  level = 0.95,
                                                  cex.main = 0.7,
                                                  n_sim = 5000,
                                                  grid_size = NULL,
                                                  seed = 271828,
                                                  band_col = "grey85") {
  U <- as.matrix(U)

  if (!is.null(file)) {
    grDevices::pdf(file, width = 8, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  p <- ncol(U)
  nr <- ceiling(sqrt(p))
  nc <- ceiling(p / nr)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nr, nc), mar = c(4, 4, 2.5, 1))

  bands <- list()

  for (j in seq_len(p)) {
    u <- U[, j]
    u <- u[is.finite(u)]
    N <- length(u)
    key <- as.character(N)

    if (is.null(bands[[key]])) {
      gs <- if (is.null(grid_size)) min(N, 400) else grid_size
      bands[[key]] <- make_ecdf_simultaneous_band(
        N = N, level = level, n_sim = n_sim,
        grid_size = gs, seed = seed
      )
    }
    band <- bands[[key]]

    Fhat <- stats::ecdf(u)(band$z)
    diff_obs <- Fhat - band$z
    ylim_use <- range(c(band$lower_diff, band$upper_diff, diff_obs), finite = TRUE)
    pad <- 0.05 * diff(ylim_use)
    if (!is.finite(pad) || pad <= 0) pad <- 0.01
    ylim_use <- ylim_use + c(-pad, pad)

    plot(
      NA,
      xlim = c(0, 1), ylim = ylim_use,
      main = paste(main_prefix, colnames(U)[j]),
      cex.main = cex.main,
      xlab = "u", ylab = expression(hat(F)(u) - u)
    )
    polygon(
      c(band$z, rev(band$z)),
      c(band$lower_diff, rev(band$upper_diff)),
      col = band_col, border = NA
    )
    abline(h = 0, lty = 2)
    lines(band$z, diff_obs, type = "s", lwd = 1.2)
  }

  invisible(bands)
}



# ------------------------------------------------------------
# 9. Helpers for the refitted RatioIS localized-SBC architecture
# ------------------------------------------------------------

.make_gllim_init_param_local <- function(fit_obj) {
  # IMPORTANT: prune_fit_components() deliberately drops singleton component
  # dimensions when K = 1.  gllim_strict_em(), however, always indexes the
  # initialization as b[, k], A[, , k], Sigma[, , k], etc.  Therefore a
  # one-component warm start must be expanded back to matrix/array form here.
  # Do not return fit_obj$optpar before this normalization.

  K <- as.integer(fit_obj$K)
  if (length(K) != 1L || !is.finite(K) || K < 1L) {
    stop("Invalid K in init_gllim_fit.")
  }

  src <- fit_obj$optpar
  if (is.null(src)) {
    src <- list(
      A = fit_obj$A,
      b = fit_obj$b,
      Sigma = fit_obj$Sigma,
      c = fit_obj$c,
      Gamma = fit_obj$Gamma,
      pi = fit_obj$prior_mix_prob
    )
  }

  A <- src$A
  b <- src$b
  Sigma <- src$Sigma
  cvec <- src$c
  Gamma <- src$Gamma
  pi <- as.numeric(src$pi)

  if (K == 1L) {
    # Infer the response/covariate dimensions from b and c, whose singleton
    # component dimension may have been dropped by prune_fit_components().
    D <- length(as.numeric(b))
    L <- length(as.numeric(cvec))

    b <- matrix(as.numeric(b), nrow = D, ncol = 1L)
    cvec <- matrix(as.numeric(cvec), nrow = L, ncol = 1L)
    A <- array(as.numeric(A), dim = c(D, L, 1L))
    Sigma <- array(as.numeric(Sigma), dim = c(D, D, 1L))
    Gamma <- array(as.numeric(Gamma), dim = c(L, L, 1L))
    pi <- 1
  }

  list(A = A, b = b, Sigma = Sigma, c = cvec, Gamma = Gamma, pi = pi)
}

fit_local_gllim_ratio_sbc <- function(theta_train,
                                      y_train,
                                      init_gllim_fit = NULL,
                                      K_used = NULL,
                                      warm_start = TRUE,
                                      maxiter = 300,
                                      verbose = 1,
                                      cov_structure = "",
                                      mixprob_thresh = 0) {
  .require_local_sbc_function("fit_gllim")

  theta_train <- as.matrix(theta_train)
  y_train <- as.matrix(y_train)
  if (ncol(theta_train) != ncol(y_train)) {
    stop("theta_train and y_train must have the same number of columns.")
  }

  if (is.null(K_used)) {
    if (!is.null(init_gllim_fit$K)) {
      K_used <- init_gllim_fit$K
    } else {
      stop("K_used must be supplied when init_gllim_fit is NULL.")
    }
  }

  init_param <- NULL
  if (isTRUE(warm_start)) {
    if (is.null(init_gllim_fit)) {
      stop("warm_start=TRUE requires init_gllim_fit.")
    }
    init_param <- .make_gllim_init_param_local(init_gllim_fit)
  }

  if (verbose >= 1) {
    cat("Fitting one local GLLiM for refitted RatioIS local SBC\n")
    cat("N_local_train =", ncol(theta_train), "\n")
    cat("K_used        =", K_used, "\n")
    cat("warm_start    =", isTRUE(warm_start), "\n")
  }

  fit <- fit_gllim(
    responses = theta_train,
    covariates = y_train,
    K = K_used,
    maxiter = maxiter,
    verbose = verbose,
    cov_structure = cov_structure,
    mixprob_thresh = mixprob_thresh,
    init_param = init_param
  )

  fit$num_param <- nrow(theta_train)
  fit$dim_data <- nrow(y_train)
  fit
}

fit_local_summary_pipeline <- function(local_fit,
                                       y_train,
                                       num_param,
                                       summary_config = list()) {
  .require_local_sbc_function("prepare_gllim_posterior_cache")
  .require_local_sbc_function("summarize_y_matrix_theta_from_cache")

  include_covariance <- summary_config$include_covariance_in_S1 %||% TRUE
  include_mixture_weights <- summary_config$include_mixture_weights_in_S1 %||% TRUE
  scale_floor <- summary_config$summary_scale_floor %||% 1e-4

  y_train <- as.matrix(y_train)
  cache <- prepare_gllim_posterior_cache(
    fit_obj = local_fit,
    num_param = num_param,
    dim_data = nrow(y_train)
  )

  S_raw <- summarize_y_matrix_theta_from_cache(
    y_mat = y_train,
    cache = cache,
    num_theta = num_param,
    include_covariance = include_covariance,
    include_mixture_weights = include_mixture_weights
  )

  center <- rowMeans(S_raw)
  scale <- apply(S_raw, 1, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  scale <- pmax(scale, scale_floor)

  S_std <- sweep(S_raw, 1, center, FUN = "-")
  S_std <- sweep(S_std, 1, scale, FUN = "/")

  S_train <- S_std

  storage.mode(S_train) <- "double"

  list(
    cache = cache,
    S_train_raw = S_raw,
    S_train = S_train,
    center = center,
    scale = scale,
    include_covariance = include_covariance,
    include_mixture_weights = include_mixture_weights,
    num_param = num_param,
    dim_summary = nrow(S_train)
  )
}

transform_local_summary_matrix <- function(y_mat, summary_pipeline) {
  .require_local_sbc_function("summarize_y_matrix_theta_from_cache")

  y_mat <- as.matrix(y_mat)
  S_raw <- summarize_y_matrix_theta_from_cache(
    y_mat = y_mat,
    cache = summary_pipeline$cache,
    num_theta = summary_pipeline$num_param,
    include_covariance = summary_pipeline$include_covariance,
    include_mixture_weights = summary_pipeline$include_mixture_weights
  )

  S_std <- sweep(S_raw, 1, summary_pipeline$center, FUN = "-")
  S_std <- sweep(S_std, 1, summary_pipeline$scale, FUN = "/")

  S <- S_std

  storage.mode(S) <- "double"
  S
}

local_refitted_ratio_pit <- function(theta_star,
                                     T_star,
                                     g_mix,
                                     local_ratio_fit,
                                     prior_pdf = NULL,
                                     M_is = 5000,
                                     ratio_lambda_choice = "lambda.1se") {
  .require_local_sbc_function("predict_log_ratio_batch")

  prop <- sample_valid_from_gmm(
    n = M_is,
    g_mix = g_mix,
    prior_pdf = prior_pdf,
    verbose = FALSE
  )
  theta_prop <- prop$theta

  log_r <- predict_log_ratio_batch(
    ratio_fit = local_ratio_fit,
    theta_mat = t(theta_prop),
    y_used = as.numeric(T_star),
    s = ratio_lambda_choice
  )

  w_raw <- normalize_log_weights(log_r)
  # Same deterministic clipping rule as Eq. (22)/(23) in the main RatioIS.
  w_clip <- clip_ratiois_weights(w_raw)

  finite_lr <- log_r[is.finite(log_r)]

  list(
    U_raw = weighted_empirical_pit(theta_star, theta_prop, w_raw),
    U_clipped = weighted_empirical_pit(theta_star, theta_prop, w_clip),
    ESS_raw = 1 / sum(w_raw^2),
    ESS_clipped = 1 / sum(w_clip^2),
    max_weight_raw = max(w_raw),
    max_weight_clipped = max(w_clip),
    proposal_acceptance_rate = prop$acceptance_rate,
    proposal_attempts = prop$attempts,
    log_ratio_sd = if (length(finite_lr) > 1) stats::sd(finite_lr) else NA_real_,
    log_ratio_min = if (length(finite_lr) > 0) min(finite_lr) else NA_real_,
    log_ratio_max = if (length(finite_lr) > 0) max(finite_lr) else NA_real_
  )
}

# ------------------------------------------------------------
# 10. Main refit-once RatioIS localized SBC driver
# ------------------------------------------------------------

run_local_sbc <- function(surrogate_out,
                          ratio_out,
                          model,
                          model_param = NULL,
                          prior_pdf = NULL,
                          observedData = NULL,
                          N_local = 2000,
                          N_local_check = N_local,
                          M_is = 5000,
                          lambda = 1,
                          ratio_lambda_choice = NULL,
                          weight_method = c("clipped", "raw"),
                          ratio_alpha = NULL,
                          ratio_feature_args = NULL,
                          ratio_use_cv = NULL,
                          ratio_nfolds = NULL,
                          K_used = NULL,
                          postprob_eps = 1e-6,
                          min_components = 1,
                          warm_start = TRUE,
                          maxiter = 300,
                          cov_structure = "",
                          mixprob_thresh = 0,
                          param_names = NULL,
                          output_dir = NULL,
                          cex.main = 0.7,
                          verbose = TRUE) {
  weight_method <- match.arg(weight_method)

  if (is.null(surrogate_out) || is.null(surrogate_out$fit1_raw)) {
    stop("surrogate_out with the already fitted final q2 is required.")
  }
  if (is.null(ratio_out) || is.null(ratio_out$ratio_fit)) {
    stop("ratio_out is required so the local classifier can inherit the main feature specification.")
  }

  if (length(N_local) != 1 || !is.finite(N_local) || N_local < 10 ||
      N_local != as.integer(N_local)) {
    stop("N_local must be an integer >= 10.")
  }
  if (length(N_local_check) != 1 || !is.finite(N_local_check) || N_local_check < 10 ||
      N_local_check != as.integer(N_local_check)) {
    stop("N_local_check must be an integer >= 10.")
  }
  if (length(M_is) != 1 || !is.finite(M_is) || M_is < 50 ||
      M_is != as.integer(M_is)) {
    stop("M_is must be an integer >= 50.")
  }
  N_local <- as.integer(N_local)
  N_local_check <- as.integer(N_local_check)
  M_is <- as.integer(M_is)

  if (is.null(observedData)) {
    observedData <- ratio_out$observedData %||% surrogate_out$observedData
  }
  if (is.null(observedData)) stop("observedData is required.")

  if (is.null(prior_pdf)) prior_pdf <- surrogate_out$prior_pdf
  if (is.null(prior_pdf)) {
    warning("prior_pdf is NULL: g_lambda draws will not be restricted to original prior support.")
  }

  num_param <- surrogate_out$num_param
  if (is.null(num_param)) stop("surrogate_out$num_param is required.")

  if (is.null(param_names)) {
    param_names <- rownames(ratio_out$theta_prop)
  }
  if (is.null(param_names)) param_names <- paste0("theta", seq_len(num_param))

  ratio_cfg <- ratio_out$ratio_config %||% list()
  if (is.null(ratio_lambda_choice)) ratio_lambda_choice <- ratio_cfg$lambda_choice %||% "lambda.1se"
  if (is.null(ratio_alpha)) ratio_alpha <- ratio_cfg$alpha %||% 1
  if (is.null(ratio_use_cv)) ratio_use_cv <- ratio_cfg$use_cv %||% TRUE
  if (is.null(ratio_nfolds)) ratio_nfolds <- ratio_cfg$nfolds %||% 10

  if (!is.finite(postprob_eps) || postprob_eps < 0) {
    stop("postprob_eps must be finite and non-negative.")
  }
  if (!is.finite(min_components) || min_components < 1) {
    stop("min_components must be >= 1.")
  }
  min_components <- as.integer(min_components)

  if (is.null(ratio_feature_args)) {
    ratio_feature_args <- ratio_out$ratio_fit$feature_args
  }
  if (is.null(ratio_feature_args)) {
    ratio_feature_args <- list(
      include_linear = TRUE,
      include_quadratic = FALSE,
      include_within_block_interactions = FALSE,
      include_cross_interactions = TRUE
    )
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Current manuscript Eq. (25): local design obtained by broadening q2.
  g_mix <- build_g_lambda_from_q2(
    surrogate_out = surrogate_out,
    ratio_out = ratio_out,
    observedData = observedData,
    lambda = lambda,
    param_names = param_names
  )

  # Match the older local-SBC warm-start convention: use only q2 components
  # with non-negligible posterior responsibility at y_obs unless K_used is
  # explicitly supplied. This affects the local GLLiM initialization/size,
  # not the definition of g_lambda itself.
  fit_K <- surrogate_out$fit1_raw$K
  if (g_mix$K == fit_K) {
    if (is.null(K_used)) {
      active_components <- which(g_mix$probs > postprob_eps)
      if (length(active_components) < min_components) {
        active_components <- order(g_mix$probs, decreasing = TRUE)[seq_len(min(min_components, fit_K))]
      }
      active_components <- sort(unique(active_components))
      K_used <- length(active_components)
    } else {
      K_used <- as.integer(K_used)
      if (K_used < 1 || K_used > fit_K) {
        stop("K_used must lie between 1 and the number of q2 components.")
      }
      active_components <- sort(order(g_mix$probs, decreasing = TRUE)[seq_len(K_used)])
    }
  } else {
    if (is.null(K_used)) K_used <- fit_K
    active_components <- seq_len(fit_K)
    warning("g_lambda component count differs from q2 fit; using all q2 components for initialization.")
  }

  init_fit <- surrogate_out$fit1_raw
  if (length(active_components) < fit_K) {
    .require_local_sbc_function("prune_fit_components")
    init_fit <- prune_fit_components(
      fit_obj = init_fit,
      active = active_components,
      stage_name = "local_sbc_ratio_refit_init"
    )
  }

  if (verbose) {
    cat("\nRefitted RatioIS localized SBC\n")
    cat("Design base             = q2(theta | y_obs)\n")
    cat("lambda                  =", lambda, "(SD inflation factor)\n")
    cat("N_local training        =", N_local, "\n")
    cat("N_local checking        =", N_local_check, "\n")
    cat("M_is per pseudo-data    =", M_is, "\n")
    cat("Local GLLiM refits      = 1\n")
    cat("Local ratio refits      = 1\n")
    cat("Inner IS proposal       = g_lambda(theta | y_obs)\n")
    cat("Primary PIT weights     =", weight_method, "\n")
    cat("Classifier alpha        =", ratio_alpha, "\n")
    cat("Active q2 component idx =", paste(active_components, collapse = ", "), "\n")
    cat("Number active q2 comps  =", length(active_components), "\n")
    cat("Local GLLiM K           =", K_used, "\n")
  }

  # ----------------------------------------------------------
  # A. Independent local TRAINING bank
  # ----------------------------------------------------------
  if (verbose) cat("\nSampling local training theta from g_lambda\n")
  tr_draw <- sample_valid_from_gmm(
    n = N_local,
    g_mix = g_mix,
    prior_pdf = prior_pdf,
    verbose = verbose
  )
  theta_train <- tr_draw$theta

  if (verbose) cat("Simulating local training data from p(y | theta)\n")
  y_train <- simulate_local_bank(
    theta_mat = theta_train,
    model = model,
    model_param = model_param,
    verbose = verbose
  )


  local_fit <- fit_local_gllim_ratio_sbc(
    theta_train = theta_train,
    y_train = y_train,
    init_gllim_fit = init_fit,
    K_used = K_used,
    warm_start = warm_start,
    maxiter = maxiter,
    verbose = if (isTRUE(verbose)) 1 else 0,
    cov_structure = cov_structure,
    mixprob_thresh = mixprob_thresh
  )

  # Fit T_lambda on the local GLLiM and use ONLY the training bank to estimate
  # standardization. The checking bank remains untouched at this stage.
  summary_pipeline <- fit_local_summary_pipeline(
    local_fit = local_fit,
    y_train = y_train,
    num_param = num_param,
    summary_config = surrogate_out$summary_config %||% list()
  )
  S_train <- summary_pipeline$S_train

  .require_local_sbc_function("fit_ratio_classifier")
  if (verbose) {
    cat("Fitting one local ratio classifier on (theta_train, T_lambda(y_train))\n")
    cat("Local summary dimension =", nrow(S_train), "\n")
  }

  t_ratio0 <- Sys.time()
  local_ratio_fit <- fit_ratio_classifier(
    theta_mat = theta_train,
    y_mat = S_train,
    alpha = ratio_alpha,
    feature_args = ratio_feature_args,
    use_cv = ratio_use_cv,
    nfolds = ratio_nfolds
  )
  local_ratio_fit_time_sec <- as.numeric(Sys.time() - t_ratio0, units = "secs")

  # ----------------------------------------------------------
  # B. Independent local CHECKING bank
  # ----------------------------------------------------------
  if (verbose) cat("\nSampling independent local checking theta from g_lambda\n")
  ck_draw <- sample_valid_from_gmm(
    n = N_local_check,
    g_mix = g_mix,
    prior_pdf = prior_pdf,
    verbose = verbose
  )
  theta_check <- ck_draw$theta

  if (verbose) cat("Simulating independent local checking data from p(y | theta)\n")
  y_check <- simulate_local_bank(
    theta_mat = theta_check,
    model = model,
    model_param = model_param,
    verbose = verbose
  )

  # Apply the already frozen local summary map to all checking pseudo-data.
  S_check <- transform_local_summary_matrix(y_check, summary_pipeline)

  p <- nrow(theta_check)
  U_raw <- matrix(NA_real_, nrow = N_local_check, ncol = p)
  U_clipped <- matrix(NA_real_, nrow = N_local_check, ncol = p)
  colnames(U_raw) <- param_names
  colnames(U_clipped) <- param_names

  is_diag <- data.frame(
    replicate = seq_len(N_local_check),
    ESS_raw = NA_real_,
    ESS_clipped = NA_real_,
    max_weight_raw = NA_real_,
    max_weight_clipped = NA_real_,
    proposal_acceptance_rate = NA_real_,
    proposal_attempts = NA_integer_,
    log_ratio_sd = NA_real_,
    log_ratio_min = NA_real_,
    log_ratio_max = NA_real_
  )

  for (i in seq_len(N_local_check)) {
    if (verbose && (i %% 25 == 0 || i == N_local_check)) {
      cat("Ratio-weighted", i, "of", N_local_check, "checking pseudo-observations\n")
    }

    one <- local_refitted_ratio_pit(
      theta_star = theta_check[, i],
      T_star = S_check[, i],
      g_mix = g_mix,
      local_ratio_fit = local_ratio_fit,
      prior_pdf = prior_pdf,
      M_is = M_is,
      ratio_lambda_choice = ratio_lambda_choice
    )

    U_raw[i, ] <- one$U_raw
    U_clipped[i, ] <- one$U_clipped
    is_diag$ESS_raw[i] <- one$ESS_raw
    is_diag$ESS_clipped[i] <- one$ESS_clipped
    is_diag$max_weight_raw[i] <- one$max_weight_raw
    is_diag$max_weight_clipped[i] <- one$max_weight_clipped
    is_diag$proposal_acceptance_rate[i] <- one$proposal_acceptance_rate
    is_diag$proposal_attempts[i] <- one$proposal_attempts
    is_diag$log_ratio_sd[i] <- one$log_ratio_sd
    is_diag$log_ratio_min[i] <- one$log_ratio_min
    is_diag$log_ratio_max[i] <- one$log_ratio_max
  }

  # The clipped empirical posterior is primary because this matches the
  # finite-M posterior approximation actually reported by RatioIS (Eq. 23).
  U <- if (weight_method == "clipped") U_clipped else U_raw
  pit_summary <- summarise_pit(U)
  pit_summary_raw <- summarise_pit(U_raw)
  pit_summary_clipped <- summarise_pit(U_clipped)

  out <- list(
    U = U,
    U_raw = U_raw,
    U_clipped = U_clipped,
    pit_summary = pit_summary,
    pit_summary_raw = pit_summary_raw,
    pit_summary_clipped = pit_summary_clipped,
    theta_train = theta_train,
    y_train = y_train,
    S_train = S_train,
    theta_check = theta_check,
    y_check = y_check,
    S_check = S_check,
    # Compatibility aliases for plotting/post-processing code.
    theta_star = theta_check,
    y_star = y_check,
    theta_sim = theta_check,
    y_sim = y_check,
    g_mix = g_mix,
    local_fit = local_fit,
    local_ratio_fit = local_ratio_fit,
    local_summary_pipeline = summary_pipeline,
    is_diagnostics = is_diag,
    settings = list(
      implementation = "refit_once_ratiois_local_sbc",
      design_base = "q2(theta | y_obs)",
      target = "p_g_lambda(theta | T_lambda(y_star); y_obs)",
      inner_proposal = "g_lambda(theta | y_obs)",
      raw_weight = "rhat_lambda(theta, T_lambda(y_star))",
      clipping_rule = "Mc=floor(sqrt(M_is)); cap=Mc-th largest normalized raw weight",
      primary_posterior = "clipped Eq.-(23)-type weighted empirical measure",
      N_local = N_local,
      N_local_train = N_local,
      N_local_check = N_local_check,
      N_model_simulations = N_local + N_local_check,
      M_is = M_is,
      lambda = lambda,
      lambda_interpretation = "within-component SD inflation factor",
      ratio_lambda_choice = ratio_lambda_choice,
      weight_method = weight_method,
      ratio_alpha = ratio_alpha,
      ratio_feature_args = ratio_feature_args,
      ratio_use_cv = ratio_use_cv,
      ratio_nfolds = ratio_nfolds,
      K_used = K_used,
      active_components = active_components,
      postprob_eps = postprob_eps,
      warm_start = warm_start,
      local_gllim_refits = 1L,
      local_ratio_refits = 1L,
      local_ratio_fit_time_sec = local_ratio_fit_time_sec,
      param_names = param_names
    )
  )

  if (!is.null(output_dir)) {
    # Save numerical output before graphics so a PDF-device problem never
    # forces the expensive simulator calls to be repeated.
    saveRDS(out, file.path(output_dir, "local_sbc_out.rds"))

    write.csv(U, file.path(output_dir, "local_sbc_pit_values.csv"), row.names = FALSE)
    write.csv(U_raw, file.path(output_dir, "local_sbc_pit_values_raw.csv"), row.names = FALSE)
    write.csv(U_clipped, file.path(output_dir, "local_sbc_pit_values_clipped.csv"), row.names = FALSE)

    write.csv(pit_summary, file.path(output_dir, "local_sbc_pit_summary.csv"), row.names = FALSE)
    write.csv(pit_summary_raw, file.path(output_dir, "local_sbc_pit_summary_raw.csv"), row.names = FALSE)
    write.csv(pit_summary_clipped, file.path(output_dir, "local_sbc_pit_summary_clipped.csv"), row.names = FALSE)
    write.csv(is_diag, file.path(output_dir, "local_sbc_is_diagnostics.csv"), row.names = FALSE)

    write.csv(t(theta_train), file.path(output_dir, "local_sbc_theta_train.csv"), row.names = FALSE)
    write.csv(t(y_train), file.path(output_dir, "local_sbc_y_train.csv"), row.names = FALSE)
    write.csv(t(theta_check), file.path(output_dir, "local_sbc_theta_check.csv"), row.names = FALSE)
    write.csv(t(y_check), file.path(output_dir, "local_sbc_y_check.csv"), row.names = FALSE)

    method_tag <- if (weight_method == "clipped") "clipped" else "raw"

    plot_pit_histograms(
      U,
      file = file.path(output_dir, "sbc_hist.pdf"),
      breaks = 20,
      main_prefix = paste0("Refit RatioIS local SBC (", method_tag, "), lambda=", lambda, ":"),
      cex.main = cex.main
    )

    plot_pit_ecdf(
      U,
      file = file.path(output_dir, "sbc_ecdf.pdf"),
      main_prefix = paste0("Refit RatioIS local SBC ECDF (", method_tag, "), lambda=", lambda, ":"),
      cex.main = cex.main
    )

    plot_pit_ecdf_simultaneous(
      U,
      file = file.path(output_dir, "sbc_ecdf_sim.pdf"),
      main_prefix = paste0("ECDF, lambda=", lambda, ":"),
      level = 0.95,
      cex.main = cex.main
    )

    plot_pit_ecdf_difference_simultaneous(
      U,
      file = file.path(output_dir, "sbc_ecdf_diff.pdf"),
      main_prefix = paste0("ECDF difference, lambda=", lambda, ":"),
      level = 0.95,
      cex.main = cex.main
    )
  }

  out
}
