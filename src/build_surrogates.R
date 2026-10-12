# ============================================================
# Two-round localized GLLiM surrogate construction
# ------------------------------------------------------------
# Workflow:
#   prior predictive data -> fit q0(theta | y)
#   -> prune q0 components by posterior responsibility at y_obs
#   -> sample first focused design and fit q1(theta | y)
#   -> sample refreshed design and fit final q2(theta | y)
#   -> construct posterior-based summaries S1(y) from q2
#   -> standardize summaries for ratio estimation.
#
# `mixprob_thresh` is an optional numerical cleanup on fitted global
# mixture proportions. `postprob_thresh` acts on conditional component
# probabilities at the observed data. Posterior conditioning is stabilized
# with a scale-adaptive ridge.
# ============================================================

posterior_mixture_from_input <- function(cache, x) {
  post_mix_prob <- compute_post_mix_prob(
    A = cache$A, b = cache$b, c = cache$c,
    Sigma = cache$Sigma, Gamma = cache$Gamma, pi = cache$pi,
    K = cache$K, observedData = x
  )
  mean_proposal <- compute_mean_proposal(
    A_star = cache$A_star, observedData = x, b_star = cache$b_star,
    num_param = cache$num_param, K = cache$K
  )
  means <- as.matrix(mean_proposal)
  if (cache$K == 1) means <- matrix(as.numeric(mean_proposal), ncol = 1)
  list(probs = as.numeric(post_mix_prob), means = means, covars = cache$Sigma_star)
}

posterior_theta_mixture_from_cache <- function(cache, y, num_theta = cache$num_param) {
  gmm <- posterior_mixture_from_input(cache, y)
  if (nrow(gmm$means) != num_theta) {
    stop("Posterior mixture dimension does not match num_theta.")
  }
  gmm
}

posterior_theta_summary_from_cache <- function(y, cache, num_theta,
                                               include_covariance = TRUE,
                                               include_mixture_weights = TRUE) {
  gmm_theta <- posterior_theta_mixture_from_cache(cache, y, num_theta)
  probs <- gmm_theta$probs
  means <- gmm_theta$means
  covars <- gmm_theta$covars
  K <- length(probs)
  post_mean <- if (K == 1) as.numeric(means[, 1]) else as.numeric(means %*% probs)
  out <- post_mean
  if (include_covariance) {
    if (K == 1) {
      post_covar <- as.matrix(covars)
    } else {
      post_covar <- matrix(0, nrow = num_theta, ncol = num_theta)
      for (k in seq_len(K)) {
        mk <- means[, k]
        Sk <- covars[, , k]
        diffk <- matrix(mk - post_mean, ncol = 1)
        post_covar <- post_covar + probs[k] * (Sk + diffk %*% t(diffk))
      }
    }
    out <- c(out, vech(post_covar))
  }
  if (include_mixture_weights) out <- c(out, as.numeric(probs))
  out
}

summarize_y_matrix_theta_from_cache <- function(y_mat, cache, num_theta,
                                                include_covariance = TRUE,
                                                include_mixture_weights = TRUE) {
  y_mat <- as.matrix(y_mat)
  out <- do.call(cbind, lapply(seq_len(ncol(y_mat)), function(i) {
    posterior_theta_summary_from_cache(
      y = y_mat[, i], cache = cache, num_theta = num_theta,
      include_covariance = include_covariance,
      include_mixture_weights = include_mixture_weights
    )
  }))
  storage.mode(out) <- "double"
  out
}

# ::::::::::::::::::::: FUNCTION DEFINITIONS ::::::::::::::::::::::::::


canonicalize_gmm_components <- function(gmm_obj) {
  probs <- as.numeric(gmm_obj$probs)
  means <- as.matrix(gmm_obj$means)
  covars <- gmm_obj$covars
  
  K <- length(probs)
  d <- nrow(means)
  
  if (K <= 1) {
    return(list(probs = probs, means = means, covars = covars))
  }
  
  if (length(dim(covars)) == 3) {
    var_score <- vapply(seq_len(K), function(k) sum(diag(covars[, , k])), numeric(1))
  } else {
    var_score <- rep(sum(diag(covars)), K)
  }
  
  if (d >= 2) {
    ord <- order(-probs, means[1, ], means[2, ], var_score)
  } else {
    ord <- order(-probs, means[1, ], var_score)
  }
  
  list(
    probs = probs[ord],
    means = means[, ord, drop = FALSE],
    covars = if (length(dim(covars)) == 3) covars[, , ord, drop = FALSE] else covars
  )
}

select_postprob_components <- function(post_mix_prob, postprob_thresh = 0, min_components = 1) {
  post_mix_prob <- as.numeric(post_mix_prob)
  K <- length(post_mix_prob)

  if (K < 1) {
    stop("select_postprob_components: empty posterior probability vector.")
  }

  if (is.null(postprob_thresh) || !is.finite(postprob_thresh)) {
    postprob_thresh <- 0
  }

  # Exact no-pruning convention used in the paper code: when the
  # threshold is zero (or negative), retain every fitted component, even if
  # its numerically evaluated responsibility at y_obs is exactly zero.
  if (postprob_thresh <= 0) {
    return(seq_len(K))
  }

  active <- which(post_mix_prob > postprob_thresh)

  if (length(active) < min_components) {
    ord <- order(post_mix_prob, decreasing = TRUE)
    active <- ord[seq_len(min(min_components, K))]
  }

  sort(unique(active))
}

prune_fit_components <- function(fit_obj, active, stage_name = "fit") {
  K_old <- fit_obj$K
  active <- sort(unique(as.integer(active)))

  if (length(active) < 1) {
    stop("prune_fit_components: no active components selected for ", stage_name, ".")
  }
  if (any(active < 1) || any(active > K_old)) {
    stop("prune_fit_components: active component index out of range for ", stage_name, ".")
  }

  if (length(active) == K_old) {
    fit_obj$K_before_postprob_pruning <- K_old
    fit_obj$postprob_active_components <- active
    return(fit_obj)
  }

  fit_obj$K_before_postprob_pruning <- K_old
  fit_obj$postprob_active_components <- active
  fit_obj$K <- length(active)

  if (K_old == 1) {
    # Nothing to subset; keep original matrix representation.
  } else {
    fit_obj$A <- fit_obj$A[, , active, drop = FALSE]
    fit_obj$b <- fit_obj$b[, active, drop = FALSE]
    fit_obj$Sigma <- fit_obj$Sigma[, , active, drop = FALSE]
    fit_obj$c <- fit_obj$c[, active, drop = FALSE]
    fit_obj$Gamma <- fit_obj$Gamma[, , active, drop = FALSE]
    fit_obj$prior_mix_prob <- fit_obj$prior_mix_prob[active]
    fit_obj$prior_mix_prob <- fit_obj$prior_mix_prob / sum(fit_obj$prior_mix_prob)
  }

  if (fit_obj$K == 1 && K_old > 1) {
    fit_obj$A <- drop(fit_obj$A)
    fit_obj$b <- drop(fit_obj$b)
    fit_obj$Sigma <- drop(fit_obj$Sigma)
    fit_obj$c <- drop(fit_obj$c)
    fit_obj$Gamma <- drop(fit_obj$Gamma)
  }

  fit_obj$optpar <- list(
    A = fit_obj$A,
    b = fit_obj$b,
    Sigma = fit_obj$Sigma,
    c = fit_obj$c,
    Gamma = fit_obj$Gamma,
    pi = fit_obj$prior_mix_prob
  )

  fit_obj
}


posterior_prune_fit_at_yobs <- function(fit_obj, observedData, num_param, dim_data,
                                         postprob_thresh = 0, stage_name = "fit",
                                         ridge_post = 1e-6, min_components = 1,
                                         verbose = 1) {
  post_mix_prob_full <- compute_post_mix_prob(
    A = fit_obj$A, b = fit_obj$b, c = fit_obj$c, Sigma = fit_obj$Sigma,
    Gamma = fit_obj$Gamma, pi = fit_obj$prior_mix_prob, K = fit_obj$K,
    observedData = observedData
  )
  active <- select_postprob_components(post_mix_prob_full, postprob_thresh, min_components)
  fit_pruned <- prune_fit_components(fit_obj, active, stage_name)
  cache_pruned <- prepare_posterior_cache_from_fit(
    fit_obj = fit_pruned, num_param = num_param, dim_data = dim_data, ridge_post = ridge_post
  )
  gmm <- posterior_mixture_from_input(cache_pruned, observedData)
  if (verbose >= 1) {
    cat("Posterior-probability pruning for", stage_name, ": K original =",
        fit_obj$K, ", K used =", fit_pruned$K,
        ", threshold =", postprob_thresh, "\n")
  }
  list(
    fit = fit_pruned, cache = cache_pruned,
    post_mix_prob = gmm$probs, mean_proposal = gmm$means,
    covar_proposal = gmm$covars,
    post_mix_prob_full = post_mix_prob_full, active_components = active
  )
}

fit_gllim <- function(responses, covariates, K, maxiter, verbose,
                      cov_structure, mixprob_thresh, init_param) {
  mod <- gllim_strict_em( # the original gllim code from xLLiM modified to have a different termination criterion for EM
    responses,
    covariates,
    in_K = K,
    maxiter = maxiter,
    in_theta = init_param,
    verb = verbose,
    cstr = list(Sigma = cov_structure)
  )
  
  A <- mod$A
  b <- mod$b
  Sigma <- mod$Sigma
  c <- mod$c
  Gamma <- mod$Gamma
  prior_mix_prob <- mod$pi
  K <- length(prior_mix_prob)
  
  id_eliminate <- which(prior_mix_prob < mixprob_thresh)
  if (length(id_eliminate) > 0) {
    K <- K - length(id_eliminate)
    warning(
      "mixprob_thresh reduced the mixture number of components from ",
      K + length(id_eliminate),
      " to ",
      K
    )
    A <- A[, , -id_eliminate, drop = FALSE]
    b <- b[, -id_eliminate, drop = FALSE]
    Sigma <- Sigma[, , -id_eliminate, drop = FALSE]
    c <- c[, -id_eliminate, drop = FALSE]
    Gamma <- Gamma[, , -id_eliminate, drop = FALSE]
    prior_mix_prob <- prior_mix_prob[-id_eliminate]
  }
  
  if (K == 1) {
    A <- drop(A)
    b <- drop(b)
    Sigma <- drop(Sigma)
    c <- drop(c)
    Gamma <- drop(Gamma)
  }
  
  optpar <- list(
    A = A,
    b = b,
    Sigma = Sigma,
    c = c,
    Gamma = Gamma,
    pi = prior_mix_prob
  )
  
  list(
    K = K,
    A = A,
    b = b,
    Sigma = Sigma,
    c = c,
    Gamma = Gamma,
    prior_mix_prob = prior_mix_prob,
    optpar = optpar
  )
}

make_valid_covariance <- function(S, eps = 1e-8) {
  S <- as.matrix(S)
  S <- 0.5 * (S + t(S))
  ee <- eigen(S, symmetric = TRUE)
  vals <- ee$values
  vecs <- ee$vectors
  vals[!is.finite(vals)] <- eps
  vals[vals < eps] <- eps
  S_valid <- vecs %*% diag(vals, length(vals)) %*% t(vecs)
  S_valid <- 0.5 * (S_valid + t(S_valid))
  S_valid
}

mixpdf <- function(theta, probs, means, covars) {
  pdf <- 0
  numcomp <- length(probs)
  
  if (numcomp > 1) {
    for (k in seq_len(numcomp)) {
      Sigma_k <- make_valid_covariance(covars[, , k])
      pdf <- pdf + probs[k] * mvtnorm::dmvnorm(
        theta,
        mean = means[, k],
        sigma = Sigma_k
      )
    }
  } else {
    Sigma_1 <- make_valid_covariance(covars)
    pdf <- mvtnorm::dmvnorm(theta, mean = means, sigma = Sigma_1)
  }
  
  pdf
}

par_inflated_mix <- function(probs, means, covars, factor_cvMH = 1) {
  numcomp <- length(probs)
  k <- which(rmultinom(1, 1, probs) == 1)
  
  if (numcomp > 1) {
    Sigma_k <- make_valid_covariance(factor_cvMH * covars[, , k])
    out <- mvtnorm::rmvnorm(1, mean = means[, k], sigma = Sigma_k)
  } else {
    Sigma_1 <- make_valid_covariance(factor_cvMH * covars)
    out <- mvtnorm::rmvnorm(1, mean = means, sigma = Sigma_1)
  }
  
  as.vector(out)
}

# ------------------------------------------------------------
# Helper: sample only prior-valid theta draws from an unconstrained
# Gaussian-mixture proposal used during surrogate building.
# The distinct name avoids collisions with ratio-estimation helpers when
# several source files are loaded in the same R session.
# ------------------------------------------------------------
sample_valid_from_mixture_build <- function(n,
                                                num_param,
                                                probs,
                                                means,
                                                covars,
                                                factor_cvMH,
                                                prior_pdf,
                                                sample_prior = NULL,
                                                eps_defensive = 0,
                                                max_tries_factor = 50,
                                                verbose = TRUE) {
  theta_mat <- matrix(NA_real_, nrow = num_param, ncol = n)

  i <- 1
  n_try <- 0
  max_try <- max_tries_factor * n

  while (i <= n && n_try < max_try) {
    n_try <- n_try + 1

    if (!is.null(sample_prior) && eps_defensive > 0 && runif(1) < eps_defensive) {
      th <- sample_prior()
    } else {
      th <- par_inflated_mix(
        probs = probs,
        means = means,
        covars = covars,
        factor_cvMH = factor_cvMH
      )
    }

    p_th <- tryCatch(prior_pdf(th), error = function(e) NA_real_)

    ok <- length(th) == num_param &&
      all(is.finite(th)) &&
      is.finite(p_th) &&
      p_th > 0

    if (ok) {
      theta_mat[, i] <- th
      i <- i + 1
    }
  }

  if (i <= n) {
    stop(
      "Could not obtain enough valid proposal draws. ",
      "Obtained ", i - 1, " out of ", n, " after ", n_try, " attempts. ",
      "Consider reducing factor_cvMH or using a transformed parametrisation."
    )
  }

  if (isTRUE(verbose)) {
    cat(
      "Valid surrogate proposal draw acceptance rate:",
      round(n / n_try, 4),
      "(", n, "accepted out of", n_try, "attempts )\n"
    )
  }

  if (n / n_try < 0.05) {
    warning(
      "Low valid proposal acceptance rate: ", round(n / n_try, 4),
      ". Consider reducing factor_cvMH or using a transformed parametrisation."
    )
  }

  theta_mat
}

safe_solve <- function(M, ridge = 1e-6, max_tries = 6) {
  M <- as.matrix(M)
  M <- 0.5 * (M + t(M))
  d <- nrow(M)

  scale_M <- mean(abs(diag(M)))
  if (!is.finite(scale_M) || scale_M <= 0) {
    scale_M <- mean(abs(M))
  }
  if (!is.finite(scale_M) || scale_M <= 0) {
    scale_M <- 1
  }

  for (j in seq_len(max_tries)) {
    ridge_abs <- ridge * scale_M * 10^(j - 1)
    out <- try(solve(M + diag(ridge_abs, d)), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
  }

  stop("safe_solve failed even after adaptive ridge up to ",
       ridge * scale_M * 10^(max_tries - 1))
}

compute_post_param <- function(A, b, c, Sigma, Gamma, num_param, K, dim_data,
                               ridge_post = 1e-6) {
  A_star <- array(NA, dim = c(num_param, dim_data, K))
  b_star <- array(NA, dim = c(num_param, K))
  Sigma_star <- array(NA, dim = c(num_param, num_param, K))

  for (k in seq_len(K)) {
    if (K == 1) {
      Gamma_inv <- safe_solve(Gamma, ridge = ridge_post)
      Sigma_inv <- safe_solve(Sigma, ridge = ridge_post)
      post_prec <- Gamma_inv + t(A) %*% Sigma_inv %*% A
      post_prec <- 0.5 * (post_prec + t(post_prec))

      Sigma_star <- safe_solve(post_prec, ridge = ridge_post)
      Sigma_star <- make_valid_covariance(Sigma_star)
      A_star <- Sigma_star %*% t(A) %*% Sigma_inv
      b_star <- Sigma_star %*% (Gamma_inv %*% c - t(A) %*% Sigma_inv %*% b)
    } else {
      Gamma_inv_k <- safe_solve(Gamma[, , k], ridge = ridge_post)
      Sigma_inv_k <- safe_solve(Sigma[, , k], ridge = ridge_post)
      post_prec_k <- Gamma_inv_k + t(A[, , k]) %*% Sigma_inv_k %*% A[, , k]
      post_prec_k <- 0.5 * (post_prec_k + t(post_prec_k))

      Sigma_star_k <- safe_solve(post_prec_k, ridge = ridge_post)
      Sigma_star_k <- make_valid_covariance(Sigma_star_k)
      Sigma_star[, , k] <- Sigma_star_k
      A_star[, , k] <- Sigma_star_k %*% t(A[, , k]) %*% Sigma_inv_k
      b_star[, k] <- Sigma_star_k %*% (
        Gamma_inv_k %*% c[, k] -
          t(A[, , k]) %*% Sigma_inv_k %*% b[, k]
      )
    }
  }

  list(A_star = A_star, b_star = b_star, Sigma_star = Sigma_star)
}

compute_post_mix_prob <- function(A, b, c, Sigma, Gamma, pi, K, observedData,
                                  cov_eps = 1e-8) {
  if (is.null(K) || length(K) != 1 || !is.finite(K)) {
    stop("compute_post_mix_prob: K is missing or invalid.")
  }
  K <- as.integer(K)
  if (K < 1) {
    stop("compute_post_mix_prob: K must be >= 1.")
  }
  if (is.null(pi) || length(pi) == 0) {
    stop("compute_post_mix_prob: pi is missing or length zero.")
  }
  pi <- as.numeric(pi)
  if (length(pi) != K) {
    stop("compute_post_mix_prob: length(pi) does not match K.")
  }

  observedData <- as.numeric(observedData)

  if (K > 1) {
    if (length(dim(c)) != 2 || ncol(c) != K) {
      stop("compute_post_mix_prob: c does not have expected dimension latent_dim x K.")
    }
    if (length(dim(b)) != 2 || ncol(b) != K) {
      stop("compute_post_mix_prob: b does not have expected dimension dim_data x K.")
    }
    if (length(dim(A)) != 3 || dim(A)[3] != K) {
      stop("compute_post_mix_prob: A does not have expected dimension dim_data x latent_dim x K.")
    }
    if (length(dim(Sigma)) != 3 || dim(Sigma)[3] != K) {
      stop("compute_post_mix_prob: Sigma does not have expected dimension dim_data x dim_data x K.")
    }
    if (length(dim(Gamma)) != 3 || dim(Gamma)[3] != K) {
      stop("compute_post_mix_prob: Gamma does not have expected dimension latent_dim x latent_dim x K.")
    }
  }

  log_post_mix_prob <- rep(-Inf, K)

  for (k in seq_len(K)) {
    if (!is.finite(pi[k]) || pi[k] <= 0) {
      log_post_mix_prob[k] <- -Inf
      next
    }

    if (K == 1) {
      c_star_k <- as.numeric(A %*% c + b)
      Gamma_star_k <- Sigma + A %*% Gamma %*% t(A)
    } else {
      c_star_k <- as.numeric(A[, , k] %*% c[, k] + b[, k])
      Gamma_star_k <- Sigma[, , k] + A[, , k] %*% Gamma[, , k] %*% t(A[, , k])
    }

    Gamma_star_k <- make_valid_covariance(Gamma_star_k, eps = cov_eps)

    log_post_mix_prob[k] <- log(pi[k]) + mvtnorm::dmvnorm(
      observedData,
      mean = c_star_k,
      sigma = Gamma_star_k,
      log = TRUE
    )
  }

  if (all(!is.finite(log_post_mix_prob))) {
    stop("compute_post_mix_prob: all log component probabilities are non-finite.")
  }

  m <- max(log_post_mix_prob)
  post_mix_prob <- exp(log_post_mix_prob - m)
  sum_probs <- sum(post_mix_prob)

  if (!is.finite(sum_probs) || sum_probs <= 0) {
    stop("compute_post_mix_prob: posterior mixture probabilities could not be normalized.")
  }

  post_mix_prob / sum_probs
}

compute_mean_proposal <- function(A_star, observedData, b_star, num_param, K) {
  mean_proposal <- array(NA, dim = c(num_param, K))
  
  for (k in seq_len(K)) {
    if (K == 1) {
      mean_proposal <- A_star %*% observedData + b_star
    } else {
      mean_proposal[, k] <- A_star[, , k] %*% observedData + b_star[, k]
    }
  }
  
  mean_proposal
}

mean_surrogate_likelihood <- function(A, b, observedData, theta, K) {
  out <- array(NA, dim = c(length(observedData), K))
  
  for (k in seq_len(K)) {
    if (K == 1) {
      out <- A %*% theta + b
    } else {
      out[, k] <- A[, , k] %*% theta + b[, k]
    }
  }
  
  out
}

likelihood_mix_prob <- function(c, Gamma, pi, K, theta) {
  out <- vector(length = K)
  
  for (k in seq_len(K)) {
    if (K == 1) {
      out <- pi * mvtnorm::dmvnorm(theta, mean = c, sigma = Gamma)
    } else {
      out[k] <- pi[k] * mvtnorm::dmvnorm(theta, mean = c[, k], sigma = Gamma[, , k])
    }
  }
  
  out / sum(out)
}

vech <- function(S) {
  S[lower.tri(S, diag = TRUE)]
}

prepare_posterior_cache_from_fit <- function(fit_obj, num_param, dim_data,
                                             ridge_post = 1e-6) {
  post_param <- compute_post_param(
    A = fit_obj$A,
    b = fit_obj$b,
    c = fit_obj$c,
    Sigma = fit_obj$Sigma,
    Gamma = fit_obj$Gamma,
    num_param = num_param,
    K = fit_obj$K,
    dim_data = dim_data,
    ridge_post = ridge_post
  )
  
  list(
    K = fit_obj$K,
    num_param = num_param,
    dim_data = dim_data,
    A = fit_obj$A,
    b = fit_obj$b,
    c = fit_obj$c,
    Sigma = fit_obj$Sigma,
    Gamma = fit_obj$Gamma,
    pi = fit_obj$prior_mix_prob,
    A_star = post_param$A_star,
    b_star = post_param$b_star,
    Sigma_star = post_param$Sigma_star
  )
}

posterior_summary_from_cache <- function(y,
                                         cache,
                                         include_covariance = TRUE,
                                         include_mixture_weights = TRUE) {
  post_mix_prob <- compute_post_mix_prob(
    A = cache$A,
    b = cache$b,
    c = cache$c,
    Sigma = cache$Sigma,
    Gamma = cache$Gamma,
    pi = cache$pi,
    K = cache$K,
    observedData = y
  )

  mean_proposal <- compute_mean_proposal(
    A_star = cache$A_star,
    observedData = y,
    b_star = cache$b_star,
    num_param = cache$num_param,
    K = cache$K
  )

  if (cache$K > 1) {
    gmm_tmp <- canonicalize_gmm_components(list(
      probs = post_mix_prob,
      means = mean_proposal,
      covars = cache$Sigma_star
    ))
    post_mix_prob <- gmm_tmp$probs
    mean_proposal <- gmm_tmp$means
    Sigma_star_use <- gmm_tmp$covars
  } else {
    mean_proposal <- matrix(as.numeric(mean_proposal), nrow = cache$num_param, ncol = 1)
    post_mix_prob <- as.numeric(post_mix_prob)
    Sigma_star_use <- array(as.matrix(cache$Sigma_star),
                            dim = c(cache$num_param, cache$num_param, 1))
  }

  out <- numeric(0)

  if (include_mixture_weights) {
    out <- c(out, as.numeric(post_mix_prob))
  }

  out <- c(out, as.numeric(mean_proposal))

  if (include_covariance) {
    log_diag_vars <- matrix(NA_real_, nrow = cache$num_param, ncol = cache$K)
    for (k in seq_len(cache$K)) {
      diag_k <- diag(Sigma_star_use[, , k])
      diag_k <- pmax(diag_k, 1e-12)
      log_diag_vars[, k] <- log(diag_k)
    }
    out <- c(out, as.numeric(log_diag_vars))
  }

  out
}


summarize_y_matrix_from_cache <- function(y_mat,
                                          cache,
                                          include_covariance = TRUE,
                                          include_mixture_weights = TRUE) {
  y_mat <- as.matrix(y_mat)
  
  out_list <- lapply(seq_len(ncol(y_mat)), function(i) {
    posterior_summary_from_cache(
      y = y_mat[, i],
      cache = cache,
      include_covariance = include_covariance,
      include_mixture_weights = include_mixture_weights
    )
  })
  
  out <- do.call(cbind, out_list)
  storage.mode(out) <- "double"
  out
}

standardize_summary_matrix <- function(S_train, S_obs, scale_floor = 1e-4) {
  S_train <- as.matrix(S_train)
  S_obs <- as.numeric(S_obs)
  
  center <- rowMeans(S_train)
  scale <- apply(S_train, 1, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  scale <- pmax(scale, scale_floor)
  
  S_train_std <- sweep(S_train, 1, center, FUN = "-")
  S_train_std <- sweep(S_train_std, 1, scale, FUN = "/")
  S_obs_std <- (S_obs - center) / scale
  
  list(
    S_train_std = S_train_std,
    S_obs_std = S_obs_std,
    center = center,
    scale = scale
  )
}

# ::::::::::::::::::::: MAIN FUNCTION :::::::::::::::::::::::::::::::::


# Hard no-pruning convention: postprob_thresh <= 0 keeps all fitted components.
# ::::::::::::::::::::: MAIN FUNCTION :::::::::::::::::::::::::::::::::

build_surrogates <- function(
    observedData,
    K_start = 10,
    cov_structure = "",
    maxiter = 300,
    prior_pdf,
    sample_prior,
    model,
    model_param = NULL,
    dim_data,
    num_priorpred_samples,
    num_focused_samples,
    mixprob_thresh = 0,
    postprob_thresh = 0,
    ridge_post = 1e-6,
    verbose = 1,
    output_dir = NULL,
    include_covariance_in_S1 = TRUE,
    include_mixture_weights_in_S1 = TRUE,
    reuse_q0_init_for_q1_raw = TRUE,
    num_refreshed_samples = NULL,
    eps_defensive_stage1 = 0.05,
    defensive_source_stage1 = c("prior", "none"),
    factor_cvMH_stage1 = 1,
    factor_cvMH_stage2 = 1,
    summary_scale_floor = 1e-4,
    eps_defensive_stage2 = 0.05,
    defensive_source_stage2 = c("prior", "none")
) {

  defensive_source_stage1 <- match.arg(defensive_source_stage1)
  defensive_source_stage2 <- match.arg(defensive_source_stage2)
  num_param <- length(sample_prior())
  if (is.null(num_refreshed_samples)) {
    num_refreshed_samples <- num_focused_samples
  }

  if (defensive_source_stage1 == "none") {
    eps_defensive_stage1 <- 0
  }
  eps_defensive_stage1 <- max(0, min(1, eps_defensive_stage1))

  if (defensive_source_stage2 == "none") {
    eps_defensive_stage2 <- 0
  }
  eps_defensive_stage2 <- max(0, min(1, eps_defensive_stage2))


  if (verbose >= 1) {
    cat("Building two-round local surrogate models...\n")
    cat("num_param =", num_param, "\n")
    cat("dim_data =", dim_data, "\n")
    cat("num_priorpred_samples =", num_priorpred_samples, "\n")
    cat("num_focused_samples (stage a) =", num_focused_samples, "\n")
    cat("num_refreshed_samples (stage b) =", num_refreshed_samples, "\n")
    cat("mixprob_thresh (global pi_k cleanup) =", mixprob_thresh, "\n")
    cat("postprob_thresh (omega_k(y_obs) pruning) =", postprob_thresh, "\n")
    cat("eps_defensive_stage1 =", eps_defensive_stage1, "\n")
    cat("defensive_source_stage1 =", defensive_source_stage1, "\n")
    cat("eps_defensive_stage2 =", eps_defensive_stage2, "\n")
    cat("defensive_source_stage2 =", defensive_source_stage2, "\n")
  }

  # --------------------------------
  # 1. Prior-predictive training data
  # --------------------------------
  if (verbose >= 1) cat("Simulating prior-predictive training data for global q0...\n")

  theta0 <- replicate(num_priorpred_samples, sample_prior())
  theta0 <- as.matrix(theta0)

  y0 <- apply(theta0, 2, model, model_param = model_param)
  y0 <- as.matrix(y0)

  if (nrow(y0) != dim_data) {
    stop("The simulated prior-predictive data y0 does not match dim_data.")
  }


  if (!is.null(output_dir)) {
    write.table(t(theta0), file = file.path(output_dir, "theta0_priorpred.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(t(y0), file = file.path(output_dir, "y0_priorpred.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
  }

  # --------------------------------
  # 2. Fit global q0_raw
  # --------------------------------
  if (verbose >= 1) cat("Fitting global surrogate q0_raw...\n")

  fit0_raw <- fit_gllim(
    responses = theta0,
    covariates = y0,
    K = K_start,
    maxiter = maxiter,
    verbose = verbose,
    cov_structure = cov_structure,
    mixprob_thresh = mixprob_thresh,
    init_param = NULL
  )

  q0_pruned <- posterior_prune_fit_at_yobs(
    fit_obj = fit0_raw,
    observedData = observedData,
    num_param = num_param,
    dim_data = dim_data,
    postprob_thresh = postprob_thresh,
    stage_name = "q0_raw",
    ridge_post = ridge_post,
    min_components = 1,
    verbose = verbose
  )

  fit0_raw_unpruned <- fit0_raw
  fit0_raw <- q0_pruned$fit
  cache0_raw <- q0_pruned$cache
  post_mix_prob0 <- q0_pruned$post_mix_prob
  mean_proposal0 <- q0_pruned$mean_proposal
  covar_proposal0 <- q0_pruned$covar_proposal
  post_mix_prob0_unpruned <- q0_pruned$post_mix_prob_full

  if (verbose >= 1) {
    cat("Global q0_raw fitted.\n")
    cat("K0_raw =", fit0_raw$K, "\n")
  }

  # --------------------------------
  # 3. First local design theta_a ~ q0(theta | y_obs), y_a ~ p(y | theta_a)
  # --------------------------------
  if (verbose >= 1) {
    cat("Sampling first local parameters theta_a from q0_raw(theta | y_obs)...\n")
    if (eps_defensive_stage1 > 0) {
      cat("Adding defensive stage-a draws from the prior with probability", eps_defensive_stage1, "per draw.\n")
    }
  }

  theta_a <- sample_valid_from_mixture_build(
    n = num_focused_samples,
    num_param = num_param,
    probs = post_mix_prob0,
    means = mean_proposal0,
    covars = covar_proposal0,
    factor_cvMH = factor_cvMH_stage1,
    prior_pdf = prior_pdf,
    sample_prior = if (defensive_source_stage1 == "prior") sample_prior else NULL,
    eps_defensive = eps_defensive_stage1,
    verbose = verbose >= 1
  )

  if (verbose >= 1) cat("Simulating first local pseudo-data y_a...\n")
  y_a <- apply(theta_a, 2, model, model_param = model_param)
  y_a <- as.matrix(y_a)

  ok_y_a <- apply(y_a, 2, function(z) all(is.finite(z)))
  if (!all(ok_y_a)) {
    stop("Non-finite first local simulations. Number bad: ", sum(!ok_y_a))
  }
  if (nrow(y_a) != dim_data) {
    stop("The first local simulated data y_a does not match dim_data.")
  }

  if (!is.null(output_dir)) {
    write.table(t(theta_a), file = file.path(output_dir, "theta_a_from_q0.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(t(y_a), file = file.path(output_dir, "y_a_from_q0.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
  }

  # --------------------------------
  # 4. Fit local q1_raw on (theta_a, y_a)
  # --------------------------------
  if (verbose >= 1) cat("Fitting first local surrogate q1_raw...\n")

  init_q1_raw <- if (isTRUE(reuse_q0_init_for_q1_raw)) fit0_raw$optpar else NULL

  fit1_raw <- fit_gllim(
    responses = theta_a,
    covariates = y_a,
    K = fit0_raw$K,
    maxiter = maxiter,
    verbose = verbose,
    cov_structure = cov_structure,
    mixprob_thresh = mixprob_thresh,
    init_param = init_q1_raw
  )

  q1_pruned <- posterior_prune_fit_at_yobs(
    fit_obj = fit1_raw,
    observedData = observedData,
    num_param = num_param,
    dim_data = dim_data,
    postprob_thresh = postprob_thresh,
    stage_name = "q1_raw",
    ridge_post = ridge_post,
    min_components = 1,
    verbose = verbose
  )

  fit1_raw_unpruned <- fit1_raw
  fit1_raw <- q1_pruned$fit
  cache1_raw <- q1_pruned$cache
  post_mix_prob1_raw <- q1_pruned$post_mix_prob
  mean_proposal1_raw <- q1_pruned$mean_proposal
  covar_proposal1_raw <- q1_pruned$covar_proposal
  post_mix_prob1_raw_unpruned <- q1_pruned$post_mix_prob_full

  if (verbose >= 1) {
    cat("First local q1_raw fitted.\n")
    cat("K1_raw =", fit1_raw$K, "\n")
  }

  # --------------------------------
  # 5. Second local design theta_b from q1 plus optional defensive draws
  # --------------------------------
  if (verbose >= 1) cat("Sampling second local parameters theta_b from q1_raw(theta | y_obs)...\n")

  n_def <- round(num_refreshed_samples * eps_defensive_stage2)
  n_main <- num_refreshed_samples - n_def

  theta_b_main <- if (n_main > 0) {
    sample_valid_from_mixture_build(
      n = n_main,
      num_param = num_param,
      probs = post_mix_prob1_raw,
      means = mean_proposal1_raw,
      covars = covar_proposal1_raw,
      factor_cvMH = factor_cvMH_stage2,
      prior_pdf = prior_pdf,
      sample_prior = sample_prior,
      eps_defensive = 0,
      verbose = verbose >= 1
    )
  } else {
    matrix(numeric(0), nrow = num_param, ncol = 0)
  }

  theta_b_def <- matrix(numeric(0), nrow = num_param, ncol = 0)
  if (n_def > 0) {
    if (defensive_source_stage2 == "prior") {
      if (verbose >= 1) cat("Adding defensive stage-b draws from the prior...\n")
      theta_b_def <- replicate(n_def, sample_prior())
      theta_b_def <- as.matrix(theta_b_def)
    } else if (defensive_source_stage2 == "none") {
      theta_b_def <- matrix(numeric(0), nrow = num_param, ncol = 0)
    } else {
      stop("Unknown defensive_source_stage2: ", defensive_source_stage2)
    }
  }

  theta_b <- cbind(theta_b_main, theta_b_def)
  if (ncol(theta_b) != num_refreshed_samples) {
    stop("Internal error: theta_b has ", ncol(theta_b), " columns, expected ", num_refreshed_samples, ".")
  }

  if (verbose >= 1) cat("Simulating second local pseudo-data y_b...\n")
  y_b <- apply(theta_b, 2, model, model_param = model_param)
  y_b <- as.matrix(y_b)

  ok_y_b <- apply(y_b, 2, function(z) all(is.finite(z)))
  if (!all(ok_y_b)) {
    stop("Non-finite second local simulations. Number bad: ", sum(!ok_y_b))
  }
  if (nrow(y_b) != dim_data) {
    stop("The second local simulated data y_b does not match dim_data.")
  }

  if (!is.null(output_dir)) {
    write.table(t(theta_b), file = file.path(output_dir, "theta_b_from_q1.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(t(y_b), file = file.path(output_dir, "y_b_from_q1.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
  }

  # --------------------------------
  # 6. Fit final q2_raw on the refreshed local design
  # --------------------------------
  if (verbose >= 1) cat("Fitting final surrogate q2_raw...\n")
  
  
  fit2_raw <- fit_gllim(
    responses = theta_b,
    covariates = y_b,
    K = fit1_raw$K,
    maxiter = maxiter,
    verbose = verbose,
    cov_structure = cov_structure,
    mixprob_thresh = mixprob_thresh,
    init_param = NULL
  )

  q2_pruned <- posterior_prune_fit_at_yobs(
    fit_obj = fit2_raw,
    observedData = observedData,
    num_param = num_param,
    dim_data = dim_data,
    postprob_thresh = postprob_thresh,
    stage_name = "q2_raw",
    ridge_post = ridge_post,
    min_components = 1,
    verbose = verbose
  )

  fit2_raw_unpruned <- fit2_raw
  fit2_raw <- q2_pruned$fit
  cache2_raw <- q2_pruned$cache
  post_mix_prob2_raw <- q2_pruned$post_mix_prob
  mean_proposal2_raw <- q2_pruned$mean_proposal
  covar_proposal2_raw <- q2_pruned$covar_proposal
  post_mix_prob2_raw_unpruned <- q2_pruned$post_mix_prob_full

  if (verbose >= 1) {
    cat("Final surrogate q2_raw fitted.\n")
    cat("K2_raw =", fit2_raw$K, "\n")
  }

  # --------------------------------
  # 7. Build posterior-based summaries S1 from final q2_raw on y_b
  # --------------------------------
  if (verbose >= 1) cat("Constructing S1 summaries from final q2_raw on second local data...\n")

  S1_y1_raw <- summarize_y_matrix_theta_from_cache(
    y_mat = y_b,
    cache = cache2_raw,
    num_theta = num_param,
    include_covariance = include_covariance_in_S1,
    include_mixture_weights = include_mixture_weights_in_S1
  )

  S1_observedData_raw <- posterior_theta_summary_from_cache(
    y = observedData,
    cache = cache2_raw,
    num_theta = num_param,
    include_covariance = include_covariance_in_S1,
    include_mixture_weights = include_mixture_weights_in_S1
  )

  S1_std <- standardize_summary_matrix(
    S_train = S1_y1_raw,
    S_obs = S1_observedData_raw,
    scale_floor = summary_scale_floor
  )

  S1_y1_std <- S1_std$S_train_std
  S1_observedData_std <- S1_std$S_obs_std

  S1_y1 <- S1_y1_std
  S1_observedData <- S1_observedData_std
  dim_summary <- length(S1_observedData)

  if (verbose >= 1) {
    cat("Raw summary dimension =", nrow(S1_y1_raw), "\n")
    cat("Range of standardized summaries: ",
        paste(range(S1_y1_std), collapse = " to "), "\n")
    cat("Smallest 10 summary scales: ",
        paste(sort(S1_std$scale)[1:min(10, length(S1_std$scale))], collapse = ", "), "\n")
    cat("Summary dimension =", dim_summary, "\n")
  }

  if (!is.null(output_dir)) {
    write.table(t(S1_y1_raw), file = file.path(output_dir, "S1_y1_raw.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(matrix(S1_observedData_raw, nrow = 1),
                file = file.path(output_dir, "S1_observedData_raw.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(t(S1_y1_std), file = file.path(output_dir, "S1_y1_standardized.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(matrix(S1_observedData_std, nrow = 1),
                file = file.path(output_dir, "S1_observedData_standardized.csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
  }

  # --------------------------------
  # 8. Optional save of training objects
  # --------------------------------
  if (!is.null(output_dir)) {
    saveRDS(
      list(
        observedData = observedData,
        theta0 = theta0,
        y0 = y0,
        theta_a = theta_a,
        y_a = y_a,
        theta_b = theta_b,
        y_b = y_b,
        eps_defensive_stage1 = eps_defensive_stage1,
        defensive_source_stage1 = defensive_source_stage1,
        eps_defensive_stage2 = eps_defensive_stage2,
        defensive_source_stage2 = defensive_source_stage2,
        S1_y1 = S1_y1,
        S1_observedData = S1_observedData
      ),
      file = file.path(output_dir, "surrogate_training_data.rds")
    )
  }

  # --------------------------------
  # 9. Return surrogate object.
  # fit1_raw/theta1/y1 are mapped to the final q2 stage for ratio correction.
  # --------------------------------
  surrogate_out <- list(
    observedData = observedData,
    num_param = num_param,
    dim_data = dim_data,
    two_round_localization = TRUE,

    fit0_raw_unpruned = fit0_raw_unpruned,
    fit1_intermediate_raw_unpruned = fit1_raw_unpruned,
    fit1_raw_unpruned = fit2_raw_unpruned,

    theta0 = theta0,
    y0 = y0,

    theta1_stage1 = theta_a,
    y1_stage1 = y_a,

    theta1 = theta_b,
    y1 = y_b,

    theta_a = theta_a,
    y_a = y_a,
    theta_b = theta_b,
    y_b = y_b,

    S1_y1_raw = S1_y1_raw,
    S1_observedData_raw = S1_observedData_raw,

    S1_y1 = S1_y1,
    S1_observedData = S1_observedData,

    S1_standardization = list(
      center = S1_std$center,
      scale = S1_std$scale
    ),

    dim_summary = dim_summary,

    summary_config = list(
      include_covariance_in_S1 = include_covariance_in_S1,
      include_mixture_weights_in_S1 = include_mixture_weights_in_S1,
      num_refreshed_samples = num_refreshed_samples,
      eps_defensive_stage1 = eps_defensive_stage1,
      defensive_source_stage1 = defensive_source_stage1,
      factor_cvMH_stage1 = factor_cvMH_stage1,
      factor_cvMH_stage2 = factor_cvMH_stage2,
      eps_defensive_stage2 = eps_defensive_stage2,
      defensive_source_stage2 = defensive_source_stage2,
      summary_scale_floor = summary_scale_floor,
      mixprob_thresh = mixprob_thresh,
      postprob_thresh = postprob_thresh,
      ridge_post = ridge_post
    ),

    fit0_raw = list(
      K = fit0_raw$K,
      A = fit0_raw$A,
      b = fit0_raw$b,
      Sigma = fit0_raw$Sigma,
      c = fit0_raw$c,
      Gamma = fit0_raw$Gamma,
      prior_mix_prob = fit0_raw$prior_mix_prob,
      post_mix_prob = post_mix_prob0,
      mean_proposal = mean_proposal0,
      covar_proposal = covar_proposal0,
      optpar = fit0_raw$optpar,
      K_before_postprob_pruning = if (!is.null(fit0_raw$K_before_postprob_pruning)) fit0_raw$K_before_postprob_pruning else fit0_raw$K,
      postprob_active_components = fit0_raw$postprob_active_components,
      post_mix_prob_unpruned = post_mix_prob0_unpruned
    ),

    fit1_intermediate_raw = list(
      K = fit1_raw$K,
      A = fit1_raw$A,
      b = fit1_raw$b,
      Sigma = fit1_raw$Sigma,
      c = fit1_raw$c,
      Gamma = fit1_raw$Gamma,
      prior_mix_prob = fit1_raw$prior_mix_prob,
      post_mix_prob = post_mix_prob1_raw,
      mean_proposal = mean_proposal1_raw,
      covar_proposal = covar_proposal1_raw,
      optpar = fit1_raw$optpar,
      K_before_postprob_pruning = if (!is.null(fit1_raw$K_before_postprob_pruning)) fit1_raw$K_before_postprob_pruning else fit1_raw$K,
      postprob_active_components = fit1_raw$postprob_active_components,
      post_mix_prob_unpruned = post_mix_prob1_raw_unpruned
    ),

    fit1_raw = list(
      K = fit2_raw$K,
      A = fit2_raw$A,
      b = fit2_raw$b,
      Sigma = fit2_raw$Sigma,
      c = fit2_raw$c,
      Gamma = fit2_raw$Gamma,
      prior_mix_prob = fit2_raw$prior_mix_prob,
      post_mix_prob = post_mix_prob2_raw,
      mean_proposal = mean_proposal2_raw,
      covar_proposal = covar_proposal2_raw,
      optpar = fit2_raw$optpar,
      K_before_postprob_pruning = if (!is.null(fit2_raw$K_before_postprob_pruning)) fit2_raw$K_before_postprob_pruning else fit2_raw$K,
      postprob_active_components = fit2_raw$postprob_active_components,
      post_mix_prob_unpruned = post_mix_prob2_raw_unpruned
    ),

    fit0 = list(
      K = fit0_raw$K,
      A = fit0_raw$A,
      b = fit0_raw$b,
      Sigma = fit0_raw$Sigma,
      c = fit0_raw$c,
      Gamma = fit0_raw$Gamma,
      prior_mix_prob = fit0_raw$prior_mix_prob,
      post_mix_prob = post_mix_prob0,
      mean_proposal = mean_proposal0,
      covar_proposal = covar_proposal0,
      optpar = fit0_raw$optpar,
      K_before_postprob_pruning = if (!is.null(fit0_raw$K_before_postprob_pruning)) fit0_raw$K_before_postprob_pruning else fit0_raw$K,
      postprob_active_components = fit0_raw$postprob_active_components,
      post_mix_prob_unpruned = post_mix_prob0_unpruned
    ),

    fit1 = list(
      K = fit2_raw$K,
      A = fit2_raw$A,
      b = fit2_raw$b,
      Sigma = fit2_raw$Sigma,
      c = fit2_raw$c,
      Gamma = fit2_raw$Gamma,
      prior_mix_prob = fit2_raw$prior_mix_prob,
      post_mix_prob = post_mix_prob2_raw,
      mean_proposal = mean_proposal2_raw,
      covar_proposal = covar_proposal2_raw,
      optpar = fit2_raw$optpar,
      K_before_postprob_pruning = if (!is.null(fit2_raw$K_before_postprob_pruning)) fit2_raw$K_before_postprob_pruning else fit2_raw$K,
      postprob_active_components = fit2_raw$postprob_active_components,
      post_mix_prob_unpruned = post_mix_prob2_raw_unpruned
    ),

    prior_pdf = prior_pdf
  )

  surrogate_out
}
