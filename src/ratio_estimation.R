# ============================================================
# Ratio estimation and ratio-corrected importance sampling
# ------------------------------------------------------------
# The final q2(theta | y_obs) GLLiM fit supplies the proposal.
# Posterior-based standardized summaries are used only by the ratio
# classifier. The default classifier uses linear terms and parameter-summary
# cross interactions with lasso regularization.
# ============================================================

library(glmnet)
library(mixtools)

# ------------------------------------------------------------
vech <- function(S) {
  S[lower.tri(S, diag = TRUE)]
}

# ------------------------------------------------------------
# 1. Feature construction
# ------------------------------------------------------------

pairwise_products <- function(X) {
  p <- ncol(X)
  if (p < 2) return(NULL)
  
  cols <- combn(p, 2, simplify = FALSE)
  out <- lapply(cols, function(idx) X[, idx[1]] * X[, idx[2]])
  out <- do.call(cbind, out)
  
  colnames(out) <- sapply(cols, function(idx) {
    paste0(colnames(X)[idx[1]], ":", colnames(X)[idx[2]])
  })
  
  out
}

cross_products <- function(A, B) {
  pA <- ncol(A)
  pB <- ncol(B)
  
  out_list <- vector("list", pA * pB)
  nm <- character(pA * pB)
  k <- 1
  
  for (i in seq_len(pA)) {
    for (j in seq_len(pB)) {
      out_list[[k]] <- A[, i] * B[, j]
      nm[k] <- paste0(colnames(A)[i], ":", colnames(B)[j])
      k <- k + 1
    }
  }
  
  out <- do.call(cbind, out_list)
  colnames(out) <- nm
  out
}

build_features <- function(theta_block,
                           y_block,
                           include_linear = TRUE,
                           include_quadratic = FALSE,
                           include_within_block_interactions = FALSE,
                           include_cross_interactions = TRUE) {
  theta_block <- as.matrix(theta_block)
  y_block     <- as.matrix(y_block)
  
  if (nrow(theta_block) != nrow(y_block)) {
    stop("theta_block and y_block must have the same number of rows.")
  }
  
  d_theta <- ncol(theta_block)
  d_y     <- ncol(y_block)
  
  colnames(theta_block) <- paste0("th", seq_len(d_theta))
  colnames(y_block)     <- paste0("y", seq_len(d_y))
  
  feat_list <- list()
  
  if (include_linear) {
    feat_list <- c(feat_list, list(theta_block, y_block))
  }
  
  if (include_quadratic) {
    theta_sq <- theta_block^2
    y_sq     <- y_block^2
    colnames(theta_sq) <- paste0(colnames(theta_block), "^2")
    colnames(y_sq)     <- paste0(colnames(y_block), "^2")
    feat_list <- c(feat_list, list(theta_sq, y_sq))
  }
  
  if (include_within_block_interactions) {
    theta_pair <- pairwise_products(theta_block)
    y_pair     <- pairwise_products(y_block)
    
    if (!is.null(theta_pair)) feat_list <- c(feat_list, list(theta_pair))
    if (!is.null(y_pair))     feat_list <- c(feat_list, list(y_pair))
  }
  
  if (include_cross_interactions) {
    cross_mat <- cross_products(theta_block, y_block)
    feat_list <- c(feat_list, list(cross_mat))
  }
  
  X <- do.call(cbind, feat_list)
  storage.mode(X) <- "double"
  X
}

build_features_one <- function(theta_vec, y_vec, ...) {
  theta_block <- matrix(theta_vec, nrow = 1)
  y_block     <- matrix(y_vec, nrow = 1)
  build_features(theta_block, y_block, ...)
}

build_features_fixed_y <- function(theta_mat,
                                   y_vec,
                                   include_linear = TRUE,
                                   include_quadratic = FALSE,
                                   include_within_block_interactions = FALSE,
                                   include_cross_interactions = TRUE) {
  theta_mat <- as.matrix(theta_mat)
  y_vec <- as.numeric(y_vec)
  
  y_block <- matrix(y_vec, nrow = nrow(theta_mat), ncol = length(y_vec), byrow = TRUE)
  
  build_features(
    theta_block = theta_mat,
    y_block = y_block,
    include_linear = include_linear,
    include_quadratic = include_quadratic,
    include_within_block_interactions = include_within_block_interactions,
    include_cross_interactions = include_cross_interactions
  )
}

# ------------------------------------------------------------
# 2. Build ratio-training dataset
# ------------------------------------------------------------

build_ratio_training_data <-
  function(theta_mat,
           y_mat,
           feature_args = list()) {
    n <- ncol(theta_mat)
    
    theta_joint <- t(theta_mat)
    y_joint     <- t(y_mat)
    
    perm <- sample.int(n)
    theta_prod <- theta_joint
    y_prod     <- y_joint[perm, , drop = FALSE]
    
    X_joint <- do.call(build_features, c(list(theta_joint, y_joint), feature_args))
    X_prod  <- do.call(build_features, c(list(theta_prod, y_prod), feature_args))
    
    X <- rbind(X_joint, X_prod)
    z <- c(rep(1, n), rep(0, n))
    
    list(X = X, z = z)
  }



# ------------------------------------------------------------
# 3. Fit ratio classifier
# ------------------------------------------------------------
fit_ratio_classifier <-
  function(theta_mat,
           y_mat,
           alpha = 1,
           feature_args = list(),
           use_cv = TRUE,
           nfolds = 10) {
    dat <- build_ratio_training_data(
      theta_mat = theta_mat,
      y_mat = y_mat,
      feature_args = feature_args
    )
    
    if (use_cv) {
      fit <- cv.glmnet(
        x = as.matrix(dat$X),
        y = dat$z,
        family = "binomial",
        alpha = alpha,
        standardize = TRUE,
        nfolds = nfolds
      #  nlambda = 100,
      #  lambda.min.ratio = 1e-4
      )
    } else {
      fit <- glmnet(
        x = as.matrix(dat$X),
        y = dat$z,
        family = "binomial",
        alpha = alpha,
        standardize = TRUE
      )
    }
    
    fit$feature_args <- feature_args
    fit$use_cv_fit <- use_cv
    fit
  }

#-------------------------------------------------------------
#3bis Resolve lambda choice
#-------------------------------------------------------------

resolve_lambda_choice <- function(ratio_fit, s) {
  if (inherits(ratio_fit, "cv.glmnet")) {
    return(s)
  }
  
  # plain glmnet object: s must be numeric
  if (is.numeric(s)) {
    return(s)
  }
  
  # fallback for temporary glmnet runs
  if (is.character(s)) {
    return(ratio_fit$lambda[length(ratio_fit$lambda)])
  }
  
  stop("Unsupported lambda choice.")
}

# ------------------------------------------------------------
# 4. Predict ratio
# ------------------------------------------------------------

predict_log_ratio <-
  function(ratio_fit,
           theta_vec,
           y_used,
           s = "lambda.1se") {
    feature_args <- ratio_fit$feature_args
    x_new <- do.call(build_features_one, c(list(theta_vec, y_used), feature_args))
    s_used <- resolve_lambda_choice(ratio_fit, s)
    log_r <- predict(ratio_fit, newx = x_new, s = s_used, type = "link")
    as.numeric(log_r)
  }

predict_log_ratio_batch <- function(ratio_fit,
                                    theta_mat,
                                    y_used,
                                    s = "lambda.1se") {
  feature_args <- ratio_fit$feature_args
  
  x_new <- do.call(
    build_features_fixed_y,
    c(list(theta_mat = theta_mat, y_vec = y_used), feature_args)
  )
  
  s_used <- resolve_lambda_choice(ratio_fit, s)
  
  log_r <- predict(ratio_fit, newx = x_new, s = s_used, type = "link")
  as.numeric(log_r)
}

predict_ratio <-
  function(ratio_fit,
           theta_vec,
           y_used,
           s = "lambda.1se") {
    exp(predict_log_ratio(
      ratio_fit = ratio_fit,
      theta_vec = theta_vec,
      y_used = y_used,
      s = s
    ))
  }

# ------------------------------------------------------------
# 5. Proposal density q2(theta | y_obs)
# ------------------------------------------------------------

proposal_density <- function(theta_vec, probs, means, covars, factor_cvMH = 1) {
  if (length(probs) > 1) {
    covars_use <- covars
    for (k in seq_len(dim(covars_use)[3])) {
      covars_use[, , k] <- factor_cvMH * covars_use[, , k]
    }
  } else {
    covars_use <- factor_cvMH * covars
  }
  
  q_val <- mixpdf(
    theta = theta_vec,
    probs = probs,
    means = means,
    covars = covars_use
  )
  as.numeric(q_val)
}

# ------------------------------------------------------------
# 6. Posterior-based pruning of proposal components
# ------------------------------------------------------------

prune_proposal_components <- function(post_mix_prob,
                                      mean_proposal,
                                      covar_proposal,
                                      postprob_thresh = 0,
                                      verbose = TRUE) {
  K <- length(post_mix_prob)
  keep <- which(post_mix_prob >= postprob_thresh)
  
  if (length(keep) == 0) {
    stop("All components were removed by pruning.")
  }
  
  if (length(keep) == K) {
    if (verbose) cat("No proposal components removed.\n")
    return(list(
      post_mix_prob = post_mix_prob,
      mean_proposal = mean_proposal,
      covar_proposal = covar_proposal,
      kept_components = keep
    ))
  }
  
  if (verbose) {
    cat("Pruning proposal components:\n")
    cat("Original K =", K, "\n")
    cat("Remaining K =", length(keep), "\n")
    cat("Removed components:", setdiff(seq_len(K), keep), "\n")
  }
  
  post_mix_prob_new <- post_mix_prob[keep]
  post_mix_prob_new <- post_mix_prob_new / sum(post_mix_prob_new)
  
  if (length(keep) == 1) {
    mean_proposal_new  <- matrix(mean_proposal[, keep], ncol = 1)
    covar_proposal_new <- covar_proposal[, , keep]
  } else {
    mean_proposal_new  <- mean_proposal[, keep, drop = FALSE]
    covar_proposal_new <- covar_proposal[, , keep, drop = FALSE]
  }
  
  list(
    post_mix_prob = post_mix_prob_new,
    mean_proposal = mean_proposal_new,
    covar_proposal = covar_proposal_new,
    kept_components = keep
  )
}


prepare_gllim_posterior_cache <- function(fit_obj, num_param, dim_data) {
  cache <- prepare_posterior_cache_from_fit(
    fit_obj = fit_obj,
    num_param = num_param,
    dim_data = dim_data
  )
  
  cache$pi <- fit_obj$prior_mix_prob
  cache$K <- fit_obj$K
  
  cache
}

posterior_mixture_from_input <- function(cache, x) {
  post_mix_prob <- compute_post_mix_prob(
    A = cache$A,
    b = cache$b,
    c = cache$c,
    Sigma = cache$Sigma,
    Gamma = cache$Gamma,
    pi = cache$pi,
    K = cache$K,
    observedData = x
  )
  
  mean_proposal <- compute_mean_proposal(
    A_star = cache$A_star,
    observedData = x,
    b_star = cache$b_star,
    num_param = cache$num_param,
    K = cache$K
  )
  
  covar_proposal <- cache$Sigma_star
  
  list(
    probs = as.numeric(post_mix_prob),
    means = mean_proposal,
    covars = covar_proposal
  )
}

select_local_ratio_subset_euclidean <- function(theta_mat,
                                                summary_mat,
                                                observed_summary,
                                                max_ratio_train,
                                                local_pool_size = max_ratio_train) {
  theta_mat <- as.matrix(theta_mat)
  summary_mat <- as.matrix(summary_mat)
  observed_summary <- as.numeric(observed_summary)
  
  if (ncol(theta_mat) != ncol(summary_mat)) {
    stop("theta_mat and summary_mat must have the same number of columns.")
  }
  
  if (length(observed_summary) != nrow(summary_mat)) {
    stop("observed_summary length must match nrow(summary_mat).")
  }
  
  # Robust scaling using median and MAD
  med_vec <- apply(summary_mat, 1, median)
  mad_vec <- apply(summary_mat, 1, mad)
  mad_vec_use <- pmax(mad_vec, 1e-4)
  
  S_train_rob <- sweep(summary_mat, 1, med_vec, "-")
  S_train_rob <- sweep(S_train_rob, 1, mad_vec_use, "/")
  
  S_obs_rob <- (observed_summary - med_vec) / mad_vec_use
  
  dists <- colSums((S_train_rob - S_obs_rob)^2)
  ord <- order(dists)
  
  local_pool_size <- min(local_pool_size, ncol(summary_mat))
  max_ratio_train <- min(max_ratio_train, local_pool_size)
  
  idx_pool <- ord[seq_len(local_pool_size)]
  
  if (max_ratio_train < local_pool_size) {
    idx_keep <- sample(idx_pool, size = max_ratio_train, replace = FALSE)
  } else {
    idx_keep <- idx_pool
  }
  
  list(
    theta_mat = theta_mat[, idx_keep, drop = FALSE],
    summary_mat = summary_mat[, idx_keep, drop = FALSE],
    idx_keep = idx_keep,
    idx_pool = idx_pool,
    dists = dists
  )
}

# ------------------------------------------------------------
# 7. Helper: sample only prior-valid final proposal draws
# ------------------------------------------------------------
# This helper samples
# in theta-space from the q2(theta | y_obs) proposal and rejects
# draws outside the prior support.  The IS denominator below is still the
# Gaussian-mixture density q2(theta | y_obs); the truncation constant
# is common to all accepted particles and cancels in normalized weights.
# ------------------------------------------------------------
sample_valid_from_mixture_ratio <- function(n,
                                                num_param,
                                                probs,
                                                means,
                                                covars,
                                                factor_cvMH,
                                                prior_pdf,
                                                max_tries_factor = 200,
                                                verbose = TRUE) {
  theta_mat <- matrix(NA_real_, nrow = num_param, ncol = n)

  i <- 1
  n_try <- 0
  max_try <- max_tries_factor * n

  while (i <= n && n_try < max_try) {
    n_try <- n_try + 1

    th <- par_inflated_mix(
      probs = probs,
      means = means,
      covars = covars,
      factor_cvMH = factor_cvMH
    )

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
      "Could not obtain enough valid final proposal draws. ",
      "Obtained ", i - 1, " out of ", n, " after ", n_try, " attempts. ",
      "Consider reducing factor_cvMH or using a transformed parametrisation."
    )
  }

  if (isTRUE(verbose)) {
    cat(
      "Valid final proposal draw acceptance rate:",
      round(n / n_try, 4),
      "(", n, "accepted out of", n_try, "attempts )\n"
    )
  }

  if (n / n_try < 0.05) {
    warning(
      "Low valid final proposal acceptance rate: ", round(n / n_try, 4),
      ". Consider reducing factor_cvMH or using a transformed parametrisation."
    )
  }

  theta_mat
}

# ------------------------------------------------------------
# 8. Main ratio-corrected sampler
# ------------------------------------------------------------

ratio_corrected_sample <-
  function(
    surrogate_out,
    num_final_samples,
    alpha = 1,
    factor_cvMH = 1,
    prune_components = FALSE,
    postprob_thresh = 0,
    feature_args = list(
      include_linear = TRUE,
      include_quadratic = FALSE,
      include_within_block_interactions = FALSE,
      include_cross_interactions = TRUE
    ),
    max_ratio_train = NULL,
    local_pool_size = NULL,
    use_local_subset = FALSE,
    use_cv = TRUE,
    nfolds = 10,
    lambda_choice = "lambda.1se",
    report_timing = TRUE
  ) {
    
    if (is.null(surrogate_out$fit1_raw) ||
        is.null(surrogate_out$theta1) ||
        is.null(surrogate_out$y1) ||
        is.null(surrogate_out$S1_y1) ||
        is.null(surrogate_out$S1_observedData) ||
        is.null(surrogate_out$observedData)) {
      stop("surrogate_out must contain fit1_raw, theta1, y1, S1_y1, S1_observedData, and observedData.")
    }
    
    # -----------------------------
    # Extract final proposal q2(theta | y_obs)
    # -----------------------------
    num_param <- surrogate_out$num_param
    dim_data <- surrogate_out$dim_data

    cache_raw <- prepare_gllim_posterior_cache(
      fit_obj = surrogate_out$fit1_raw,
      num_param = num_param,
      dim_data = dim_data
    )
    gmm_raw_obs <- posterior_mixture_from_input(cache_raw, surrogate_out$observedData)

    post_mix_prob  <- gmm_raw_obs$probs
    mean_proposal  <- gmm_raw_obs$means
    covar_proposal <- gmm_raw_obs$covars
    
    K_original <- length(post_mix_prob)
    kept_components <- seq_len(K_original)
    
    if (prune_components) {
      pruned <- prune_proposal_components(
        post_mix_prob = post_mix_prob,
        mean_proposal = mean_proposal,
        covar_proposal = covar_proposal,
        postprob_thresh = postprob_thresh,
        verbose = FALSE
      )
      
      post_mix_prob  <- pruned$post_mix_prob
      mean_proposal  <- pruned$mean_proposal
      covar_proposal <- pruned$covar_proposal
      
      K_used <- length(post_mix_prob)
      kept_components <- pruned$kept_components
      
      cat(
        "Proposal pruning applied: ",
        "K original = ", K_original,
        ", K used = ", K_used,
        ", threshold = ", postprob_thresh, "
",
        sep = ""
      )
    } else {
      K_used <- K_original
      cat("No proposal pruning: K = ", K_used, "
", sep = "")
    }
    
    prior_pdf <- surrogate_out$prior_pdf
    observed_summary <- surrogate_out$S1_observedData
    
    # -----------------------------
    # Focused classifier training set: summaries only for discriminator
    # -----------------------------
    theta_train <- surrogate_out$theta1
    y_train <- surrogate_out$S1_y1
    observed_summary <- surrogate_out$S1_observedData
    
    local_subset <- NULL
    
    if (!is.null(max_ratio_train) && ncol(theta_train) > max_ratio_train) {
      if (isTRUE(use_local_subset)) {
        if (is.null(local_pool_size)) {
          local_pool_size <- min(2000, ncol(theta_train))
        }
        
        local_subset <- select_local_ratio_subset_euclidean(
          theta_mat = theta_train,
          summary_mat = y_train,
          observed_summary = observed_summary,
          max_ratio_train = max_ratio_train,
          local_pool_size = local_pool_size
        )
        
        theta_train <- local_subset$theta_mat
        y_train <- local_subset$summary_mat
        
        cat("Summary-space robust (MAD-scaled) local subset selection applied: ",
            "pool size = ", local_pool_size,
            ", classifier training size = ", ncol(theta_train), "
",
            sep = "")
      } else {
        idx_keep <- sample.int(ncol(theta_train), max_ratio_train)
        theta_train <- theta_train[, idx_keep, drop = FALSE]
        y_train <- y_train[, idx_keep, drop = FALSE]
        
        cat("Random subset selection applied: ",
            "classifier training size = ", ncol(theta_train), "
",
            sep = "")
      }
    }
    
    # -----------------------------
    # Fit ratio classifier
    # -----------------------------
    t0 <- Sys.time()
    
    ratio_fit <- fit_ratio_classifier(
      theta_mat = theta_train,
      y_mat = y_train,
      alpha = alpha,
      feature_args = feature_args,
      use_cv = use_cv,
      nfolds = nfolds
    )
    
    t1 <- Sys.time()
    fit_time_sec <- as.numeric(t1 - t0, units = "secs")
    
    if (report_timing) {
      cat("Time for fit_ratio_classifier:", fit_time_sec, "seconds
")
    }
    
    # -----------------------------
    # Importance sampling from q2
    # -----------------------------
    d <- nrow(theta_train)
    theta_prop <- sample_valid_from_mixture_ratio(
      n = num_final_samples,
      num_param = d,
      probs = post_mix_prob,
      means = mean_proposal,
      covars = covar_proposal,
      factor_cvMH = factor_cvMH,
      prior_pdf = prior_pdf,
      verbose = report_timing
    )
    weights     <- numeric(num_final_samples)
    log_weights <- numeric(num_final_samples)
    
    theta_prop_t <- t(theta_prop)
    
    log_r_hat_all <- predict_log_ratio_batch(
      ratio_fit = ratio_fit,
      theta_mat = theta_prop_t,
      y_used = observed_summary,
      s = lambda_choice
    )
    
    for (j in seq_len(num_final_samples)) {
      theta_j <- theta_prop[, j]
      
      q_hat <- proposal_density(
        theta_vec = theta_j,
        probs = post_mix_prob,
        means = mean_proposal,
        covars = covar_proposal,
        factor_cvMH = factor_cvMH
      )
      
      p_hat <- prior_pdf(theta_j)
      
      if (!is.finite(p_hat) || p_hat <= 0) {
        log_weights[j] <- -Inf
      } else if (!is.finite(q_hat) || q_hat <= 0) {
        log_weights[j] <- -Inf
      } else if (!is.finite(log_r_hat_all[j])) {
        log_weights[j] <- -Inf
      } else {
        log_weights[j] <- log(p_hat) + log_r_hat_all[j] - log(q_hat)
      }
    }
    
    t2 <- Sys.time()
    is_time_sec <- as.numeric(t2 - t1, units = "secs")
    
    if (report_timing) {
      cat("Time for IS loop:", is_time_sec, "seconds
")
    }
    
    # -----------------------------
    # Stabilize raw log weights
    # -----------------------------
    log_weights_raw <- log_weights
    max_log_w <- max(log_weights_raw)
    
    if (!is.finite(max_log_w)) {
      stop("All importance weights are zero or non-finite. The proposal produced no samples with positive prior support and finite ratio/proposal density.")
    }
    
    weights <- exp(log_weights_raw - max_log_w)
    weights[!is.finite(weights)] <- 0
    
    if (sum(weights) <= 0) {
      stop("All importance weights are zero or non-finite after stabilization.")
    }
    
    weights <- weights / sum(weights)
    ESS <- 1 / sum(weights^2)
    
    # -----------------------------
    # Pareto-smoothed importance sampling
    # -----------------------------
    psis_obj <- loo::psis(log_weights_raw)
    weights_psis <- weights(psis_obj, log = FALSE, normalize = TRUE)
    ESS_psis <- 1 / sum(weights_psis^2)
    pareto_k <- loo::pareto_k_values(psis_obj)
    
    # -----------------------------
    # Truncated weights (diagnostic)
    # -----------------------------
    Nc <- floor(sqrt(length(weights)))
    ord <- order(weights, decreasing = TRUE)
    cap <- weights[ord[Nc]]
    weights_clip <- pmin(weights, cap)
    weights_trunc <- weights_clip / sum(weights_clip)
    ESS_trunc <- 1 / sum(weights_trunc^2)
    
    # -----------------------------
    # Weighted resampling
    # -----------------------------
    idx_psis <- sample.int(
      num_final_samples,
      size = num_final_samples,
      replace = TRUE,
      prob = weights_psis
    )
    theta_resampled <- theta_prop[, idx_psis, drop = FALSE]
    
    idx_trunc <- sample.int(
      num_final_samples,
      size = num_final_samples,
      replace = TRUE,
      prob = weights_trunc
    )
    theta_resampled_trunc <- theta_prop[, idx_trunc, drop = FALSE]
    
    list(
      ratio_fit = ratio_fit,
      
      theta_prop = theta_prop,
      weights = weights,
      weights_psis = weights_psis,
      theta_resampled = theta_resampled,
      
      ESS = ESS,
      ESS_psis = ESS_psis,
      pareto_k = pareto_k,
      
      K_original = K_original,
      K_used = K_used,
      kept_components = kept_components,
      
      post_mix_prob_used = post_mix_prob,
      mean_proposal_used = mean_proposal,
      covar_proposal_used = covar_proposal,
      
      local_subset = local_subset,
      
      fit_time_sec = fit_time_sec,
      is_time_sec = is_time_sec,
      
      log_weights = log_weights,
      weights_trunc = weights_trunc,
      theta_resampled_trunc = theta_resampled_trunc,
      ESS_trunc = ESS_trunc,
      
      # ----------------------------------------------------------
      # Added for posterior-local SBC
      # ----------------------------------------------------------
      # This is the final GLLiM fit used to build q(theta | y_obs)
      # before ratio correction. It is needed to condition on new
      # pseudo-observations y_i^* in local SBC.
      last_gllim_fit = surrogate_out$fit1_raw,
      
      # Cache for conditioning the same GLLiM model at arbitrary y.
      # This avoids recomputing the posterior matrices each time.
      last_gllim_cache = cache_raw,
      dim_data = dim_data,

      ratio_config = list(
        alpha = alpha, feature_args = feature_args, use_cv = use_cv,
        nfolds = nfolds, lambda_choice = lambda_choice
      ),
      
      # The original observed data and summary used in the ratio correction.
      observedData = surrogate_out$observedData,
      observed_summary = observed_summary
    )
  }

