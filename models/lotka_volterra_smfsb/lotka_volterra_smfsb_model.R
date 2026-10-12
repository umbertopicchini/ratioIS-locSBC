# ============================================================
# models/lotka_volterra_smfsb/lotka_volterra_smfsb_model.R
# ------------------------------------------------------------
# Stochastic Lotka--Volterra model based on the smfsb package.
#
# Model:
#   X1 -> 2 X1          rate th1 * X1       prey reproduction
#   X1 + X2 -> 2 X2     rate th2 * X1 * X2  predation/predator reproduction
#   X2 -> empty         rate th3 * X2       predator death
#
# Inference parameterisation:
#   phi = (log_th1, log_th2, log_th3)
#   th  = exp(phi)
#
# Default prior:
#   independent Normal priors on phi, centred at the canonical smfsb values
#   th = (1, 0.005, 0.6), with sd 0.5 on each log-rate.
#
# Observation model:
#   simulate a stochastic LV trajectory on the observed time grid and add
#   independent Gaussian observation noise to each species/time point.
#   The default noise_sd = 10 matches the smfsb LVnoise10 data.
# ============================================================

if (!requireNamespace("smfsb", quietly = TRUE)) {
  stop("Package 'smfsb' is required. Install it with install.packages('smfsb').")
}

# -----------------------------
# Global dimensions and labels
# -----------------------------
dim_param <- 3
num_param <- 3

theta_labels <- c("th1", "th2", "th3")
phi_labels <- paste0("log_", theta_labels)

# A default dimension; the run script should reset dim_data after loading data.
dim_data <- NA_integer_

# -----------------------------
# Default model settings
# -----------------------------
lv_smfsb_settings <- list(
  x0 = c(x1 = 50, x2 = 100),
  t0 = 0,
  times = seq(0, 20, by = 2),
  noise_sd = c(x1 = 10, x2 = 10),
  prior_logmean = log(c(th1 = 1, th2 = 0.005, th3 = 0.6)),
  prior_logsd = c(th1 = 0.5, th2 = 0.5, th3 = 0.5),
  use_compiled_step = TRUE,
  max_attempts = 5,
  min_count_for_likelihood = 0
)

# -----------------------------
# Helper functions
# -----------------------------
inference_to_output_scale <- function(phi_mat) {
  exp(phi_mat)
}

theta_for_feature_diagnostics <- function(phi_mat) {
  inference_to_output_scale(phi_mat)
}

lv_extract_times <- function(observed_object) {
  if (inherits(observed_object, "ts")) {
    return(as.numeric(stats::time(observed_object)))
  }
  if (is.matrix(observed_object) || is.data.frame(observed_object)) {
    if ("time" %in% colnames(observed_object)) {
      return(as.numeric(observed_object[, "time"]))
    }
    if ("t" %in% colnames(observed_object)) {
      return(as.numeric(observed_object[, "t"]))
    }
  }
  stop("Could not infer observation times. Provide times explicitly to make_lv_smfsb_model_param().")
}

lv_observed_matrix <- function(observed_object) {
  mat <- as.matrix(observed_object)

  # smfsb irregular timedData matrices may include a time column.
  if ("time" %in% colnames(mat)) {
    mat <- mat[, setdiff(colnames(mat), "time"), drop = FALSE]
  }
  if ("t" %in% colnames(mat)) {
    mat <- mat[, setdiff(colnames(mat), "t"), drop = FALSE]
  }

  if (ncol(mat) < 2) {
    stop("Expected a two-species LV data object with columns x1 and x2.")
  }
  mat <- mat[, seq_len(2), drop = FALSE]
  colnames(mat) <- c("x1", "x2")
  storage.mode(mat) <- "double"
  mat
}

lv_flatten_observation <- function(y_mat) {
  y_mat <- as.matrix(y_mat)
  if (ncol(y_mat) != 2 && nrow(y_mat) == 2) {
    y_mat <- t(y_mat)
  }
  if (ncol(y_mat) != 2) {
    stop("Expected matrix with two columns: x1 and x2.")
  }
  c(as.numeric(y_mat[, 1]), as.numeric(y_mat[, 2]))
}

lv_unflatten_observation <- function(y_vec, n_times) {
  y_vec <- as.numeric(y_vec)
  if (length(y_vec) != 2 * n_times) {
    stop("y_vec length is incompatible with n_times.")
  }
  cbind(
    x1 = y_vec[seq_len(n_times)],
    x2 = y_vec[n_times + seq_len(n_times)]
  )
}

make_lv_smfsb_model_param <- function(times = NULL,
                                      observed_object = NULL,
                                      x0 = lv_smfsb_settings$x0,
                                      t0 = NULL,
                                      noise_sd = lv_smfsb_settings$noise_sd,
                                      prior_logmean = lv_smfsb_settings$prior_logmean,
                                      prior_logsd = lv_smfsb_settings$prior_logsd,
                                      use_compiled_step = TRUE,
                                      max_attempts = lv_smfsb_settings$max_attempts) {
  if (is.null(times)) {
    if (!is.null(observed_object)) {
      times <- lv_extract_times(observed_object)
    } else {
      times <- lv_smfsb_settings$times
    }
  }
  times <- as.numeric(times)
  if (length(times) < 2 || any(!is.finite(times))) {
    stop("times must contain at least two finite observation times.")
  }

  if (is.null(t0)) {
    t0 <- min(times)
  }

  if (length(noise_sd) == 1) {
    noise_sd <- rep(as.numeric(noise_sd), 2)
    names(noise_sd) <- c("x1", "x2")
  }
  noise_sd <- as.numeric(noise_sd)
  if (length(noise_sd) != 2 || any(!is.finite(noise_sd)) || any(noise_sd <= 0)) {
    stop("noise_sd must be a positive scalar or positive length-2 vector.")
  }
  names(noise_sd) <- c("x1", "x2")

  prior_logmean <- as.numeric(prior_logmean)
  prior_logsd <- as.numeric(prior_logsd)
  if (length(prior_logmean) != dim_param || length(prior_logsd) != dim_param) {
    stop("prior_logmean and prior_logsd must have length 3.")
  }
  names(prior_logmean) <- theta_labels
  names(prior_logsd) <- theta_labels

  list(
    x0 = as.numeric(x0),
    t0 = as.numeric(t0),
    times = times,
    n_times = length(times),
    noise_sd = noise_sd,
    prior_logmean = prior_logmean,
    prior_logsd = prior_logsd,
    use_compiled_step = isTRUE(use_compiled_step),
    max_attempts = as.integer(max_attempts)
  )
}

# -----------------------------
# Prior on inference scale phi
# -----------------------------
sample_prior <- function(n = NULL, model_param = lv_smfsb_settings) {
  prior_logmean <- as.numeric(model_param$prior_logmean)
  prior_logsd <- as.numeric(model_param$prior_logsd)

  if (is.null(n)) {
    phi <- stats::rnorm(dim_param, mean = prior_logmean, sd = prior_logsd)
    names(phi) <- phi_labels
    return(phi)
  }

  n <- as.integer(n)
  phi <- matrix(NA_real_, nrow = dim_param, ncol = n)
  for (j in seq_len(dim_param)) {
    phi[j, ] <- stats::rnorm(n, mean = prior_logmean[j], sd = prior_logsd[j])
  }
  rownames(phi) <- phi_labels
  phi
}

prior_pdf <- function(phi, model_param = lv_smfsb_settings) {
  phi <- as.matrix(phi)
  if (nrow(phi) != dim_param && ncol(phi) == dim_param) {
    phi <- t(phi)
  }
  if (nrow(phi) != dim_param) {
    stop("phi must be a 3 x n matrix, n x 3 matrix, or length-3 vector.")
  }

  prior_logmean <- as.numeric(model_param$prior_logmean)
  prior_logsd <- as.numeric(model_param$prior_logsd)

  dens <- rep(1, ncol(phi))
  for (j in seq_len(dim_param)) {
    dens <- dens * stats::dnorm(phi[j, ], mean = prior_logmean[j], sd = prior_logsd[j])
  }
  dens
}

# -----------------------------
# Stochastic LV simulator
# -----------------------------
simulate_lv_smfsb_latent <- function(th, model_param) {
  th <- as.numeric(th)
  if (length(th) != dim_param || any(!is.finite(th)) || any(th <= 0)) {
    return(NULL)
  }

  x0 <- model_param$x0
  names(x0) <- c("x1", "x2")
  times <- as.numeric(model_param$times)
  t0 <- as.numeric(model_param$t0)

  # Prefer the compiled smfsb LV transition kernel when requested.
  step_fun <- if (isTRUE(model_param$use_compiled_step)) {
    smfsb::stepLVc
  } else {
    data(spnModels, package = "smfsb", envir = environment())
    smfsb::StepGillespie(LV)
  }

  out <- tryCatch(
    smfsb::simTimes(x0, t0, times, step_fun, th = th),
    error = function(e) NULL,
    warning = function(w) invokeRestart("muffleWarning")
  )

  if (is.null(out)) {
    return(NULL)
  }

  out <- as.matrix(out)
  if (ncol(out) < 2) {
    return(NULL)
  }
  out <- out[, seq_len(2), drop = FALSE]
  colnames(out) <- c("x1", "x2")

  if (!all(dim(out) == c(length(times), 2)) || any(!is.finite(out))) {
    return(NULL)
  }

  out
}

add_lv_observation_noise <- function(latent_mat, model_param) {
  if (is.null(latent_mat) || any(!is.finite(latent_mat))) {
    return(NULL)
  }

  noise_sd <- as.numeric(model_param$noise_sd)
  y <- latent_mat
  y[, 1] <- y[, 1] + stats::rnorm(nrow(y), mean = 0, sd = noise_sd[1])
  y[, 2] <- y[, 2] + stats::rnorm(nrow(y), mean = 0, sd = noise_sd[2])

  y
}

simulate_lv_smfsb_one <- function(phi, model_param) {
  th <- exp(as.numeric(phi))
  latent <- simulate_lv_smfsb_latent(th, model_param)
  y_mat <- add_lv_observation_noise(latent, model_param)
  if (is.null(y_mat)) {
    return(rep(NaN, 2 * model_param$n_times))
  }
  lv_flatten_observation(y_mat)
}

# -----------------------------
# Main simulator
# -----------------------------
model <- function(theta, model_param = NULL) {
  if (is.null(model_param)) {
    model_param <- make_lv_smfsb_model_param()
  }

  theta <- as.matrix(theta)
  if (nrow(theta) != dim_param && ncol(theta) == dim_param) {
    theta <- t(theta)
  }
  if (nrow(theta) != dim_param) {
    stop("theta must be a 3 x n matrix, an n x 3 matrix, or a length-3 vector.")
  }

  n <- ncol(theta)
  y <- matrix(NA_real_, nrow = 2 * model_param$n_times, ncol = n)
  max_attempts <- as.integer(model_param$max_attempts)
  retry_count <- 0L

  for (i in seq_len(n)) {
    
    phi_i <- theta[, i]
    success <- FALSE
    
    for (attempt in seq_len(max_attempts)) {
      
      yi <- simulate_lv_smfsb_one(
        phi = phi_i,
        model_param = model_param
      )
      
      if (all(is.finite(yi))) {
        y[, i] <- yi
        success <- TRUE
        break
      }
      
      retry_count <- retry_count + 1L
    }
    
    if (!success) {
      stop(
        "model(): failed at phi = (",
        paste(signif(phi_i, 6), collapse = ", "),
        ") after ",
        max_attempts,
        " attempts."
      )
    }
  }
  
  cat(
    "\nLV simulator diagnostic: total failed attempts =",
    retry_count,
    "\n"
  )

  rownames(y) <- c(
    paste0("x1_t", seq_len(model_param$n_times)),
    paste0("x2_t", seq_len(model_param$n_times))
  )
  y
}

# Noise-free latent observation vector, useful for plotting/debugging.
model_latent_observation <- function(theta, model_param = NULL) {
  if (is.null(model_param)) {
    model_param <- make_lv_smfsb_model_param()
  }

  theta <- as.matrix(theta)
  if (nrow(theta) != dim_param && ncol(theta) == dim_param) {
    theta <- t(theta)
  }
  if (nrow(theta) != dim_param) {
    stop("theta must be a 3 x n matrix, an n x 3 matrix, or a length-3 vector.")
  }

  n <- ncol(theta)
  y <- matrix(NA_real_, nrow = 2 * model_param$n_times, ncol = n)
  for (i in seq_len(n)) {
    latent <- simulate_lv_smfsb_latent(exp(theta[, i]), model_param)
    y[, i] <- if (is.null(latent)) rep(NaN, nrow(y)) else lv_flatten_observation(latent)
  }
  rownames(y) <- c(
    paste0("x1_t", seq_len(model_param$n_times)),
    paste0("x2_t", seq_len(model_param$n_times))
  )
  y
}
