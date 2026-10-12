# models/ecoli_st131/ecoli_st131_model_nature_param.R
# ST131 SIR-colonisation / BSI observation model.
# This file contains only model constants, priors, and the generative model.
# External data are read in the run scripts and passed through model_param.
#
# Nature-paper parameterisation:
#   theta = (tau, R0, Dt)
# where
#   tau = beta - gamma is the net transmission / early growth rate,
#   R0  = beta / gamma is the basic reproduction number,
#   Dt  is the observation delay in weeks.
#
# Priors:
#   tau ~ Uniform(0, 0.1),
#   R0  ~ Normal(0, 1) truncated to [1.01, 3],
#   Dt  ~ Uniform(0, additional_years * time_period).
#
# Transformation used by the simulator:
#   gamma = tau / (R0 - 1),
#   beta  = tau * R0 / (R0 - 1).

if (!requireNamespace("deSolve", quietly = TRUE)) {
  stop("Package 'deSolve' is required. Install it with install.packages('deSolve').")
}

# -----------------------------
# Default model settings
# -----------------------------

ecoli_settings <- list(
  clade = "A",
  simulator_model = "SIR",
  theta_bsi = 0.00001249048,
  theta_c = 1,
  pop_size = 1000000,
  n_years = 14,
  time_period = 52,
  additional_years = 0.5,
  I0 = 1,
  R_init = 0,
  return_abc_summaries = FALSE
)

num_param <- 3
dim_data <- if (isTRUE(ecoli_settings$return_abc_summaries)) 3 else ecoli_settings$n_years

# -----------------------------
# Model-parameter constructor
# -----------------------------

make_ecoli_model_param <- function(
    clade = ecoli_settings$clade,
    mu_OR,
    return_abc_summaries = ecoli_settings$return_abc_summaries,
    theta_bsi = ecoli_settings$theta_bsi,
    theta_c = ecoli_settings$theta_c,
    pop_size = ecoli_settings$pop_size,
    n_years = ecoli_settings$n_years,
    time_period = ecoli_settings$time_period,
    additional_years = NULL,
    I0 = ecoli_settings$I0,
    R_init = ecoli_settings$R_init) {

  if (missing(mu_OR) || is.null(mu_OR) || !is.finite(mu_OR)) {
    stop("make_ecoli_model_param() requires a finite clade-specific mu_OR value.")
  }

  if (is.null(additional_years)) {
    additional_years <- if (identical(clade, "B")) 6.5 else ecoli_settings$additional_years
  }

  list(
    clade = clade,
    mu_OR = as.numeric(mu_OR),
    return_abc_summaries = return_abc_summaries,
    theta_bsi = theta_bsi,
    theta_c = theta_c,
    pop_size = pop_size,
    n_years = n_years,
    time_period = time_period,
    additional_years = additional_years,
    I0 = I0,
    R_init = R_init
  )
}

get_ecoli_param <- function(model_param = NULL, name) {
  if (!is.null(model_param) && !is.null(model_param[[name]])) return(model_param[[name]])
  ecoli_settings[[name]]
}

# -----------------------------
# Prior on theta = (tau, R0, Dt)
# -----------------------------

truncnorm01_pdf <- function(x, lower = 1.01, upper = 3.0) {
  if (!is.finite(x) || x < lower || x > upper) return(0)
  dnorm(x, mean = 0, sd = 1) / (pnorm(upper, 0, 1) - pnorm(lower, 0, 1))
}

get_Dt_max <- function(model_param = NULL) {
  get_ecoli_param(model_param, "additional_years") * get_ecoli_param(model_param, "time_period")
}

sample_prior <- function() {
  tau <- runif(1, min = 0, max = 0.1)

  repeat {
    R0 <- rnorm(1, mean = 0, sd = 1)
    if (R0 >= 1.01 && R0 <= 3.0) break
  }

  Dt_max <- ecoli_settings$additional_years * ecoli_settings$time_period
  Dt <- runif(1, min = 0, max = Dt_max)

  c(tau = tau, R0 = R0, Dt = Dt)
}

prior_pdf <- function(theta) {
  theta <- as.numeric(theta)
  if (length(theta) != 3) return(0)

  tau <- theta[1]
  R0 <- theta[2]
  Dt <- theta[3]

  if (!all(is.finite(theta))) return(0)

  Dt_max <- ecoli_settings$additional_years * ecoli_settings$time_period

  if (tau <= 0 || tau > 0.1) return(0)
  if (R0 < 1.01 || R0 > 3.0) return(0)
  if (Dt < 0 || Dt > Dt_max) return(0)

  d_tau <- 1 / 0.1
  d_R0 <- truncnorm01_pdf(R0, lower = 1.01, upper = 3.0)
  d_Dt <- 1 / Dt_max

  d_tau * d_R0 * d_Dt
}

prior_logpdf <- function(theta) {
  p <- prior_pdf(theta)
  if (p <= 0) return(-Inf)
  log(p)
}

# -----------------------------
# Parameter transformations
# -----------------------------

ecoli_eta_to_bg <- function(theta) {
  theta <- as.numeric(theta)
  if (length(theta) != 3) stop("theta must be c(tau, R0, Dt).")

  tau <- theta[1]
  R0 <- theta[2]
  Dt <- theta[3]

  if (!all(is.finite(theta)) || tau <= 0 || R0 <= 1 || Dt < 0) {
    return(c(beta = NA_real_, gamma = NA_real_, Dt = Dt))
  }

  gamma <- tau / (R0 - 1)
  beta <- tau * R0 / (R0 - 1)

  c(beta = beta, gamma = gamma, Dt = Dt)
}

ecoli_bg_to_eta <- function(theta_bg) {
  theta_bg <- as.numeric(theta_bg)
  if (length(theta_bg) != 3) stop("theta_bg must be c(beta, gamma, Dt).")

  beta <- theta_bg[1]
  gamma <- theta_bg[2]
  Dt <- theta_bg[3]

  if (!all(is.finite(theta_bg)) || beta <= 0 || gamma <= 0) {
    return(c(tau = NA_real_, R0 = NA_real_, Dt = Dt))
  }

  c(tau = beta - gamma, R0 = beta / gamma, Dt = Dt)
}

# -----------------------------
# Latent SIR colonisation model
# -----------------------------

ecoli_sir_rhs <- function(t, state, parms) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  beta <- parms[1]
  gamma <- parms[2]
  N <- parms[3]

  dS <- -beta * S * I / N
  dI <- beta * S * I / N - gamma * I
  dR <- gamma * I

  list(c(dS, dI, dR))
}

simulate_ecoli_sir_ode <- function(beta, gamma, model_param = NULL) {
  N <- get_ecoli_param(model_param, "pop_size")
  I0 <- get_ecoli_param(model_param, "I0")
  R_init <- get_ecoli_param(model_param, "R_init")
  n_years <- get_ecoli_param(model_param, "n_years")
  time_period <- get_ecoli_param(model_param, "time_period")
  additional_years <- get_ecoli_param(model_param, "additional_years")

  n_weeks_total <- as.integer((n_years + additional_years) * time_period)
  times <- seq(0, n_weeks_total - 1, by = 1)

  state0 <- c(S = N - I0 - R_init, I = I0, R = R_init)
  out <- deSolve::ode(
    y = state0,
    times = times,
    func = ecoli_sir_rhs,
    parms = c(beta, gamma, N),
    method = "lsoda"
  )

  as.data.frame(out)
}

# -----------------------------
# Observation model and annual aggregation
# -----------------------------

colonisation_to_bsi <- function(I_weekly, mu_OR, model_param = NULL, poisson = TRUE) {
  N <- get_ecoli_param(model_param, "pop_size")
  theta_c <- get_ecoli_param(model_param, "theta_c")
  theta_bsi <- get_ecoli_param(model_param, "theta_bsi")

  theta_c_a_counts <- I_weekly
  theta_c_counts <- theta_c * N
  theta_bsi_counts <- N * theta_bsi
  theta_c_0_counts <- theta_c_counts - theta_c_a_counts
  rho_c_C_c <- mu_OR * theta_c_a_counts

  lambda <- 0.1 + theta_bsi_counts * rho_c_C_c / (theta_c_0_counts + rho_c_C_c)
  lambda <- pmax(lambda, 0)

  if (isTRUE(poisson)) {
    rpois(length(lambda), lambda = lambda)
  } else {
    lambda
  }
}

sync_timewindow_R <- function(bsi_weekly, Dt, model_param = NULL) {
  n_years <- get_ecoli_param(model_param, "n_years")
  time_period <- get_ecoli_param(model_param, "time_period")
  n_obs_weeks <- n_years * time_period

  start <- as.integer(floor(Dt)) + 1L
  idx <- start + seq_len(n_obs_weeks) - 1L

  if (max(idx) > length(bsi_weekly)) {
    stop("Delay Dt requests weeks beyond the simulated BSI trajectory.")
  }

  bsi_weekly[idx]
}

sum_over_bsi_R <- function(bsi_weekly, model_param = NULL) {
  n_years <- get_ecoli_param(model_param, "n_years")
  time_period <- get_ecoli_param(model_param, "time_period")

  out <- numeric(n_years)
  for (i in seq_len(n_years)) {
    start <- (i - 1L) * time_period + 1L
    end <- i * time_period
    out[i] <- sum(bsi_weekly[start:end])
  }

  out
}

ecoli_abc_summaries <- function(y) {
  y <- as.numeric(y)
  log_y <- log(y + 0.5)
  c(
    log_BSI_max = max(log_y),
    BSI_max_t = which.max(y) - 1L,
    log_BSI_t0 = log_y[1]
  )
}

simulate_ecoli_st131 <- function(theta, model_param = NULL, poisson = TRUE) {
  theta <- as.numeric(theta)
  if (length(theta) != 3) stop("theta must be c(tau, R0, Dt).")

  tau <- theta[1]
  R0 <- theta[2]
  Dt <- theta[3]

  out_len <- if (isTRUE(get_ecoli_param(model_param, "return_abc_summaries"))) 3 else get_ecoli_param(model_param, "n_years")
  if (!all(is.finite(theta)) || tau <= 0 || R0 <= 1 || Dt < 0) {
    return(rep(NA_real_, out_len))
  }

  bg <- ecoli_eta_to_bg(theta)
  beta <- bg["beta"]
  gamma <- bg["gamma"]

  if (!all(is.finite(c(beta, gamma)))) {
    return(rep(NA_real_, out_len))
  }

  mu_OR <- get_ecoli_param(model_param, "mu_OR")
  if (is.null(mu_OR) || !is.finite(mu_OR)) {
    stop("model_param$mu_OR must be provided by the run script.")
  }

  return_abc_summaries <- get_ecoli_param(model_param, "return_abc_summaries")

  ode_out <- simulate_ecoli_sir_ode(beta, gamma, model_param)
  bsi_weekly <- colonisation_to_bsi(ode_out$I, mu_OR = mu_OR, model_param = model_param, poisson = poisson)
  bsi_window <- sync_timewindow_R(bsi_weekly, Dt = Dt, model_param = model_param)
  y_yearly <- sum_over_bsi_R(bsi_window, model_param = model_param)

  if (isTRUE(return_abc_summaries)) return(ecoli_abc_summaries(y_yearly))

  as.numeric(y_yearly)
}

model <- function(theta, model_param = NULL) {
  simulate_ecoli_st131(theta, model_param = model_param, poisson = TRUE)
}

# -----------------------------
# Optional candidate features / augmentations in Nature coordinates
# -----------------------------

augment_ecoli_log_duration <- function(theta_mat) {
  theta_mat <- as.matrix(theta_mat)
  if (nrow(theta_mat) != 3 && ncol(theta_mat) == 3) theta_mat <- t(theta_mat)

  tau <- theta_mat[1, ]
  R0 <- theta_mat[2, ]
  log_duration <- log(R0 - 1) - log(tau)

  matrix(log_duration, nrow = 1, dimnames = list("log_duration", NULL))
}

augment_ecoli_duration <- function(theta_mat) {
  theta_mat <- as.matrix(theta_mat)
  if (nrow(theta_mat) != 3 && ncol(theta_mat) == 3) theta_mat <- t(theta_mat)

  tau <- theta_mat[1, ]
  R0 <- theta_mat[2, ]
  duration <- (R0 - 1) / tau

  matrix(duration, nrow = 1, dimnames = list("duration", NULL))
}

candidate_features_ecoli_nature <- function(theta_mat) {
  theta_mat <- as.matrix(theta_mat)
  if (nrow(theta_mat) != 3 && ncol(theta_mat) == 3) theta_mat <- t(theta_mat)

  tau <- theta_mat[1, ]
  R0 <- theta_mat[2, ]
  Dt <- theta_mat[3, ]
  gamma <- tau / (R0 - 1)
  beta <- tau * R0 / (R0 - 1)
  duration <- 1 / gamma

  rbind(
    tau = tau,
    R0 = R0,
    Dt = Dt,
    beta = beta,
    gamma = gamma,
    duration = duration,
    log_duration = log(duration),
    beta_sq = beta^2,
    gamma_sq = gamma^2,
    tau_sq = tau^2,
    R0_minus_1 = R0 - 1,
    log_R0 = log(R0),
    R0_times_Dt = R0 * Dt
  )
}
