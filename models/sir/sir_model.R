# models/sir/sir_model.R

sir_settings <- list(
  N = 1000000,
  I0 = 1,
  R_init = 0,
  days = 160,
  saveat = 1,
  total_count = 1000,
  times_full = seq(0, 160, by = 1),
  obs_indices = seq(1, 161, by = 17)  # R indices corresponding to Python ::17
)

dim_data <- 10
num_param <- 2

sample_prior <- function() {
  beta <- rlnorm(1, meanlog = log(0.4), sdlog = 0.5)
  gamma <- rlnorm(1, meanlog = log(0.125), sdlog = 0.2)
  
  c(beta, gamma)
}

prior_pdf <- function(theta) {
  theta <- as.numeric(theta)
  
  if (length(theta) != 2) return(0)
  
  beta <- theta[1]
  gamma <- theta[2]
  
  if (!is.finite(beta) || !is.finite(gamma) || beta <= 0 || gamma <= 0) {
    return(0)
  }
  
  dlnorm(beta, meanlog = log(0.4), sdlog = 0.5) *
    dlnorm(gamma, meanlog = log(0.125), sdlog = 0.2)
}

sir_rhs <- function(t, state, parms) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parms[1]
  gamma <- parms[2]
  N <- sir_settings$N
  
  dS <- -beta * S * I / N
  dI <- beta * S * I / N - gamma * I
  dR <- gamma * I
  
  list(c(dS, dI, dR))
}

simulate_sir_ode <- function(theta) {
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("Package 'deSolve' is required. Install it with install.packages('deSolve').")
  }
  
  theta <- as.numeric(theta)
  
  beta <- theta[1]
  gamma <- theta[2]
  
  N <- sir_settings$N
  I0 <- sir_settings$I0
  R_init <- sir_settings$R_init
  S0 <- N - I0 - R_init
  
  state0 <- c(S = S0, I = I0, R = R_init)
  
  out <- deSolve::ode(
    y = state0,
    times = sir_settings$times_full,
    func = sir_rhs,
    parms = c(beta, gamma),
    method = "lsoda"
  )
  
  out <- as.data.frame(out)
  
  I_full <- out$I
  
  I_sub <- I_full[sir_settings$obs_indices]
  
  p <- pmin(pmax(I_sub / N, 0), 1)
  
  y <- rbinom(
    n = length(p),
    size = sir_settings$total_count,
    prob = p
  )
  
  as.numeric(y)
}

model <- function(theta, model_param = NULL) {
  simulate_sir_ode(theta)
}