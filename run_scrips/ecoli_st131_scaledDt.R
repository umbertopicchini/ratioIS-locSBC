# ============================================================
# ST131 epidemiological application
# RatioIS inference + localized SBC + paper plotting/PPC
# ============================================================

library(xLLiM)
library(deSolve)
library(readxl)

source("src/gllim_strict_em.R")
source("src/build_surrogates.R")
source("src/ratio_estimation.R")
source("src/local_sbc_corrected_posterior.R")
source("models/ecoli_st131/ecoli_st131_model.R")

# ------------------------------------------------------------
# Data and model
# ------------------------------------------------------------

model_name <- "ecoli_st131"
clade_to_fit <- "A"  # main paper case; optional: "C1", "C2"
subfolder <- paste0(
  "ratio_correction_ScaledDt_Kstart10_100k_clade_",
  clade_to_fit
)

# The paper analysis uses the complete annual BSI vector.
return_abc_summaries <- FALSE

bsi_data_file <- file.path(
  "models", model_name, "data", "NORM_incidence_updated.csv"
)
or_data_file <- file.path(
  "models", model_name, "data",
  "ST131_clades_OR_E_coli_carriage_disease_collapsed_wBSAC2.xlsx"
)
OR_type <- "mean"

ecoli_settings$clade <- clade_to_fit
ecoli_settings$return_abc_summaries <- return_abc_summaries

# Read external data once, in the run script.
bsi_dat <- read.csv(bsi_data_file, check.names = FALSE)
if (!("Year" %in% names(bsi_dat))) {
  stop("Column 'Year' not found in BSI data file.")
}
if (!(clade_to_fit %in% names(bsi_dat))) {
  stop("Clade ", clade_to_fit, " not found in BSI data file.")
}

or_dat <- readxl::read_excel(or_data_file)
required_or_cols <- c("Label", "OR", "lower", "upper")
missing_or_cols <- setdiff(required_or_cols, names(or_dat))
if (length(missing_or_cols) > 0) {
  stop(
    "Missing columns in OR data file: ",
    paste(missing_or_cols, collapse = ", ")
  )
}

or_label <- paste0(clade_to_fit, " (BSAC2)")
or_row <- or_dat[or_dat$Label == or_label, , drop = FALSE]
if (nrow(or_row) != 1) {
  stop("Could not find a unique OR row for ", or_label)
}

mu_OR <- switch(
  OR_type,
  mean = as.numeric(or_row$OR),
  lower = as.numeric(or_row$lower),
  upper = as.numeric(or_row$upper),
  stop("OR_type must be one of 'mean', 'lower', or 'upper'.")
)
if (!is.finite(mu_OR)) stop("Extracted mu_OR is not finite.")

model_param <- make_ecoli_model_param(
  clade = clade_to_fit,
  mu_OR = mu_OR,
  return_abc_summaries = return_abc_summaries
)

# Observed annual BSI incidence series.
observedData_full <- as.numeric(bsi_dat[[clade_to_fit]])
observedData_full <- observedData_full[seq_len(model_param$n_years)]
years_full <- bsi_dat$Year[seq_len(model_param$n_years)]

keep_idx <- seq_along(observedData_full)
if (clade_to_fit == "C2") {
  keep_idx <- which(years_full >= 2007)
}

observedData <- observedData_full[keep_idx]
years_obs <- years_full[keep_idx]

# Restrict simulator output to exactly the observation window used for fitting.
model_full <- model
model <- function(theta, model_param) {
  y_full <- model_full(theta, model_param)
  y_full[keep_idx]
}

# ------------------------------------------------------------
# Ojala parameter scaling used for surrogate fitting
# ------------------------------------------------------------
# Inference coordinates are eta = (tau, R0, Dt_scaled), with
# Dt_scaled = 0.001 * Dt and Dt measured in weeks.

Dt_fit_factor <- 0.001
Dt_unscale_factor <- 1 / Dt_fit_factor
Dt_max <- get_Dt_max(model_param)

to_raw_theta <- function(eta) {
  eta <- as.numeric(eta)
  if (length(eta) != 3) {
    stop("eta must be c(tau, R0, Dt_scaled).")
  }

  c(
    tau = eta[1],
    R0 = eta[2],
    Dt = eta[3] * Dt_unscale_factor
  )
}

from_raw_theta <- function(theta) {
  theta <- as.numeric(theta)
  if (length(theta) != 3) {
    stop("theta must be c(tau, R0, Dt).")
  }

  c(
    tau = theta[1],
    R0 = theta[2],
    Dt_scaled = theta[3] * Dt_fit_factor
  )
}

sample_prior_eta <- function() {
  from_raw_theta(sample_prior())
}

prior_pdf_eta <- function(eta) {
  theta_raw <- to_raw_theta(eta)

  # Change of variables: Dt = 1000 * Dt_scaled.
  prior_pdf(theta_raw) * Dt_unscale_factor
}

model_eta <- function(eta, model_param) {
  model(to_raw_theta(eta), model_param)
}

param_names_eta <- c("tau", "R0", "Dt_scaled")
dim_data <- length(observedData)

# ------------------------------------------------------------
# Simulation budget used for the paper analysis
# ------------------------------------------------------------
# total number of model simulations is: num_priorpred_samples+num_focused_samples+num_refreshed_samples
num_priorpred_samples <- 50000 # N0
num_focused_samples <- 30000   # Na
num_refreshed_samples <- 20000 # Nb

# Fresh q2 proposal draws used for RatioIS correction.
num_final_samples <- 10000  # this does not require simulations from the generative model

# ------------------------------------------------------------
# User-exposed inference settings
# ------------------------------------------------------------
# Initial number of GLLiM mixture components and covariance structure
K_start <- 10
cov_structure <- ""

# Prior-defensive fractions in the two localization stages.
# These correspond to epsilon_1 and epsilon_2 in the paper.
eps_defensive_stage1 <- 0.05
eps_defensive_stage2 <- 0.05

# Ratio-classifier feature map. The baseline specification contains
# linear terms and all parameter-summary cross interactions.
feature_args <- list(
  include_linear = TRUE,
  include_quadratic = FALSE,
  include_within_block_interactions = FALSE,
  include_cross_interactions = TRUE
)

# Localized-SBC settings.
run_local_sbc_diagnostic <- TRUE
local_sbc_lambdas <- c(1, 1.05, 1.1, 1.25, 1.5, 2)
local_sbc_N <- 2000

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

output_dir <- file.path(
  getwd(),
  "results",
  model_name,
  "ratio_correction",
  subfolder
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(t(observedData)),
  file = file.path(output_dir, "observedData.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# q0 -> q1 -> q2 surrogate construction
# ------------------------------------------------------------

set.seed(123)

surrogate_out <- build_surrogates(
  observedData = observedData,
  prior_pdf = prior_pdf_eta,
  sample_prior = sample_prior_eta,
  model = model_eta,
  model_param = model_param,
  dim_data = dim_data,
  num_priorpred_samples = num_priorpred_samples,
  num_focused_samples = num_focused_samples,
  num_refreshed_samples = num_refreshed_samples,
  K_start = K_start,
  cov_structure = cov_structure,
  eps_defensive_stage1 = eps_defensive_stage1,
  eps_defensive_stage2 = eps_defensive_stage2,
  output_dir = output_dir
)

saveRDS(
  surrogate_out,
  file = file.path(output_dir, "surrogate_out.rds")
)

# ------------------------------------------------------------
# RatioIS correction of fresh q2 proposal draws
# ------------------------------------------------------------

set.seed(123)

ratio_out <- ratio_corrected_sample(
  surrogate_out = surrogate_out,
  num_final_samples = num_final_samples,
  feature_args = feature_args
)

saveRDS(
  ratio_out,
  file = file.path(output_dir, "ratio_out.rds")
)

write_theta_eta_csv <- function(theta_eta_mat, file) {
  out <- as.data.frame(t(theta_eta_mat))
  names(out) <- c("tau", "R0", "Dt_scaled")
  out$Dt <- out$Dt_scaled * Dt_unscale_factor
  write.csv(out, file = file, row.names = FALSE)
}

write_theta_eta_csv(
  ratio_out$theta_prop,
  file.path(output_dir, "proposal_draws_q2.csv")
)

write_theta_eta_csv(
  ratio_out$theta_resampled,
  file.path(output_dir, "posterior_draws_corrected_psis.csv")
)

# Clipped RatioIS is the primary posterior approximation.
write_theta_eta_csv(
  ratio_out$theta_resampled_trunc,
  file.path(output_dir, "posterior_draws_corrected_clipped.csv")
)

diagnostics <- data.frame(
  ESS_raw = ratio_out$ESS,
  pareto_k = max(ratio_out$pareto_k),
  ESS_psis = ratio_out$ESS_psis,
  ESS_clipped = ratio_out$ESS_trunc,
  K_original = ratio_out$K_original,
  K_used = ratio_out$K_used
)

write.csv(
  diagnostics,
  file.path(output_dir, "diagnostics.csv"),
  row.names = FALSE
)

settings <- list(
  model_name = model_name,
  clade = clade_to_fit,
  OR_type = OR_type,
  parameterisation =
    "(tau, R0, Dt_scaled), with Dt_scaled = 0.001 * Dt",
  simulation_budget = list(
    N0 = num_priorpred_samples,
    Na = num_focused_samples,
    Nb = num_refreshed_samples,
    M = num_final_samples
  ),
  inference_settings = list(
    K_start = K_start,
    cov_structure = cov_structure,
    eps_defensive_stage1 = eps_defensive_stage1,
    eps_defensive_stage2 = eps_defensive_stage2,
    feature_args = feature_args
  ),
  localized_sbc = list(
    lambdas = local_sbc_lambdas,
    N_local = local_sbc_N,
    N_check = local_sbc_N,
    M_is = 5000
  ),
  surrogate_config = surrogate_out$summary_config,
  ratio_config = ratio_out$ratio_config
)

write(
  rjson::toJSON(settings),
  file = file.path(output_dir, "settings.json")
)

cat("\nST131 inference completed for clade", clade_to_fit, "\n")
cat("ESS raw     =", ratio_out$ESS, "\n")
cat("Pareto-k    =", max(ratio_out$pareto_k), "\n")
cat("ESS PSIS    =", ratio_out$ESS_psis, "\n")
cat("ESS clipped =", ratio_out$ESS_trunc, "\n\n")

# ------------------------------------------------------------
# Localized SBC
# ------------------------------------------------------------

if (isTRUE(run_local_sbc_diagnostic)) {
  seed_clade <- switch(
    clade_to_fit,
    A = 1,
    C1 = 2,
    C2 = 3,
    99
  )

  for (lambda_local in local_sbc_lambdas) {

    # Seed convention used for the reported ST131 localized-SBC runs.
    set.seed(
      100000 +
        10000 * seed_clade +
        round(1000 * lambda_local)
    )

    local_sbc_dir <- file.path(
      output_dir,
      paste0(
        "local_sbc_refit_ratiois_lambda",
        gsub("\\.", "p", as.character(lambda_local))
      )
    )

    local_sbc_out <- run_local_sbc(
      surrogate_out = surrogate_out,
      ratio_out = ratio_out,
      model = model_eta,
      model_param = model_param,
      prior_pdf = prior_pdf_eta,
      observedData = observedData,
      N_local = local_sbc_N,
      N_local_check = local_sbc_N,
      M_is = 5000,
      lambda = lambda_local,
      weight_method = "clipped",
      postprob_eps = 1e-6,
      cov_structure = cov_structure,
      param_names = param_names_eta,
      output_dir = local_sbc_dir
    )

    cat(
      "\nLocalized SBC completed for clade",
      clade_to_fit,
      "with lambda =",
      lambda_local,
      "\n"
    )
    print(local_sbc_out$pit_summary_clipped)
  }
}

#::::::: PLOTTING + PPC : OJALA PARAMETRIZATION :::::::::::::::::::::::::

output_dir <- file.path(
  getwd(),
  "results",
  model_name,
  "ratio_correction",
  subfolder
)

data_dir <- file.path(
  getwd(),
  "models",
  model_name,
  "data"
)

# ------------------------------------------------------------
# Load saved objects
# ------------------------------------------------------------

surrogate_out <- readRDS(
  file.path(output_dir, "surrogate_out.rds")
)

ratio_out <- readRDS(
  file.path(output_dir, "ratio_out.rds")
)

diagnostics <- read.csv(
  file.path(output_dir, "diagnostics.csv")
)


theta_resampled_weightstruncated <- read.csv(
  file.path(output_dir, "posterior_draws_corrected_clipped.csv")
)

q2_samp <- read.csv(
  file.path(output_dir, "proposal_draws_q2.csv")
)


# ------------------------------------------------------------
# Parameter names
#
# Second parametrization:
#   theta = (tau, R0, Dt_scaled)
#   Dt_scaled = 0.001 * Dt
#   Dt = 1000 * Dt_scaled
# ------------------------------------------------------------


names(theta_resampled_weightstruncated)[1:3] <-
  c("tau", "R0", "Dt_scaled")

names(q2_samp)[1:3] <-
  c("tau", "R0", "Dt_scaled")


# ------------------------------------------------------------
# Convert delay back to weeks for interpretation and plotting
# ------------------------------------------------------------

Dt_unscale_factor <- 1000

theta_resampled_weightstruncated$Dt <-
  theta_resampled_weightstruncated$Dt_scaled * Dt_unscale_factor

q2_samp$Dt <-
  q2_samp$Dt_scaled * Dt_unscale_factor


# ============================================================
# 1. MARGINAL POSTERIOR DENSITY OF R0
# ============================================================

R0_q2 <- q2_samp[, "R0"]

R0_post <-
  theta_resampled_weightstruncated[, "R0"]

dens_R0_q2 <- density(R0_q2)

dens_R0_post <- density(R0_post)

ylim_R0_density <- range(
  c(
    dens_R0_q2$y,
    dens_R0_post$y
  ),
  finite = TRUE
)

textmain <- bquote(
  "Posterior density of " *
    R[0] *
    " for ST131-" *
    .(clade_to_fit)
)

plot(
  dens_R0_q2,
  xlab = expression(R[0]),
  ylab = "Density",
  main = textmain,
  lwd = 2,
  lty = 2,
  ylim = ylim_R0_density
)

lines(
  dens_R0_post,
  lwd = 2,
  lty = 1
)

abline(
  v = median(R0_q2),
  lty = 2
)

abline(
  v = median(R0_post),
  lty = 1
)

legend(
  "topright",
  legend = c(
    "Final q2 surrogate",
    "Clipped ratio-corrected"
  ),
  lty = c(2, 1),
  lwd = 2,
  bty = "n"
)


# ============================================================
# 2. THREE PAIRWISE POSTERIOR CONTOUR PANELS
#
# Display Dt on the original scientific scale in weeks,
# even though inference was performed using Dt_scaled.
# ============================================================

tau_q2 <- q2_samp[, "tau"]
R0_q2  <- q2_samp[, "R0"]
Dt_q2  <- q2_samp[, "Dt"]

tau_post <-
  theta_resampled_weightstruncated[, "tau"]

R0_post <-
  theta_resampled_weightstruncated[, "R0"]

Dt_post <-
  theta_resampled_weightstruncated[, "Dt"]


# ------------------------------------------------------------
# Helper function for one pairwise panel
# ------------------------------------------------------------

plot_pairwise_contour <- function(
    x_q2,
    y_q2,
    x_post,
    y_post,
    xlab_expr,
    ylab_expr,
    panel_title = NULL,
    n_grid = 200,
    legend_position = NULL,
    pad_fraction = 0.3,
    xlim_fixed = NULL,
    ylim_fixed = NULL
) {
  
  if (is.null(xlim_fixed)) {
    xlim <- range(c(x_q2, x_post), finite = TRUE)
    
    xpad <- pad_fraction * diff(xlim)
    
    if (!is.finite(xpad) || xpad == 0) {
      xpad <- 0.05 * max(abs(xlim), 1)
    }
    
    xlim <- xlim + c(-xpad, xpad)
    
  } else {
    xlim <- xlim_fixed
  }
  
  
  if (is.null(ylim_fixed)) {
    ylim <- range(c(y_q2, y_post), finite = TRUE)
    
    ypad <- pad_fraction * diff(ylim)
    
    if (!is.finite(ypad) || ypad == 0) {
      ypad <- 0.05 * max(abs(ylim), 1)
    }
    
    ylim <- ylim + c(-ypad, ypad)
    
  } else {
    ylim <- ylim_fixed
  }
  
  # q2 density
  dens_q2 <- MASS::kde2d(
    x = x_q2,
    y = y_q2,
    n = n_grid,
    lims = c(
      xlim,
      ylim
    )
  )
  
  # clipped ratio-corrected density
  dens_post <- MASS::kde2d(
    x = x_post,
    y = y_post,
    n = n_grid,
    lims = c(
      xlim,
      ylim
    )
  )
  
  plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab_expr,
    ylab = ylab_expr,
    main = panel_title,
    cex.axis = 1.3,
    cex.lab = 1.5,
    cex.main = 1.4
  )
  
  contour(
    dens_q2$x,
    dens_q2$y,
    dens_q2$z,
    add = TRUE,
    drawlabels = FALSE,
    col = "darkorange2",
    lwd = 2
  )
  
  contour(
    dens_post$x,
    dens_post$y,
    dens_post$z,
    add = TRUE,
    drawlabels = FALSE,
    col = "dodgerblue3",
    lwd = 2
  )
  
  if (!is.null(legend_position)) {
    
    legend(
      legend_position,
      legend = c(
        "q2 surrogate",
        "Clipped ratio-corrected"
      ),
      col = c(
        "darkorange2",
        "dodgerblue3"
      ),
      lty = c(1, 1),
      lwd = c(2, 2),
      bty = "n"
    )
  }
}


# ------------------------------------------------------------
# Draw three panels
# ------------------------------------------------------------

old_par <- par(no.readonly = TRUE)

par(
  mfrow = c(1, 3),
  mar = c(5, 5, 4, 2.5),
  oma = c(0, 0, 3.5, 0)
)

# tau versus R0
plot_pairwise_contour(
  x_q2 = tau_q2,
  y_q2 = R0_q2,
  x_post = tau_post,
  y_post = R0_post,
  xlab_expr = expression(tau),
  ylab_expr = expression(R[0]),
  panel_title = expression(
    tau ~ "vs" ~ R[0]
  )
)

# tau versus Dt
plot_pairwise_contour(
  x_q2 = tau_q2,
  y_q2 = Dt_q2,
  x_post = tau_post,
  y_post = Dt_post,
  xlab_expr = expression(tau),
  ylab_expr = expression(
    D[t] ~ "(weeks)"
  ),
  panel_title = expression(
    tau ~ "vs" ~ D[t]
  ),
  ylim_fixed = c(-3, 30)
)

# R0 versus Dt
plot_pairwise_contour(
  x_q2 = R0_q2,
  y_q2 = Dt_q2,
  x_post = R0_post,
  y_post = Dt_post,
  xlab_expr = expression(R[0]),
  ylab_expr = expression(
    D[t] ~ "(weeks)"
  ),
  panel_title = expression(
    R[0] ~ "vs" ~ D[t]
  ),
  legend_position = "topright",
  ylim_fixed = c(-3, 30)
)

mtext(
  paste0(
    "ST131-",
    clade_to_fit,
    ": q2 surrogate versus clipped ratio-corrected posterior"
  ),
  outer = TRUE,
  side = 3,
  line = 1,
  font = 2
)

par(old_par)


# ============================================================
# 3. POSTERIOR PREDICTIVE CHECK
#
# IMPORTANT:
# inference was performed using Dt_scaled,
# but the original model() expects Dt in weeks.
# ============================================================

theta_ppc_df <-
  theta_resampled_weightstruncated

# Explicit back-transformation
theta_ppc_df$Dt <-
  theta_ppc_df$Dt_scaled *
  Dt_unscale_factor

# Original simulator receives:
# theta = (tau, R0, Dt),
# with Dt measured in weeks
theta_post <- as.matrix(
  theta_ppc_df[
    ,
    c(
      "tau",
      "R0",
      "Dt"
    )
  ]
)

theta_post <- t(theta_post)


# ------------------------------------------------------------
# Choose posterior draws for PPC
# ------------------------------------------------------------

n_ppc <- min(
  10000,
  ncol(theta_post)
)

set.seed(1)

theta_ppc <- theta_post[
  ,
  sample(
    seq_len(ncol(theta_post)),
    n_ppc
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# Generate posterior predictive datasets
# ------------------------------------------------------------

y_rep <- apply(
  theta_ppc,
  2,
  model,
  model_param = model_param
)

y_rep <- as.matrix(y_rep)

# Ensure:
# rows    = observation times
# columns = posterior predictive draws
if (
  nrow(y_rep) != length(observedData) &&
  ncol(y_rep) == length(observedData)
) {
  y_rep <- t(y_rep)
}


# ------------------------------------------------------------
# Pointwise predictive quantiles
# ------------------------------------------------------------

ppc_q <- apply(
  y_rep,
  1,
  quantile,
  probs = c(
    0.025,
    0.25,
    0.5,
    0.75,
    0.975
  ),
  na.rm = TRUE
)

years <- years_obs

textmain <- paste0(
  "Posterior predictive check: ST131-",
  clade_to_fit
)


# ------------------------------------------------------------
# Plot PPC
# ------------------------------------------------------------

plot(
  years,
  observedData,
  type = "b",
  pch = 19,
  ylim = range(
    c(
      observedData,
      ppc_q
    ),
    finite = TRUE
  ),
  xlab = "Year",
  ylab = "Annual BSI incidence",
  main = textmain
)

# 95% predictive interval
polygon(
  c(
    years,
    rev(years)
  ),
  c(
    ppc_q[1, ],
    rev(ppc_q[5, ])
  ),
  col = adjustcolor(
    "grey70",
    alpha.f = 0.4
  ),
  border = NA
)

# 50% predictive interval
polygon(
  c(
    years,
    rev(years)
  ),
  c(
    ppc_q[2, ],
    rev(ppc_q[4, ])
  ),
  col = adjustcolor(
    "grey50",
    alpha.f = 0.4
  ),
  border = NA
)

# Posterior predictive median
lines(
  years,
  ppc_q[3, ],
  lwd = 2
)

# Observed data
points(
  years,
  observedData,
  pch = 19
)

lines(
  years,
  observedData,
  lty = 2
)

legend(
  "topleft",
  legend = c(
    "Observed",
    "Posterior median",
    "50% band",
    "95% band"
  ),
  lty = c(
    2,
    1,
    NA,
    NA
  ),
  lwd = c(
    1,
    2,
    NA,
    NA
  ),
  pch = c(
    19,
    NA,
    15,
    15
  ),
  pt.cex = c(
    1,
    NA,
    2,
    2
  ),
  col = c(
    "black",
    "black",
    "grey45",
    "grey80"
  ),
  bty = "n"
)
