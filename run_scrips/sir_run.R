# ============================================================
# SIR case study: RatioIS inference + localized SBC + plotting/PPC
# ============================================================

library(xLLiM)
library(mixtools)
library(mvtnorm)
library(rjson)
library(loo)
library(deSolve)

source("src/gllim_strict_em.R")
source("src/build_surrogates.R")
source("src/ratio_estimation.R")
source("src/local_sbc_corrected_posterior.R")
source("src/plot_ratio_results.R")
source("models/sir/sir_model.R")

# ------------------------------------------------------------
# Case-study settings
# ------------------------------------------------------------

model_name <- "sir"
model_param <- NULL
subfolder <- "ratio_correction_tworounds_Kstart10_50k"

# Simulator budget used for the SIR results in the paper:
# total number of model simulations is: num_priorpred_samples+num_focused_samples+num_refreshed_samples
num_priorpred_samples <- 20000  # N0
num_focused_samples <- 20000    # Na
num_refreshed_samples <- 10000  # Nb

# Fresh q2 proposal draws used for RatioIS correction.
num_final_samples <- 10000 # this does not require simulations from the generative model

# ------------------------------------------------------------
# User-exposed inference settings
# ------------------------------------------------------------

# Initial number of GLLiM mixture components and covariance structure.
K_start <- 10
cov_structure <- ""  # see the options in the xLLiM package

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

observation_indices <- 1:10  # specifying 1:10 will run inferences independently on each of the ten datasets

# ------------------------------------------------------------
# Run inference and localized SBC for each sbibm observation
# ------------------------------------------------------------

for (num_observation in observation_indices) {
  set.seed(num_observation)

  output_dir <- file.path(
    getwd(),
    "results",
    model_name,
    "ratio_correction",
    subfolder,
    paste0("obs", num_observation)
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  observedData <- unlist(
    read.csv(
      file.path(
        "models",
        model_name,
        paste0("num_observation_", num_observation),
        "observation.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # q0 -> q1 -> q2 surrogate construction
  # ----------------------------------------------------------

  surrogate_out <- build_surrogates(
    observedData = observedData,
    prior_pdf = prior_pdf,
    sample_prior = sample_prior,
    model = model,
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

  # ----------------------------------------------------------
  # RatioIS correction of fresh q2 proposal draws
  # ----------------------------------------------------------

  set.seed(num_observation)

  ratio_out <- ratio_corrected_sample(
    surrogate_out = surrogate_out,
    num_final_samples = num_final_samples,
    feature_args = feature_args
  )

  saveRDS(
    ratio_out,
    file = file.path(output_dir, "ratio_out.rds")
  )

  q2_csv <- as.data.frame(t(ratio_out$theta_prop))
  names(q2_csv) <- c("beta", "gamma")
  write.csv(
    q2_csv,
    file = file.path(output_dir, "proposal_draws_q2.csv"),
    row.names = FALSE
  )

  psis_csv <- as.data.frame(t(ratio_out$theta_resampled))
  names(psis_csv) <- c("beta", "gamma")
  write.csv(
    psis_csv,
    file = file.path(output_dir, "posterior_draws_corrected_psis.csv"),
    row.names = FALSE
  )

  # The clipped RatioIS posterior is the primary approximation
  # reported in the paper.
  clipped_csv <- as.data.frame(t(ratio_out$theta_resampled_trunc))
  names(clipped_csv) <- c("beta", "gamma")
  write.csv(
    clipped_csv,
    file = file.path(output_dir, "posterior_draws_corrected_clipped.csv"),
    row.names = FALSE
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
    file = file.path(output_dir, "diagnostics.csv"),
    row.names = FALSE
  )

  settings <- list(
    model_name = model_name,
    observation = num_observation,
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
      eps_defensive_stage2 = eps_defensive_stage2
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

  cat("\nSIR observation", num_observation, "completed.\n")
  cat("ESS raw     =", ratio_out$ESS, "\n")
  cat("Pareto-k    =", max(ratio_out$pareto_k), "\n")
  cat("ESS PSIS    =", ratio_out$ESS_psis, "\n")
  cat("ESS clipped =", ratio_out$ESS_trunc, "\n\n")

  # ----------------------------------------------------------
  # Localized SBC
  # ----------------------------------------------------------

  if (isTRUE(run_local_sbc_diagnostic)) {
    for (lambda_local in local_sbc_lambdas) {
      set.seed(
        100000 +
          1000 * num_observation +
          round(100 * lambda_local)
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
        model = model,
        model_param = model_param,
        prior_pdf = prior_pdf,
        observedData = observedData,
        N_local = local_sbc_N,
        N_local_check = local_sbc_N,
        M_is = 5000,
        lambda = lambda_local,
        cov_structure = cov_structure,
        param_names = c("beta", "gamma"),
        output_dir = local_sbc_dir
      )

      cat(
        "\nLocalized SBC completed for observation",
        num_observation,
        "with lambda =",
        lambda_local,
        "\n"
      )
      print(local_sbc_out$pit_summary_clipped)
    }
  }
}

# ============================================================
# Posterior comparison plots for one selected observation
# (plots are saved in the results folder)
# ============================================================

# Observation ID for the posterior-comparison and PPC plots below.
plot_observation <- 1  # can be any number in 1:10

#**********
#*YOU WILL HAVE TO ADJUST MANUALLY THE NUMBERS IN param_limits INSIDE plot_ratio_results DEPENDING ON THE OBSERVATION YOU ARE FITTING
#**********

num_observation <- plot_observation

output_dir <- file.path(
  getwd(),
  "results",
  model_name,
  "ratio_correction",
  subfolder,
  paste0("obs", num_observation)
)

data_dir <- file.path(
  getwd(),
  "models",
  model_name,
  paste0("num_observation_", num_observation)
)

surrogate_out <- readRDS(file.path(output_dir, "surrogate_out.rds"))
ratio_out <- readRDS(file.path(output_dir, "ratio_out.rds"))

reference_samples <- as.matrix(
  read.csv(file.path(data_dir, "reference_posterior_samples.csv"))
)

observedData <- unlist(
  read.csv(file.path(data_dir, "observation.csv"))
)

theta_true <- unlist(
  read.csv(file.path(data_dir, "true_parameters.csv"))
)

diagnostics <- read.csv(file.path(output_dir, "diagnostics.csv"))
q2_samp <- read.csv(file.path(output_dir, "proposal_draws_q2.csv"))
theta_resampled_clipped <- read.csv(
  file.path(output_dir, "posterior_draws_corrected_clipped.csv")
)

# Save the same posterior-comparison plots produced by the earlier
# SIR run script. Because plot_ratio_results() produces several plot
# types by default, this is a multi-page PDF.
pdf(
  file.path(output_dir, "posterior_comparison_reference_q2_ratiois.pdf"),
  width = 9,
  height = 7
)

#**********
#*YOU WILL HAVE TO ADJUST MANUALLY THE NUMBERS IN param_limits INSIDE plot_ratio_results DEPENDING ON THE OBSERVATION YOU ARE FITTING
# (plots are saved in the results folder)
#**********
plot_ratio_results(
  reference_samples = reference_samples,
  q2_samp = q2_samp,
  theta_resampled = theta_resampled_clipped,
  theta_true = theta_true,
  theta_labels = c("beta", "gamma"),
  plot_2d_density = TRUE,
  plot_q2 = TRUE,
  reference_label = "Reference",
  q2_label = "q2",
  corrected_label = "clipped ratioIS",
  print_labels = TRUE,
  param_limits = list(
    beta = c(0.55, 0.75),
    gamma = c(0.05, 0.25)
  ),
  title_cex = 1.2,
  axis_label_cex = 1.5,
  axis_tick_cex = 1.2,
  legend_cex = 1.5
)

dev.off()

# ============================================================
# Posterior predictive checks for SIR
# (plots are saved in the results folder)
# ============================================================
simulate_sir_ppc <- function(theta_post,
                             observedData,
                             model,
                             model_param = NULL,
                             n_rep = 1000,
                             seed = 1,
                             output_dir = NULL,
                             verbose = TRUE) {
  set.seed(seed)
  
  theta_post <- as.matrix(theta_post)
  
  # Accept either 2 x N or N x 2.
  if (nrow(theta_post) != 2 && ncol(theta_post) == 2) {
    theta_post <- t(theta_post)
  }
  if (nrow(theta_post) != 2) {
    stop("theta_post must be 2 x N or N x 2, with rows/columns beta,gamma.")
  }
  
  observedData <- as.numeric(observedData)
  n_post <- ncol(theta_post)
  
  idx <- sample(
    seq_len(n_post),
    size = n_rep,
    replace = n_rep > n_post
  )
  
  y_rep <- matrix(
    NA_real_,
    nrow = n_rep,
    ncol = length(observedData)
  )
  
  for (i in seq_len(n_rep)) {
    if (verbose && (i %% 100 == 0 || i == n_rep)) {
      cat("PPC simulation", i, "of", n_rep, "\n")
    }
    
    y_rep[i, ] <- model(
      theta = theta_post[, idx[i]],
      model_param = model_param
    )
  }
  
  if (exists("sir_settings", inherits = TRUE)) {
    obs_times <- sir_settings$times_full[sir_settings$obs_indices]
  } else {
    obs_times <- seq_along(observedData)
  }
  
  envelope <- data.frame(
    time = obs_times,
    observed = observedData,
    pred_median = apply(y_rep, 2, median),
    pred_q025 = apply(y_rep, 2, quantile, probs = 0.025),
    pred_q05 = apply(y_rep, 2, quantile, probs = 0.05),
    pred_q95 = apply(y_rep, 2, quantile, probs = 0.95),
    pred_q975 = apply(y_rep, 2, quantile, probs = 0.975)
  )
  
  out <- list(
    y_rep = y_rep,
    observedData = observedData,
    obs_times = obs_times,
    envelope = envelope
  )
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    write.csv(
      envelope,
      file = file.path(output_dir, "ppc_time_envelope.csv"),
      row.names = FALSE
    )
    
    saveRDS(
      out,
      file = file.path(output_dir, "ppc_out.rds")
    )
  }
  
  out
}


plot_sir_ppc <- function(ppc_out,
                         output_file = NULL,
                         main = "SIR posterior predictive check") {
  envelope <- ppc_out$envelope
  
  if (!is.null(output_file)) {
    pdf(output_file, width = 8, height = 5)
    on.exit(dev.off())
  }
  
  ylim <- range(
    c(
      envelope$observed,
      envelope$pred_q025,
      envelope$pred_q975
    ),
    finite = TRUE
  )
  
  plot(
    envelope$time,
    envelope$pred_median,
    type = "n",
    ylim = ylim,
    xlab = "time",
    ylab = "infected counts",
    main = main
  )
  
  polygon(
    c(envelope$time, rev(envelope$time)),
    c(envelope$pred_q975, rev(envelope$pred_q025)),
    border = NA,
    col = adjustcolor("gray80", alpha.f = 0.7)
  )
  
  polygon(
    c(envelope$time, rev(envelope$time)),
    c(envelope$pred_q95, rev(envelope$pred_q05)),
    border = NA,
    col = adjustcolor("gray60", alpha.f = 0.7)
  )
  
  lines(
    envelope$time,
    envelope$pred_median,
    lwd = 2
  )
  
  points(
    envelope$time,
    envelope$observed,
    pch = 19,
    cex = 1.1
  )
  
  legend(
    "topright",
    legend = c(
      "Observed",
      "Predictive median",
      "90% interval",
      "95% interval"
    ),
    pch = c(19, NA, 15, 15),
    lty = c(NA, 1, NA, NA),
    col = c(
      "black",
      "black",
      adjustcolor("gray60", alpha.f = 0.7),
      adjustcolor("gray80", alpha.f = 0.7)
    ),
    pt.cex = c(1.1, NA, 2, 2),
    bty = "n"
  )
}


# Clipped RatioIS posterior used for the PPC.
theta_post <- ratio_out$theta_resampled_trunc

ppc_out <- simulate_sir_ppc(
  theta_post = theta_post,
  observedData = observedData,
  model = model,
  model_param = model_param,
  n_rep = 1000,
  seed = 1,
  output_dir = file.path(output_dir, "ppc"),
  verbose = TRUE
)

plot_sir_ppc(
  ppc_out,
  output_file = file.path(
    output_dir,
    "ppc",
    "ppc_time_envelope.pdf"
  ),
  main = "SIR PPC: clipped RatioIS posterior"
)

