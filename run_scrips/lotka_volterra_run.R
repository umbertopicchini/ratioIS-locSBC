# ============================================================
# Stochastic Lotka--Volterra case study
# RatioIS inference + localized SBC + posterior plots + PPC
# ============================================================

library(xLLiM)
library(smfsb)

source("src/gllim_strict_em.R")
source("src/build_surrogates.R")
source("src/ratio_estimation.R")
source("src/local_sbc_corrected_posterior.R")
source("models/lotka_volterra_smfsb/lotka_volterra_smfsb_model.R")

# ------------------------------------------------------------
# Data and model
# ------------------------------------------------------------

model_name <- "lotka_volterra_smfsb"
subfolder <- "ratio_correction_LVnoise10_Kstart10_50k"

# Built-in smfsb dataset used in the paper.
data(LVdata, package = "smfsb")
observed_object_name <- "LVnoise10"
observed_object <- get(observed_object_name)

observed_mat <- lv_observed_matrix(observed_object)
observed_times <- lv_extract_times(observed_object)
observedData <- lv_flatten_observation(observed_mat)

model_param <- make_lv_smfsb_model_param(
  observed_object = observed_object,
  x0 = c(x1 = 50, x2 = 100),
  noise_sd = c(x1 = 10, x2 = 10),
  prior_logmean = log(c(th1 = 1, th2 = 0.005, th3 = 0.6)),
  prior_logsd = c(th1 = 0.5, th2 = 0.5, th3 = 0.5),
  use_compiled_step = TRUE,
  max_attempts = 5
)

dim_data <- length(observedData)
stopifnot(dim_data == 2 * model_param$n_times)

prior_pdf_run <- function(phi) {
  prior_pdf(phi, model_param = model_param)
}

sample_prior_run <- function(n = NULL) {
  sample_prior(n = n, model_param = model_param)
}

# ------------------------------------------------------------
# Simulation budget used for the LV results in the paper
# ------------------------------------------------------------
# total number of model simulations is: num_priorpred_samples+num_focused_samples+num_refreshed_samples
num_priorpred_samples <- 20000  # N0
num_focused_samples <- 20000    # Na
num_refreshed_samples <- 10000  # Nb

# Fresh q2 proposal draws used for RatioIS correction.
num_final_samples <- 10000 # this does not require simulations from the generative model

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
  data.frame(time = observed_times, observed_mat),
  file = file.path(output_dir, "observedData_timeseries.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# q0 -> q1 -> q2 surrogate construction
# ------------------------------------------------------------

set.seed(123)

surrogate_out <- build_surrogates(
  observedData = observedData,
  prior_pdf = prior_pdf_run,
  sample_prior = sample_prior_run,
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

# Save inference-scale and natural-scale posterior draws.
write_theta_csv <- function(phi_mat, file) {
  out <- as.data.frame(t(phi_mat))
  names(out) <- phi_labels

  out$th1 <- exp(out$log_th1)
  out$th2 <- exp(out$log_th2)
  out$th3 <- exp(out$log_th3)

  out$x1_star <- out$th3 / out$th2
  out$x2_star <- out$th1 / out$th2

  write.csv(out, file = file, row.names = FALSE)
}

write_theta_csv(
  ratio_out$theta_prop,
  file.path(output_dir, "proposal_draws_q2.csv")
)

write_theta_csv(
  ratio_out$theta_resampled,
  file.path(output_dir, "posterior_draws_corrected_psis.csv")
)

# Clipped RatioIS is the primary posterior approximation.
write_theta_csv(
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
  file = file.path(output_dir, "diagnostics.csv"),
  row.names = FALSE
)

settings <- list(
  model_name = model_name,
  observed_object = observed_object_name,
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

cat("\nStochastic LV inference completed.\n")
cat("ESS raw     =", ratio_out$ESS, "\n")
cat("Pareto-k    =", max(ratio_out$pareto_k), "\n")
cat("ESS PSIS    =", ratio_out$ESS_psis, "\n")
cat("ESS clipped =", ratio_out$ESS_trunc, "\n\n")

# ------------------------------------------------------------
# Localized SBC
# ------------------------------------------------------------

if (isTRUE(run_local_sbc_diagnostic)) {
  for (lambda_local in local_sbc_lambdas) {

    # Seed convention used for the reported LV localized-SBC runs.
    set.seed(
      100000 +
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
      model = model,
      model_param = model_param,
      prior_pdf = prior_pdf_run,
      observedData = observedData,
      N_local = local_sbc_N,
      N_local_check = local_sbc_N,
      M_is = 5000,
      lambda = lambda_local,
      cov_structure = cov_structure,
      param_names = phi_labels,
      output_dir = local_sbc_dir
    )

    cat(
      "\nLocalized SBC completed for LV with lambda =",
      lambda_local,
      "\n"
    )
    print(local_sbc_out$pit_summary_clipped)
  }
}

# ============================================================
# Posterior marginal comparisons with PMCMC
# ============================================================

pmmc_file <- file.path(
  getwd(),
  "results",
  model_name,
  "PMCMC",
  "pmmc_lvnoise10_fixed_sd10.rds"
)

if (file.exists(pmmc_file)) {
  pmmc <- as.data.frame(readRDS(pmmc_file))
  pmmc <- pmmc[(20000 + 1):nrow(pmmc), c("th1", "th2", "th3"), drop = FALSE]
  
  q2_draws <- read.csv(
    file.path(output_dir, "proposal_draws_q2.csv")
  )
  
  ratiois_draws <- read.csv(
    file.path(output_dir, "posterior_draws_corrected_clipped.csv")
  )
  
  xvar <- "th1"
  yvar <- "th2"
  
  xlim <- range(
    pmmc[[xvar]],
    q2_draws[[xvar]],
    ratiois_draws[[xvar]],
    finite = TRUE
  )
  
  ylim <- range(
    pmmc[[yvar]],
    q2_draws[[yvar]],
    ratiois_draws[[yvar]],
    finite = TRUE
  )
  
  ref_kde <- MASS::kde2d(
    pmmc[[xvar]],
    pmmc[[yvar]],
    n = 100,
    lims = c(xlim, ylim)
  )
  
  plot_density_shaded <- function(x,
                                  y,
                                  fill_col,
                                  panel_title,
                                  ref_kde = NULL) {
    
    k <- MASS::kde2d(
      x,
      y,
      n = 100,
      lims = c(xlim, ylim)
    )
    
    z <- k$z
    
    lev <- seq(
      min(z),
      max(z),
      length.out = 9
    )
    
    plot(
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = expression(theta[1]),
      ylab = expression(theta[2]),
      main = panel_title
    )
    
    image(
      k$x,
      k$y,
      z,
      col = colorRampPalette(
        c("white", fill_col)
      )(length(lev) - 1),
      breaks = lev,
      add = TRUE,
      useRaster = TRUE
    )
    
    if (!is.null(ref_kde)) {
      contour(
        ref_kde,
        add = TRUE,
        nlevels = 5,
        drawlabels = FALSE,
        lwd = 2
      )
    }
  }
  
  pdf(
    file.path(output_dir, "lv_posterior_2d.pdf"),
    width = 12,
    height = 4.5
  )
  
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.5, 1))
  
  plot_density_shaded(
    pmmc[[xvar]],
    pmmc[[yvar]],
    fill_col = "grey70",
    panel_title = "PMCMC"
  )
  
  plot_density_shaded(
    q2_draws[[xvar]],
    q2_draws[[yvar]],
    fill_col = "gold",
    panel_title = expression(q[2]),
    ref_kde = ref_kde
  )
  
  plot_density_shaded(
    ratiois_draws[[xvar]],
    ratiois_draws[[yvar]],
    fill_col = "skyblue",
    panel_title = "clipped RatioIS",
    ref_kde = ref_kde
  )
  
  par(old_par)
  dev.off()
} else {
  warning("PMCMC reference file not found.")
}

# ============================================================
# Posterior predictive time envelopes
# ============================================================

simulate_lv_ppc <- function(theta_post,
                            observed_mat,
                            observed_times,
                            model,
                            model_param,
                            n_rep = 1000,
                            seed = 1,
                            output_dir = NULL,
                            verbose = TRUE) {
  set.seed(seed)
  
  theta_post <- as.matrix(theta_post)
  if (nrow(theta_post) != 3 && ncol(theta_post) == 3) {
    theta_post <- t(theta_post)
  }
  if (nrow(theta_post) != 3) {
    stop("theta_post must be 3 x N or N x 3 on the log-rate scale.")
  }
  
  n_post <- ncol(theta_post)
  idx <- sample(
    seq_len(n_post),
    size = n_rep,
    replace = n_rep > n_post
  )
  
  y_rep <- model(
    theta = theta_post[, idx, drop = FALSE],
    model_param = model_param
  )
  
  n_times <- model_param$n_times
  x1_rep <- y_rep[seq_len(n_times), , drop = FALSE]
  x2_rep <- y_rep[n_times + seq_len(n_times), , drop = FALSE]
  
  if (verbose) {
    cat("Generated", n_rep, "LV posterior predictive trajectories.\n")
  }
  
  make_envelope <- function(mat, observed) {
    data.frame(
      time = observed_times,
      observed = observed,
      pred_median = apply(mat, 1, median),
      pred_q025 = apply(mat, 1, quantile, probs = 0.025),
      pred_q25 = apply(mat, 1, quantile, probs = 0.25),
      pred_q75 = apply(mat, 1, quantile, probs = 0.75),
      pred_q975 = apply(mat, 1, quantile, probs = 0.975)
    )
  }
  
  out <- list(
    x1 = make_envelope(x1_rep, observed_mat[, "x1"]),
    x2 = make_envelope(x2_rep, observed_mat[, "x2"])
  )
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(out$x1, file.path(output_dir, "ppc_x1_envelope.csv"), row.names = FALSE)
    write.csv(out$x2, file.path(output_dir, "ppc_x2_envelope.csv"), row.names = FALSE)
    saveRDS(out, file.path(output_dir, "ppc_out.rds"))
  }
  
  out
}

plot_lv_ppc <- function(ppc_out,
                        output_file = NULL,
                        main = "Lotka--Volterra posterior predictive check") {
  if (!is.null(output_file)) {
    pdf(output_file, width = 9, height = 7)
    on.exit(dev.off())
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))
  
  plot_one <- function(envelope, species_label) {
    ylim <- range(
      c(envelope$observed, envelope$pred_q025, envelope$pred_q975),
      finite = TRUE
    )
    
    plot(
      envelope$time,
      envelope$pred_median,
      type = "n",
      ylim = ylim,
      xlab = "time",
      ylab = species_label,
      main = paste(main, "-", species_label)
    )
    
    polygon(
      c(envelope$time, rev(envelope$time)),
      c(envelope$pred_q975, rev(envelope$pred_q025)),
      border = NA,
      col = adjustcolor("gray80", alpha.f = 0.7)
    )
    
    polygon(
      c(envelope$time, rev(envelope$time)),
      c(envelope$pred_q75, rev(envelope$pred_q25)),
      border = NA,
      col = adjustcolor("gray60", alpha.f = 0.7)
    )
    
    lines(envelope$time, envelope$pred_median, lwd = 2)
    points(envelope$time, envelope$observed, pch = 19)
    lines(envelope$time, envelope$observed, lty = 2)
    
    legend(
      "topright",
      legend = c("Observed", "Predictive median", "50% interval", "95% interval"),
      pch = c(19, NA, 15, 15),
      lty = c(2, 1, NA, NA),
      col = c(
        "black",
        "black",
        adjustcolor("gray60", alpha.f = 0.7),
        adjustcolor("gray80", alpha.f = 0.7)
      ),
      pt.cex = c(1, NA, 2, 2),
      bty = "n"
    )
  }
  
  plot_one(ppc_out$x1, "prey x1")
  plot_one(ppc_out$x2, "predator x2")
}

# The clipped RatioIS posterior is used for the PPC.
theta_post <- ratio_out$theta_resampled_trunc

ppc_out <- simulate_lv_ppc(
  theta_post = theta_post,
  observed_mat = observed_mat,
  observed_times = observed_times,
  model = model,
  model_param = model_param,
  n_rep = 1000,
  seed = 1,
  output_dir = file.path(output_dir, "ppc")
)

plot_lv_ppc(
  ppc_out,
  output_file = file.path(output_dir, "ppc", "ppc_time_envelope.pdf")
)
