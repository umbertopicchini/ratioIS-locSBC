# ============================================================
# plot_ratio_results_v15_q2_interface.R
# ------------------------------------------------------------
# Generic plotting helpers for ratio-correction experiments.
# Works for arbitrary parameter dimension d >= 1.
#
# Summary of changes relative to earlier versions:
#   - Supports q1_samp = NULL and plot_q1 = TRUE/FALSE.
#   - If plot_q1 = FALSE, q1/proposal samples are ignored and
#     no q1/proposal curves, points, contours, panels, or legend
#     entries are plotted.
#   - plot_2d_density = TRUE produces filled 2D KDE panels.
#   - method_labels can be supplied to relabel reference, q1,
#     corrected, and true entries in legends and panel titles.
#     Allowed names are: reference, q1, corrected, true.
#   - Removed remaining hard-coded plotting labels from the KDE
#     branch so custom corrected labels are not silently ignored.
#   - param_limits can be supplied by theta label, by theta1/theta2
#     style names, by param1/param2 style names, or positionally.
#   - v9: added explicit label arguments reference_label, q1_label,
#     corrected_label, and true_label.
#   - v10: label overrides are applied with a shared helper and every
#     low-level plotting function accepts the explicit label arguments,
#     so custom labels work whether functions are called directly or
#     through plot_ratio_results().
#   - v11: the old hard-coded corrected-method title has been removed
#     from the file entirely; the fallback label is now "corrected".
#     print_labels = TRUE can be used to confirm the active labels.
#   - v12: unified the visual encoding across scatter, KDE, and filled
#     density plots:
#       * reference = black/grey
#       * q1/proposal = orange
#       * corrected method = blue
#       * true parameter = red cross
#     Filled 2D density panels now use method-specific palettes instead
#     of the same orange/yellow heat palette for all methods.
#   - v13: added font-size controls:
#       title_cex for plot titles,
#       axis_label_cex for x/y axis labels,
#       axis_tick_cex for axis tick labels,
#       legend_cex for legends.
#   - v14: added show_legend = TRUE/FALSE and propagated it to
#     all low-level plotting functions. print_labels only controls
#     the console diagnostic message; show_legend controls whether
#     legends are drawn.
#   - v15: adds q2-facing aliases in the main wrapper:
#       q2_samp, plot_q2, q2_label.
#     Internally the older q1_samp/plot_q1 machinery is retained for
#     backward compatibility, but user-facing labels/messages use q2
#     when the q2 aliases are supplied.
# ============================================================

# -----------------------------
# Default colors
# -----------------------------
default_ratio_plot_colors <- function() {
  list(
    reference_line = "black",
    proposal_line  = "darkorange2",
    corrected_line = "dodgerblue3",

    reference_fill = rgb(0, 0, 0, 0.08),
    proposal_fill  = rgb(1, 0.55, 0, 0.08),
    corrected_fill = rgb(0.12, 0.47, 0.95, 0.08),

    reference_density_low  = "white",
    reference_density_high = "grey20",
    proposal_density_low   = "white",
    proposal_density_high  = "darkorange2",
    corrected_density_low  = "white",
    corrected_density_high = "dodgerblue3",

    true_point = "red3"
  )
}

# -----------------------------
# Default font sizes
# -----------------------------
default_ratio_plot_fonts <- function() {
  list(
    title_cex = 1.15,
    axis_label_cex = 1.35,
    axis_tick_cex = 1.15,
    legend_cex = 1.15
  )
}

normalize_ratio_plot_fonts <- function(title_cex = NULL,
                                       axis_label_cex = NULL,
                                       axis_tick_cex = NULL,
                                       legend_cex = NULL) {
  fonts <- default_ratio_plot_fonts()

  if (!is.null(title_cex))       fonts$title_cex <- title_cex
  if (!is.null(axis_label_cex))  fonts$axis_label_cex <- axis_label_cex
  if (!is.null(axis_tick_cex))   fonts$axis_tick_cex <- axis_tick_cex
  if (!is.null(legend_cex))      fonts$legend_cex <- legend_cex

  fonts <- lapply(fonts, as.numeric)

  if (any(!is.finite(unlist(fonts))) || any(unlist(fonts) <= 0)) {
    stop("Font-size controls must be positive finite numbers.")
  }

  fonts
}

# -----------------------------
# Default method labels
# -----------------------------
default_ratio_plot_labels <- function() {
  c(
    reference = "Reference",
    q1 = "q1",
    corrected = "corrected",
    true = "True"
  )
}

plot_ratio_results_version <- function() {
  "plot_ratio_results_v15_q2_interface"
}

normalize_ratio_plot_labels <- function(method_labels = NULL) {
  defaults <- default_ratio_plot_labels()

  if (is.null(method_labels)) {
    return(defaults)
  }

  if (is.list(method_labels)) {
    method_labels <- unlist(method_labels, use.names = TRUE)
  }

  method_labels <- as.character(method_labels)

  if (is.null(names(method_labels))) {
    if (length(method_labels) > length(defaults)) {
      stop("method_labels without names can have at most ", length(defaults), " entries.")
    }
    names(method_labels) <- names(defaults)[seq_along(method_labels)]
  }

  bad <- setdiff(names(method_labels), names(defaults))
  if (length(bad) > 0) {
    stop(
      "Unknown method_labels names: ", paste(bad, collapse = ", "),
      ". Allowed names are: ", paste(names(defaults), collapse = ", "), "."
    )
  }

  defaults[names(method_labels)] <- method_labels
  defaults
}


apply_ratio_plot_label_overrides <- function(method_labels = NULL,
                                             reference_label = NULL,
                                             q1_label = NULL,
                                             corrected_label = NULL,
                                             true_label = NULL) {
  labels <- normalize_ratio_plot_labels(method_labels)

  if (!is.null(reference_label)) labels["reference"] <- as.character(reference_label)
  if (!is.null(q1_label))        labels["q1"]        <- as.character(q1_label)
  if (!is.null(corrected_label)) labels["corrected"] <- as.character(corrected_label)
  if (!is.null(true_label))      labels["true"]      <- as.character(true_label)

  labels
}

# -----------------------------
# Input checking / harmonizing
# -----------------------------
coerce_sample_matrix <- function(x, name = "sample") {
  if (is.null(x)) {
    return(NULL)
  }

  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (nrow(x) < 2 && ncol(x) >= 2) {
    return(x)
  }

  if (ncol(x) < 2 && nrow(x) >= 2) {
    return(t(x))
  }

  x
}

validate_plot_inputs <- function(reference_samples,
                                 q1_samp = NULL,
                                 theta_resampled,
                                 theta_true = NULL,
                                 theta_labels = NULL,
                                 plot_q1 = TRUE) {
  reference_samples <- coerce_sample_matrix(reference_samples, "reference_samples")
  theta_resampled   <- coerce_sample_matrix(theta_resampled, "theta_resampled")

  if (!isTRUE(plot_q1)) {
    q1_samp <- NULL
  } else {
    q1_samp <- coerce_sample_matrix(q1_samp, "q1_samp")
  }

  d_ref <- ncol(reference_samples)
  d_cor <- ncol(theta_resampled)

  if (d_ref != d_cor) {
    stop("reference_samples and theta_resampled must have the same number of columns.")
  }

  if (!is.null(q1_samp)) {
    d_q1 <- ncol(q1_samp)
    if (d_q1 != d_ref) {
      stop("q1_samp must have the same number of columns as reference_samples and theta_resampled.")
    }
  }

  d <- d_ref

  if (!is.null(theta_true)) {
    theta_true <- as.numeric(theta_true)
    if (length(theta_true) != d) {
      stop("theta_true must have length equal to the number of inferred parameters.")
    }
  }

  if (is.null(theta_labels)) {
    theta_labels <- paste0("theta", seq_len(d))
  } else {
    if (length(theta_labels) != d) {
      stop("theta_labels must have length equal to the number of inferred parameters.")
    }
  }

  list(
    reference_samples = reference_samples,
    q1_samp = q1_samp,
    theta_resampled = theta_resampled,
    theta_true = theta_true,
    theta_labels = theta_labels,
    d = d,
    plot_q1 = !is.null(q1_samp)
  )
}

# -----------------------------
# Helper: all 2D pairs
# -----------------------------
all_param_pairs <- function(d) {
  if (d < 2) {
    return(matrix(integer(0), nrow = 2))
  }
  combn(d, 2)
}

# -----------------------------
# Helper: pretty plot layout
# -----------------------------
square_layout <- function(n_panels) {
  if (n_panels <= 0) {
    return(c(1, 1))
  }
  ncol <- ceiling(sqrt(n_panels))
  nrow <- ceiling(n_panels / ncol)
  c(nrow, ncol)
}

# -----------------------------
# Helper: optional parameter axis limits
# -----------------------------
get_param_limit <- function(param_limits, theta_labels, j) {
  if (is.null(param_limits)) {
    return(NULL)
  }

  lim <- NULL
  nm <- names(param_limits)

  # Accept several naming conventions:
  #   1. current plotting label, e.g. "beta", "gamma"
  #   2. generic theta names, e.g. "theta1", "theta2"
  #   3. generic parameter names, e.g. "param1", "param2"
  #   4. numeric names, e.g. "1", "2"
  #   5. unnamed list entries by position
  key_candidates <- c(
    theta_labels[j],
    paste0("theta", j),
    paste0("theta_", j),
    paste0("param", j),
    paste0("param_", j),
    as.character(j)
  )

  if (!is.null(nm) && length(nm) > 0) {
    hit <- key_candidates[key_candidates %in% nm]
    if (length(hit) > 0) {
      lim <- param_limits[[hit[1]]]
    }
  }

  if (is.null(lim) && length(param_limits) >= j) {
    lim <- param_limits[[j]]
  }

  if (is.null(lim)) {
    return(NULL)
  }

  lim <- as.numeric(lim)
  if (length(lim) != 2 || any(!is.finite(lim)) || lim[1] >= lim[2]) {
    stop(
      "Each entry of param_limits must be a finite increasing vector of length 2. ",
      "Problem with parameter index ", j, " / label ", theta_labels[j], "."
    )
  }

  lim
}

get_axis_range <- function(values, param_limits, theta_labels, j) {
  lim <- get_param_limit(param_limits, theta_labels, j)
  if (!is.null(lim)) {
    return(lim)
  }
  range(values, finite = TRUE)
}

# -----------------------------
# 1D marginal overlays
# -----------------------------
plot_1d_marginals <- function(reference_samples,
                              q1_samp = NULL,
                              theta_resampled,
                              theta_true = NULL,
                              theta_labels = NULL,
                              dims = NULL,
                              colors = default_ratio_plot_colors(),
                              lwd = 2,
                              main_prefix = "Marginal",
                              legend_position = "topright",
                              param_limits = NULL,
                              plot_q1 = TRUE,
                              method_labels = NULL,
                              reference_label = NULL,
                              q1_label = NULL,
                              corrected_label = NULL,
                              true_label = NULL,
                              title_cex = NULL,
                              axis_label_cex = NULL,
                              axis_tick_cex = NULL,
                              legend_cex = NULL,
                              show_legend = TRUE) {
  inp <- validate_plot_inputs(
    reference_samples = reference_samples,
    q1_samp = q1_samp,
    theta_resampled = theta_resampled,
    theta_true = theta_true,
    theta_labels = theta_labels,
    plot_q1 = plot_q1
  )

  reference_samples <- inp$reference_samples
  q1_samp           <- inp$q1_samp
  theta_resampled   <- inp$theta_resampled
  theta_true        <- inp$theta_true
  theta_labels      <- inp$theta_labels
  d                 <- inp$d
  plot_q1           <- inp$plot_q1
  method_labels    <- apply_ratio_plot_label_overrides(
    method_labels = method_labels,
    reference_label = reference_label,
    q1_label = q1_label,
    corrected_label = corrected_label,
    true_label = true_label
  )

  plot_fonts <- normalize_ratio_plot_fonts(
    title_cex = title_cex,
    axis_label_cex = axis_label_cex,
    axis_tick_cex = axis_tick_cex,
    legend_cex = legend_cex
  )

  if (is.null(dims)) {
    dims <- seq_len(d)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  lay <- square_layout(length(dims))
  par(mfrow = lay)

  for (j in dims) {
    dens_ref <- density(reference_samples[, j])
    dens_cor <- density(theta_resampled[, j])
    dens_q1 <- NULL
    if (plot_q1) {
      dens_q1 <- density(q1_samp[, j])
    }

    x_values <- c(dens_ref$x, dens_cor$x, if (plot_q1) dens_q1$x else NULL,
                  if (!is.null(theta_true)) theta_true[j] else NULL)
    y_values <- c(0, dens_ref$y, dens_cor$y, if (plot_q1) dens_q1$y else NULL)

    xlim <- get_axis_range(
      values = x_values,
      param_limits = param_limits,
      theta_labels = theta_labels,
      j = j
    )
    ylim <- range(y_values)

    plot(dens_ref,
         xlim = xlim,
         ylim = ylim,
         lwd = lwd,
         col = colors$reference_line,
         main = paste(main_prefix, theta_labels[j]),
         xlab = theta_labels[j],
         ylab = "Density",
         cex.main = plot_fonts$title_cex,
         cex.lab = plot_fonts$axis_label_cex,
         cex.axis = plot_fonts$axis_tick_cex)

    if (plot_q1) {
      lines(dens_q1, lwd = lwd, col = colors$proposal_line)
    }
    lines(dens_cor, lwd = lwd, col = colors$corrected_line)

    if (!is.null(theta_true)) {
      abline(v = theta_true[j], col = colors$true_point, lwd = 2, lty = 2)
    }

    if (isTRUE(show_legend)) {
      legend_items <- c(unname(method_labels["reference"]), if (plot_q1) unname(method_labels["q1"]) else NULL, unname(method_labels["corrected"]) )
      legend_cols  <- c(colors$reference_line, if (plot_q1) colors$proposal_line else NULL, colors$corrected_line)
      legend_lty   <- rep(1, length(legend_items))

      if (!is.null(theta_true)) {
        legend_items <- c(legend_items, unname(method_labels["true"]))
        legend_cols  <- c(legend_cols, colors$true_point)
        legend_lty   <- c(legend_lty, 2)
      }

      legend(legend_position,
             legend = legend_items,
             col = legend_cols,
             lty = legend_lty,
             lwd = 2,
             bty = "n",
             cex = plot_fonts$legend_cex)
    }
  }
}

# -----------------------------
# 2D scatter overlays
# -----------------------------
plot_2d_scatter_overlays <- function(reference_samples,
                                     q1_samp = NULL,
                                     theta_resampled,
                                     theta_true = NULL,
                                     theta_labels = NULL,
                                     pairs = NULL,
                                     colors = default_ratio_plot_colors(),
                                     point_cex = 0.45,
                                     legend_position = "topright",
                                     param_limits = NULL,
                                     plot_q1 = TRUE,
                                     method_labels = NULL,
                                     reference_label = NULL,
                                     q1_label = NULL,
                                     corrected_label = NULL,
                                     true_label = NULL,
                                     title_cex = NULL,
                                     axis_label_cex = NULL,
                                     axis_tick_cex = NULL,
                                     legend_cex = NULL,
                                     show_legend = TRUE) {
  inp <- validate_plot_inputs(
    reference_samples = reference_samples,
    q1_samp = q1_samp,
    theta_resampled = theta_resampled,
    theta_true = theta_true,
    theta_labels = theta_labels,
    plot_q1 = plot_q1
  )

  reference_samples <- inp$reference_samples
  q1_samp           <- inp$q1_samp
  theta_resampled   <- inp$theta_resampled
  theta_true        <- inp$theta_true
  theta_labels      <- inp$theta_labels
  d                 <- inp$d
  plot_q1           <- inp$plot_q1
  method_labels    <- apply_ratio_plot_label_overrides(
    method_labels = method_labels,
    reference_label = reference_label,
    q1_label = q1_label,
    corrected_label = corrected_label,
    true_label = true_label
  )

  plot_fonts <- normalize_ratio_plot_fonts(
    title_cex = title_cex,
    axis_label_cex = axis_label_cex,
    axis_tick_cex = axis_tick_cex,
    legend_cex = legend_cex
  )

  if (d < 2) {
    message("plot_2d_scatter_overlays(): fewer than 2 parameters, nothing to plot.")
    return(invisible(NULL))
  }

  if (is.null(pairs)) {
    pairs <- all_param_pairs(d)
  }

  n_panels <- ncol(pairs)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  lay <- square_layout(n_panels)
  par(mfrow = lay)

  for (k in seq_len(n_panels)) {
    i <- pairs[1, k]
    j <- pairs[2, k]

    x_all <- c(reference_samples[, i], theta_resampled[, i],
               if (plot_q1) q1_samp[, i] else NULL,
               if (!is.null(theta_true)) theta_true[i] else NULL)
    y_all <- c(reference_samples[, j], theta_resampled[, j],
               if (plot_q1) q1_samp[, j] else NULL,
               if (!is.null(theta_true)) theta_true[j] else NULL)

    plot(reference_samples[, i], reference_samples[, j],
         pch = 16, cex = point_cex,
         col = colors$reference_fill,
         xlab = theta_labels[i],
         ylab = theta_labels[j],
         main = paste(theta_labels[i], "vs", theta_labels[j]),
         xlim = get_axis_range(x_all, param_limits, theta_labels, i),
         ylim = get_axis_range(y_all, param_limits, theta_labels, j),
         cex.main = plot_fonts$title_cex,
         cex.lab = plot_fonts$axis_label_cex,
         cex.axis = plot_fonts$axis_tick_cex)

    if (plot_q1) {
      points(q1_samp[, i], q1_samp[, j],
             pch = 16, cex = point_cex,
             col = colors$proposal_fill)
    }

    points(theta_resampled[, i], theta_resampled[, j],
           pch = 16, cex = point_cex,
           col = colors$corrected_fill)

    if (!is.null(theta_true)) {
      points(theta_true[i], theta_true[j],
             pch = 4, cex = 1.5, lwd = 2,
             col = colors$true_point)
    }

    legend_items <- c(unname(method_labels["reference"]), if (plot_q1) unname(method_labels["q1"]) else NULL, unname(method_labels["corrected"]) )
    legend_cols  <- c(colors$reference_line, if (plot_q1) colors$proposal_line else NULL, colors$corrected_line)
    legend_pch   <- rep(16, length(legend_items))

    if (!is.null(theta_true)) {
      legend_items <- c(legend_items, unname(method_labels["true"]))
      legend_cols  <- c(legend_cols, colors$true_point)
      legend_pch   <- c(legend_pch, 4)
    }
    if (isTRUE(show_legend)) {
    legend(legend_position,
           legend = legend_items,
           col = legend_cols,
           pch = legend_pch,
           bty = "n",
           cex = plot_fonts$legend_cex)
    }
  }
}

# -----------------------------
# 2D KDE contour overlays
# -----------------------------
plot_2d_kde_overlays <- function(reference_samples,
                                 q1_samp = NULL,
                                 theta_resampled,
                                 theta_true = NULL,
                                 theta_labels = NULL,
                                 pairs = NULL,
                                 colors = default_ratio_plot_colors(),
                                 kde_n = 120,
                                 lwd = 2,
                                 n_levels = 4,
                                 legend_position = "topright",
                                 param_limits = NULL,
                                 plot_q1 = TRUE,
                                 method_labels = NULL,
                                 reference_label = NULL,
                                 q1_label = NULL,
                                 corrected_label = NULL,
                                 true_label = NULL,
                                 title_cex = NULL,
                                 axis_label_cex = NULL,
                                 axis_tick_cex = NULL,
                                 legend_cex = NULL,
                                 show_legend = TRUE) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for 2D KDE plots.")
  }

  inp <- validate_plot_inputs(
    reference_samples = reference_samples,
    q1_samp = q1_samp,
    theta_resampled = theta_resampled,
    theta_true = theta_true,
    theta_labels = theta_labels,
    plot_q1 = plot_q1
  )

  reference_samples <- inp$reference_samples
  q1_samp           <- inp$q1_samp
  theta_resampled   <- inp$theta_resampled
  theta_true        <- inp$theta_true
  theta_labels      <- inp$theta_labels
  d                 <- inp$d
  plot_q1           <- inp$plot_q1
  method_labels    <- apply_ratio_plot_label_overrides(
    method_labels = method_labels,
    reference_label = reference_label,
    q1_label = q1_label,
    corrected_label = corrected_label,
    true_label = true_label
  )

  plot_fonts <- normalize_ratio_plot_fonts(
    title_cex = title_cex,
    axis_label_cex = axis_label_cex,
    axis_tick_cex = axis_tick_cex,
    legend_cex = legend_cex
  )

  if (d < 2) {
    message("plot_2d_kde_overlays(): fewer than 2 parameters, nothing to plot.")
    return(invisible(NULL))
  }

  if (is.null(pairs)) {
    pairs <- all_param_pairs(d)
  }

  n_panels <- ncol(pairs)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  lay <- square_layout(n_panels)
  par(mfrow = lay)

  for (k in seq_len(n_panels)) {
    i <- pairs[1, k]
    j <- pairs[2, k]

    x_all <- c(reference_samples[, i], theta_resampled[, i],
               if (plot_q1) q1_samp[, i] else NULL,
               if (!is.null(theta_true)) theta_true[i] else NULL)
    y_all <- c(reference_samples[, j], theta_resampled[, j],
               if (plot_q1) q1_samp[, j] else NULL,
               if (!is.null(theta_true)) theta_true[j] else NULL)

    xlim <- get_axis_range(x_all, param_limits, theta_labels, i)
    ylim <- get_axis_range(y_all, param_limits, theta_labels, j)
    kde_lims <- c(xlim, ylim)

    kde_list <- list(
      reference = MASS::kde2d(reference_samples[, i], reference_samples[, j], n = kde_n,
                              lims = kde_lims),
      corrected = MASS::kde2d(theta_resampled[, i], theta_resampled[, j], n = kde_n,
                              lims = kde_lims)
    )
    line_cols <- c(colors$reference_line, colors$corrected_line)

    if (plot_q1) {
      kde_list <- append(
        kde_list,
        list(q1 = MASS::kde2d(q1_samp[, i], q1_samp[, j], n = kde_n, lims = kde_lims)),
        after = 1
      )
      line_cols <- c(colors$reference_line,
                     colors$proposal_line,
                     colors$corrected_line)
    }

    z_all <- unlist(lapply(kde_list, function(obj) obj$z))
    levs <- pretty(range(z_all), n = n_levels)
    levs <- levs[levs > 0]

    plot(NA,
         xlim = xlim,
         ylim = ylim,
         xlab = theta_labels[i],
         ylab = theta_labels[j],
         main = paste("KDE:", theta_labels[i], "vs", theta_labels[j]),
         cex.main = plot_fonts$title_cex,
         cex.lab = plot_fonts$axis_label_cex,
         cex.axis = plot_fonts$axis_tick_cex)

    for (m in seq_along(kde_list)) {
      contour(kde_list[[m]],
              add = TRUE,
              drawlabels = FALSE,
              col = line_cols[m],
              lwd = lwd,
              levels = levs)
    }

    if (!is.null(theta_true)) {
      points(theta_true[i], theta_true[j],
             pch = 4, cex = 1.5, lwd = 2,
             col = colors$true_point)
    }

    legend_items <- c(unname(method_labels["reference"]), if (plot_q1) unname(method_labels["q1"]) else NULL, unname(method_labels["corrected"]) )
    legend_cols  <- c(colors$reference_line, if (plot_q1) colors$proposal_line else NULL, colors$corrected_line)
    legend_lty   <- rep(1, length(legend_items))
    legend_pch   <- rep(NA, length(legend_items))
    legend_lwd   <- rep(2, length(legend_items))

    if (!is.null(theta_true)) {
      legend_items <- c(legend_items, unname(method_labels["true"]))
      legend_cols  <- c(legend_cols, colors$true_point)
      legend_lty   <- c(legend_lty, NA)
      legend_pch   <- c(legend_pch, 4)
      legend_lwd   <- c(legend_lwd, NA)
    }
    if (isTRUE(show_legend)) {
    legend(legend_position,
           legend = legend_items,
           col = legend_cols,
           lty = legend_lty,
           pch = legend_pch,
           lwd = legend_lwd,
           bty = "n",
           cex = plot_fonts$legend_cex)
    }
  }
}



# -----------------------------
# Helper: method-specific density palettes
# -----------------------------
make_ratio_density_palette <- function(method = c("reference", "q1", "corrected"),
                                       colors = default_ratio_plot_colors(),
                                       n = 64) {
  method <- match.arg(method)

  if (method == "reference") {
    return(grDevices::colorRampPalette(
      c(colors$reference_density_low, colors$reference_density_high)
    )(n))
  }

  if (method == "q1") {
    return(grDevices::colorRampPalette(
      c(colors$proposal_density_low, colors$proposal_density_high)
    )(n))
  }

  grDevices::colorRampPalette(
    c(colors$corrected_density_low, colors$corrected_density_high)
  )(n)
}

# -----------------------------
# 2D filled KDE density panels
# -----------------------------
plot_2d_density_panels <- function(reference_samples,
                                   q1_samp = NULL,
                                   theta_resampled,
                                   theta_true = NULL,
                                   theta_labels = NULL,
                                   pairs = NULL,
                                   colors = default_ratio_plot_colors(),
                                   kde_n = 120,
                                   n_levels = 4,
                                   param_limits = NULL,
                                   plot_q1 = TRUE,
                                   add_reference_contour = TRUE,
                                   method_labels = NULL,
                                   reference_label = NULL,
                                   q1_label = NULL,
                                   corrected_label = NULL,
                                   true_label = NULL,
                                   title_cex = NULL,
                                   axis_label_cex = NULL,
                                   axis_tick_cex = NULL,
                                   legend_cex = NULL,
                                   show_legend = TRUE) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for 2D density plots.")
  }

  inp <- validate_plot_inputs(
    reference_samples = reference_samples,
    q1_samp = q1_samp,
    theta_resampled = theta_resampled,
    theta_true = theta_true,
    theta_labels = theta_labels,
    plot_q1 = plot_q1
  )

  reference_samples <- inp$reference_samples
  q1_samp           <- inp$q1_samp
  theta_resampled   <- inp$theta_resampled
  theta_true        <- inp$theta_true
  theta_labels      <- inp$theta_labels
  d                 <- inp$d
  plot_q1           <- inp$plot_q1
  method_labels    <- apply_ratio_plot_label_overrides(
    method_labels = method_labels,
    reference_label = reference_label,
    q1_label = q1_label,
    corrected_label = corrected_label,
    true_label = true_label
  )

  plot_fonts <- normalize_ratio_plot_fonts(
    title_cex = title_cex,
    axis_label_cex = axis_label_cex,
    axis_tick_cex = axis_tick_cex,
    legend_cex = legend_cex
  )

  if (d < 2) {
    message("plot_2d_density_panels(): fewer than 2 parameters, nothing to plot.")
    return(invisible(NULL))
  }

  if (is.null(pairs)) {
    pairs <- all_param_pairs(d)
  }

  # Helper used for each method panel.
  draw_density_panel <- function(kde_obj,
                                 ref_kde_obj,
                                 panel_title,
                                 xlab,
                                 ylab,
                                 xlim,
                                 ylim,
                                 fill_col,
                                 method_contour_col = NULL,
                                 theta_true_pair = NULL) {
    zlim <- range(kde_obj$z, finite = TRUE)

    image(kde_obj$x, kde_obj$y, kde_obj$z,
          xlim = xlim,
          ylim = ylim,
          zlim = zlim,
          col = fill_col,
          xlab = xlab,
          ylab = ylab,
          main = panel_title,
          cex.main = plot_fonts$title_cex,
          cex.lab = plot_fonts$axis_label_cex,
          cex.axis = plot_fonts$axis_tick_cex)

    if (isTRUE(add_reference_contour) && !is.null(ref_kde_obj)) {
      ref_levels <- pretty(range(ref_kde_obj$z, finite = TRUE), n = n_levels)
      ref_levels <- ref_levels[ref_levels > 0]
      if (length(ref_levels) > 0) {
        contour(ref_kde_obj,
                add = TRUE,
                drawlabels = FALSE,
                levels = ref_levels,
                col = colors$reference_line,
                lwd = 1)
      }
    }

    if (!is.null(method_contour_col)) {
      method_levels <- pretty(range(kde_obj$z, finite = TRUE), n = n_levels)
      method_levels <- method_levels[method_levels > 0]
      if (length(method_levels) > 0) {
        contour(kde_obj,
                add = TRUE,
                drawlabels = FALSE,
                levels = method_levels,
                col = method_contour_col,
                lwd = 1)
      }
    }

    if (!is.null(theta_true_pair)) {
      points(theta_true_pair[1], theta_true_pair[2],
             pch = 4, cex = 1.5, lwd = 2,
             col = colors$true_point)
    }
  }

  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]
    j <- pairs[2, k]

    x_all <- c(reference_samples[, i], theta_resampled[, i],
               if (plot_q1) q1_samp[, i] else NULL,
               if (!is.null(theta_true)) theta_true[i] else NULL)
    y_all <- c(reference_samples[, j], theta_resampled[, j],
               if (plot_q1) q1_samp[, j] else NULL,
               if (!is.null(theta_true)) theta_true[j] else NULL)

    xlim <- get_axis_range(x_all, param_limits, theta_labels, i)
    ylim <- get_axis_range(y_all, param_limits, theta_labels, j)
    kde_lims <- c(xlim, ylim)

    kde_ref <- MASS::kde2d(reference_samples[, i], reference_samples[, j],
                           n = kde_n, lims = kde_lims)
    kde_cor <- MASS::kde2d(theta_resampled[, i], theta_resampled[, j],
                           n = kde_n, lims = kde_lims)

    panel_names <- c(unname(method_labels["reference"]), if (plot_q1) unname(method_labels["q1"]) else NULL, unname(method_labels["corrected"]))
    n_panels <- length(panel_names)
    lay <- square_layout(n_panels)

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = lay)

    true_pair <- if (!is.null(theta_true)) c(theta_true[i], theta_true[j]) else NULL

    draw_density_panel(
      kde_obj = kde_ref,
      ref_kde_obj = kde_ref,
      panel_title = paste0(unname(method_labels["reference"]), ": ", theta_labels[i], " vs ", theta_labels[j]),
      xlab = theta_labels[i],
      ylab = theta_labels[j],
      xlim = xlim,
      ylim = ylim,
      fill_col = make_ratio_density_palette("reference", colors = colors),
      method_contour_col = colors$reference_line,
      theta_true_pair = true_pair
    )

    if (plot_q1) {
      kde_q1 <- MASS::kde2d(q1_samp[, i], q1_samp[, j],
                            n = kde_n, lims = kde_lims)
      draw_density_panel(
        kde_obj = kde_q1,
        ref_kde_obj = kde_ref,
        panel_title = paste0(unname(method_labels["q1"]), ": ", theta_labels[i], " vs ", theta_labels[j]),
        xlab = theta_labels[i],
        ylab = theta_labels[j],
        xlim = xlim,
        ylim = ylim,
        fill_col = make_ratio_density_palette("q1", colors = colors),
        method_contour_col = colors$proposal_line,
        theta_true_pair = true_pair
      )
    }

    draw_density_panel(
      kde_obj = kde_cor,
      ref_kde_obj = kde_ref,
      panel_title = paste0(unname(method_labels["corrected"]), ": ", theta_labels[i], " vs ", theta_labels[j]),
      xlab = theta_labels[i],
      ylab = theta_labels[j],
      xlim = xlim,
      ylim = ylim,
      fill_col = make_ratio_density_palette("corrected", colors = colors),
      method_contour_col = colors$corrected_line,
      theta_true_pair = true_pair
    )

    par(old_par)
  }

  invisible(NULL)
}

# -----------------------------
# Main wrapper
# -----------------------------
plot_ratio_results <- function(reference_samples,
                               q1_samp = NULL,
                               q2_samp = NULL,
                               theta_resampled,
                               theta_true = NULL,
                               theta_labels = NULL,
                               plot_1d = TRUE,
                               plot_2d_scatter = TRUE,
                               plot_2d_kde = TRUE,
                               plot_2d_density = FALSE,
                               plot_q1 = TRUE,
                               plot_q2 = NULL,
                               kde_n = 120,
                               n_levels = 4,
                               colors = default_ratio_plot_colors(),
                               param_limits = NULL,
                               method_labels = NULL,
                               reference_label = NULL,
                               q1_label = NULL,
                               q2_label = NULL,
                               corrected_label = NULL,
                               true_label = NULL,
                               print_labels = FALSE,
                               show_legend = TRUE,
                               title_cex = NULL,
                               axis_label_cex = NULL,
                               axis_tick_cex = NULL,
                               legend_cex = NULL) {
  using_q2_alias <- !is.null(q2_samp) || !is.null(plot_q2) || !is.null(q2_label)

  if (!is.null(q2_samp)) {
    if (!is.null(q1_samp)) {
      stop("Supply only one of q1_samp or q2_samp, not both.")
    }
    q1_samp <- q2_samp
  }

  if (!is.null(plot_q2)) {
    if (!identical(plot_q1, TRUE)) {
      stop("Supply only one of plot_q1 or plot_q2 when changing the default.")
    }
    plot_q1 <- plot_q2
  }

  if (!is.null(q2_label)) {
    if (!is.null(q1_label)) {
      stop("Supply only one of q1_label or q2_label, not both.")
    }
    q1_label <- q2_label
  } else if (isTRUE(using_q2_alias) && is.null(q1_label)) {
    q1_label <- "q2"
  }

  inp <- validate_plot_inputs(
    reference_samples = reference_samples,
    q1_samp = q1_samp,
    theta_resampled = theta_resampled,
    theta_true = theta_true,
    theta_labels = theta_labels,
    plot_q1 = plot_q1
  )

  reference_samples <- inp$reference_samples
  q1_samp           <- inp$q1_samp
  theta_resampled   <- inp$theta_resampled
  theta_true        <- inp$theta_true
  theta_labels      <- inp$theta_labels
  plot_q1           <- inp$plot_q1
  method_labels    <- apply_ratio_plot_label_overrides(
    method_labels = method_labels,
    reference_label = reference_label,
    q1_label = q1_label,
    corrected_label = corrected_label,
    true_label = true_label
  )

  plot_fonts <- normalize_ratio_plot_fonts(
    title_cex = title_cex,
    axis_label_cex = axis_label_cex,
    axis_tick_cex = axis_tick_cex,
    legend_cex = legend_cex
  )

  if (isTRUE(print_labels)) {
    proposal_name <- if (isTRUE(using_q2_alias)) "q2" else "q1"
    message(
      "plot_ratio_results active labels: reference='", method_labels[["reference"]],
      "', ", proposal_name, "='", method_labels[["q1"]],
      "', corrected='", method_labels[["corrected"]],
      "', true='", method_labels[["true"]], "'"
    )
  }


  if (plot_1d) {
    plot_1d_marginals(
      reference_samples = reference_samples,
      q1_samp = q1_samp,
      theta_resampled = theta_resampled,
      theta_true = theta_true,
      theta_labels = theta_labels,
      colors = colors,
      param_limits = param_limits,
      plot_q1 = plot_q1,
      method_labels = method_labels,
      reference_label = reference_label,
      q1_label = q1_label,
      corrected_label = corrected_label,
      true_label = true_label,
      title_cex = title_cex,
      axis_label_cex = axis_label_cex,
      axis_tick_cex = axis_tick_cex,
      legend_cex = legend_cex,
      show_legend = show_legend
    )
  }

  if (plot_2d_scatter) {
    plot_2d_scatter_overlays(
      reference_samples = reference_samples,
      q1_samp = q1_samp,
      theta_resampled = theta_resampled,
      theta_true = theta_true,
      theta_labels = theta_labels,
      colors = colors,
      param_limits = param_limits,
      plot_q1 = plot_q1,
      method_labels = method_labels,
      reference_label = reference_label,
      q1_label = q1_label,
      corrected_label = corrected_label,
      true_label = true_label,
      title_cex = title_cex,
      axis_label_cex = axis_label_cex,
      axis_tick_cex = axis_tick_cex,
      legend_cex = legend_cex,
      show_legend = show_legend
    )
  }

  if (plot_2d_kde) {
    plot_2d_kde_overlays(
      reference_samples = reference_samples,
      q1_samp = q1_samp,
      theta_resampled = theta_resampled,
      theta_true = theta_true,
      theta_labels = theta_labels,
      colors = colors,
      kde_n = kde_n,
      n_levels = n_levels,
      param_limits = param_limits,
      plot_q1 = plot_q1,
      method_labels = method_labels,
      reference_label = reference_label,
      q1_label = q1_label,
      corrected_label = corrected_label,
      true_label = true_label,
      title_cex = title_cex,
      axis_label_cex = axis_label_cex,
      axis_tick_cex = axis_tick_cex,
      legend_cex = legend_cex,
      show_legend = show_legend
    )
  }

  if (isTRUE(plot_2d_density)) {
    plot_2d_density_panels(
      reference_samples = reference_samples,
      q1_samp = q1_samp,
      theta_resampled = theta_resampled,
      theta_true = theta_true,
      theta_labels = theta_labels,
      colors = colors,
      kde_n = kde_n,
      n_levels = n_levels,
      param_limits = param_limits,
      plot_q1 = plot_q1,
      method_labels = method_labels,
      reference_label = reference_label,
      q1_label = q1_label,
      corrected_label = corrected_label,
      true_label = true_label,
      title_cex = title_cex,
      axis_label_cex = axis_label_cex,
      axis_tick_cex = axis_tick_cex,
      legend_cex = legend_cex,
      show_legend = show_legend
    )
  }

  invisible(NULL)
}
