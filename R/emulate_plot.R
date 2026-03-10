#' Plot results from target trial emulation
#'
#' Creates \pkg{ggplot2} visualizations of TTE analysis results. Four plot
#' types are available, each designed for a different stage of the analysis.
#'
#' @section Plot types:
#' \describe{
#'   \item{\strong{\code{"km"}} -- Kaplan-Meier curves}{
#'     Plots weighted (if weights exist) or unweighted Kaplan-Meier survival
#'     curves for the treatment and control arms using the expanded trial
#'     data. This gives a non-parametric view of the survival experience in
#'     each arm. Requires the data to have been expanded
#'     (\code{\link{emulate_expand}}). Use this for a quick visual check before
#'     fitting a model.
#'   }
#'   \item{\strong{\code{"cumhaz"}} -- Marginal prediction curves}{
#'     Plots the G-formula standardized marginal predictions (cumulative
#'     incidence or survival, depending on what was computed) for both arms,
#'     with Monte Carlo confidence intervals. This is the primary results
#'     plot. Requires \code{\link{emulate_predict}} to have been run first.
#'     Despite the name "cumhaz," this plots whatever prediction type
#'     (cumulative incidence or survival) was specified in
#'     \code{\link{emulate_predict}}.
#'   }
#'   \item{\strong{\code{"weights"}} -- Weight distribution}{
#'     Plots kernel density curves of the inverse probability weights,
#'     separately for the treatment and control arms. Use this to visually
#'     inspect the weight distribution for extreme values or heavy tails.
#'     Requires \code{\link{emulate_weight}} to have been run.
#'   }
#'   \item{\strong{\code{"balance"}} -- Love plot (covariate balance)}{
#'     Plots a Love plot showing the absolute standardized mean differences
#'     (SMDs) for each covariate, comparing unweighted (open circles) vs.
#'     weighted (filled circles) balance. A dashed vertical line at
#'     \code{SMD = 0.1} indicates the common balance threshold. Requires
#'     \code{\link{emulate_diagnose}} to have been run with
#'     \code{balance_covariates}.
#'   }
#' }
#'
#' @param obj A \code{emulate} object. The required pipeline state depends on the
#'   plot type: \code{"km"} and \code{"weights"} require expansion,
#'   \code{"cumhaz"} requires predictions, and \code{"balance"} requires
#'   diagnostics with balance covariates.
#' @param type A character string specifying the plot type. One of
#'   \code{"km"}, \code{"cumhaz"} (default), \code{"weights"}, or
#'   \code{"balance"}.
#' @param ci A logical value. If \code{TRUE} (the default), confidence
#'   intervals are shown as shaded ribbons on \code{"km"} and \code{"cumhaz"}
#'   plots. Ignored for \code{"weights"} and \code{"balance"} plots.
#' @param title An optional character string for the plot title. If
#'   \code{NULL} (the default), a sensible default title is used.
#'
#' @return A \code{ggplot2} object that can be further customized, saved
#'   with \code{ggsave()}, or printed to display.
#'
#' @seealso \code{\link{emulate_predict}} for computing predictions (needed for
#'   \code{type = "cumhaz"}), \code{\link{emulate_diagnose}} for balance
#'   assessment (needed for \code{type = "balance"}).
#'
#' @examples
#' \donttest{
#' # Full pipeline to produce a cumulative incidence plot
#' set.seed(42)
#' n_ids <- 30
#' n_per <- 8
#' dat <- data.frame(
#'   id = rep(seq_len(n_ids), each = n_per),
#'   period = rep(0:(n_per - 1), times = n_ids),
#'   treatment = 0L, outcome = 0L, eligible = 1L,
#'   age = rep(rnorm(n_ids, 50, 10), each = n_per)
#' )
#' dat$treatment <- ifelse(dat$period >= sample(2:6, nrow(dat), replace = TRUE), 1L, 0L)
#' dat$outcome[sample(which(dat$period > 4), 4)] <- 1L
#'
#' obj <- emulate_prepare(dat, id = "id", period = "period",
#'                    treatment = "treatment", outcome = "outcome",
#'                    eligible = "eligible", covariates = "age",
#'                    estimand = "ITT")
#' obj <- emulate_expand(obj, maxfollowup = 5)
#' obj <- emulate_weight(obj)
#' obj <- emulate_fit(obj, outcome_cov = "age")
#' obj <- emulate_predict(obj, times = 1:5, samples = 50, seed = 123)
#'
#' # Cumulative incidence plot (default)
#' emulate_plot(obj, type = "cumhaz")
#'
#' # Kaplan-Meier plot
#' emulate_plot(obj, type = "km")
#'
#' # Weight distribution plot
#' emulate_plot(obj, type = "weights")
#' }
#'
#' @export
emulate_plot <- function(obj, type = "cumhaz", ci = TRUE, title = NULL) {
  .check_emulate(obj)

  if (!type %in% c("km", "cumhaz", "weights", "balance")) {
    stop("type must be 'km', 'cumhaz', 'weights', or 'balance'", call. = FALSE)
  }

  switch(type,
    km = .plot_km(obj, ci, title),
    cumhaz = .plot_cumhaz(obj, ci, title),
    weights = .plot_weights(obj, title),
    balance = .plot_balance(obj, title)
  )
}

#' @keywords internal
.plot_km <- function(obj, ci, title) {
  .check_expanded(obj)
  dt <- obj$data
  s <- obj$settings
  prefix <- s$prefix

  fu_col <- paste0(prefix, "followup")
  arm_col <- paste0(prefix, "arm")
  outobs_col <- paste0(prefix, "outcome_obs")
  trial_col <- paste0(prefix, "trial")
  weight_col <- paste0(prefix, "weight")

  # One row per person-trial-arm: max followup and max outcome
  km_data <- dt[, .(time = max(get(fu_col)) + 1L,
                     event = max(get(outobs_col))),
                by = c(s$id, trial_col, arm_col)]

  has_weights <- weight_col %in% names(dt)

  if (has_weights) {
    # Merge exit-time weight (last follow-up per person-trial-arm)
    w_exit <- dt[, .SD[which.max(get(fu_col))],
                 by = c(s$id, trial_col, arm_col),
                 .SDcols = weight_col]
    km_data <- merge(km_data, w_exit, by = c(s$id, trial_col, arm_col),
                     all.x = TRUE)

    sf <- survival::survfit(
      survival::Surv(time, event) ~ get(arm_col),
      data = km_data,
      weights = km_data[[weight_col]]
    )
  } else {
    sf <- survival::survfit(
      survival::Surv(time, event) ~ get(arm_col),
      data = km_data
    )
  }

  # Convert survfit to data.frame
  sf_df <- .survfit_to_df(sf)

  p <- ggplot2::ggplot(sf_df, ggplot2::aes(x = time, y = surv,
                                             color = group)) +
    ggplot2::geom_step(linewidth = 0.8)

  if (ci && "lower" %in% names(sf_df)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = group),
      alpha = 0.2, color = NA
    )
  }

  p <- p +
    ggplot2::scale_color_manual(values = c("navy", "firebrick"),
                                 labels = c("Control", "Treated")) +
    ggplot2::labs(
      x = "Follow-up time", y = "Survival probability",
      title = if (!is.null(title)) title else "Kaplan-Meier Curves",
      color = "Arm"
    ) +
    ggplot2::theme_minimal()

  p
}

#' @keywords internal
.plot_cumhaz <- function(obj, ci, title) {
  if (is.null(obj$predictions)) {
    stop("No predictions found. Run emulate_predict() first.", call. = FALSE)
  }

  pred <- obj$predictions

  # Build plot data
  plot_df <- data.frame(
    time = c(pred$time, pred$time),
    estimate = c(pred$est_0, pred$est_1),
    arm = rep(c("Control", "Treated"), each = nrow(pred))
  )

  if (ci) {
    plot_df$lower <- c(pred$ci_lo_0, pred$ci_lo_1)
    plot_df$upper <- c(pred$ci_hi_0, pred$ci_hi_1)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = estimate,
                                               color = arm))

  if (ci && "lower" %in% names(plot_df)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = arm),
      alpha = 0.15, color = NA
    )
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("Control" = "navy",
                                            "Treated" = "firebrick")) +
    ggplot2::labs(
      x = "Follow-up time",
      y = if (is.null(obj$prediction_meta)) "Estimate"
          else if (obj$prediction_meta$type == "cum_inc") "Cumulative Incidence"
          else "Survival Probability",
      title = if (!is.null(title)) title else "Marginal Predictions",
      color = "Arm"
    ) +
    ggplot2::theme_minimal()

  p
}

#' @keywords internal
.plot_weights <- function(obj, title) {
  .check_expanded(obj)
  dt <- obj$data
  prefix <- obj$settings$prefix

  weight_col <- paste0(prefix, "weight")
  arm_col <- paste0(prefix, "arm")

  if (!weight_col %in% names(dt)) {
    stop("No weight variable found", call. = FALSE)
  }

  plot_df <- data.frame(
    weight = dt[[weight_col]],
    arm = factor(dt[[arm_col]], labels = c("Control", "Treated"))
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = weight, color = arm)) +
    ggplot2::geom_density(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("Control" = "navy",
                                            "Treated" = "firebrick")) +
    ggplot2::labs(
      x = "Weight", y = "Density",
      title = if (!is.null(title)) title else "IP Weight Distribution",
      color = "Arm"
    ) +
    ggplot2::theme_minimal()

  p
}

#' @keywords internal
.plot_balance <- function(obj, title) {
  if (is.null(obj$diagnostics$balance)) {
    stop("No balance data. Run emulate_diagnose() with balance_covariates first.",
         call. = FALSE)
  }

  bal <- obj$diagnostics$balance
  n_covs <- nrow(bal)

  plot_df <- data.frame(
    covariate = rep(bal$covariate, 2),
    smd = c(abs(bal$smd_unwt), abs(bal$smd_wt)),
    type = rep(c("Unweighted", "Weighted"), each = n_covs),
    y = rep(seq_len(n_covs), 2)
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = smd, y = y,
                                               color = type, shape = type)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey50") +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n_covs),
      labels = bal$covariate
    ) +
    ggplot2::scale_color_manual(values = c("Unweighted" = "navy",
                                            "Weighted" = "firebrick")) +
    ggplot2::scale_shape_manual(values = c("Unweighted" = 1,
                                            "Weighted" = 16)) +
    ggplot2::labs(
      x = "Absolute Standardized Mean Difference",
      y = "",
      title = if (!is.null(title)) title else "Love Plot: Covariate Balance",
      color = "Type", shape = "Type"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Convert survfit object to data.frame
#' @keywords internal
.survfit_to_df <- function(sf) {
  if (is.null(sf$strata)) {
    data.frame(time = sf$time, surv = sf$surv,
               lower = sf$lower, upper = sf$upper,
               group = "All")
  } else {
    strata_names <- names(sf$strata)
    strata_lengths <- sf$strata
    group <- rep(strata_names, strata_lengths)

    data.frame(time = sf$time, surv = sf$surv,
               lower = sf$lower, upper = sf$upper,
               group = group)
  }
}
