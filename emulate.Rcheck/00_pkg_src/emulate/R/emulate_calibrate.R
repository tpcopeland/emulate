#' Negative control calibration for target trial emulation
#'
#' Adjusts the primary treatment effect estimate using negative control outcomes
#' (NCOs) to correct for residual systematic error. This implements the
#' "debiasing" approach from the empirical calibration framework.
#'
#' @section How calibration works:
#' \enumerate{
#'   \item Run the TTE pipeline on each negative control outcome (outcomes
#'     with a known true effect of zero) to get estimated treatment effects.
#'   \item Pass these NCO results to \code{emulate_calibrate()}.
#'   \item The function fits a null distribution to the NCO estimates using
#'     \code{EmpiricalCalibration::fitNull()}.
#'   \item The primary estimate's confidence interval is widened to account
#'     for the estimated systematic error.
#' }
#'
#' @section Negative control outcomes:
#' NCOs are outcomes that are \emph{known} (or strongly believed) to not be
#' causally affected by the treatment. Any nonzero estimated effect on an NCO
#' must therefore reflect bias (confounding, selection, measurement error).
#' By characterizing this bias across multiple NCOs, we can adjust the primary
#' estimate.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_fit}}.
#'   Must have \code{state$fitted == TRUE}.
#' @param nco_results A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{outcome_name}}{Character: name of the negative control outcome.}
#'     \item{\code{log_estimate}}{Numeric: log-scale treatment effect estimate
#'       (log-OR, log-HR, or log-RR depending on model type).}
#'     \item{\code{se_log_estimate}}{Numeric: standard error of the log-scale
#'       estimate.}
#'   }
#'   At least 3 NCOs are recommended; minimum is 2.
#' @param level Confidence level as a percentage. Default \code{95}.
#' @param quiet Suppress messages. Default \code{FALSE}.
#'
#' @return The input \code{emulate} object with:
#'   \itemize{
#'     \item \code{state$calibrated} set to \code{TRUE}
#'     \item \code{calibration} slot populated with calibrated estimates
#'   }
#'
#' @seealso \code{\link{emulate_fit}} for the previous step,
#'   \code{\link{emulate_plot}} with \code{type = "calibration"} for
#'   visualization.
#'
#' @references
#' Schuemie MJ, Hripcsak G, Ryan PB, Madigan D, Suchard MA (2018).
#' "Empirical confidence interval calibration for population-level effect
#' estimation studies in observational healthcare data."
#' \emph{PNAS}, 115(11), 2571-2577.
#'
#' @export
emulate_calibrate <- function(obj, nco_results, level = 95, quiet = FALSE) {
  .check_fitted(obj)

  if (!requireNamespace("EmpiricalCalibration", quietly = TRUE)) {
    stop("Package 'EmpiricalCalibration' is required for calibration. ",
         "Install with: install.packages('EmpiricalCalibration')",
         call. = FALSE)
  }

  if (missing(nco_results) || is.null(nco_results)) {
    stop("nco_results is required for calibration", call. = FALSE)
  }

  required_cols <- c("outcome_name", "log_estimate", "se_log_estimate")
  missing_cols <- setdiff(required_cols, names(nco_results))
  if (length(missing_cols) > 0) {
    stop("nco_results must have columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (nrow(nco_results) < 2) {
    stop("At least 2 negative control outcomes are required for calibration",
         call. = FALSE)
  }

  message("emulate_calibrate - Negative Control Calibration")
  message(strrep("-", 50))
  message("NCOs:          ", nrow(nco_results))

  # Get primary treatment effect
  b_treat <- obj$model$b_treat
  se_treat <- obj$model$se_treat

  # Fit systematic error model from NCOs (all have true effect = 0)
  error_model <- EmpiricalCalibration::fitSystematicErrorModel(
    logRr = nco_results$log_estimate,
    seLogRr = nco_results$se_log_estimate,
    trueLogRr = rep(0, nrow(nco_results))
  )

  # Also fit null distribution for reporting
  null_dist <- EmpiricalCalibration::fitNull(
    logRr = nco_results$log_estimate,
    seLogRr = nco_results$se_log_estimate
  )

  message("Null distribution:")
  message("  Mean:        ", sprintf("%.4f", null_dist["mean"]))
  message("  SD:          ", sprintf("%.4f", null_dist["sd"]))

  # Calibrate the primary estimate using the systematic error model
  calibrated <- EmpiricalCalibration::calibrateConfidenceInterval(
    logRr = b_treat,
    seLogRr = se_treat,
    model = error_model,
    ciWidth = level / 100
  )

  cal_log_est <- calibrated$logRr
  cal_lo <- calibrated$logLb95Rr
  cal_hi <- calibrated$logUb95Rr
  cal_se <- calibrated$seLogRr

  # Display results
  z_crit <- qnorm((100 + level) / 200)

  message("\nPrimary estimate (uncalibrated):")
  message("  Log-scale:   ", sprintf("%.4f", b_treat),
          " (SE: ", sprintf("%.4f", se_treat), ")")
  message("  Exp:         ", sprintf("%.4f", exp(b_treat)),
          " (", level, "% CI: ",
          sprintf("%.4f", exp(b_treat - z_crit * se_treat)),
          " - ",
          sprintf("%.4f", exp(b_treat + z_crit * se_treat)), ")")

  message("\nCalibrated estimate:")
  message("  Log-scale:   ", sprintf("%.4f", cal_log_est),
          " (SE: ", sprintf("%.4f", cal_se), ")")
  message("  Exp:         ", sprintf("%.4f", exp(cal_log_est)),
          " (", level, "% CI: ",
          sprintf("%.4f", exp(cal_lo)),
          " - ",
          sprintf("%.4f", exp(cal_hi)), ")")

  message(strrep("-", 50))

  obj$state$calibrated <- TRUE
  obj$calibration <- list(
    nco_results = nco_results,
    null_distribution = null_dist,
    error_model = error_model,
    uncalibrated = list(
      log_estimate = b_treat,
      se_log_estimate = se_treat,
      ci_lo = b_treat - z_crit * se_treat,
      ci_hi = b_treat + z_crit * se_treat
    ),
    calibrated = list(
      log_estimate = cal_log_est,
      se_log_estimate = cal_se,
      ci_lo = cal_lo,
      ci_hi = cal_hi
    ),
    level = level
  )

  invisible(obj)
}
