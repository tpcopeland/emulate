#' Publication-quality results table
#'
#' Displays or exports a comprehensive analysis summary including the
#' estimand, model type, sample sizes, weight summary, coefficient table
#' with confidence intervals, and optionally the prediction table from
#' \code{\link{emulate_predict}}.
#'
#' @section Output formats:
#' \describe{
#'   \item{\code{"display"}}{(Default) Prints the results to the console
#'     via \code{message()}. This is useful for interactive analysis and for
#'     reviewing results before exporting.}
#'   \item{\code{"csv"}}{Exports the coefficient table to a CSV file.
#'     Requires the \code{export} argument. Only the coefficient table is
#'     exported (not the analysis summary header).}
#'   \item{\code{"excel"}}{Exports results to an Excel workbook with
#'     separate sheets for the analysis summary, coefficient table, and
#'     (optionally) the prediction table. Requires the \pkg{openxlsx}
#'     package and the \code{export} argument.}
#' }
#'
#' @section The eform option:
#' When \code{eform = TRUE}, coefficients are exponentiated before display.
#' For a logistic model, this converts log-odds ratios to \strong{odds
#' ratios (OR)}. For a Cox model, this converts log-hazard ratios to
#' \strong{hazard ratios (HR)}. The confidence intervals are also
#' exponentiated. The column header automatically changes to "OR" or "HR"
#' depending on the model type.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_fit}}. The
#'   \code{fitted} state flag must be \code{TRUE}.
#' @param format A character string specifying the output format. One of
#'   \code{"display"} (default), \code{"csv"}, or \code{"excel"}.
#' @param export An optional character string giving the file path for CSV
#'   or Excel export. Required when \code{format} is \code{"csv"} or
#'   \code{"excel"}; ignored when \code{format = "display"}.
#' @param decimals An integer specifying the number of decimal places for
#'   coefficient estimates, standard errors, and confidence intervals.
#'   Default is \code{3}. P-values use \code{decimals + 1} places.
#' @param eform A logical value. If \code{TRUE}, coefficients and confidence
#'   intervals are exponentiated to produce odds ratios (logistic model) or
#'   hazard ratios (Cox model). Default is \code{FALSE} (show log-scale
#'   coefficients).
#' @param ci_separator A character string used to separate the lower and
#'   upper confidence interval bounds in display output. Default is
#'   \code{" - "}.
#' @param title An optional character string used as the table title in
#'   display and Excel output. If \code{NULL}, no title is shown.
#' @param predictions A logical value. If \code{TRUE} and predictions have
#'   been computed via \code{\link{emulate_predict}}, the prediction table is
#'   included in the output. For \code{"display"} format, it is printed
#'   below the coefficient table. For \code{"excel"} format, it is added
#'   as a separate worksheet. Default is \code{FALSE}.
#'
#' @return The input \code{emulate} object, returned invisibly.
#'
#' @seealso \code{\link{emulate_fit}} for fitting the model,
#'   \code{\link{emulate_predict}} for predictions,
#'   \code{\link{emulate_plot}} for visualization.
#'
#' @examples
#' \donttest{
#' # Full pipeline
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
#'
#' # Display results with exponentiated coefficients (odds ratios)
#' emulate_report(obj, eform = TRUE, title = "TTE Analysis Results")
#'
#' # Include predictions
#' obj <- emulate_predict(obj, times = 1:5, samples = 50, seed = 123)
#' emulate_report(obj, eform = TRUE, predictions = TRUE)
#'
#' # Export to CSV
#' tmp <- tempfile(fileext = ".csv")
#' emulate_report(obj, format = "csv", export = tmp)
#' }
#'
#' @export
emulate_report <- function(obj, format = "display", export = NULL,
                       decimals = 3, eform = FALSE,
                       ci_separator = " - ", title = NULL,
                       predictions = FALSE) {
  .check_fitted(obj)

  s <- obj$settings
  dt <- obj$data
  m <- obj$model
  prefix <- s$prefix
  level <- m$level

  arm_col <- paste0(prefix, "arm")
  cens_col <- paste0(prefix, "censored")
  outobs_col <- paste0(prefix, "outcome_obs")
  trial_col <- paste0(prefix, "trial")
  weight_col <- paste0(prefix, "weight")
  esample_col <- paste0(prefix, "esample")

  # Summary statistics
  est_dt <- if (esample_col %in% names(dt)) dt[get(esample_col) == 1] else dt
  n_obs <- nrow(est_dt)
  n_treat <- sum(est_dt[[arm_col]] == 1)
  n_control <- sum(est_dt[[arm_col]] == 0)
  n_events <- sum(est_dt[[outobs_col]] == 1, na.rm = TRUE)
  n_trials <- data.table::uniqueN(est_dt[[trial_col]])

  # Coefficients
  b <- coef(m$object)
  V <- m$vcov
  se <- sqrt(diag(V))
  z_crit <- qnorm((100 + level) / 200)

  n_coefs <- length(b)
  coef_table <- data.frame(
    Variable = names(b),
    stringsAsFactors = FALSE
  )

  if (eform) {
    effect_label <- switch(m$type,
      logistic = "OR",
      poisson  = "RR",
      cox      = "HR",
      "HR"
    )
    coef_table[[effect_label]] <- round(exp(b), decimals)
    coef_table$CI_lower <- round(exp(b - z_crit * se), decimals)
    coef_table$CI_upper <- round(exp(b + z_crit * se), decimals)
  } else {
    coef_table$Estimate <- round(b, decimals)
    coef_table$SE <- round(se, decimals)
    coef_table$CI_lower <- round(b - z_crit * se, decimals)
    coef_table$CI_upper <- round(b + z_crit * se, decimals)
  }

  z_vals <- b / se
  coef_table$p_value <- round(2 * pnorm(-abs(z_vals)), decimals + 1)

  if (format == "display") {
    if (!is.null(title)) message(title)
    message(strrep("=", 60))
    message("Analysis Summary")
    message(strrep("-", 60))
    message("  Estimand:       ", s$estimand)
    message("  Model:          ", m$type)
    message("  Observations:   ", format(n_obs, big.mark = ","))
    message("  Treatment arm:  ", format(n_treat, big.mark = ","))
    message("  Control arm:    ", format(n_control, big.mark = ","))
    message("  Events:         ", format(n_events, big.mark = ","))
    message("  Trials:         ", n_trials)

    # Weight summary (if available)
    if (weight_col %in% names(dt) && s$estimand != "ITT") {
      w <- dt[[weight_col]]
      ess <- sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)
      message("\n  Weight mean:    ", sprintf("%.4f", mean(w, na.rm = TRUE)))
      message("  Weight ESS:     ", sprintf("%.1f", ess))
    }

    message(strrep("-", 60))
    message("\nCoefficients:")
    print(coef_table, row.names = FALSE)
    message(strrep("=", 60))

    if (isTRUE(obj$state$calibrated) && !is.null(obj$calibration)) {
      cal <- obj$calibration
      message("\nCalibrated Results (via Negative Control Calibration):")
      message("  NCOs used:      ", nrow(cal$nco_results))
      if (eform) {
        message("  Calibrated ", effect_label, ": ",
                sprintf(paste0("%.", decimals, "f"), exp(cal$calibrated$log_estimate)),
                " (", level, "% CI: ",
                sprintf(paste0("%.", decimals, "f"), exp(cal$calibrated$ci_lo)),
                ci_separator,
                sprintf(paste0("%.", decimals, "f"), exp(cal$calibrated$ci_hi)), ")")
      } else {
        message("  Calibrated est:  ",
                sprintf(paste0("%.", decimals, "f"), cal$calibrated$log_estimate),
                " (", level, "% CI: ",
                sprintf(paste0("%.", decimals, "f"), cal$calibrated$ci_lo),
                ci_separator,
                sprintf(paste0("%.", decimals, "f"), cal$calibrated$ci_hi), ")")
      }
    }

    if (predictions && !is.null(obj$predictions)) {
      message("\nPredictions:")
      print(obj$predictions, row.names = FALSE)
    }
  } else if (format == "csv") {
    if (is.null(export)) stop("export path required for csv format", call. = FALSE)
    utils::write.csv(coef_table, export, row.names = FALSE)
    message("Coefficients exported to: ", export)
  } else if (format == "excel") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("openxlsx package required for Excel export", call. = FALSE)
    }
    if (is.null(export)) stop("export path required for excel format", call. = FALSE)

    wb <- openxlsx::createWorkbook()

    # Summary sheet
    openxlsx::addWorksheet(wb, "Summary")
    summary_df <- data.frame(
      Item = c("Estimand", "Model", "Observations", "Treatment arm",
               "Control arm", "Events", "Trials"),
      Value = c(s$estimand, m$type, n_obs, n_treat, n_control,
                n_events, n_trials)
    )
    openxlsx::writeData(wb, "Summary", summary_df)

    # Coefficients sheet
    openxlsx::addWorksheet(wb, "Coefficients")
    openxlsx::writeData(wb, "Coefficients", coef_table)

    # Predictions sheet (optional)
    if (predictions && !is.null(obj$predictions)) {
      openxlsx::addWorksheet(wb, "Predictions")
      openxlsx::writeData(wb, "Predictions", obj$predictions)
    }

    openxlsx::saveWorkbook(wb, export, overwrite = TRUE)
    message("Results exported to: ", export)
  }

  invisible(obj)
}
