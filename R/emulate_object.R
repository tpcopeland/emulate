#' Create a new emulate object (internal constructor)
#'
#' This is the internal constructor for the \code{emulate} S3 class. It is called
#' by \code{\link{emulate_prepare}} and is not intended to be called directly by
#' users.
#'
#' A \code{emulate} object is a named list with the following components:
#' \describe{
#'   \item{\code{data}}{A \code{data.table} containing the analysis dataset.
#'     This is updated at each pipeline step (prepare, expand, weight, fit).}
#'   \item{\code{settings}}{A list of variable mappings and configuration
#'     options set during \code{\link{emulate_prepare}}, including \code{id},
#'     \code{period}, \code{treatment}, \code{outcome}, \code{eligible},
#'     \code{censor}, \code{covariates}, \code{baseline_covariates},
#'     \code{estimand}, and \code{prefix}.}
#'   \item{\code{state}}{A list of logical flags tracking which pipeline
#'     steps have been completed: \code{prepared}, \code{expanded},
#'     \code{weighted}, \code{fitted}.}
#'   \item{\code{expansion}}{A list of expansion metadata populated by
#'     \code{\link{emulate_expand}}: number of trials, trial periods, counts
#'     of expanded observations, treatment/control arms, censored, events,
#'     and the expansion ratio.}
#'   \item{\code{weights}}{A list of weighting metadata populated by
#'     \code{\link{emulate_weight}}: weight variable name, truncation bounds,
#'     mean, SD, min, max, effective sample size (ESS), and count of
#'     truncated observations.}
#'   \item{\code{model}}{A list of model metadata populated by
#'     \code{\link{emulate_fit}}: model type, the fitted model object, the
#'     cluster-robust variance-covariance matrix, specification details,
#'     and the treatment effect estimate.}
#'   \item{\code{predictions}}{A \code{data.frame} of marginal predictions
#'     populated by \code{\link{emulate_predict}}, or \code{NULL} if predictions
#'     have not been computed.}
#'   \item{\code{diagnostics}}{A list of diagnostic results populated by
#'     \code{\link{emulate_diagnose}}: weight distribution statistics, ESS by
#'     arm, covariate balance, and per-trial weight summaries.}
#' }
#'
#' @param data A \code{data.table} containing the analysis dataset. This
#'   should be person-period data that has already been validated.
#' @param settings A named list of variable mappings and configuration. Must
#'   include elements: \code{id}, \code{period}, \code{treatment},
#'   \code{outcome}, \code{eligible}, \code{censor}, \code{covariates},
#'   \code{baseline_covariates}, \code{estimand}, and \code{prefix}.
#'
#' @return A \code{emulate} object (an S3 list of class \code{"emulate"}) with all
#'   pipeline state flags initialized to \code{FALSE} and all result slots
#'   set to empty lists or \code{NULL}.
#'
#' @keywords internal
new_emulate <- function(data, settings) {
  obj <- list(
    data         = data,
    settings     = settings,
    state        = list(prepared = FALSE, expanded = FALSE,
                        weighted = FALSE, matched = FALSE,
                        stratified = FALSE, fitted = FALSE,
                        calibrated = FALSE),
    expansion    = list(),
    weights      = list(),
    model        = list(),
    predictions  = NULL,
    diagnostics  = list(),
    matching     = list(),
    stratification = list(),
    calibration  = list()
  )
  class(obj) <- "emulate"
  obj
}

#' Print method for emulate objects
#'
#' Displays a concise summary of a \code{emulate} object, including the estimand,
#' variable mappings, the current pipeline state (which steps have been
#' completed), observation and individual counts, and model information if a
#' model has been fitted.
#'
#' @param x A \code{emulate} object created by \code{\link{emulate_prepare}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The \code{emulate} object \code{x}, returned invisibly.
#'
#' @examples
#' \donttest{
#' # Create a small test dataset
#' dat <- data.frame(
#'   id = rep(1:10, each = 5),
#'   period = rep(0:4, times = 10),
#'   treatment = sample(0:1, 50, replace = TRUE),
#'   outcome = sample(c(0, 0, 0, 0, 1), 50, replace = TRUE),
#'   eligible = 1L
#' )
#' obj <- emulate_prepare(dat, id = "id", period = "period",
#'                    treatment = "treatment", outcome = "outcome",
#'                    eligible = "eligible", estimand = "ITT")
#' print(obj)
#' }
#'
#' @export
print.emulate <- function(x, ...) {
  cat("Target Trial Emulation Object\n")
  cat(strrep("-", 50), "\n")

  s <- x$settings
  cat("Estimand:    ", s$estimand, "\n")
  cat("ID:          ", s$id, "\n")
  cat("Period:      ", s$period, "\n")
  cat("Treatment:   ", s$treatment, "\n")
  cat("Outcome:     ", s$outcome, "\n")

  # Pipeline state
  states <- c(
    if (x$state$prepared)    "prepared",
    if (x$state$expanded)    "expanded",
    if (x$state$weighted)    "weighted",
    if (isTRUE(x$state$matched))     "matched",
    if (isTRUE(x$state$stratified))  "stratified",
    if (x$state$fitted)      "fitted",
    if (isTRUE(x$state$calibrated))  "calibrated"
  )
  cat("Pipeline:    ", if (length(states)) paste(states, collapse = " > ") else "not started", "\n")

  cat("Observations:", nrow(x$data), "\n")
  cat("Individuals: ", data.table::uniqueN(x$data[[s$id]]), "\n")

  if (x$state$expanded) {
    cat("Trials:      ", x$expansion$n_trials, "\n")
  }
  if (x$state$fitted) {
    cat("Model:       ", x$model$type, "\n")
  }
  invisible(x)
}

#' Summary method for emulate objects
#'
#' Displays everything shown by \code{\link{print.emulate}}, plus additional
#' detail about the weight distribution (mean, SD, range, and effective sample
#' size) if weights have been computed, and the full coefficient table if a
#' model has been fitted.
#'
#' @param object A \code{emulate} object created by \code{\link{emulate_prepare}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The \code{emulate} object \code{object}, returned invisibly.
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   id = rep(1:10, each = 5),
#'   period = rep(0:4, times = 10),
#'   treatment = sample(0:1, 50, replace = TRUE),
#'   outcome = sample(c(0, 0, 0, 0, 1), 50, replace = TRUE),
#'   eligible = 1L
#' )
#' obj <- emulate_prepare(dat, id = "id", period = "period",
#'                    treatment = "treatment", outcome = "outcome",
#'                    eligible = "eligible", estimand = "ITT")
#' summary(obj)
#' }
#'
#' @export
summary.emulate <- function(object, ...) {
  print(object)

  if (object$state$weighted) {
    cat("\nWeight Summary:\n")
    w <- object$data[[object$weights$weight_var]]
    cat("  Mean:   ", sprintf("%.4f", mean(w, na.rm = TRUE)), "\n")
    cat("  SD:     ", sprintf("%.4f", sd(w, na.rm = TRUE)), "\n")
    cat("  Range:  ", sprintf("%.4f - %.4f", min(w, na.rm = TRUE), max(w, na.rm = TRUE)), "\n")
    ess <- sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)
    cat("  ESS:    ", sprintf("%.1f", ess), "\n")
  }

  if (object$state$fitted) {
    cat("\nModel Coefficients:\n")
    print(coef(object$model$object))
  }

  invisible(object)
}
