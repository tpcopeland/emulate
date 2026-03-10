#' Propensity score matching for target trial emulation
#'
#' An alternative to inverse probability weighting that matches treated and
#' control observations on propensity scores. Matching subsets the data to
#' matched pairs (or sets) and fits the outcome model on the matched sample
#' without IP weights.
#'
#' @section How matching works:
#' \enumerate{
#'   \item Propensity scores are estimated on the expanded data using logistic
#'     regression or lasso (via \code{ps_method}).
#'   \item \code{MatchIt::matchit()} performs nearest-neighbor matching on the
#'     propensity score (logit scale by default).
#'   \item The matched data is extracted and the outcome model is fitted on it.
#' }
#'
#' @section Estimand:
#' By default, matching estimates the ATT (Average Treatment Effect on the
#' Treated). The \code{estimand_match} argument controls this. Note that this
#' is the matching estimand, which may differ from the trial estimand (ITT/PP/AT)
#' set in \code{\link{emulate_prepare}}.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_expand}}.
#'   The \code{expanded} state flag must be \code{TRUE}.
#' @param match_cov A character vector of covariate names for propensity score
#'   estimation. Required.
#' @param match_method Matching method passed to \code{MatchIt::matchit()}.
#'   Default is \code{"nearest"} (nearest-neighbor).
#' @param match_ratio Number of control matches per treated unit. Default is
#'   \code{5L} (5:1 matching).
#' @param caliper Caliper width in standard deviations of the propensity score
#'   (logit scale). Default is \code{0.2}.
#' @param ps_method Propensity score estimation method: \code{"glm"} (default)
#'   or \code{"lasso"}.
#' @param estimand_match Estimand for matching: \code{"ATT"} (default) or
#'   \code{"ATE"}.
#' @param outcome_cov Optional character vector of covariates for the outcome
#'   model.
#' @param model Outcome model type: \code{"logistic"} (default), \code{"cox"},
#'   or \code{"poisson"}.
#' @param followup_spec How follow-up time enters the model. Default
#'   \code{"quadratic"}.
#' @param trial_period_spec How trial period enters the model. Default
#'   \code{"quadratic"}.
#' @param cluster Clustering variable for robust SEs. Default: ID variable.
#' @param level Confidence level as percentage. Default \code{95}.
#' @param quiet Suppress messages. Default \code{FALSE}.
#'
#' @return The input \code{emulate} object with:
#'   \itemize{
#'     \item \code{state$matched} set to \code{TRUE}
#'     \item \code{state$fitted} set to \code{TRUE}
#'     \item \code{matching} slot populated with MatchIt object and PS values
#'     \item \code{model} slot populated with fitted outcome model
#'   }
#'
#' @seealso \code{\link{emulate_weight}} for the IPTW alternative,
#'   \code{\link{emulate_stratify}} for stratification.
#'
#' @export
emulate_match <- function(obj, match_cov, match_method = "nearest",
                          match_ratio = 5L, caliper = 0.2,
                          ps_method = "glm", estimand_match = "ATT",
                          outcome_cov = NULL, model = "logistic",
                          followup_spec = "quadratic",
                          trial_period_spec = "quadratic",
                          cluster = NULL, level = 95, quiet = FALSE) {
  .check_expanded(obj)

  if (missing(match_cov) || is.null(match_cov) || length(match_cov) == 0) {
    stop("match_cov is required for propensity score matching", call. = FALSE)
  }

  if (!ps_method %in% c("glm", "lasso")) {
    stop("ps_method must be 'glm' or 'lasso'", call. = FALSE)
  }

  if (!estimand_match %in% c("ATT", "ATE")) {
    stop("estimand_match must be 'ATT' or 'ATE'", call. = FALSE)
  }

  s <- obj$settings
  dt <- data.table::copy(obj$data)
  prefix <- s$prefix
  id_col <- s$id
  treat_col <- s$treatment

  arm_col <- paste0(prefix, "arm")
  fu_col <- paste0(prefix, "followup")
  trial_col <- paste0(prefix, "trial")
  cens_col <- paste0(prefix, "censored")

  if (is.null(cluster)) cluster <- id_col

  message("emulate_match - Propensity Score Matching")
  message(strrep("-", 50))
  message("Method:        ", match_method)
  message("Ratio:         ", match_ratio, ":1")
  message("Caliper:       ", caliper)
  message("PS method:     ", ps_method)
  message("Estimand:      ", estimand_match)

  # Build PS formula
  bq <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))
  ps_formula <- as.formula(paste(bq(arm_col), "~",
                                  paste(bq(match_cov), collapse = " + ")))

  # Fit PS
  if (ps_method == "lasso") {
    message("Fitting lasso propensity score model...")
    ps_result <- .fit_ps_lasso(ps_formula, data = as.data.frame(dt))
    dt[, .ps_match := ps_result$ps]
  } else {
    message("Fitting GLM propensity score model...")
    ps_fit <- tryCatch(
      suppressWarnings(glm(ps_formula, data = dt, family = binomial())),
      error = function(e) NULL
    )
    if (!is.null(ps_fit)) {
      dt[, .ps_match := predict(ps_fit, type = "response")]
    } else {
      stop("GLM propensity score model failed to converge", call. = FALSE)
    }
  }

  # Run MatchIt
  message("Running matching...")
  match_f <- as.formula(paste(bq(arm_col), "~",
                               paste(bq(match_cov), collapse = " + ")))

  # MatchIt needs a data.frame
  dt_df <- as.data.frame(dt)

  m_out <- tryCatch(
    MatchIt::matchit(match_f, data = dt_df,
                     method = match_method,
                     distance = dt$.ps_match,
                     ratio = match_ratio,
                     caliper = caliper,
                     estimand = estimand_match),
    error = function(e) {
      stop("Matching failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  # Extract matched data
  matched_df <- MatchIt::match.data(m_out)
  matched_dt <- data.table::as.data.table(matched_df)

  n_orig <- nrow(dt)
  n_matched <- nrow(matched_dt)
  message("Matched:       ", n_matched, " / ", n_orig, " observations")

  # Store PS values before cleanup
  ps_values <- dt$.ps_match

  # Clean up temp column
  if (".ps_match" %in% names(matched_dt)) {
    matched_dt[, .ps_match := NULL]
  }

  # Fit outcome model on matched data (no weights)
  obj$data <- matched_dt
  obj$state$matched <- TRUE
  obj$matching <- list(
    matchit_object = m_out,
    ps_values = ps_values,
    n_original = n_orig,
    n_matched = n_matched,
    match_method = match_method,
    match_ratio = match_ratio,
    caliper = caliper,
    estimand_match = estimand_match,
    ps_method = ps_method
  )

  # Now fit the outcome model on matched data
  message("Fitting outcome model on matched data...")
  obj <- emulate_fit(obj, outcome_cov = outcome_cov, model = model,
                     followup_spec = followup_spec,
                     trial_period_spec = trial_period_spec,
                     cluster = cluster, level = level)

  obj$model$adjustment_method <- "matching"

  message(strrep("-", 50))
  invisible(obj)
}
