#' Validate prepared data for target trial emulation
#'
#' Runs 10 data quality checks on the prepared dataset and reports the
#' results. This is the second step in the TTE pipeline, used between
#' \code{\link{emulate_prepare}} and \code{\link{emulate_expand}} to catch data
#' problems before they cause errors or biased estimates downstream.
#'
#' @section The 10 validation checks:
#' Each check is reported as PASS, WARNING, or FAIL:
#' \enumerate{
#'   \item \strong{Person-period format}: Verifies that there are no
#'     duplicate \code{(id, period)} rows. Duplicates would mean the data is
#'     not in proper person-period format.
#'   \item \strong{No gaps in period sequences}: Checks that each
#'     individual's period values are consecutive (e.g., 0, 1, 2, 3 with no
#'     jumps like 0, 1, 3). Gaps may indicate missing data or data errors.
#'   \item \strong{Outcome is terminal}: Verifies that no rows exist after
#'     an outcome event for any individual. If the outcome is death or a
#'     first event, no further follow-up should be recorded.
#'   \item \strong{Treatment consistency}: For PP and AT estimands, reports
#'     the number of treatment discontinuations (switches from 1 to 0). For
#'     ITT, this check is informational only since treatment switching is
#'     expected.
#'   \item \strong{Missing data}: Checks for \code{NA} values in all core
#'     variables (id, period, treatment, outcome, eligible) and covariates.
#'     Missing data in these variables can cause model fitting failures.
#'   \item \strong{Eligibility consistency}: Checks that no individual is
#'     marked as eligible (\code{eligible = 1}) in a period after they have
#'     already experienced the outcome event.
#'   \item \strong{Sufficient eligible observations per period}: Warns if
#'     any time period has fewer than 10 eligible individuals, which can
#'     lead to unstable estimates.
#'   \item \strong{Positivity}: Checks that among eligible observations,
#'     there is variation in treatment (i.e., both treated and untreated
#'     individuals exist). Lack of positivity violates a key causal
#'     inference assumption.
#'   \item \strong{Period numbering}: Checks whether the period variable
#'     starts at 0 or 1 (the expected convention) and reports the range.
#'   \item \strong{Event rate}: Reports the number and rate of outcome
#'     events. Warns if there are fewer than 5 events, which may lead to
#'     unreliable model estimates.
#' }
#'
#' @section Strict mode:
#' By default (\code{strict = FALSE}), issues in checks 2, 3, 5, 6, and 7
#' are reported as warnings but do not stop execution. When
#' \code{strict = TRUE}, these issues are treated as errors and the function
#' will call \code{stop()} if any errors are found, preventing you from
#' proceeding with problematic data.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_prepare}}.
#'   The object must have its \code{prepared} state flag set to \code{TRUE}.
#' @param strict A logical value. If \code{TRUE}, warnings are promoted to
#'   errors and the function will stop with an error message if any check
#'   fails. If \code{FALSE} (the default), issues are reported as warnings
#'   and the function continues, allowing you to decide whether to proceed.
#' @param verbose A logical value. If \code{TRUE}, additional detail is
#'   printed for certain checks (e.g., per-variable missing value counts for
#'   check 5). Default is \code{FALSE}.
#'
#' @return The input \code{emulate} object, returned invisibly, with a
#'   \code{validation} list appended containing:
#'   \describe{
#'     \item{\code{n_checks}}{Integer: total number of checks run (always 10).}
#'     \item{\code{n_errors}}{Integer: number of checks that failed.}
#'     \item{\code{n_warnings}}{Integer: number of checks that produced
#'       warnings.}
#'     \item{\code{n_events}}{Integer: total number of outcome events.}
#'     \item{\code{event_rate}}{Numeric: percentage of observations with
#'       an event.}
#'   }
#'
#' @seealso \code{\link{emulate_prepare}} for the previous pipeline step,
#'   \code{\link{emulate_expand}} for the next step.
#'
#' @examples
#' \donttest{
#' # Create a small test dataset
#' set.seed(42)
#' dat <- data.frame(
#'   id = rep(1:20, each = 5),
#'   period = rep(0:4, times = 20),
#'   treatment = sample(0:1, 100, replace = TRUE),
#'   outcome = sample(c(rep(0, 19), 1), 100, replace = TRUE),
#'   eligible = 1L
#' )
#'
#' obj <- emulate_prepare(dat, id = "id", period = "period",
#'                    treatment = "treatment", outcome = "outcome",
#'                    eligible = "eligible", estimand = "ITT")
#'
#' # Run validation (default: warnings only)
#' obj <- emulate_validate(obj)
#'
#' # Run validation in strict mode (stops on any issue)
#' # obj <- emulate_validate(obj, strict = TRUE)
#'
#' # Verbose mode shows per-variable missing counts
#' obj <- emulate_validate(obj, verbose = TRUE)
#' }
#'
#' @export
emulate_validate <- function(obj, strict = FALSE, verbose = FALSE) {
  .check_prepared(obj)

  dt <- obj$data
  s <- obj$settings
  id_col <- s$id
  period_col <- s$period
  treat_col <- s$treatment
  out_col <- s$outcome
  elig_col <- s$eligible
  censor_col <- s$censor

  n_warnings <- 0L
  n_errors <- 0L
  n_checks <- 0L
  results <- list()

  .msg_check <- function(num, name) {
    message("Check ", num, ": ", name)
  }
  .pass <- function(detail = NULL) {
    if (is.null(detail)) message("  PASS") else message("  PASS (", detail, ")")
  }
  .warn <- function(msg) {
    message("  WARNING: ", msg)
    n_warnings <<- n_warnings + 1L
  }
  .fail <- function(msg) {
    message("  FAIL: ", msg)
    n_errors <<- n_errors + 1L
  }

  message("\nemulate_validate - Data Quality Checks")
  message(strrep("-", 50))

  # Check 1: Person-period format
  n_checks <- n_checks + 1L
  .msg_check(1, "Person-period format")
  dup_n <- dt[, .N, by = c(id_col, period_col)][N > 1, .N]
  if (dup_n > 0) .fail(paste(dup_n, "duplicate (id, period) rows")) else .pass()

  # Check 2: No gaps in period sequences
  n_checks <- n_checks + 1L
  .msg_check(2, "No gaps in period sequences")
  dt_sorted <- data.table::copy(dt)
  data.table::setorderv(dt_sorted, c(id_col, period_col))
  dt_sorted[, pdiff := get(period_col) - shift(get(period_col), 1L),
            by = c(id_col)]
  n_gaps <- sum(dt_sorted$pdiff > 1, na.rm = TRUE)
  dt_sorted[, pdiff := NULL]
  if (n_gaps > 0) {
    if (strict) .fail(paste(n_gaps, "gaps in period sequences"))
    else .warn(paste(n_gaps, "gaps in period sequences"))
  } else {
    .pass()
  }

  # Check 3: Outcome is terminal
  n_checks <- n_checks + 1L
  .msg_check(3, "Outcome is terminal (no rows after event)")
  dt_sorted2 <- data.table::copy(dt)
  data.table::setorderv(dt_sorted2, c(id_col, period_col))
  dt_sorted2[, cum_out := cumsum(get(out_col)), by = c(id_col)]
  dt_sorted2[, prior_out := shift(cum_out, 1L, fill = 0L), by = c(id_col)]
  n_post <- sum(dt_sorted2$prior_out >= 1, na.rm = TRUE)
  if (n_post > 0) {
    if (strict) .fail(paste(n_post, "rows found after outcome event"))
    else .warn(paste(n_post, "rows found after outcome event"))
  } else {
    .pass()
  }

  # Check 4: Treatment consistency
  n_checks <- n_checks + 1L
  .msg_check(4, "Treatment consistency")
  if (s$estimand %in% c("PP", "AT")) {
    dt_sorted3 <- data.table::copy(dt)
    data.table::setorderv(dt_sorted3, c(id_col, period_col))
    dt_sorted3[, prev_treat := shift(get(treat_col), 1L), by = c(id_col)]
    n_switchoff <- sum(dt_sorted3$prev_treat == 1 & dt_sorted3[[treat_col]] == 0,
                       na.rm = TRUE)
    if (n_switchoff > 0) {
      message("  NOTE: ", n_switchoff, " treatment discontinuations found")
    } else {
      .pass("treatment is absorbing")
    }
  } else {
    .pass("ITT: treatment switching is allowed")
  }

  # Check 5: Missing data
  n_checks <- n_checks + 1L
  .msg_check(5, "Missing data")
  core_vars <- c(id_col, period_col, treat_col, out_col, elig_col)
  any_missing <- FALSE
  for (v in c(core_vars, s$covariates, s$baseline_covariates)) {
    n_miss <- sum(is.na(dt[[v]]))
    if (n_miss > 0) {
      any_missing <- TRUE
      if (verbose) message("    ", v, ": ", n_miss, " missing values")
    }
  }
  if (any_missing) {
    if (strict) .fail("missing values in core variables")
    else .warn("missing values found")
  } else {
    .pass("no missing values in core variables")
  }

  # Check 6: Eligibility consistency
  n_checks <- n_checks + 1L
  .msg_check(6, "Eligibility consistency")
  dt_sorted4 <- data.table::copy(dt)
  data.table::setorderv(dt_sorted4, c(id_col, period_col))
  dt_sorted4[, prior_out2 := shift(cumsum(get(out_col)), 1L, fill = 0L),
             by = c(id_col)]
  n_elig_post <- sum(dt_sorted4[[elig_col]] == 1 & dt_sorted4$prior_out2 > 0,
                     na.rm = TRUE)
  if (n_elig_post > 0) {
    if (strict) .fail(paste(n_elig_post, "eligible obs after prior outcome"))
    else .warn(paste(n_elig_post, "eligible observations after prior outcome"))
  } else {
    .pass()
  }

  # Check 7: Sufficient eligible observations per period
  n_checks <- n_checks + 1L
  .msg_check(7, "Sufficient eligible observations per period")
  total_eligible <- sum(dt[[elig_col]] == 1, na.rm = TRUE)
  if (total_eligible == 0) {
    .fail("no eligible observations")
  } else {
    elig_per_period <- dt[get(elig_col) == 1, .N, by = c(period_col)]
    min_elig <- min(elig_per_period$N)
    if (min_elig < 10) {
      .warn(paste0("some periods have fewer than 10 eligible individuals (min: ",
                    min_elig, ")"))
    } else {
      .pass(paste("min eligible per period:", min_elig))
    }
  }

  # Check 8: Positivity
  n_checks <- n_checks + 1L
  .msg_check(8, "Positivity (treatment variation among eligible)")
  elig_rows <- dt[get(elig_col) == 1]
  n_elig_treat <- sum(elig_rows[[treat_col]] == 1, na.rm = TRUE)
  n_elig_untreat <- sum(elig_rows[[treat_col]] == 0, na.rm = TRUE)
  if (n_elig_treat == 0 || n_elig_untreat == 0) {
    .fail(paste0("no treatment variation (treated: ", n_elig_treat,
                 ", untreated: ", n_elig_untreat, ")"))
  } else {
    pct <- 100 * n_elig_treat / (n_elig_treat + n_elig_untreat)
    .pass(sprintf("treatment prevalence: %.1f%%", pct))
  }

  # Check 9: Period numbering
  n_checks <- n_checks + 1L
  .msg_check(9, "Period numbering")
  p_min <- min(dt[[period_col]], na.rm = TRUE)
  p_max <- max(dt[[period_col]], na.rm = TRUE)
  if (!p_min %in% c(0, 1)) {
    message("  NOTE: period starts at ", p_min, " (expected 0 or 1)")
  } else {
    .pass(paste("period range:", p_min, "to", p_max))
  }

  # Check 10: Event rate
  n_checks <- n_checks + 1L
  .msg_check(10, "Event rate")
  n_events <- sum(dt[[out_col]] == 1, na.rm = TRUE)
  event_rate <- 100 * n_events / nrow(dt)
  if (n_events < 5) {
    .warn(paste0("very few events (", n_events, "); estimates may be unreliable"))
  } else {
    .pass(sprintf("%d events, rate: %.2f%%", n_events, event_rate))
  }

  # Summary
  message(strrep("-", 50))
  message("Checks: ", n_checks, "  Errors: ", n_errors,
          "  Warnings: ", n_warnings)
  if (n_errors > 0) {
    msg <- "Data validation failed. Fix errors before proceeding."
    if (strict) stop(msg, call. = FALSE)
    message(msg)
  } else if (n_warnings > 0) {
    message("Passed with warnings.")
  } else {
    message("All checks passed.")
  }
  message(strrep("-", 50))

  obj$validation <- list(
    n_checks = n_checks,
    n_errors = n_errors,
    n_warnings = n_warnings,
    n_events = n_events,
    event_rate = event_rate
  )

  invisible(obj)
}
