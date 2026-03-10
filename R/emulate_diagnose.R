#' Weight diagnostics and covariate balance assessment
#'
#' Reports weight distribution statistics, effective sample size (ESS), and
#' optionally computes standardized mean differences (SMDs) for covariate
#' balance between treatment arms. This function helps you assess whether
#' the inverse probability weights from \code{\link{emulate_weight}} are
#' well-behaved and whether weighting has achieved adequate covariate balance.
#'
#' @section Weight distribution:
#' Reports the mean, standard deviation, minimum, maximum, and selected
#' percentiles (1st, 5th, 25th, 50th, 75th, 95th, 99th) of the weight
#' distribution. Ideally, stabilized weights should have a mean close to 1.
#' Very large maximum values or high variability (large SD) suggest
#' near-violations of the positivity assumption and may indicate that
#' truncation is needed.
#'
#' @section Effective sample size (ESS):
#' The ESS measures the "effective" number of independent observations after
#' weighting:
#' \deqn{ESS = \frac{(\sum w_i)^2}{\sum w_i^2}}{ESS = (sum(w))^2 / sum(w^2)}
#' If all weights are equal (e.g., all 1.0), the ESS equals the actual
#' sample size. Highly variable weights reduce the ESS, sometimes
#' dramatically. An ESS much smaller than the actual sample size (e.g.,
#' ESS/N < 0.5) indicates that a few observations dominate the analysis,
#' which can lead to unstable estimates. The ESS is reported overall and
#' separately for each treatment arm.
#'
#' @section Covariate balance (SMD):
#' When \code{balance_covariates} is provided, the function computes the
#' \strong{standardized mean difference} (SMD) for each covariate, both
#' before and after weighting:
#' \deqn{SMD = \frac{\bar{x}_1 - \bar{x}_0}{\sqrt{(s_1^2 + s_0^2)/2}}}{
#'   SMD = (mean_1 - mean_0) / sqrt((var_1 + var_0) / 2)}
#' A common threshold for adequate balance is \code{|SMD| < 0.1}. The
#' function reports both unweighted and weighted SMDs so you can see whether
#' weighting improved balance. A "Love plot" of these SMDs can be produced
#' using \code{\link{emulate_plot}(obj, type = "balance")}.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_weight}}. The
#'   \code{expanded} state flag must be \code{TRUE} and a weight variable
#'   must exist in the data.
#' @param balance_covariates An optional character vector of covariate names
#'   for which to compute standardized mean differences. These should be
#'   numeric variables (binary or continuous). If \code{NULL} (the default),
#'   no balance assessment is performed (only weight distribution and ESS
#'   are reported).
#' @param by_trial A logical value. If \code{TRUE}, the weight distribution
#'   (mean, SD, min, max) is reported separately for each emulated trial
#'   period. This can help identify specific trials with problematic weights.
#'   Default is \code{FALSE}.
#'
#' @return The input \code{emulate} object, returned invisibly, with the
#'   \code{diagnostics} slot populated with:
#'   \describe{
#'     \item{weight_mean / weight_sd / weight_min / weight_max}{Summary
#'       statistics of the weight distribution.}
#'     \item{weight_quantiles}{Named numeric vector of weight
#'       percentiles (1 pct, 5 pct, 25 pct, 50 pct, 75 pct, 95 pct, 99 pct).}
#'     \item{ess}{Overall effective sample size.}
#'     \item{ess_by_arm}{Numeric vector of length 2 giving the ESS
#'       for arm 0 (control) and arm 1 (treatment).}
#'     \item{balance}{A \code{data.frame} with columns
#'       \code{covariate}, \code{smd_unwt} (unweighted SMD), \code{smd_wt}
#'       (weighted SMD), and \code{threshold} (0.1). \code{NULL} if
#'       \code{balance_covariates} was not provided.}
#'     \item{trial_weights}{A \code{data.frame} with per-trial weight
#'       statistics. \code{NULL} if \code{by_trial = FALSE}.}
#'   }
#'
#' @seealso \code{\link{emulate_weight}} for computing weights,
#'   \code{\link{emulate_plot}} with \code{type = "weights"} or
#'   \code{type = "balance"} for visualization.
#'
#' @references
#' Austin PC, Stuart EA (2015). "Moving towards best practice when using
#' inverse probability of treatment weighting (IPTW) using the propensity
#' score to estimate causal treatment effects in observational studies."
#' \emph{Statistics in Medicine}, 34(28), 3661-3679.
#'
#' @examples
#' \donttest{
#' # Build a weighted dataset
#' set.seed(42)
#' n_ids <- 30
#' n_per <- 8
#' dat <- data.frame(
#'   id = rep(seq_len(n_ids), each = n_per),
#'   period = rep(0:(n_per - 1), times = n_ids),
#'   treatment = 0L, outcome = 0L, eligible = 1L,
#'   age = rep(rnorm(n_ids, 50, 10), each = n_per),
#'   sex = rep(sample(0:1, n_ids, replace = TRUE), each = n_per)
#' )
#' dat$treatment <- ifelse(dat$period >= sample(2:6, nrow(dat), replace = TRUE), 1L, 0L)
#' dat$outcome[sample(which(dat$period > 4), 4)] <- 1L
#'
#' obj <- emulate_prepare(dat, id = "id", period = "period",
#'                    treatment = "treatment", outcome = "outcome",
#'                    eligible = "eligible", covariates = c("age", "sex"),
#'                    estimand = "PP")
#' obj <- emulate_expand(obj, maxfollowup = 5)
#' obj <- emulate_weight(obj, switch_d_cov = c("age", "sex"))
#'
#' # Weight diagnostics only
#' obj <- emulate_diagnose(obj)
#'
#' # Weight diagnostics + covariate balance
#' obj <- emulate_diagnose(obj, balance_covariates = c("age", "sex"))
#'
#' # Per-trial weight distribution
#' obj <- emulate_diagnose(obj, by_trial = TRUE)
#' }
#'
#' @export
emulate_diagnose <- function(obj, balance_covariates = NULL, by_trial = FALSE) {
  .check_expanded(obj)

  s <- obj$settings
  dt <- obj$data
  prefix <- s$prefix
  id_col <- s$id

  arm_col <- paste0(prefix, "arm")
  trial_col <- paste0(prefix, "trial")
  cens_col <- paste0(prefix, "censored")
  weight_col <- paste0(prefix, "weight")

  has_weights <- weight_col %in% names(dt)
  if (!has_weights) {
    stop("No weight variable found. Run emulate_weight() first.", call. = FALSE)
  }

  w <- dt[[weight_col]]

  message("emulate_diagnose - Weight Diagnostics")
  message(strrep("-", 50))

  # Weight distribution
  w_mean <- mean(w, na.rm = TRUE)
  w_sd <- sd(w, na.rm = TRUE)
  w_min <- min(w, na.rm = TRUE)
  w_max <- max(w, na.rm = TRUE)
  w_q <- quantile(w, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99),
                  na.rm = TRUE)

  message("Weight distribution:")
  message("  Mean:   ", sprintf("%.4f", w_mean))
  message("  SD:     ", sprintf("%.4f", w_sd))
  message("  Min:    ", sprintf("%.4f", w_min))
  message("  P1:     ", sprintf("%.4f", w_q["1%"]))
  message("  P25:    ", sprintf("%.4f", w_q["25%"]))
  message("  Median: ", sprintf("%.4f", w_q["50%"]))
  message("  P75:    ", sprintf("%.4f", w_q["75%"]))
  message("  P99:    ", sprintf("%.4f", w_q["99%"]))
  message("  Max:    ", sprintf("%.4f", w_max))

  # ESS overall
  sum_w <- sum(w, na.rm = TRUE)
  sum_w2 <- sum(w^2, na.rm = TRUE)
  ess <- sum_w^2 / sum_w2

  # ESS by arm
  ess_by_arm <- numeric(2)
  for (a in 0:1) {
    wa <- w[dt[[arm_col]] == a]
    sw <- sum(wa, na.rm = TRUE)
    sw2 <- sum(wa^2, na.rm = TRUE)
    ess_by_arm[a + 1] <- sw^2 / sw2
  }

  message("\nEffective Sample Size:")
  message("  Overall: ", sprintf("%.1f", ess), " / ", length(w))
  message("  Arm 0:   ", sprintf("%.1f", ess_by_arm[1]))
  message("  Arm 1:   ", sprintf("%.1f", ess_by_arm[2]))

  # By trial (optional)
  trial_weights <- NULL
  if (by_trial) {
    message("\nWeight distribution by trial:")
    trials <- sort(unique(dt[[trial_col]]))
    trial_weights <- data.frame(
      trial = trials,
      mean = NA_real_, sd = NA_real_,
      min = NA_real_, max = NA_real_
    )
    for (i in seq_along(trials)) {
      wt <- w[dt[[trial_col]] == trials[i]]
      trial_weights$mean[i] <- mean(wt, na.rm = TRUE)
      trial_weights$sd[i] <- sd(wt, na.rm = TRUE)
      trial_weights$min[i] <- min(wt, na.rm = TRUE)
      trial_weights$max[i] <- max(wt, na.rm = TRUE)
    }
    print(trial_weights, row.names = FALSE)
  }

  # Covariate balance
  balance <- NULL
  if (!is.null(balance_covariates) && length(balance_covariates) > 0) {
    message("\nCovariate Balance (SMD):")
    message(sprintf("  %-20s %10s %10s %10s", "Covariate",
                    "Unweighted", "Weighted", "Threshold"))
    message(strrep("-", 55))

    balance <- data.frame(
      covariate = balance_covariates,
      smd_unwt = NA_real_,
      smd_wt = NA_real_,
      threshold = 0.1,
      stringsAsFactors = FALSE
    )

    arm0 <- dt[[arm_col]] == 0
    arm1 <- dt[[arm_col]] == 1

    for (i in seq_along(balance_covariates)) {
      v <- balance_covariates[i]
      x <- dt[[v]]

      # Unweighted SMD
      m0 <- mean(x[arm0], na.rm = TRUE)
      m1 <- mean(x[arm1], na.rm = TRUE)
      v0 <- var(x[arm0], na.rm = TRUE)
      v1 <- var(x[arm1], na.rm = TRUE)
      pooled_sd <- sqrt((v0 + v1) / 2)
      smd_unwt <- if (pooled_sd > 0) (m1 - m0) / pooled_sd else 0

      # Weighted SMD
      wm0 <- weighted.mean(x[arm0], w[arm0], na.rm = TRUE)
      wm1 <- weighted.mean(x[arm1], w[arm1], na.rm = TRUE)
      # Weighted variance
      wv0 <- .weighted_var(x[arm0], w[arm0])
      wv1 <- .weighted_var(x[arm1], w[arm1])
      wpooled_sd <- sqrt((wv0 + wv1) / 2)
      smd_wt <- if (wpooled_sd > 0) (wm1 - wm0) / wpooled_sd else 0

      balance$smd_unwt[i] <- smd_unwt
      balance$smd_wt[i] <- smd_wt

      message(sprintf("  %-20s %10.4f %10.4f %10.1f",
                      v, abs(smd_unwt), abs(smd_wt), 0.1))
    }

    max_smd <- max(abs(balance$smd_wt))
    balanced <- max_smd < 0.1
    message(strrep("-", 55))
    message("Max weighted SMD: ", sprintf("%.4f", max_smd))
    message("Balance achieved: ", if (balanced) "Yes" else "No")
  }

  message(strrep("-", 50))

  obj$diagnostics <- list(
    weight_mean = w_mean,
    weight_sd = w_sd,
    weight_min = w_min,
    weight_max = w_max,
    weight_quantiles = w_q,
    ess = ess,
    ess_by_arm = ess_by_arm,
    balance = balance,
    trial_weights = trial_weights
  )

  invisible(obj)
}

#' Compute weighted variance (internal helper)
#'
#' Computes the weighted variance of a numeric vector using the formula:
#' \code{sum(w * (x - weighted_mean)^2) / sum(w)}. This is a population
#' variance (denominator is \code{sum(w)}, not \code{sum(w) - 1}).
#'
#' @param x A numeric vector of values.
#' @param w A numeric vector of non-negative weights, same length as \code{x}.
#'
#' @return A single numeric value: the weighted variance. Returns 0 if
#'   fewer than 2 non-missing observations.
#'
#' @keywords internal
.weighted_var <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  if (length(x) < 2) return(0)
  wm <- weighted.mean(x, w)
  sum(w * (x - wm)^2) / sum(w)
}
