#' Marginal predictions with Monte Carlo confidence intervals
#'
#' Computes G-formula standardized marginal predictions for both treatment
#' arms at specified follow-up times, with Monte Carlo confidence intervals.
#' This is the step that translates the fitted model into clinically
#' interpretable survival or cumulative incidence curves.
#'
#' @section G-formula standardization:
#' The G-formula (also called the parametric g-formula or g-computation)
#' works as follows:
#' \enumerate{
#'   \item Use the fitted outcome model to predict each individual's
#'     probability of the outcome at each follow-up time, under each
#'     treatment strategy (arm = 0 and arm = 1).
#'   \item Convert these per-period probabilities into cumulative survival
#'     curves: \code{S(t) = product of (1 - p_s) for s = 0 to t}.
#'   \item \strong{Average} across all individuals in the reference
#'     population to get the marginal (population-average) survival curve
#'     for each arm.
#' }
#' This averaging step is what makes the predictions "standardized" -- they
#' represent the expected outcome in the entire study population, not just
#' for a specific covariate profile. The reference population consists of
#' individuals at follow-up time 0 in the estimation sample.
#'
#' @section Cumulative incidence vs. survival:
#' The \code{type} argument controls which scale the results are reported on:
#' \describe{
#'   \item{\code{"cum_inc"}}{(Default) Cumulative incidence = \code{1 - S(t)}.
#'     This is the probability of having experienced the event by time t.
#'     Also known as the "risk" or "failure probability.""}
#'   \item{\code{"survival"}}{Survival probability = \code{S(t)}. This is the
#'     probability of being event-free at time t.}
#' }
#'
#' @section Monte Carlo confidence intervals:
#' Confidence intervals are computed by simulation:
#' \enumerate{
#'   \item Draw a set of coefficients from the multivariate normal
#'     distribution centered at the estimated coefficients with the
#'     cluster-robust variance-covariance matrix.
#'   \item Recompute the marginal predictions using the drawn coefficients.
#'   \item Repeat \code{samples} times.
#'   \item Take the appropriate percentiles of the simulated predictions as
#'     the confidence interval bounds.
#' }
#' The default of 100 samples is a reasonable starting point; use 500-1000
#' for final results.
#'
#' @section Risk differences:
#' When \code{difference = TRUE}, the function also computes the
#' \strong{risk difference} (or survival difference) at each time point:
#' \code{estimate_arm1 - estimate_arm0}. The confidence interval for the
#' difference is computed from the Monte Carlo samples of the difference.
#' This provides a direct measure of the treatment effect on the absolute
#' scale.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_fit}}. The
#'   \code{fitted} state flag must be \code{TRUE}. Currently, predictions are
#'   only supported for the pooled logistic model (\code{model = "logistic"}).
#' @param times An integer vector of follow-up times at which to compute
#'   predictions. For example, \code{times = c(1, 2, 3, 4, 5)} computes
#'   predictions at follow-up times 1 through 5. Must be non-negative
#'   integers.
#' @param type A character string specifying the prediction scale. One of:
#'   \describe{
#'     \item{\code{"cum_inc"}}{(Default) Cumulative incidence (1 - survival).}
#'     \item{\code{"survival"}}{Survival probability.}
#'   }
#' @param samples An integer specifying the number of Monte Carlo simulation
#'   draws for confidence interval estimation. Default is \code{100}. Must be
#'   at least 10. Higher values give more precise CIs but take longer. Use
#'   500-1000 for publication-quality results.
#' @param seed An optional integer specifying the random seed for
#'   reproducibility of Monte Carlo simulations. If \code{NULL} (the
#'   default), no seed is set.
#' @param level A numeric value specifying the confidence level as a
#'   percentage. Default is \code{95} (for 95 percent confidence intervals).
#' @param difference A logical value. If \code{TRUE}, the risk (or survival)
#'   difference between arm 1 and arm 0 is computed at each time point,
#'   along with its Monte Carlo confidence interval. Default is \code{FALSE}.
#'
#' @return The input \code{emulate} object with the \code{predictions} slot set
#'   to a \code{data.frame} with the following columns:
#'   \describe{
#'     \item{\code{time}}{The follow-up time.}
#'     \item{\code{est_0}}{Point estimate for arm 0 (control).}
#'     \item{\code{ci_lo_0}}{Lower confidence bound for arm 0.}
#'     \item{\code{ci_hi_0}}{Upper confidence bound for arm 0.}
#'     \item{\code{est_1}}{Point estimate for arm 1 (treatment).}
#'     \item{\code{ci_lo_1}}{Lower confidence bound for arm 1.}
#'     \item{\code{ci_hi_1}}{Upper confidence bound for arm 1.}
#'     \item{\code{diff}}{(Only if \code{difference = TRUE}) Risk/survival
#'       difference (arm 1 minus arm 0).}
#'     \item{\code{diff_lo}}{(Only if \code{difference = TRUE}) Lower CI
#'       bound for the difference.}
#'     \item{\code{diff_hi}}{(Only if \code{difference = TRUE}) Upper CI
#'       bound for the difference.}
#'   }
#'   The \code{prediction_meta} slot is also populated with the prediction
#'   type, number of MC samples, confidence level, and the raw MC simulation
#'   matrices (\code{mc_0} and \code{mc_1}) for further analysis if needed.
#'   The object is returned invisibly.
#'
#' @seealso \code{\link{emulate_fit}} for the previous step,
#'   \code{\link{emulate_plot}} for visualization,
#'   \code{\link{emulate_report}} for formatted results.
#'
#' @references
#' Hernan MA, Robins JM (2020). \emph{Causal Inference: What If}.
#' Chapman & Hall/CRC. Chapter 13 (standardization / g-formula).
#'
#' @examples
#' \donttest{
#' # Full pipeline: prepare -> expand -> weight -> fit -> predict
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
#' # Predict cumulative incidence at times 1-5
#' obj <- emulate_predict(obj, times = 1:5, type = "cum_inc",
#'                    samples = 50, seed = 123)
#' print(obj$predictions)
#'
#' # With risk differences
#' obj <- emulate_predict(obj, times = 1:5, type = "cum_inc",
#'                    samples = 50, seed = 123, difference = TRUE)
#' }
#'
#' @export
emulate_predict <- function(obj, times, type = "cum_inc", samples = 100L,
                        seed = NULL, level = 95, difference = FALSE) {
  .check_fitted(obj)

  if (obj$model$type != "logistic") {
    stop("emulate_predict currently only supports logistic model", call. = FALSE)
  }

  if (!type %in% c("cum_inc", "survival")) {
    stop("type must be 'cum_inc' or 'survival'", call. = FALSE)
  }
  if (samples < 10) stop("samples must be at least 10", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  s <- obj$settings
  dt <- obj$data
  prefix <- s$prefix
  id_col <- s$id
  m <- obj$model

  fu_col <- paste0(prefix, "followup")
  trial_col <- paste0(prefix, "trial")
  esample_col <- paste0(prefix, "esample")

  message("emulate_predict - Marginal Predictions")
  message(strrep("-", 50))
  message("Type:          ", type)
  message("MC samples:    ", samples)
  message("Level:         ", level, "%")

  # Get coefficient vector and variance matrix
  b_hat <- coef(m$object)
  V_hat <- m$vcov
  n_coefs <- length(b_hat)

  # Reference population: unique individuals at followup=0, in estimation sample
  ref <- data.table::copy(dt[get(fu_col) == 0])
  if (esample_col %in% names(ref)) {
    ref <- ref[get(esample_col) == 1]
  }
  # One row per id-trial
  ref <- ref[, .SD[1], by = c(id_col, trial_col)]

  times <- sort(as.integer(times))
  n_times <- length(times)
  last_time <- max(times)

  # CI percentiles
  alpha <- (100 - level) / 2
  lo_pct <- alpha
  hi_pct <- 100 - alpha

  # Initialize results
  n_cols <- if (difference) 10L else 7L
  results <- matrix(NA_real_, nrow = n_times, ncol = n_cols)
  results[, 1] <- times

  # Point estimates
  message("Computing predictions...")
  for (arm in 0:1) {
    cum_surv <- rep(1, nrow(ref))

    for (ss in 0:last_time) {
      prob <- .predict_prob(ref, b_hat, ss, arm,
                            m$model_var, prefix,
                            m$followup_spec, m$trial_period_spec,
                            m$outcome_cov, m$fu_ns_info, m$tr_ns_info)
      cum_surv <- cum_surv * (1 - prob)

      if (ss %in% times) {
        mean_surv <- mean(cum_surv)
        est <- if (type == "cum_inc") 1 - mean_surv else mean_surv
        tidx <- which(times == ss)
        col_offset <- 2L + arm * 3L  # arm 0: col 2, arm 1: col 5
        results[tidx, col_offset] <- est
      }
    }
  }

  # MC confidence intervals
  message("Running ", samples, " Monte Carlo simulations...")

  mc_0 <- matrix(NA_real_, nrow = samples, ncol = n_times)
  mc_1 <- matrix(NA_real_, nrow = samples, ncol = n_times)

  # Cholesky decomposition
  # Handle dropped coefficients (zero variance)
  keep_idx <- which(diag(V_hat) > 0)
  n_keep <- length(keep_idx)
  chol_ok <- TRUE

  if (n_keep < n_coefs && n_keep > 0) {
    V_use <- V_hat[keep_idx, keep_idx, drop = FALSE]
    b_use <- b_hat[keep_idx]
    L_chol <- tryCatch(chol(V_use), error = function(e) NULL)
    if (is.null(L_chol)) chol_ok <- FALSE
  } else if (n_keep == n_coefs) {
    V_use <- V_hat
    b_use <- b_hat
    L_chol <- tryCatch(chol(V_use), error = function(e) NULL)
    if (is.null(L_chol)) chol_ok <- FALSE
  } else {
    chol_ok <- FALSE
  }

  if (!chol_ok) {
    warning("Variance matrix is not positive definite. ",
            "Using diagonal approximation for MC CIs.", call. = FALSE)
  }

  for (sim in seq_len(samples)) {
    # Draw from MVN(b_hat, V_hat)
    if (chol_ok) {
      z <- rnorm(n_keep)
      b_reduced <- b_use + as.numeric(z %*% L_chol)
      # Reconstruct full coefficient vector
      b_draw <- b_hat
      b_draw[keep_idx] <- b_reduced
    } else {
      b_draw <- b_hat
      for (j in seq_len(n_coefs)) {
        se_j <- sqrt(V_hat[j, j])
        if (se_j > 0) b_draw[j] <- b_hat[j] + rnorm(1) * se_j
      }
    }

    for (arm in 0:1) {
      cum_surv <- rep(1, nrow(ref))

      for (ss in 0:last_time) {
        prob <- .predict_prob(ref, b_draw, ss, arm,
                              m$model_var, prefix,
                              m$followup_spec, m$trial_period_spec,
                              m$outcome_cov, m$fu_ns_info, m$tr_ns_info)
        cum_surv <- cum_surv * (1 - prob)

        if (ss %in% times) {
          mean_surv <- mean(cum_surv)
          pred <- if (type == "cum_inc") 1 - mean_surv else mean_surv
          tidx <- which(times == ss)
          if (arm == 0) mc_0[sim, tidx] <- pred
          else mc_1[sim, tidx] <- pred
        }
      }
    }

    if (sim %% 50 == 0) {
      message("  ... ", sim, " of ", samples, " samples completed")
    }
  }

  # Compute CIs from MC samples
  for (tidx in seq_len(n_times)) {
    # Arm 0 CIs
    results[tidx, 3] <- .emulate_pctile(mc_0[, tidx], lo_pct)
    results[tidx, 4] <- .emulate_pctile(mc_0[, tidx], hi_pct)

    # Arm 1 CIs
    results[tidx, 6] <- .emulate_pctile(mc_1[, tidx], lo_pct)
    results[tidx, 7] <- .emulate_pctile(mc_1[, tidx], hi_pct)

    # Risk difference
    if (difference) {
      results[tidx, 8] <- results[tidx, 5] - results[tidx, 2]
      diffs <- mc_1[, tidx] - mc_0[, tidx]
      results[tidx, 9] <- .emulate_pctile(diffs, lo_pct)
      results[tidx, 10] <- .emulate_pctile(diffs, hi_pct)
    }
  }

  # Column names
  if (difference) {
    colnames(results) <- c("time", "est_0", "ci_lo_0", "ci_hi_0",
                           "est_1", "ci_lo_1", "ci_hi_1",
                           "diff", "diff_lo", "diff_hi")
  } else {
    colnames(results) <- c("time", "est_0", "ci_lo_0", "ci_hi_0",
                           "est_1", "ci_lo_1", "ci_hi_1")
  }

  # Display
  message("")
  type_label <- if (type == "cum_inc") "Cumulative Incidence" else "Survival"
  message(type_label, " Estimates (", level, "% CI)")
  message(strrep("-", 50))

  for (tidx in seq_len(n_times)) {
    t <- results[tidx, 1]
    e0 <- results[tidx, 2]
    e1 <- results[tidx, 5]
    msg <- sprintf("  t=%d: Control=%.4f  Treated=%.4f", t, e0, e1)
    if (difference) {
      rd <- results[tidx, 8]
      msg <- paste0(msg, sprintf("  Diff=%.4f", rd))
    }
    message(msg)
  }
  message(strrep("-", 50))

  obj$predictions <- as.data.frame(results)
  obj$prediction_meta <- list(type = type, samples = samples,
                               level = level, mc_0 = mc_0, mc_1 = mc_1)

  invisible(obj)
}
