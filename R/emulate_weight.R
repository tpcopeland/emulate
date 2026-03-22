#' Compute inverse probability weights for target trial emulation
#'
#' Calculates stabilized inverse probability of treatment weights (IPTW) and
#' optionally inverse probability of censoring weights (IPCW) for the
#' expanded sequential trials dataset. These weights correct for the
#' selection bias introduced by artificial censoring in the
#' clone-censor-weight framework.
#'
#' @section Why weights are needed:
#' When the PP or AT estimand is chosen, \code{\link{emulate_expand}} artificially
#' censors individuals who deviate from their assigned treatment strategy.
#' This censoring is \emph{informative} -- individuals who deviate are likely
#' different from those who adhere. Inverse probability weights re-weight the
#' remaining (uncensored) observations so that they represent the full
#' population, as if no one had deviated. This is the "weight" part of
#' clone-censor-weight.
#'
#' @section IPTW (treatment switching weights):
#' The treatment switching weights model the probability of remaining on the
#' assigned treatment at each time point, conditional on covariates. Two
#' logistic regression models are fitted:
#' \describe{
#'   \item{\strong{Denominator model}}{Models \code{P(treatment_t | lagged
#'     treatment, covariates, follow-up time)}. This is the "full" model that
#'     conditions on all predictors of treatment switching. Covariates are
#'     specified via \code{switch_d_cov}.}
#'   \item{\strong{Numerator model}}{Models \code{P(treatment_t | lagged
#'     treatment [, baseline covariates])}. This is a simpler model used for
#'     stabilization. Covariates (if any) are specified via
#'     \code{switch_n_cov}. If \code{switch_n_cov = NULL}, only lagged
#'     treatment is used.}
#' }
#' The stabilized weight at each time point is
#' \code{numerator_probability / denominator_probability} for the observed
#' treatment value. The cumulative product over time gives the final weight.
#'
#' @section IPCW (censoring weights):
#' If your data has an external censoring indicator (e.g., loss to follow-up)
#' in addition to the artificial censoring from the clone-censor-weight step,
#' you can model this censoring similarly using \code{censor_d_cov} and
#' \code{censor_n_cov}. The censoring weights are multiplied with the
#' treatment switching weights to produce a combined weight.
#'
#' @section Stabilization:
#' Stabilized weights (the default, \code{stabilized = TRUE}) divide the
#' numerator model probability by the denominator model probability at each
#' time point, rather than using \code{1 / denominator}. Stabilized weights
#' have mean closer to 1, smaller variance, and better finite-sample
#' properties.
#'
#' @section Weight truncation:
#' Extreme weights can destabilize estimates. The \code{truncate} argument
#' allows you to cap weights at specified percentiles. For example,
#' \code{truncate = c(1, 99)} replaces any weight below the 1st percentile
#' with the 1st percentile value and any weight above the 99th percentile
#' with the 99th percentile value.
#'
#' @section Pooled vs. stratified models:
#' By default, separate weight models are fitted for each arm
#' (\code{pool_switch = FALSE}). Setting \code{pool_switch = TRUE} fits a
#' single model across both arms (with arm as a covariate), which can be more
#' efficient when sample sizes are small. The same applies to censoring
#' models via \code{pool_censor}.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_expand}}. The
#'   \code{expanded} state flag must be \code{TRUE}.
#' @param switch_d_cov A character vector of covariate names for the treatment
#'   switch \strong{denominator} model. These should be time-varying and
#'   baseline covariates that predict treatment switching. Required for PP
#'   and AT estimands. Follow-up time is automatically included.
#' @param switch_n_cov An optional character vector of covariate names for
#'   the treatment switch \strong{numerator} (stabilization) model. Typically
#'   a subset of \code{switch_d_cov}, often just baseline covariates. If
#'   \code{NULL} (the default), the numerator model uses only lagged
#'   treatment as a predictor, which is the simplest stabilization.
#' @param censor_d_cov An optional character vector of covariate names for
#'   the censoring \strong{denominator} model. If \code{NULL} (the default),
#'   no censoring weights are computed. Specify this when your data has an
#'   external censoring variable (set in \code{\link{emulate_prepare}}).
#' @param censor_n_cov An optional character vector of covariate names for
#'   the censoring \strong{numerator} (stabilization) model. If \code{NULL},
#'   the numerator model uses only an intercept.
#' @param pool_switch A logical value. If \code{TRUE}, a single pooled
#'   treatment switch model is fitted across both arms (with arm as a
#'   covariate). If \code{FALSE} (the default), separate models are fitted
#'   for each arm.
#' @param pool_censor A logical value. If \code{TRUE}, a single pooled
#'   censoring model is fitted across both arms. If \code{FALSE} (the
#'   default), separate models are fitted for each arm.
#' @param truncate An optional numeric vector of length 2 giving the lower
#'   and upper truncation percentiles for weight trimming. For example,
#'   \code{c(1, 99)} truncates at the 1st and 99th percentiles. The lower
#'   value must be less than the upper value. If \code{NULL} (the default),
#'   no truncation is applied.
#' @param stabilized A logical value indicating whether to use stabilized
#'   weights. Default is \code{TRUE}. Set to \code{FALSE} to use unstabilized
#'   weights (\code{1 / denominator}), though this is rarely recommended.
#' @param weight_name An optional character string specifying the column name
#'   for the weight variable. If \code{NULL} (the default), the name is
#'   constructed as \code{paste0(prefix, "weight")} (e.g.,
#'   \code{"_emulate_weight"}).
#' @param ps_method A character string specifying the propensity score
#'   estimation method. One of:
#'   \describe{
#'     \item{\code{"glm"}}{(Default) Standard logistic regression, the same
#'       method used in previous versions.}
#'     \item{\code{"lasso"}}{Cross-validated L1 (lasso) regularized logistic
#'       regression via \code{glmnet::cv.glmnet}. Useful when there are many
#'       covariates or potential collinearity.}
#'   }
#' @param ps_trim An optional numeric vector of length 2 specifying propensity
#'   score trimming quantiles (e.g., \code{c(0.05, 0.95)}). Observations with
#'   propensity scores outside these quantiles are \strong{removed} before
#'   weight computation. This is distinct from \code{truncate}, which caps
#'   weight \emph{values} after computation. Values must be between 0 and 1.
#'   Default is \code{NULL} (no trimming).
#' @param quiet A logical value. If \code{TRUE}, suppresses progress
#'   messages during model fitting. Default is \code{FALSE}.
#'
#' @return The input \code{emulate} object with the weight variable added to the
#'   data, \code{state$weighted} set to \code{TRUE}, and the \code{weights}
#'   slot populated with:
#'   \describe{
#'     \item{\code{weight_var}}{Character: name of the weight variable.}
#'     \item{\code{truncation}}{Numeric vector of truncation bounds, or
#'       \code{NULL}.}
#'     \item{\code{mean}, \code{sd}, \code{min}, \code{max}}{Summary
#'       statistics of the final weight distribution.}
#'     \item{\code{ess}}{Effective sample size, computed as
#'       \code{(sum(w))^2 / sum(w^2)}. Values much smaller than the actual
#'       sample size indicate highly variable weights.}
#'     \item{\code{n_truncated}}{Number of observations whose weights were
#'       truncated.}
#'   }
#'   The object is returned invisibly. For the ITT estimand, all weights are
#'   set to 1.0 and no models are fitted.
#'
#' @seealso \code{\link{emulate_expand}} for the previous step,
#'   \code{\link{emulate_fit}} for the next step,
#'   \code{\link{emulate_diagnose}} for weight diagnostics.
#'
#' @references
#' Hernan MA, Robins JM (2020). \emph{Causal Inference: What If}.
#' Chapman & Hall/CRC. Chapters 12 and 16.
#'
#' Danaei G, Rodriguez LAG, Cantero OF, Logan R, Hernan MA (2013).
#' "Observational data for comparative effectiveness research."
#' \emph{Statistical Methods in Medical Research}, 22(1), 70-96.
#'
#' @examples
#' \donttest{
#' # Create and prepare test data
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
#'
#' # Compute stabilized weights with truncation
#' obj <- emulate_weight(obj,
#'   switch_d_cov = c("age", "sex"),
#'   switch_n_cov = c("age"),
#'   truncate = c(1, 99)
#' )
#'
#' # ITT: weights are all 1 (no modeling needed)
#' obj_itt <- emulate_prepare(dat, id = "id", period = "period",
#'                        treatment = "treatment", outcome = "outcome",
#'                        eligible = "eligible", estimand = "ITT")
#' obj_itt <- emulate_expand(obj_itt, maxfollowup = 5)
#' obj_itt <- emulate_weight(obj_itt)  # All weights = 1
#' }
#'
#' @export
emulate_weight <- function(obj, switch_d_cov = NULL, switch_n_cov = NULL,
                       censor_d_cov = NULL, censor_n_cov = NULL,
                       pool_switch = FALSE, pool_censor = FALSE,
                       truncate = NULL, stabilized = TRUE,
                       weight_name = NULL, ps_method = "glm",
                       ps_trim = NULL, quiet = FALSE) {
  .check_expanded(obj)

  s <- obj$settings
  dt <- obj$data
  estimand <- s$estimand
  prefix <- s$prefix
  id_col <- s$id
  treat_col <- s$treatment

  fu_col <- paste0(prefix, "followup")
  trial_col <- paste0(prefix, "trial")
  arm_col <- paste0(prefix, "arm")
  cens_col <- paste0(prefix, "censored")

  if (is.null(weight_name)) weight_name <- paste0(prefix, "weight")

  message("emulate_weight - Inverse Probability Weights")
  message(strrep("-", 50))
  message("Estimand:      ", estimand)

  # ITT without censoring covariates: all weights = 1
  if (estimand == "ITT" && (is.null(censor_d_cov) || length(censor_d_cov) == 0)) {
    message("ITT estimand - all weights set to 1")
    dt[, (weight_name) := 1.0]
    obj$data <- dt
    obj$state$weighted <- TRUE
    obj$weights <- list(weight_var = weight_name, truncation = NULL,
                        mean = 1, sd = 0, min = 1, max = 1,
                        ess = nrow(dt))
    return(invisible(obj))
  }

  # ITT with censoring covariates: compute IPCW only (no switch weights)
  if (estimand == "ITT" && !is.null(censor_d_cov) && length(censor_d_cov) > 0) {
    message("ITT estimand - computing IPCW for censoring only (no switch weights)")
  }

  # Validate ps_method
  if (!ps_method %in% c("glm", "lasso")) {
    stop("ps_method must be 'glm' or 'lasso'", call. = FALSE)
  }
  if (ps_method == "lasso") {
    message("PS method:     lasso (cross-validated L1)")
  }

  # PP/AT require switch covariates
  if ((is.null(switch_d_cov) || length(switch_d_cov) == 0) && estimand != "ITT") {
    stop("switch_d_cov required for ", estimand, " estimand", call. = FALSE)
  }

  # Map covariates to time-varying versions for weight models.
  # emulate_expand saves original (unfrozen) values as {prefix}tv_{varname}.
  # Weight models should condition on time-varying L_t, not frozen L_0.
  .map_to_tv <- function(cov_names) {
    vapply(cov_names, function(v) {
      tv_name <- paste0(prefix, "tv_", v)
      if (tv_name %in% names(dt)) tv_name else v
    }, character(1), USE.NAMES = FALSE)
  }

  switch_d_cov_w <- .map_to_tv(switch_d_cov)
  switch_n_cov_w <- if (!is.null(switch_n_cov)) .map_to_tv(switch_n_cov) else NULL

  n_mapped <- sum(switch_d_cov_w != switch_d_cov)
  if (n_mapped > 0 && !quiet) {
    message("  Using time-varying values for ", n_mapped, " weight covariate(s)")
  }

  # PS trimming (before weight computation)
  n_ps_trimmed <- 0L
  if (!is.null(ps_trim)) {
    message("Computing propensity scores for trimming...")

    # Fit PS model on uncensored observations (use time-varying covariates)
    bq <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))
    ps_formula <- as.formula(paste(bq(treat_col), "~",
                                    paste(bq(switch_d_cov_w), collapse = " + ")))

    if (ps_method == "lasso") {
      ps_result <- .fit_ps_lasso(ps_formula, data = as.data.frame(dt))
      dt[, .ps_trim_val := ps_result$ps]
    } else {
      ps_fit <- tryCatch(
        suppressWarnings(glm(ps_formula, data = dt, family = binomial())),
        error = function(e) NULL
      )
      if (!is.null(ps_fit)) {
        dt[, .ps_trim_val := predict(ps_fit, type = "response")]
      } else {
        dt[, .ps_trim_val := 0.5]
      }
    }

    extreme <- !.trim_ps(dt$.ps_trim_val, ps_trim)
    n_extreme_obs <- sum(extreme, na.rm = TRUE)

    if (n_extreme_obs > 0) {
      # Propagate to full id-trial-arm trajectories
      dt[, .traj_has_extreme := as.integer(any(extreme[.I])),
         by = c(id_col, trial_col, arm_col)]
      # Workaround: compute per-group max of extreme flag
      dt[, .extreme_flag := as.integer(extreme)]
      dt[, .traj_has_extreme := max(.extreme_flag),
         by = c(id_col, trial_col, arm_col)]

      n_traj_dropped <- uniqueN(dt[.traj_has_extreme == 1L,
                                    .SD, .SDcols = c(id_col, trial_col, arm_col)])
      n_before <- nrow(dt)
      dt <- dt[.traj_has_extreme == 0L]
      n_ps_trimmed <- n_before - nrow(dt)
      dt[, c(".traj_has_extreme", ".extreme_flag") := NULL]

      message("  PS trimming: removed ", n_traj_dropped,
              " id-trial-arm trajectories (", n_ps_trimmed, " observations)")
      message("  Note: entire trajectories dropped to preserve longitudinal structure.")
    } else {
      n_ps_trimmed <- 0L
      message("  PS trimming: no observations removed")
    }

    # Store PS values for diagnostics
    obj$weights$ps_values <- dt$.ps_trim_val
    dt[, .ps_trim_val := NULL]
  }

  # Initialize weight to 1
  dt[, (weight_name) := 1.0]

  # Treatment switch weights (skip for ITT+IPCW)
  if (!is.null(switch_d_cov) && length(switch_d_cov) > 0) {
  message("Fitting treatment switch models...")
  data.table::setorderv(dt, c(id_col, trial_col, arm_col, fu_col))

  if (pool_switch) {
    .weight_switch_pooled(dt, id_col, treat_col, arm_col,
                          fu_col, trial_col, cens_col,
                          switch_d_cov_w, switch_n_cov_w, weight_name, quiet,
                          ps_method)
  } else {
    for (a in 0:1) {
      if (!quiet) message("  Switch model for arm ", a, "...")
      .weight_switch_arm(dt, id_col, treat_col, arm_col, a,
                         fu_col, trial_col, cens_col,
                         switch_d_cov_w, switch_n_cov_w, weight_name, quiet,
                         ps_method)
    }
  }

  } # end switch_d_cov guard

  # Censoring weights (optional)
  if (!is.null(censor_d_cov) && length(censor_d_cov) > 0) {
    message("Fitting censoring models...")

    # P5: Combined natural + artificial censoring for non-ITT
    nat_censor <- s$censor
    art_censor <- cens_col
    if (!is.null(nat_censor) && nat_censor %in% names(dt)) {
      if (estimand == "ITT") {
        censor_indicator <- nat_censor
      } else {
        # Combine: censored if either natural or artificial censoring occurred
        dt[, .combined_censor := pmax(get(nat_censor), get(art_censor), na.rm = TRUE)]
        censor_indicator <- ".combined_censor"
        message("  Using combined natural + artificial censoring indicator")
      }
    } else {
      censor_indicator <- art_censor
    }

    # Map censor covariates to time-varying versions
    censor_d_cov_w <- .map_to_tv(censor_d_cov)
    censor_n_cov_w <- if (!is.null(censor_n_cov)) .map_to_tv(censor_n_cov) else NULL

    n_cens_mapped <- sum(censor_d_cov_w != censor_d_cov)
    if (n_cens_mapped > 0 && !quiet) {
      message("  Using time-varying values for ", n_cens_mapped, " censor covariate(s)")
    }

    if (pool_censor) {
      .weight_censor_pooled(dt, id_col, censor_indicator, arm_col,
                             fu_col, trial_col,
                             censor_d_cov_w, censor_n_cov_w, weight_name, quiet)
    } else {
      for (a in 0:1) {
        if (!quiet) message("  Censor model for arm ", a, "...")
        .weight_censor_arm(dt, id_col, censor_indicator, arm_col, a,
                            fu_col, trial_col,
                            censor_d_cov_w, censor_n_cov_w, weight_name, quiet)
      }
    }
  }

  # Clean up combined censor indicator if created
  if (".combined_censor" %in% names(dt)) {
    dt[, .combined_censor := NULL]
  }

  # Truncation
  n_truncated <- 0L
  if (!is.null(truncate)) {
    if (length(truncate) != 2 || truncate[1] >= truncate[2]) {
      stop("truncate must be a length-2 vector with lower < upper", call. = FALSE)
    }
    message("Truncating weights at ", truncate[1], "th and ", truncate[2],
            "th percentiles...")

    w <- dt[[weight_name]]
    w_nona <- w[!is.na(w)]
    lo_val <- quantile(w_nona, truncate[1] / 100, type = 2)
    hi_val <- quantile(w_nona, truncate[2] / 100, type = 2)

    n_lo <- sum(w < lo_val, na.rm = TRUE)
    n_hi <- sum(w > hi_val, na.rm = TRUE)
    n_truncated <- n_lo + n_hi

    dt[get(weight_name) < lo_val & !is.na(get(weight_name)),
       (weight_name) := lo_val]
    dt[get(weight_name) > hi_val & !is.na(get(weight_name)),
       (weight_name) := hi_val]

    message("  Truncated ", n_truncated, " observations (",
            n_lo, " low, ", n_hi, " high)")
  }

  # Diagnostics
  w <- dt[[weight_name]]
  w_stats <- summary(w[!is.na(w)])
  w_mean <- mean(w, na.rm = TRUE)
  w_sd <- sd(w, na.rm = TRUE)
  w_min <- min(w, na.rm = TRUE)
  w_max <- max(w, na.rm = TRUE)
  sum_w <- sum(w, na.rm = TRUE)
  sum_w2 <- sum(w^2, na.rm = TRUE)
  ess <- sum_w^2 / sum_w2

  message("\nWeight distribution:")
  message("  Mean:   ", sprintf("%9.4f", w_mean))
  message("  SD:     ", sprintf("%9.4f", w_sd))
  message("  Min:    ", sprintf("%9.4f", w_min))
  message("  Max:    ", sprintf("%9.4f", w_max))
  message("  ESS:    ", sprintf("%9.1f", ess))
  message(strrep("-", 50))

  obj$data <- dt
  obj$state$weighted <- TRUE
  # Preserve ps_values from trimming step if present
  ps_vals <- obj$weights$ps_values
  obj$weights <- list(
    weight_var = weight_name,
    truncation = truncate,
    ps_method = ps_method,
    ps_trim = ps_trim,
    mean = w_mean, sd = w_sd,
    min = w_min, max = w_max,
    ess = ess,
    n_truncated = n_truncated,
    n_ps_trimmed = n_ps_trimmed,
    ps_values = ps_vals
  )

  invisible(obj)
}
