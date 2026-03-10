# Weight model helpers for emulate_weight
# Matches Stata's _emulate_weight_switch_arm, _emulate_weight_switch_pooled,
# _emulate_weight_censor_arm, _emulate_weight_censor_pooled

# Helper to backtick-quote non-syntactic names
.bq <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))

#' Fit switch weight model for one arm (stratified)
#' @keywords internal
.weight_switch_arm <- function(dt, id_col, treat_col, arm_col, arm_val,
                               fu_col, trial_col, cens_col,
                               d_cov, n_cov, weight_col, quiet = FALSE,
                               ps_method = "glm") {
  in_arm <- dt[[arm_col]] == arm_val

  # Ensure ordering for correct lag computation
  data.table::setorderv(dt, c(id_col, trial_col, arm_col, fu_col))

  # Lagged treatment within person-trial-arm
  dt[, lag_treat := shift(get(treat_col), 1L, type = "lag"),
     by = c(id_col, trial_col, arm_col)]

  has_lag <- in_arm & !is.na(dt$lag_treat)

  # Denominator model: P(A_t | A_{t-1}, L, followup)
  d_formula <- paste(.bq(treat_col), "~ lag_treat +",
                     paste(.bq(c(d_cov, fu_col)), collapse = " + "))

  if (ps_method == "lasso") {
    lasso_d <- tryCatch(
      .fit_ps_lasso(as.formula(d_formula), data = as.data.frame(dt[has_lag])),
      error = function(e) NULL
    )
    if (is.null(lasso_d) || is.null(lasso_d$fit)) {
      dt[has_lag, denom_pr := 0.5]
    } else {
      dt[has_lag, denom_pr := lasso_d$ps]
    }
  } else {
    denom_fit <- tryCatch(
      suppressWarnings(glm(as.formula(d_formula), data = dt[has_lag],
                           family = binomial())),
      error = function(e) NULL
    )
    if (is.null(denom_fit)) {
      dt[has_lag, denom_pr := 0.5]
    } else {
      dt[has_lag, denom_pr := predict(denom_fit, newdata = dt[has_lag],
                                      type = "response")]
    }
  }

  # Numerator model: P(A_t | A_{t-1} [, baseline covs])
  if (!is.null(n_cov) && length(n_cov) > 0) {
    n_formula <- paste(.bq(treat_col), "~ lag_treat +",
                       paste(.bq(n_cov), collapse = " + "))
  } else {
    n_formula <- paste(.bq(treat_col), "~ lag_treat")
  }

  if (ps_method == "lasso") {
    lasso_n <- tryCatch(
      .fit_ps_lasso(as.formula(n_formula), data = as.data.frame(dt[has_lag])),
      error = function(e) NULL
    )
    if (is.null(lasso_n) || is.null(lasso_n$fit)) {
      dt[has_lag, numer_pr := 0.5]
    } else {
      dt[has_lag, numer_pr := lasso_n$ps]
    }
  } else {
    numer_fit <- tryCatch(
      suppressWarnings(glm(as.formula(n_formula), data = dt[has_lag],
                           family = binomial())),
      error = function(e) NULL
    )
    if (is.null(numer_fit)) {
      dt[has_lag, numer_pr := 0.5]
    } else {
      dt[has_lag, numer_pr := predict(numer_fit, newdata = dt[has_lag],
                                      type = "response")]
    }
  }

  # Stabilized weight contribution at each time
  dt[in_arm, sw_t := 1.0]

  dt[in_arm & get(treat_col) == 1 & !is.na(denom_pr) & denom_pr > 0.001,
     sw_t := numer_pr / denom_pr]
  dt[in_arm & get(treat_col) == 0 & !is.na(denom_pr) & denom_pr < 0.999,
     sw_t := (1 - numer_pr) / (1 - denom_pr)]

  # Cumulative product via log-sum
  dt[in_arm & !is.na(sw_t) & sw_t > 0, log_sw := log(sw_t)]
  dt[in_arm & is.na(log_sw) & is.na(lag_treat), log_sw := 0]

  dt[in_arm, cum_log_sw := cumsum(fifelse(is.na(log_sw), 0, log_sw)),
     by = c(id_col, trial_col, arm_col)]

  # Update weight
  dt[in_arm & !is.na(cum_log_sw),
     (weight_col) := get(weight_col) * exp(cum_log_sw)]

  # Clean up
  dt[, c("lag_treat", "denom_pr", "numer_pr", "sw_t",
         "log_sw", "cum_log_sw") := NULL]

  invisible(dt)
}

#' Fit pooled switch weight model (across arms)
#' @keywords internal
.weight_switch_pooled <- function(dt, id_col, treat_col, arm_col,
                                  fu_col, trial_col, cens_col,
                                  d_cov, n_cov, weight_col, quiet = FALSE,
                                  ps_method = "glm") {
  data.table::setorderv(dt, c(id_col, trial_col, arm_col, fu_col))

  # Lagged treatment
  dt[, lag_treat := shift(get(treat_col), 1L, type = "lag"),
     by = c(id_col, trial_col, arm_col)]

  has_lag <- !is.na(dt$lag_treat)

  # Denominator: P(A_t | A_{t-1}, arm, L, followup)
  d_formula <- paste(.bq(treat_col), "~ lag_treat +", .bq(arm_col), "+",
                     paste(.bq(c(d_cov, fu_col)), collapse = " + "))

  if (ps_method == "lasso") {
    lasso_d <- tryCatch(
      .fit_ps_lasso(as.formula(d_formula), data = as.data.frame(dt[has_lag])),
      error = function(e) NULL
    )
    if (is.null(lasso_d) || is.null(lasso_d$fit)) {
      dt[has_lag, denom_pr := 0.5]
    } else {
      dt[has_lag, denom_pr := lasso_d$ps]
    }
  } else {
    denom_fit <- tryCatch(
      suppressWarnings(glm(as.formula(d_formula), data = dt[has_lag],
                           family = binomial())),
      error = function(e) NULL
    )
    if (is.null(denom_fit)) {
      dt[has_lag, denom_pr := 0.5]
    } else {
      dt[has_lag, denom_pr := predict(denom_fit, newdata = dt[has_lag],
                                      type = "response")]
    }
  }

  # Numerator: P(A_t | A_{t-1}, arm [, baseline])
  if (!is.null(n_cov) && length(n_cov) > 0) {
    n_formula <- paste(.bq(treat_col), "~ lag_treat +", .bq(arm_col), "+",
                       paste(.bq(n_cov), collapse = " + "))
  } else {
    n_formula <- paste(.bq(treat_col), "~ lag_treat +", .bq(arm_col))
  }

  if (ps_method == "lasso") {
    lasso_n <- tryCatch(
      .fit_ps_lasso(as.formula(n_formula), data = as.data.frame(dt[has_lag])),
      error = function(e) NULL
    )
    if (is.null(lasso_n) || is.null(lasso_n$fit)) {
      dt[has_lag, numer_pr := 0.5]
    } else {
      dt[has_lag, numer_pr := lasso_n$ps]
    }
  } else {
    numer_fit <- tryCatch(
      suppressWarnings(glm(as.formula(n_formula), data = dt[has_lag],
                           family = binomial())),
      error = function(e) NULL
    )
    if (is.null(numer_fit)) {
      dt[has_lag, numer_pr := 0.5]
    } else {
      dt[has_lag, numer_pr := predict(numer_fit, newdata = dt[has_lag],
                                      type = "response")]
    }
  }

  dt[, sw_t := 1.0]
  dt[get(treat_col) == 1 & !is.na(denom_pr) & denom_pr > 0.001,
     sw_t := numer_pr / denom_pr]
  dt[get(treat_col) == 0 & !is.na(denom_pr) & denom_pr < 0.999,
     sw_t := (1 - numer_pr) / (1 - denom_pr)]

  dt[!is.na(sw_t) & sw_t > 0, log_sw := log(sw_t)]
  dt[is.na(log_sw) & is.na(lag_treat), log_sw := 0]

  dt[, cum_log_sw := cumsum(fifelse(is.na(log_sw), 0, log_sw)),
     by = c(id_col, trial_col, arm_col)]

  dt[!is.na(cum_log_sw),
     (weight_col) := get(weight_col) * exp(cum_log_sw)]

  dt[, c("lag_treat", "denom_pr", "numer_pr", "sw_t",
         "log_sw", "cum_log_sw") := NULL]

  invisible(dt)
}

#' Fit censor weight model for one arm (stratified)
#' @keywords internal
.weight_censor_arm <- function(dt, id_col, censor_col, arm_col, arm_val,
                                fu_col, trial_col,
                                d_cov, n_cov, weight_col, quiet = FALSE) {
  in_arm <- dt[[arm_col]] == arm_val

  # Denominator: P(C_t=1 | L)
  d_formula <- paste(.bq(censor_col), "~",
                     paste(.bq(c(d_cov, fu_col)), collapse = " + "))
  denom_fit <- tryCatch(
    suppressWarnings(glm(as.formula(d_formula), data = dt[in_arm],
                         family = binomial())),
    error = function(e) NULL
  )
  if (is.null(denom_fit)) {
    dt[in_arm, denom_pr := 0.05]
  } else {
    dt[in_arm, denom_pr := predict(denom_fit, newdata = dt[in_arm],
                                    type = "response")]
  }

  # Numerator
  if (!is.null(n_cov) && length(n_cov) > 0) {
    n_formula <- paste(.bq(censor_col), "~", paste(.bq(n_cov), collapse = " + "))
  } else {
    n_formula <- paste(.bq(censor_col), "~ 1")
  }
  numer_fit <- tryCatch(
    suppressWarnings(glm(as.formula(n_formula), data = dt[in_arm],
                         family = binomial())),
    error = function(e) NULL
  )
  if (is.null(numer_fit)) {
    dt[in_arm, numer_pr := 0.05]
  } else {
    dt[in_arm, numer_pr := predict(numer_fit, newdata = dt[in_arm],
                                    type = "response")]
  }

  # Weight: (1-numer)/(1-denom) for uncensored
  dt[in_arm & denom_pr < 0.999,
     cw_t := (1 - numer_pr) / (1 - denom_pr)]
  dt[in_arm & denom_pr >= 0.999 & !is.na(denom_pr), cw_t := 1.0]

  dt[in_arm & !is.na(cw_t), log_cw := log(cw_t)]
  dt[in_arm & is.na(log_cw) & is.na(denom_pr), log_cw := 0]

  dt[in_arm, cum_log_cw := cumsum(fifelse(is.na(log_cw), 0, log_cw)),
     by = c(id_col, trial_col, arm_col)]

  dt[in_arm & !is.na(cum_log_cw),
     (weight_col) := get(weight_col) * exp(cum_log_cw)]

  dt[, c("denom_pr", "numer_pr", "cw_t", "log_cw", "cum_log_cw") := NULL]

  invisible(dt)
}

#' Fit pooled censor weight model
#' @keywords internal
.weight_censor_pooled <- function(dt, id_col, censor_col, arm_col,
                                   fu_col, trial_col,
                                   d_cov, n_cov, weight_col, quiet = FALSE) {
  d_formula <- paste(.bq(censor_col), "~", .bq(arm_col), "+",
                     paste(.bq(c(d_cov, fu_col)), collapse = " + "))
  denom_fit <- tryCatch(
    suppressWarnings(glm(as.formula(d_formula), data = dt,
                         family = binomial())),
    error = function(e) NULL
  )
  if (is.null(denom_fit)) {
    dt[, denom_pr := 0.05]
  } else {
    dt[, denom_pr := predict(denom_fit, newdata = dt, type = "response")]
  }

  if (!is.null(n_cov) && length(n_cov) > 0) {
    n_formula <- paste(.bq(censor_col), "~", .bq(arm_col), "+",
                       paste(.bq(n_cov), collapse = " + "))
  } else {
    n_formula <- paste(.bq(censor_col), "~", .bq(arm_col))
  }
  numer_fit <- tryCatch(
    suppressWarnings(glm(as.formula(n_formula), data = dt,
                         family = binomial())),
    error = function(e) NULL
  )
  if (is.null(numer_fit)) {
    dt[, numer_pr := 0.05]
  } else {
    dt[, numer_pr := predict(numer_fit, newdata = dt, type = "response")]
  }

  dt[denom_pr < 0.999, cw_t := (1 - numer_pr) / (1 - denom_pr)]
  dt[denom_pr >= 0.999 & !is.na(denom_pr), cw_t := 1.0]

  dt[!is.na(cw_t), log_cw := log(cw_t)]
  dt[is.na(log_cw) & is.na(denom_pr), log_cw := 0]

  dt[, cum_log_cw := cumsum(fifelse(is.na(log_cw), 0, log_cw)),
     by = c(id_col, trial_col, arm_col)]

  dt[!is.na(cum_log_cw),
     (weight_col) := get(weight_col) * exp(cum_log_cw)]

  dt[, c("denom_pr", "numer_pr", "cw_t", "log_cw", "cum_log_cw") := NULL]

  invisible(dt)
}
