# Prediction helpers for emulate_predict
# Matches Stata's _emulate_predict_xb and _emulate_pctile

#' Compute linear predictor at given (time, arm) for each individual
#'
#' Manually applies coefficient vector to build xb, matching Stata's
#' _emulate_predict_xb exactly.
#'
#' @param ref_data data.table of reference population (at followup=0)
#' @param b_vec Named numeric vector of coefficients
#' @param time Integer: follow-up time point
#' @param arm Integer: treatment arm (0 or 1)
#' @param model_var Character: name of treatment variable in model
#' @param prefix Character: variable prefix
#' @param followup_spec Character: follow-up time specification
#' @param trial_period_spec Character: trial period specification
#' @param outcome_cov Character vector: outcome covariate names
#' @param fu_ns_info List with knots, df for follow-up NS (or NULL)
#' @param tr_ns_info List with knots, df for trial NS (or NULL)
#' @return Numeric vector of predicted probabilities
#' @keywords internal
.predict_prob <- function(ref_data, b_vec, time, arm, model_var, prefix,
                          followup_spec, trial_period_spec, outcome_cov,
                          fu_ns_info, tr_ns_info) {
  # Strip backticks from coefficient names so lookups match unquoted variable names
  names(b_vec) <- gsub("^`|`$", "", names(b_vec))
  # Replace NA coefficients (dropped by GLM due to collinearity) with 0
  b_vec[is.na(b_vec)] <- 0
  coef_names <- names(b_vec)
  n <- nrow(ref_data)
  xb <- rep(0, n)

  fu_col <- paste0(prefix, "followup")
  trial_col <- paste0(prefix, "trial")

  # Constant
  if ("(Intercept)" %in% coef_names) {
    xb <- xb + b_vec["(Intercept)"]
  }

  # Treatment variable (scalar arm value)
  if (model_var %in% coef_names) {
    xb <- xb + b_vec[model_var] * arm
  }

  # Follow-up time terms
  if (!is.null(fu_ns_info)) {
    # Natural spline basis at prediction time
    basis <- .emulate_rcs_basis(time, fu_ns_info$knots, fu_ns_info$df)
    for (j in seq_len(fu_ns_info$df)) {
      vname <- fu_ns_info$varnames[j]
      if (vname %in% coef_names) {
        xb <- xb + b_vec[vname] * basis[1, j]
      }
    }
  } else {
    # Polynomial terms
    fu_name <- fu_col
    if (fu_name %in% coef_names) {
      xb <- xb + b_vec[fu_name] * time
    }
    sq_name <- paste0(prefix, "followup_sq")
    if (sq_name %in% coef_names) {
      xb <- xb + b_vec[sq_name] * time^2
    }
    cu_name <- paste0(prefix, "followup_cu")
    if (cu_name %in% coef_names) {
      xb <- xb + b_vec[cu_name] * time^3
    }
  }

  # Trial period terms (use each observation's actual trial period)
  if (!is.null(tr_ns_info)) {
    # NS basis: use existing data variables
    for (j in seq_len(tr_ns_info$df)) {
      vname <- tr_ns_info$varnames[j]
      if (vname %in% coef_names && vname %in% names(ref_data)) {
        xb <- xb + b_vec[vname] * ref_data[[vname]]
      }
    }
  } else {
    if (trial_col %in% coef_names) {
      xb <- xb + b_vec[trial_col] * ref_data[[trial_col]]
    }
    sq_name <- paste0(prefix, "trial_sq")
    if (sq_name %in% coef_names && sq_name %in% names(ref_data)) {
      xb <- xb + b_vec[sq_name] * ref_data[[trial_col]]^2
    }
    cu_name <- paste0(prefix, "trial_cu")
    if (cu_name %in% coef_names && cu_name %in% names(ref_data)) {
      xb <- xb + b_vec[cu_name] * ref_data[[trial_col]]^3
    }
  }

  # Outcome covariates (use each observation's actual values)
  if (!is.null(outcome_cov)) {
    for (v in outcome_cov) {
      if (v %in% coef_names && v %in% names(ref_data)) {
        xb <- xb + b_vec[v] * ref_data[[v]]
      }
    }
  }

  # Convert to probability via logit link
  1 / (1 + exp(-xb))
}

#' Compute percentile matching Stata's ceiling-index method
#' @keywords internal
.emulate_pctile <- function(x, pct) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0L) return(NA_real_)
  x <- sort(x)
  idx <- max(1L, ceiling(n * pct / 100))
  idx <- min(idx, n)
  x[idx]
}
