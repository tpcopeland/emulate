# Propensity score helpers shared by emulate_weight, emulate_match, emulate_stratify
# Provides lasso PS estimation, PS trimming, and preference scores

#' Fit lasso-regularized propensity score model
#'
#' Uses cross-validated L1 logistic regression via \pkg{glmnet} to estimate
#' propensity scores. Automatically selects penalty via 10-fold CV.
#'
#' @param formula A formula with treatment on the left and covariates on the right.
#' @param data A data.frame or data.table containing the variables in \code{formula}.
#' @param nfolds Number of cross-validation folds (default 10).
#' @return A list with components:
#'   \describe{
#'     \item{\code{fit}}{The \code{cv.glmnet} object.}
#'     \item{\code{ps}}{Numeric vector of fitted propensity scores.}
#'   }
#' @keywords internal
.fit_ps_lasso <- function(formula, data, nfolds = 10L) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for ps_method = 'lasso'. ",
         "Install with: install.packages('glmnet')", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  x <- stats::model.matrix(formula, data = data)[, -1, drop = FALSE]

  # Remove rows with NA
  complete <- stats::complete.cases(x, y)
  x_cc <- x[complete, , drop = FALSE]
  y_cc <- y[complete]

  if (length(unique(y_cc)) < 2L) {
    warning("Treatment variable has fewer than 2 levels in lasso PS. ",
            "Returning constant PS.", call. = FALSE)
    ps <- rep(mean(y_cc), nrow(data))
    return(list(fit = NULL, ps = ps))
  }

  # glmnet requires >= 2 columns; fall back to GLM for single predictor
  if (ncol(x_cc) < 2L) {
    glm_fit <- tryCatch(
      suppressWarnings(stats::glm(formula, data = data, family = stats::binomial())),
      error = function(e) NULL
    )
    if (!is.null(glm_fit)) {
      ps <- stats::predict(glm_fit, type = "response")
      return(list(fit = glm_fit, ps = as.numeric(ps)))
    } else {
      ps <- rep(mean(y_cc), nrow(data))
      return(list(fit = NULL, ps = ps))
    }
  }

  cv_fit <- glmnet::cv.glmnet(x_cc, y_cc, family = "binomial",
                                alpha = 1, nfolds = nfolds)

  # Predict on full data
  x_full <- stats::model.matrix(formula, data = data)[, -1, drop = FALSE]
  ps <- rep(NA_real_, nrow(data))
  ps[complete] <- as.numeric(
    stats::predict(cv_fit, newx = x_cc, s = "lambda.min", type = "response")
  )
  # For incomplete rows, use marginal prevalence
  if (any(!complete)) {
    ps[!complete] <- mean(y_cc)
  }

  list(fit = cv_fit, ps = ps)
}

#' Compute preference scores from propensity scores
#'
#' Transforms PS to preference score using the formula from Walker et al. (2013):
#' \code{preference = log(ps / (1 - ps)) - log(prop / (1 - prop))}, then
#' applies the logistic transformation to return values between 0 and 1.
#'
#' @param ps Numeric vector of propensity scores in the open interval (0, 1).
#' @param prop_treated Scalar: proportion of treated in the population.
#' @return Numeric vector of preference scores in the closed interval from 0 to 1.
#' @keywords internal
.compute_preference_score <- function(ps, prop_treated) {
  # Clip to avoid log(0) or division by zero
  ps <- pmax(pmin(ps, 1 - 1e-8), 1e-8)
  prop_treated <- max(min(prop_treated, 1 - 1e-8), 1e-8)

  # Log-odds difference
  log_odds_ps <- log(ps / (1 - ps))
  log_odds_prop <- log(prop_treated / (1 - prop_treated))

  # Preference score = logistic(log_odds_ps - log_odds_prop)
  plogis(log_odds_ps - log_odds_prop)
}

#' Trim propensity scores by quantile
#'
#' Returns a logical vector indicating which observations to keep based on
#' PS quantile bounds.
#'
#' @param ps Numeric vector of propensity scores.
#' @param quantiles Numeric vector of length 2 giving lower and upper quantile
#'   bounds (e.g., \code{c(0.05, 0.95)}).
#' @return A logical vector of the same length as \code{ps}. \code{TRUE} means
#'   the observation should be kept.
#' @keywords internal
.trim_ps <- function(ps, quantiles) {
  if (length(quantiles) != 2 || quantiles[1] >= quantiles[2]) {
    stop("ps_trim must be a length-2 vector with lower < upper (e.g., c(0.05, 0.95))",
         call. = FALSE)
  }
  if (any(quantiles < 0 | quantiles > 1)) {
    stop("ps_trim values must be between 0 and 1", call. = FALSE)
  }

  lo_val <- stats::quantile(ps, quantiles[1], na.rm = TRUE)
  hi_val <- stats::quantile(ps, quantiles[2], na.rm = TRUE)

  !is.na(ps) & ps >= lo_val & ps <= hi_val
}
