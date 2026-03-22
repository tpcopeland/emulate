#' Propensity score stratification for target trial emulation
#'
#' An alternative to inverse probability weighting that stratifies observations
#' into propensity score quantile bins and fits the outcome model with stratum
#' as a fixed effect (or \code{strata()} for Cox models).
#'
#' @section How stratification works:
#' \enumerate{
#'   \item Propensity scores are estimated on the expanded data.
#'   \item Observations are divided into \code{n_strata} quantile-based strata.
#'   \item The outcome model is fitted with stratum as a factor covariate
#'     (logistic/Poisson) or using \code{survival::strata()} (Cox).
#' }
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_expand}}.
#'   The \code{expanded} state flag must be \code{TRUE}.
#' @param strat_cov A character vector of covariate names for propensity score
#'   estimation. Required.
#' @param n_strata Number of propensity score strata. Default \code{5L}
#'   (quintiles).
#' @param ps_method Propensity score estimation method: \code{"glm"} (default)
#'   or \code{"lasso"}.
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
#'     \item \code{state$stratified} set to \code{TRUE}
#'     \item \code{state$fitted} set to \code{TRUE}
#'     \item \code{stratification} slot populated with strata info and PS values
#'     \item \code{model} slot populated with fitted outcome model
#'   }
#'
#' @seealso \code{\link{emulate_weight}} for the IPTW alternative,
#'   \code{\link{emulate_match}} for matching.
#'
#' @export
emulate_stratify <- function(obj, strat_cov, n_strata = 5L,
                              ps_method = "glm",
                              outcome_cov = NULL, model = "logistic",
                              followup_spec = "quadratic",
                              trial_period_spec = "quadratic",
                              cluster = NULL, level = 95, quiet = FALSE) {
  .check_expanded(obj)

  if (missing(strat_cov) || is.null(strat_cov) || length(strat_cov) == 0) {
    stop("strat_cov is required for propensity score stratification",
         call. = FALSE)
  }

  if (!ps_method %in% c("glm", "lasso")) {
    stop("ps_method must be 'glm' or 'lasso'", call. = FALSE)
  }

  if (n_strata < 2L) {
    stop("n_strata must be at least 2", call. = FALSE)
  }

  s <- obj$settings
  dt <- data.table::copy(obj$data)
  prefix <- s$prefix
  id_col <- s$id

  arm_col <- paste0(prefix, "arm")
  fu_col <- paste0(prefix, "followup")
  trial_col <- paste0(prefix, "trial")
  cens_col <- paste0(prefix, "censored")
  outobs_col <- paste0(prefix, "outcome_obs")

  if (is.null(cluster)) cluster <- id_col
  tr_ns_info <- NULL

  message("emulate_stratify - Propensity Score Stratification")
  message(strrep("-", 50))
  message("Strata:        ", n_strata)
  message("PS method:     ", ps_method)

  # Build PS formula
  bq <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))
  ps_formula <- as.formula(paste(bq(arm_col), "~",
                                  paste(bq(strat_cov), collapse = " + ")))

  # Fit PS
  if (ps_method == "lasso") {
    message("Fitting lasso propensity score model...")
    ps_result <- .fit_ps_lasso(ps_formula, data = as.data.frame(dt))
    ps_values <- ps_result$ps
  } else {
    message("Fitting GLM propensity score model...")
    ps_fit <- tryCatch(
      suppressWarnings(glm(ps_formula, data = dt, family = binomial())),
      error = function(e) NULL
    )
    if (!is.null(ps_fit)) {
      ps_values <- predict(ps_fit, type = "response")
    } else {
      stop("GLM propensity score model failed to converge", call. = FALSE)
    }
  }

  # Create strata by PS quantiles
  breaks <- stats::quantile(ps_values, probs = seq(0, 1, length.out = n_strata + 1L),
                             na.rm = TRUE)
  # Ensure unique breaks
  breaks <- unique(breaks)
  actual_strata <- length(breaks) - 1L

  if (actual_strata < 2L) {
    warning("Propensity scores have too little variation for ", n_strata,
            " strata. Using ", actual_strata, " strata.", call. = FALSE)
  }

  stratum_ids <- cut(ps_values, breaks = breaks, labels = FALSE,
                     include.lowest = TRUE)
  stratum_ids[is.na(stratum_ids)] <- 1L

  dt[, ps_stratum := factor(stratum_ids)]

  message("Created:       ", actual_strata, " strata")

  # Report stratum sizes
  strat_tab <- table(dt$ps_stratum, dt[[arm_col]])
  for (st in seq_len(actual_strata)) {
    n0 <- if (as.character(st) %in% rownames(strat_tab)) strat_tab[as.character(st), "0"] else 0
    n1 <- if (as.character(st) %in% rownames(strat_tab)) strat_tab[as.character(st), "1"] else 0
    if (!quiet) message("  Stratum ", st, ": control=", n0, ", treated=", n1)
  }

  # Store in object before fitting
  obj$data <- dt
  obj$state$stratified <- TRUE
  obj$stratification <- list(
    ps_values = ps_values,
    stratum_ids = stratum_ids,
    n_strata = actual_strata,
    breaks = breaks,
    ps_method = ps_method
  )

  # Add ps_stratum to outcome covariates
  strat_outcome_cov <- c(outcome_cov, "ps_stratum")

  # For Cox models, use strata() term instead of fixed effect
  if (model == "cox") {
    # Build a custom formula with strata() for Cox
    message("Fitting Cox model with strata()...")

    .parse_ns <- function(spec) {
      m <- regmatches(spec, regexec("^ns\\(([0-9]+)\\)$", spec))[[1]]
      if (length(m) == 2) return(as.integer(m[2]))
      return(NULL)
    }

    model_var <- arm_col
    if (is.null(cluster)) cluster <- id_col

    # Build time vars
    time_vars <- character(0)
    # Trial period vars for Cox
    if (trial_period_spec != "none") {
      ns_df <- .parse_ns(trial_period_spec)
      if (!is.null(ns_df)) {
        tr_ns_info <- .emulate_natural_spline(dt, trial_col, ns_df,
                                               paste0(prefix, "tr_ns"))
        time_vars <- c(time_vars, tr_ns_info$varnames)
      } else {
        time_vars <- c(time_vars, trial_col)
        if (trial_period_spec %in% c("quadratic", "cubic")) {
          sq_col <- paste0(prefix, "trial_sq")
          dt[, (sq_col) := get(trial_col)^2]
          time_vars <- c(time_vars, sq_col)
        }
        if (trial_period_spec == "cubic") {
          cu_col <- paste0(prefix, "trial_cu")
          dt[, (cu_col) := get(trial_col)^3]
          time_vars <- c(time_vars, cu_col)
        }
      }
    }

    cox_covars <- c(model_var, time_vars)
    if (!is.null(outcome_cov)) cox_covars <- c(cox_covars, outcome_cov)

    # Build Surv formula with strata
    dt[, ..time_enter := get(fu_col)]
    dt[, ..time_exit := get(fu_col) + 1L]

    est_sample <- dt[[cens_col]] == 0
    rhs <- paste(c(bq(cox_covars), "survival::strata(ps_stratum)"),
                 collapse = " + ")
    surv_f <- as.formula(paste0("survival::Surv(..time_enter, ..time_exit, ",
                                 bq(outobs_col), ") ~ ", rhs))

    est_dt <- dt[est_sample]
    fit <- survival::coxph(surv_f, data = est_dt,
                            cluster = est_dt[[cluster]])
    vcov_cl <- vcov(fit)

    dt[, c("..time_enter", "..time_exit") := NULL]

    # Mark estimation sample
    esample_col <- paste0(prefix, "esample")
    dt[, (esample_col) := 0L]
    dt[est_sample, (esample_col) := 1L]

    # Treatment effect
    model_var_coef <- bq(model_var)
    b_treat <- coef(fit)[model_var_coef]
    se_treat <- sqrt(vcov_cl[model_var_coef, model_var_coef])

    obj$data <- dt
    obj$state$fitted <- TRUE
    obj$model <- list(
      type = model,
      object = fit,
      vcov = vcov_cl,
      model_var = model_var,
      followup_spec = followup_spec,
      trial_period_spec = trial_period_spec,
      outcome_cov = c(outcome_cov, "ps_stratum"),
      time_vars = time_vars,
      cluster = cluster,
      level = level,
      fu_ns_info = NULL,
      tr_ns_info = tr_ns_info,
      b_treat = unname(b_treat),
      se_treat = unname(se_treat),
      adjustment_method = "stratification"
    )

    z_crit <- qnorm((100 + level) / 200)
    hr <- exp(b_treat)
    hr_lo <- exp(b_treat - z_crit * se_treat)
    hr_hi <- exp(b_treat + z_crit * se_treat)
    message("\nTreatment effect (Cox with strata):")
    message("  Log-HR:      ", sprintf("%.4f", b_treat))
    message("  Hazard ratio:", sprintf("%.4f", hr),
            " (", level, "% CI: ", sprintf("%.4f", hr_lo),
            " - ", sprintf("%.4f", hr_hi), ")")
  } else {
    # For logistic/poisson, add ps_stratum as outcome covariate
    obj$data <- dt
    obj <- emulate_fit(obj, outcome_cov = strat_outcome_cov, model = model,
                       followup_spec = followup_spec,
                       trial_period_spec = trial_period_spec,
                       cluster = cluster, level = level)
    obj$model$adjustment_method <- "stratification"
  }

  message(strrep("-", 50))
  invisible(obj)
}
