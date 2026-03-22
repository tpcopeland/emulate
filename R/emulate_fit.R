#' Fit outcome model for target trial emulation
#'
#' Fits the outcome model on the expanded (and optionally weighted) data.
#' This is the core estimation step that produces the treatment effect
#' estimate. The function supports two model types: pooled logistic
#' regression (the default and most common choice for TTE) and Cox
#' proportional hazards. Standard errors are clustered at the individual
#' level to account for the correlation between repeated observations of
#' the same person across emulated trials.
#'
#' @section Pooled logistic regression:
#' The default model (\code{model = "logistic"}) is a logistic regression
#' fitted to the \strong{pooled} expanded dataset, where the dependent
#' variable is the binary outcome indicator at each follow-up time. This is
#' the discrete-time analogue of a survival model. Pooled logistic regression
#' is the standard approach recommended by Hernan and Robins (2020) for
#' target trial emulation because it naturally accommodates the
#' clone-censor-weight structure and supports G-formula predictions via
#' \code{\link{emulate_predict}}.
#'
#' The coefficient on the treatment variable (the arm indicator) is a
#' log-odds ratio. To convert to an odds ratio, exponentiate the coefficient
#' (or use \code{eform = TRUE} in \code{\link{emulate_report}}).
#'
#' @section Cox proportional hazards:
#' Alternatively, \code{model = "cox"} fits a Cox PH model using counting
#' process notation (\code{Surv(time_enter, time_exit, event)}). The
#' coefficient on the treatment variable is a log-hazard ratio. Note that
#' G-formula predictions via \code{\link{emulate_predict}} are currently only
#' supported for the logistic model.
#'
#' @section Time specifications:
#' Both follow-up time and trial period can be modeled flexibly:
#' \describe{
#'   \item{\code{"linear"}}{A single linear term.}
#'   \item{\code{"quadratic"}}{Linear + squared terms (the default). This
#'     allows the hazard to change nonlinearly over time, which is common
#'     in practice.}
#'   \item{\code{"cubic"}}{Linear + squared + cubed terms.}
#'   \item{\code{"ns(k)"}}{A natural (restricted cubic) spline with \code{k}
#'     degrees of freedom (e.g., \code{"ns(3)"}). Uses Harrell's restricted
#'     cubic spline basis, matching the Stata implementation exactly. This
#'     is useful for highly nonlinear time trends.}
#'   \item{\code{"none"}}{No time terms included (rarely appropriate).}
#' }
#' The follow-up time specification controls how the baseline hazard varies
#' over time within each trial. The trial period specification controls how
#' the baseline hazard varies across different trial start times (capturing
#' secular trends).
#'
#' @section Clustered standard errors:
#' Because the expanded dataset contains multiple rows per individual (from
#' different emulated trials and follow-up times), observations are not
#' independent. The function computes cluster-robust standard errors using
#' the \pkg{sandwich} package (\code{vcovCL} with \code{type = "HC1"} and
#' \code{cadjust = TRUE}), clustered by default on the individual ID
#' variable. This matches Stata's \code{vce(cluster)} behavior.
#'
#' @section Estimation sample:
#' The model is fitted only on \strong{uncensored} observations (rows where
#' \code{_emulate_censored == 0}). Censored observations carry the information
#' about whether someone deviated but do not contribute to the outcome model.
#' An estimation sample indicator (\code{_emulate_esample}) is added to the data.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_expand}} or
#'   \code{\link{emulate_weight}}. The \code{expanded} state flag must be
#'   \code{TRUE}. If weights have been computed, they are automatically used.
#' @param outcome_cov An optional character vector of covariate names to
#'   include as additional predictors in the outcome model (beyond the
#'   treatment indicator and time terms). These are typically baseline or
#'   time-varying covariates that predict the outcome.
#' @param model A character string specifying the model type. One of:
#'   \describe{
#'     \item{\code{"logistic"}}{(Default) Pooled logistic regression
#'       (\code{glm} with \code{family = binomial()}).}
#'     \item{\code{"cox"}}{Cox proportional hazards model
#'       (\code{survival::coxph}) in counting process format.}
#'   }
#' @param model_var An optional character string specifying the name of the
#'   treatment variable in the model. Defaults to the arm indicator variable
#'   (e.g., \code{"_emulate_arm"}). You generally do not need to change this.
#' @param trial_period_spec A character string specifying how the trial
#'   period (calendar time) enters the model. One of \code{"linear"},
#'   \code{"quadratic"} (default), \code{"cubic"}, \code{"ns(k)"} (where
#'   \code{k} is an integer, e.g., \code{"ns(3)"}), or \code{"none"}.
#' @param followup_spec A character string specifying how follow-up time
#'   enters the model. Same options as \code{trial_period_spec}. Default is
#'   \code{"quadratic"}.
#' @param cluster An optional character string giving the name of the
#'   clustering variable for robust standard errors. Defaults to the
#'   individual ID variable set in \code{\link{emulate_prepare}}.
#' @param level A numeric value specifying the confidence level as a
#'   percentage. Default is \code{95} (for 95 percent confidence intervals).
#'
#' @return The input \code{emulate} object with \code{state$fitted} set to
#'   \code{TRUE} and the \code{model} slot populated with:
#'   \describe{
#'     \item{\code{type}}{Character: \code{"logistic"} or \code{"cox"}.}
#'     \item{\code{object}}{The fitted model object (\code{glm} or
#'       \code{coxph}).}
#'     \item{\code{vcov}}{The cluster-robust variance-covariance matrix.}
#'     \item{\code{model_var}}{Character: the treatment variable name.}
#'     \item{followup_spec}{The follow-up time specification used.}
#'     \item{trial_period_spec}{The trial period specification used.}
#'     \item{\code{outcome_cov}}{Character vector of outcome covariates.}
#'     \item{\code{time_vars}}{Character vector of all time-related variables
#'       included in the model.}
#'     \item{\code{cluster}}{Character: the clustering variable.}
#'     \item{\code{level}}{Numeric: the confidence level.}
#'     \item{fu_ns_info}{Follow-up spline knot information (if natural
#'       splines were used), or \code{NULL}.}
#'     \item{tr_ns_info}{Trial period spline knot information (if natural
#'       splines were used), or \code{NULL}.}
#'     \item{\code{b_treat}}{Numeric: the treatment effect coefficient.}
#'     \item{\code{se_treat}}{Numeric: the cluster-robust SE of the treatment
#'       effect.}
#'   }
#'   The data is also updated with an estimation sample indicator and any
#'   generated time variables (polynomial or spline terms). The object is
#'   returned invisibly.
#'
#' @seealso \code{\link{emulate_weight}} for the previous step,
#'   \code{\link{emulate_predict}} for marginal predictions,
#'   \code{\link{emulate_report}} for formatted results.
#'
#' @references
#' Hernan MA, Robins JM (2020). \emph{Causal Inference: What If}.
#' Chapman & Hall/CRC. Chapter 17.
#'
#' @examples
#' \donttest{
#' # Minimal pipeline: prepare -> expand -> fit
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
#'
#' # Fit pooled logistic with quadratic follow-up time
#' obj <- emulate_fit(obj, outcome_cov = "age", followup_spec = "quadratic")
#'
#' # Fit with natural splines for follow-up
#' # obj <- emulate_fit(obj, outcome_cov = "age", followup_spec = "ns(3)")
#'
#' # Fit Cox PH model
#' # obj <- emulate_fit(obj, outcome_cov = "age", model = "cox")
#' }
#'
#' @export
emulate_fit <- function(obj, outcome_cov = NULL, model = "logistic",
                    model_var = NULL, trial_period_spec = "quadratic",
                    followup_spec = "quadratic", cluster = NULL,
                    level = 95) {
  .check_expanded(obj)

  s <- obj$settings
  dt <- obj$data
  prefix <- s$prefix
  id_col <- s$id
  estimand <- s$estimand

  fu_col <- paste0(prefix, "followup")
  trial_col <- paste0(prefix, "trial")
  arm_col <- paste0(prefix, "arm")
  cens_col <- paste0(prefix, "censored")
  outobs_col <- paste0(prefix, "outcome_obs")

  if (is.null(model_var)) model_var <- arm_col
  if (is.null(cluster)) cluster <- id_col

  if (!model %in% c("logistic", "cox", "poisson")) {
    stop("model must be 'logistic', 'cox', or 'poisson'", call. = FALSE)
  }

  # Validate specs
  valid_specs <- c("linear", "quadratic", "cubic", "none")
  .parse_ns <- function(spec) {
    m <- regmatches(spec, regexec("^ns\\(([0-9]+)\\)$", spec))[[1]]
    if (length(m) == 2) return(as.integer(m[2]))
    return(NULL)
  }

  for (spec_name in c("followup_spec", "trial_period_spec")) {
    spec_val <- get(spec_name)
    if (!spec_val %in% valid_specs && is.null(.parse_ns(spec_val))) {
      stop(spec_name, " must be linear, quadratic, cubic, ns(#), or none",
           call. = FALSE)
    }
  }

  # Check weight variable
  weight_col <- paste0(prefix, "weight")
  has_weights <- weight_col %in% names(dt)
  if (!has_weights && estimand != "ITT") {
    warning("No weight variable found for ", estimand, " estimand. ",
            "Unweighted analysis is generally biased.", call. = FALSE)
  }

  message("emulate_fit - Outcome Model")
  message(strrep("-", 50))
  message("Model type:     ", model)
  message("Estimand:       ", estimand)
  message("Follow-up spec: ", followup_spec)
  message("Trial spec:     ", trial_period_spec)
  if (!is.null(outcome_cov)) message("Covariates:     ", paste(outcome_cov, collapse = " "))

  # Build time specification variables
  time_vars <- character(0)
  fu_ns_info <- NULL
  tr_ns_info <- NULL

  # Follow-up time
  if (followup_spec != "none") {
    ns_df <- .parse_ns(followup_spec)
    if (!is.null(ns_df)) {
      fu_ns_info <- .emulate_natural_spline(dt, fu_col, ns_df,
                                         paste0(prefix, "fu_ns"))
      time_vars <- c(time_vars, fu_ns_info$varnames)
    } else {
      time_vars <- c(time_vars, fu_col)
      if (followup_spec %in% c("quadratic", "cubic")) {
        sq_col <- paste0(prefix, "followup_sq")
        dt[, (sq_col) := get(fu_col)^2]
        time_vars <- c(time_vars, sq_col)
      }
      if (followup_spec == "cubic") {
        cu_col <- paste0(prefix, "followup_cu")
        dt[, (cu_col) := get(fu_col)^3]
        time_vars <- c(time_vars, cu_col)
      }
    }
  }

  # Trial period
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

  # Build covariate list
  all_covars <- c(model_var, time_vars)
  if (!is.null(outcome_cov)) all_covars <- c(all_covars, outcome_cov)

  # Estimation sample: exclude censored observations
  est_sample <- dt[[cens_col]] == 0
  n_est <- sum(est_sample, na.rm = TRUE)
  if (n_est == 0) {
    stop("Estimation sample is empty (all observations censored). ",
         "Check expansion settings or grace period.", call. = FALSE)
  }
  n_events_est <- sum(dt[[outobs_col]][est_sample] == 1, na.rm = TRUE)
  if (n_events_est == 0) {
    stop("No outcome events in estimation sample. ",
         "Cannot fit the model.", call. = FALSE)
  }

  # Helper to backtick-quote non-syntactic names (starting with _ etc.)
  bq <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))

  if (model == "logistic") {
    message("Fitting pooled logistic regression...")

    rhs <- paste(bq(all_covars), collapse = " + ")
    f <- as.formula(paste(bq(outobs_col), "~", rhs))

    if (has_weights) {
      fit <- glm(f, data = dt[est_sample], family = binomial(),
                 weights = dt[[weight_col]][est_sample])
    } else {
      fit <- glm(f, data = dt[est_sample], family = binomial())
    }

    # Clustered SEs matching Stata: HC1 with cadjust
    cluster_var <- dt[[cluster]][est_sample]
    vcov_cl <- sandwich::vcovCL(fit, cluster = cluster_var,
                                 type = "HC1", cadjust = TRUE)
  } else if (model == "poisson") {
    message("Fitting Poisson regression (log-linear risk model)...")

    rhs <- paste(bq(all_covars), collapse = " + ")
    f <- as.formula(paste(bq(outobs_col), "~", rhs))

    if (has_weights) {
      fit <- glm(f, data = dt[est_sample], family = poisson(link = "log"),
                 weights = dt[[weight_col]][est_sample])
    } else {
      fit <- glm(f, data = dt[est_sample], family = poisson(link = "log"))
    }

    # Clustered SEs with sandwich
    cluster_var <- dt[[cluster]][est_sample]
    vcov_cl <- sandwich::vcovCL(fit, cluster = cluster_var,
                                 type = "HC1", cadjust = TRUE)
  } else {
    # Cox model
    message("Fitting Cox proportional hazards model...")

    # Build Cox covariates (exclude followup time vars, keep trial period vars)
    cox_covars <- model_var
    if (trial_period_spec != "none") {
      if (!is.null(tr_ns_info)) {
        cox_covars <- c(cox_covars, tr_ns_info$varnames)
      } else {
        cox_covars <- c(cox_covars, trial_col)
        if (trial_period_spec %in% c("quadratic", "cubic"))
          cox_covars <- c(cox_covars, paste0(prefix, "trial_sq"))
        if (trial_period_spec == "cubic")
          cox_covars <- c(cox_covars, paste0(prefix, "trial_cu"))
      }
    }
    if (!is.null(outcome_cov)) cox_covars <- c(cox_covars, outcome_cov)

    # Counting process format
    dt[, ..time_enter := get(fu_col)]
    dt[, ..time_exit := get(fu_col) + 1L]

    rhs <- paste(bq(cox_covars), collapse = " + ")
    surv_f <- as.formula(paste0("survival::Surv(..time_enter, ..time_exit, ",
                                 bq(outobs_col), ") ~ ", rhs))

    est_dt <- dt[est_sample]
    if (has_weights) {
      fit <- survival::coxph(surv_f, data = est_dt,
                              weights = est_dt[[weight_col]],
                              cluster = est_dt[[cluster]])
    } else {
      fit <- survival::coxph(surv_f, data = est_dt,
                              cluster = est_dt[[cluster]])
    }

    vcov_cl <- vcov(fit)
    dt[, c("..time_enter", "..time_exit") := NULL]
  }

  # Mark estimation sample
  esample_col <- paste0(prefix, "esample")
  dt[, (esample_col) := 0L]
  dt[est_sample, (esample_col) := 1L]

  # Treatment effect — coefficient names may be backtick-quoted
  model_var_coef <- bq(model_var)
  b_treat <- coef(fit)[model_var_coef]
  if (is.na(b_treat)) {
    stop("Treatment variable '", model_var, "' was dropped from the model ",
         "(coefficient is NA). This typically means the treatment variable is ",
         "collinear with other predictors, or all observations have the same ",
         "treatment value. Check the expanded data.", call. = FALSE)
  }
  if (!model_var_coef %in% rownames(vcov_cl)) {
    stop("Treatment variable '", model_var, "' not found in the variance matrix. ",
         "The coefficient may have been aliased (dropped due to collinearity).",
         call. = FALSE)
  }
  se_treat <- sqrt(vcov_cl[model_var_coef, model_var_coef])
  z <- b_treat / se_treat
  p_val <- 2 * pnorm(-abs(z))
  z_crit <- qnorm((100 + level) / 200)

  message("\nTreatment effect:")
  if (model == "logistic") {
    or <- exp(b_treat)
    or_lo <- exp(b_treat - z_crit * se_treat)
    or_hi <- exp(b_treat + z_crit * se_treat)
    message("  Log-odds:   ", sprintf("%.4f", b_treat),
            " (SE: ", sprintf("%.4f", se_treat), ")")
    message("  Odds ratio: ", sprintf("%.4f", or),
            " (", level, "% CI: ", sprintf("%.4f", or_lo),
            " - ", sprintf("%.4f", or_hi), ")")
  } else if (model == "poisson") {
    rr <- exp(b_treat)
    rr_lo <- exp(b_treat - z_crit * se_treat)
    rr_hi <- exp(b_treat + z_crit * se_treat)
    message("  Log-RR:     ", sprintf("%.4f", b_treat),
            " (SE: ", sprintf("%.4f", se_treat), ")")
    message("  Risk ratio: ", sprintf("%.4f", rr),
            " (", level, "% CI: ", sprintf("%.4f", rr_lo),
            " - ", sprintf("%.4f", rr_hi), ")")
  } else {
    hr <- exp(b_treat)
    hr_lo <- exp(b_treat - z_crit * se_treat)
    hr_hi <- exp(b_treat + z_crit * se_treat)
    message("  Log-HR:      ", sprintf("%.4f", b_treat),
            " (SE: ", sprintf("%.4f", se_treat), ")")
    message("  Hazard ratio:", sprintf("%.4f", hr),
            " (", level, "% CI: ", sprintf("%.4f", hr_lo),
            " - ", sprintf("%.4f", hr_hi), ")")
  }
  message("  p-value:    ", sprintf("%.4f", p_val))
  message(strrep("-", 50))

  # Store in object
  obj$data <- dt
  obj$state$fitted <- TRUE
  obj$model <- list(
    type = model,
    object = fit,
    vcov = vcov_cl,
    model_var = model_var,
    followup_spec = followup_spec,
    trial_period_spec = trial_period_spec,
    outcome_cov = outcome_cov,
    time_vars = time_vars,
    cluster = cluster,
    level = level,
    fu_ns_info = fu_ns_info,
    tr_ns_info = tr_ns_info,
    b_treat = unname(b_treat),
    se_treat = unname(se_treat)
  )

  invisible(obj)
}
