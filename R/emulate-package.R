#' emulate: Target Trial Emulation via Sequential Trials
#'
#' @description
#' The \pkg{emulate} package implements the sequential trials framework for
#' **target trial emulation** (TTE) from observational data. Target trial
#' emulation is a structured approach to causal inference that designs an
#' observational study to mimic a hypothetical randomized controlled trial.
#' The method was formalized by Hernan and Robins (2016).
#'
#' The core idea is: when a randomized trial is not feasible, you can still
#' estimate causal treatment effects from longitudinal observational data by
#' (1) defining the protocol of the trial you \emph{wish} you could run, and
#' (2) using the observational data to emulate that trial as closely as
#' possible. The sequential trials design emulates a new trial at every
#' eligible time point, then pools the results.
#'
#' @section Pipeline workflow:
#' The package is organized as a step-by-step pipeline. Each function takes a
#' \code{emulate} object produced by the previous step and returns an updated
#' \code{emulate} object:
#'
#' \enumerate{
#'   \item \strong{\code{\link{emulate_prepare}}}: Map variable names, validate
#'     person-period structure, choose the estimand (ITT, PP, or AT).
#'   \item \strong{\code{\link{emulate_validate}}}: Run 10 data quality checks
#'     (duplicates, gaps, missing data, positivity, etc.).
#'   \item \strong{\code{\link{emulate_expand}}}: Clone-censor-weight expansion
#'     into sequential emulated trials.
#'   \item \strong{\code{\link{emulate_weight}}}: Compute inverse probability of
#'     treatment (IPTW) and censoring (IPCW) weights.
#'   \item \strong{\code{\link{emulate_fit}}}: Fit the outcome model (pooled
#'     logistic regression or Cox PH) with clustered standard errors.
#'   \item \strong{\code{\link{emulate_predict}}}: G-formula standardized
#'     marginal predictions with Monte Carlo confidence intervals.
#'   \item \strong{\code{\link{emulate_diagnose}}}: Weight diagnostics, effective
#'     sample size (ESS), and covariate balance via standardized mean
#'     differences (SMD).
#'   \item \strong{\code{\link{emulate_plot}}}: Visualization (Kaplan-Meier,
#'     cumulative incidence, weight distributions, Love plots).
#'   \item \strong{\code{\link{emulate_report}}}: Publication-quality results
#'     tables in display, CSV, or Excel format.
#' }
#'
#' A standalone helper, \code{\link{emulate_protocol}}, generates the 7-component
#' target trial protocol table (Hernan & Robins framework) and can be used
#' independently of the pipeline.
#'
#' @section Estimands:
#' The package supports three estimands:
#' \describe{
#'   \item{\strong{ITT} (Intention-to-Treat)}{Compares groups as originally
#'     assigned, regardless of whether they adhered to treatment. No cloning or
#'     weighting is needed.}
#'   \item{\strong{PP} (Per-Protocol)}{Compares groups who adhered to their
#'     assigned treatment strategy. Uses clone-censor-weight to handle
#'     treatment switching.}
#'   \item{\strong{AT} (As-Treated)}{Compares groups based on the treatment
#'     actually received. Also uses clone-censor-weight.}
#' }
#'
#' @section Input data requirements:
#' Data must be in \strong{person-period} (long) format with one row per
#' individual per time period. Required variables are: a patient identifier,
#' an integer-valued time period, binary treatment (0/1), binary outcome (0/1),
#' and binary eligibility (0/1). An optional binary censoring indicator and
#' time-varying or baseline covariates may also be provided.
#'
#' @section Key references:
#' \itemize{
#'   \item Hernan MA, Robins JM (2016). "Using Big Data to Emulate a Target
#'     Trial When a Randomized Trial Is Not Available." \emph{American Journal
#'     of Epidemiology}, 183(8), 758-764.
#'   \item Hernan MA, Robins JM (2020). \emph{Causal Inference: What If}.
#'     Chapman & Hall/CRC.
#'   \item Danaei G, Rodriguez LAG, Cantero OF, Logan R, Hernan MA (2013).
#'     "Observational data for comparative effectiveness research: An emulation
#'     of randomised trials of statins and primary prevention of coronary heart
#'     disease." \emph{Statistical Methods in Medical Research}, 22(1), 70-96.
#' }
#'
#' @name emulate-package
#' @aliases emulate
#' @import data.table
#' @importFrom stats glm binomial predict vcov coef pnorm qnorm rnorm
#'   quantile sd var weighted.mean complete.cases model.matrix formula
#'   setNames as.formula time
#' @importFrom survival coxph Surv survfit
#' @importFrom sandwich vcovCL
#' @importFrom MASS mvrnorm
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point
#'   geom_vline geom_hline geom_segment labs theme_minimal
#'   scale_color_manual scale_linetype_manual facet_wrap coord_flip
#'   element_text theme
"_PACKAGE"

NULL

# Suppress R CMD check NOTEs for data.table column references
utils::globalVariables(c(
  ".", ".N", ".I", ".GRP", ".SD", ":=",
  # Common column names used in data.table operations
  "followup", "trial", "arm", "censored", "outcome_obs", "weight",
  "esample", "lag_treat", "deviated", "first_dev", "cens_time",
  "min_cens", "denom_pr", "numer_pr", "sw_t", "log_sw", "cum_log_sw",
  "cw_t", "log_cw", "cum_log_cw", "w2", "cum_surv", "prob",
  "followup_sq", "followup_cu", "trial_sq", "trial_cu",
  "assigned", "baseline_treat", "is_elig_id", "elig_at_t",
  "..keep_cols", "V1",
  # R CMD check: ggplot aes variables
  "smd", "y", "type", "estimate", "group", "surv", "lower", "upper",
  # R CMD check: data.table computed columns
  "..time_enter", "..time_exit", "N", "pdiff", "cum_out", "prior_out",
  "prev_treat", "prior_out2", "bl_treat"
))
