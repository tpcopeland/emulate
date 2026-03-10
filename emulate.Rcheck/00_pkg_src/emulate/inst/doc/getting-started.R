## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-data----------------------------------------------------------------
library(emulate)

# Load the bundled trial_example.csv dataset
# This ships with emulate in inst/extdata/
trial_data <- read.csv(
  system.file("extdata", "trial_example.csv", package = "emulate")
)

## ----str-data-----------------------------------------------------------------
str(trial_data)

## ----summary-data-------------------------------------------------------------
summary(trial_data)

## ----dim-data-----------------------------------------------------------------
cat("Rows:", nrow(trial_data), "\n")
cat("Unique patients:", length(unique(trial_data$id)), "\n")
cat("Period range:", min(trial_data$period), "-", max(trial_data$period), "\n")
cat("Events:", sum(trial_data$outcome == 1), "\n")
cat("Eligible observations:", sum(trial_data$eligible == 1), "\n")

## ----protocol, eval = FALSE---------------------------------------------------
# protocol <- emulate_protocol(
#   eligibility    = "Adults meeting clinical criteria, not currently on treatment,
#                     no prior outcome event",
#   treatment      = "Strategy A: Initiate treatment at the start of each eligible
#                     period. Strategy B: Do not initiate treatment.",
#   assignment     = "Patients are assigned to treatment strategies at each period
#                     when they meet eligibility criteria.",
#   followup_start = "Start of the period when eligibility is first met for each
#                     emulated trial.",
#   outcome        = "First occurrence of the binary outcome event during follow-up.
#                     Maximum follow-up: 8 periods.",
#   causal_contrast = "Intention-to-treat effect and per-protocol effect of
#                      treatment initiation vs. no treatment.",
#   analysis       = "Sequential trial emulation with pooled logistic regression.
#                     PP analysis uses stabilized IPTW with truncation at 1st/99th
#                     percentiles. G-formula marginal predictions with 100 MC
#                     samples for confidence intervals.",
#   format = "display"
# )

## ----itt-prepare--------------------------------------------------------------
itt <- emulate_prepare(
  trial_data,
  id        = "id",
  period    = "period",
  treatment = "treatment",
  outcome   = "outcome",
  eligible  = "eligible",
  covariates = c("catvarA", "catvarB", "catvarC",
                 "nvarA", "nvarB", "nvarC"),
  estimand  = "ITT"
)

## ----itt-validate-------------------------------------------------------------
itt <- emulate_validate(itt)

## ----itt-expand---------------------------------------------------------------
itt <- emulate_expand(itt, maxfollowup = 8)

## ----itt-fit------------------------------------------------------------------
itt <- emulate_fit(
  itt,
  outcome_cov = c("catvarA", "catvarB", "catvarC",
                   "nvarA", "nvarB", "nvarC"),
  followup_spec      = "quadratic",
  trial_period_spec  = "quadratic"
)

## ----itt-predict--------------------------------------------------------------
itt <- emulate_predict(
  itt,
  times      = 0:8,
  type       = "cum_inc",
  difference = TRUE,
  samples    = 100,
  seed       = 12345
)

## ----itt-plot, eval = FALSE---------------------------------------------------
# emulate_plot(itt, type = "cumhaz", title = "ITT: Cumulative Incidence")

## ----itt-report---------------------------------------------------------------
emulate_report(itt, format = "display", eform = TRUE, predictions = TRUE)

## ----itt-excel, eval = FALSE--------------------------------------------------
# emulate_report(itt, format = "excel", export = "itt_results.xlsx",
#            eform = TRUE, predictions = TRUE)

## ----pp-prepare---------------------------------------------------------------
pp <- emulate_prepare(
  trial_data,
  id        = "id",
  period    = "period",
  treatment = "treatment",
  outcome   = "outcome",
  eligible  = "eligible",
  covariates = c("catvarA", "catvarB", "catvarC",
                 "nvarA", "nvarB", "nvarC"),
  estimand  = "PP"
)

## ----pp-validate--------------------------------------------------------------
pp <- emulate_validate(pp)

## ----pp-expand----------------------------------------------------------------
pp <- emulate_expand(pp, maxfollowup = 8)

## ----pp-weight----------------------------------------------------------------
pp <- emulate_weight(
  pp,
  switch_d_cov = c("catvarA", "catvarB", "catvarC",
                    "nvarA", "nvarB", "nvarC"),
  switch_n_cov = c("catvarA", "catvarB"),
  truncate     = c(1, 99),
  stabilized   = TRUE
)

## ----pp-fit-------------------------------------------------------------------
pp <- emulate_fit(
  pp,
  outcome_cov = c("catvarA", "catvarB", "catvarC",
                   "nvarA", "nvarB", "nvarC"),
  followup_spec      = "quadratic",
  trial_period_spec  = "quadratic"
)

## ----pp-predict---------------------------------------------------------------
pp <- emulate_predict(
  pp,
  times      = 0:8,
  type       = "cum_inc",
  difference = TRUE,
  samples    = 100,
  seed       = 12345
)

## ----pp-plot-ci, eval = FALSE-------------------------------------------------
# emulate_plot(pp, type = "cumhaz", title = "PP: Cumulative Incidence")

## ----pp-plot-weights, eval = FALSE--------------------------------------------
# emulate_plot(pp, type = "weights", title = "PP: IP Weight Distribution")

## ----pp-diagnose--------------------------------------------------------------
pp <- emulate_diagnose(
  pp,
  balance_covariates = c("catvarA", "catvarB", "catvarC",
                         "nvarA", "nvarB", "nvarC")
)

## ----pp-plot-balance, eval = FALSE--------------------------------------------
# emulate_plot(pp, type = "balance", title = "PP: Covariate Balance (Love Plot)")

## ----pp-report----------------------------------------------------------------
emulate_report(pp, format = "display", eform = TRUE, predictions = TRUE)

## ----compare------------------------------------------------------------------
cat("ITT treatment log-OR:", round(itt$model$b_treat, 4), "\n")
cat("PP  treatment log-OR:", round(pp$model$b_treat, 4), "\n")
cat("\n")
cat("ITT treatment OR:", round(exp(itt$model$b_treat), 4), "\n")
cat("PP  treatment OR:", round(exp(pp$model$b_treat), 4), "\n")

## ----cox-example, eval = FALSE------------------------------------------------
# # You can use a Cox PH model instead of pooled logistic regression
# obj <- emulate_fit(obj, model = "cox",
#                outcome_cov = c("x1", "x2"),
#                trial_period_spec = "quadratic")
# # Note: emulate_predict() currently only supports the logistic model.
# # After a Cox fit, predictions are not available.

## ----sensitivity, eval = FALSE------------------------------------------------
# # Grace period analysis
# pp_grace <- emulate_expand(pp, maxfollowup = 8, grace = 2)

## ----at-example, eval = FALSE-------------------------------------------------
# at <- emulate_prepare(trial_data, ..., estimand = "AT")
# at <- emulate_expand(at, maxfollowup = 8)
# at <- emulate_weight(at, switch_d_cov = c("x1", "x2"), truncate = c(1, 99))
# at <- emulate_fit(at, outcome_cov = c("x1", "x2"))

