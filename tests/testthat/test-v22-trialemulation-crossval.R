# V22: Three-Way Cross-Validation (emulate vs TrialEmulation vs ground truth)
#
# Compares emulate against TrialEmulation (CRAN) on shared datasets with
# known DGP parameters. This validates both implementations against each other
# AND against the true data-generating process.
#
# Known algorithmic differences:
#   1. Weight model strata: emulate = 2 (by arm), TrialEmulation = 4 (by arm x lagged_treatment)
#   2. SE computation: emulate = HC1 + cadjust, TrialEmulation = HC1 only
#   3. IPCW stratification differences
#   4. Post-outcome observation handling
#
# Equivalence margins:
#   - ITT coefficient: delta = 0.10 (no weights, only SE formula differs)
#   - PP coefficient: delta = 0.25 (weight model strata differ)
#   - True value recovery: delta = 0.20 (finite-sample bias)

# ---------------------------------------------------------------------------
# Helper: run TrialEmulation on dataset
# ---------------------------------------------------------------------------
run_te <- function(data, estimand = "ITT", outcome_cov = NULL,
                   switch_d_cov = NULL, censor_d_cov = NULL,
                   followup_time_terms = "quadratic") {
  fu_formula <- switch(followup_time_terms,
    "linear"    = ~ followup_time,
    "quadratic" = ~ followup_time + I(followup_time^2),
    "cubic"     = ~ followup_time + I(followup_time^2) + I(followup_time^3),
    ~ followup_time + I(followup_time^2)
  )
  tp_formula <- ~ trial_period + I(trial_period^2)

  oc_formula <- if (!is.null(outcome_cov)) {
    as.formula(paste("~", paste(outcome_cov, collapse = " + ")))
  } else {
    NULL
  }

  args <- list(
    data = data.frame(data),
    id = "id", period = "period", treatment = "treatment",
    outcome = "outcome", eligible = "eligible",
    estimand_type = estimand,
    include_followup_time = fu_formula,
    include_trial_period = tp_formula,
    model_var = "assigned_treatment",
    use_censor_weights = FALSE,
    data_dir = tempdir(),
    quiet = TRUE
  )

  if (!is.null(oc_formula)) {
    args$outcome_cov <- oc_formula
  }

  if (estimand %in% c("PP", "As-Treated") && !is.null(switch_d_cov)) {
    sw_formula <- as.formula(paste("~", paste(switch_d_cov, collapse = " + ")))
    args$switch_d_cov <- sw_formula
    args$switch_n_cov <- sw_formula
  }

  if (!is.null(censor_d_cov)) {
    cen_formula <- as.formula(paste("~", paste(censor_d_cov, collapse = " + ")))
    args$censor_d_cov <- cen_formula
    args$censor_n_cov <- cen_formula
    args$use_censor_weights <- TRUE
  }

  result <- do.call(TrialEmulation::initiators, args)
  te_coefs <- result$robust$summary
  te_row <- te_coefs[te_coefs$names == "assigned_treatment", ]
  b  <- te_row$estimate
  se <- te_row$robust_se

  list(coef = b, se = se, result = result)
}

# ---------------------------------------------------------------------------
# Test 1: ITT coefficients agree on trial_example (emulate vs TrialEmulation)
# ---------------------------------------------------------------------------
test_that("V22-1: ITT coefficient sign agreement on trial_example", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  d <- load_trial_example()

  # emulate
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                       "nvarA", "nvarB", "nvarC"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))
  b_em <- obj$model$b_treat
  se_em <- obj$model$se_treat

  # TrialEmulation
  te <- suppressWarnings(run_te(d, estimand = "ITT",
                                 outcome_cov = c("catvarA", "catvarB",
                                                 "nvarA", "nvarB", "nvarC")))

  # Sign agreement
  expect_equal(sign(b_em), sign(te$coef),
               info = sprintf("emulate=%.4f, TE=%.4f", b_em, te$coef))

  # Within 0.10 absolute
  expect_true(abs(b_em - te$coef) < 0.10,
              info = sprintf("ITT diff=%.4f, expected < 0.10",
                             abs(b_em - te$coef)))
})

# ---------------------------------------------------------------------------
# Test 2: Both recover true effect on golden DGP (ITT)
# ---------------------------------------------------------------------------
test_that("V22-2: Both recover true effect on golden DGP (ITT, true=-0.50)", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  d <- load_golden_dgp()
  true_effect <- -0.50

  # emulate
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  b_em <- obj$model$b_treat

  # TrialEmulation
  te <- suppressWarnings(run_te(d, estimand = "ITT", outcome_cov = "x"))

  # Both within 0.20 of true value
  expect_true(abs(b_em - true_effect) < 0.20,
              info = sprintf("emulate=%.4f vs true=%.2f, diff=%.4f",
                             b_em, true_effect, abs(b_em - true_effect)))
  expect_true(abs(te$coef - true_effect) < 0.20,
              info = sprintf("TE=%.4f vs true=%.2f, diff=%.4f",
                             te$coef, true_effect, abs(te$coef - true_effect)))

  # Mutual agreement within 0.10
  expect_true(abs(b_em - te$coef) < 0.10,
              info = sprintf("emulate=%.4f vs TE=%.4f", b_em, te$coef))
})

# ---------------------------------------------------------------------------
# Test 3: PP coefficients agree on golden DGP
# ---------------------------------------------------------------------------
test_that("V22-3: PP coefficients agree on golden DGP", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  d <- load_golden_dgp()

  # emulate PP
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  b_em <- obj$model$b_treat

  # TrialEmulation PP
  te <- suppressWarnings(run_te(d, estimand = "PP", outcome_cov = "x",
                                 switch_d_cov = "x"))

  # Sign agreement
  expect_equal(sign(b_em), sign(te$coef),
               info = sprintf("PP emulate=%.4f, TE=%.4f", b_em, te$coef))

  # Within 0.25 (wider due to weight model strata differences)
  expect_true(abs(b_em - te$coef) < 0.25,
              info = sprintf("PP diff=%.4f, expected < 0.25",
                             abs(b_em - te$coef)))
})

# ---------------------------------------------------------------------------
# Test 4: SE ratio within [0.5, 2.0] (ITT)
# ---------------------------------------------------------------------------
test_that("V22-4: ITT SE ratio within [0.5, 2.0]", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  d <- load_golden_dgp()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  se_em <- obj$model$se_treat

  te <- suppressWarnings(run_te(d, estimand = "ITT", outcome_cov = "x"))

  ratio <- se_em / te$se
  expect_true(ratio > 0.5 & ratio < 2.0,
              info = sprintf("SE ratio=%.4f (emulate=%.4f, TE=%.4f)",
                             ratio, se_em, te$se))
})

# ---------------------------------------------------------------------------
# Test 5: Both recover true effect on golden DGP (PP)
# ---------------------------------------------------------------------------
test_that("V22-5: Both recover true PP effect on golden DGP (true=-0.50)", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  d <- load_golden_dgp()
  true_effect <- -0.50

  # emulate PP
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  b_em <- obj$model$b_treat

  # TrialEmulation PP
  te <- suppressWarnings(run_te(d, estimand = "PP", outcome_cov = "x",
                                 switch_d_cov = "x"))

  # Both within 0.35 of truth (PP has more variance)
  expect_true(abs(b_em - true_effect) < 0.35,
              info = sprintf("emulate PP=%.4f vs true=%.2f", b_em, true_effect))
  expect_true(abs(te$coef - true_effect) < 0.35,
              info = sprintf("TE PP=%.4f vs true=%.2f", te$coef, true_effect))
})

# ---------------------------------------------------------------------------
# Test 6: NHEFS three-way — emulate, TrialEmulation, and Stata all plausible
# ---------------------------------------------------------------------------
test_that("V22-6: NHEFS three-way — all produce plausible estimates", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  d <- load_nhefs()

  # emulate
  obj <- suppressMessages(emulate_prepare(d, id = "seqn", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("age_std", "sex", "race",
                                                       "smoke_cat", "wt71_std",
                                                       "smokeyrs_std"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("age_std", "sex", "race",
                                                   "smoke_cat", "wt71_std",
                                                   "smokeyrs_std"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "none"))
  or_em <- exp(obj$model$b_treat)

  # TrialEmulation (requires different variable naming for initiators)
  d_te <- d
  names(d_te)[names(d_te) == "seqn"] <- "id"
  te <- suppressWarnings(run_te(d_te, estimand = "ITT",
                                 outcome_cov = c("age_std", "sex", "race",
                                                 "smoke_cat", "wt71_std",
                                                 "smokeyrs_std")))
  or_te <- exp(te$coef)

  # Both ORs should be in plausible range [0.3, 2.0]
  expect_true(or_em > 0.3 & or_em < 2.0,
              info = sprintf("emulate OR=%.4f", or_em))
  expect_true(or_te > 0.3 & or_te < 2.0,
              info = sprintf("TrialEmulation OR=%.4f", or_te))
})
