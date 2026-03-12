# V23: Formal Equivalence Testing (TOST)
#
# Uses Two One-Sided Tests (TOST) where appropriate, and absolute agreement
# checks where TOST is not applicable (same-data comparisons with correlated
# estimates inflate pooled SE, making TOST underpowered).
#
# Equivalence margins (delta) and rationale:
#   emulate vs TrialEmulation ITT:  |diff| < 0.10 (absolute, known algorithmic diffs)
#   emulate vs TrialEmulation PP:   |diff| < 0.25 (absolute, weight model strata differ)
#   emulate vs DGP truth ITT:       delta = 0.20 (TOST, one-sample test vs constant)
#   MC emulate vs TE across seeds:  delta = 0.10 (TOST on paired differences)
#   emulate vs TE on NHEFS:         |diff| < 0.15 (absolute, real data)

# ---------------------------------------------------------------------------
# Test 1: emulate vs TrialEmulation ITT — absolute agreement
# ---------------------------------------------------------------------------
test_that("V23-1: emulate vs TE ITT agree within 0.10 on golden DGP", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  d <- load_golden_dgp()

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
  te <- suppressWarnings(TrialEmulation::initiators(
    data.frame(d), id = "id", period = "period", treatment = "treatment",
    outcome = "outcome", eligible = "eligible", estimand_type = "ITT",
    outcome_cov = ~ x,
    include_followup_time = ~ followup_time + I(followup_time^2),
    include_trial_period = ~ trial_period + I(trial_period^2),
    model_var = "assigned_treatment", use_censor_weights = FALSE,
    data_dir = tempdir(), quiet = TRUE
  ))
  te_coefs <- te$robust$summary
  te_row <- te_coefs[te_coefs$names == "assigned_treatment", ]
  b_te <- te_row$estimate

  # Absolute agreement within 0.10
  expect_true(abs(b_em - b_te) < 0.10,
              info = sprintf("emulate=%.4f, TE=%.4f, diff=%.4f",
                             b_em, b_te, b_em - b_te))
})

# ---------------------------------------------------------------------------
# Test 2: emulate vs TrialEmulation PP — absolute agreement
# ---------------------------------------------------------------------------
test_that("V23-2: emulate vs TE PP agree within 0.25 on golden DGP", {
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
  te <- suppressWarnings(TrialEmulation::initiators(
    data.frame(d), id = "id", period = "period", treatment = "treatment",
    outcome = "outcome", eligible = "eligible", estimand_type = "PP",
    outcome_cov = ~ x, switch_d_cov = ~ x, switch_n_cov = ~ x,
    include_followup_time = ~ followup_time + I(followup_time^2),
    include_trial_period = ~ trial_period + I(trial_period^2),
    model_var = "assigned_treatment", use_censor_weights = FALSE,
    data_dir = tempdir(), quiet = TRUE
  ))
  te_coefs <- te$robust$summary
  te_row <- te_coefs[te_coefs$names == "assigned_treatment", ]
  b_te <- te_row$estimate

  # Absolute agreement within 0.25 (wider for PP due to weight model differences)
  expect_true(abs(b_em - b_te) < 0.25,
              info = sprintf("PP emulate=%.4f, TE=%.4f, diff=%.4f",
                             b_em, b_te, b_em - b_te))
})

# ---------------------------------------------------------------------------
# Test 3: TOST — emulate ITT recovers true DGP effect
# ---------------------------------------------------------------------------
test_that("V23-3: TOST equivalence — emulate ITT vs true effect (delta=0.20)", {
  skip_on_cran()

  d <- load_golden_dgp()
  true_effect <- -0.50

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  b_em  <- obj$model$b_treat
  se_em <- obj$model$se_treat

  # TOST against known truth (se_truth = 0 since it's a constant)
  result <- tost_z(b_em, true_effect, se_em, 0, delta = 0.20)

  expect_true(result$equivalent,
              info = sprintf("TOST failed: emulate=%.4f vs true=%.2f, p=%.4f",
                             b_em, true_effect, result$p_tost))
})

# ---------------------------------------------------------------------------
# Test 4: MC-based TOST — emulate vs TrialEmulation across 20 DGP seeds
# ---------------------------------------------------------------------------
test_that("V23-4: MC TOST — emulate vs TE across 20 seeds (delta=0.10)", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  diffs <- numeric(20)

  for (k in seq_len(20)) {
    seed_k <- 20260300 + k
    d <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = seed_k,
                     outcome_intercept = -3.5)

    # emulate
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment", outcome = "outcome",
                                         eligible = "eligible", covariates = "x",
                                         estimand = "ITT"))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
    obj <- suppressMessages(emulate_weight(obj))
    obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
    b_em <- obj$model$b_treat

    # TrialEmulation
    te <- tryCatch(
      suppressWarnings(TrialEmulation::initiators(
        data.frame(d), id = "id", period = "period", treatment = "treatment",
        outcome = "outcome", eligible = "eligible", estimand_type = "ITT",
        outcome_cov = ~ x,
        include_followup_time = ~ followup_time + I(followup_time^2),
        include_trial_period = ~ trial_period + I(trial_period^2),
        model_var = "assigned_treatment", use_censor_weights = FALSE,
        data_dir = tempdir(), quiet = TRUE
      )),
      error = function(e) NULL
    )

    if (!is.null(te)) {
      te_row <- te$robust$summary[te$robust$summary$names == "assigned_treatment", ]
      b_te <- te_row$estimate
      diffs[k] <- b_em - b_te
    } else {
      diffs[k] <- NA
    }
  }

  # MC TOST on paired differences
  result <- tost_mc(diffs, delta = 0.10)

  expect_true(result$equivalent,
              info = sprintf("MC TOST failed: mean_diff=%.4f, p=%.4f, n=%d",
                             result$mean_diff, result$p_tost, result$n))
})

# ---------------------------------------------------------------------------
# Test 5: emulate vs TE on NHEFS (absolute agreement)
# ---------------------------------------------------------------------------
test_that("V23-5: emulate vs TE agree within 0.15 on NHEFS", {
  skip_on_cran()
  skip_if_not_installed("TrialEmulation")

  d <- load_nhefs()
  covs <- c("age_std", "sex", "race", "smoke_cat", "wt71_std", "smokeyrs_std")

  # emulate
  obj <- suppressMessages(emulate_prepare(d, id = "seqn", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = covs,
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = covs,
                                   followup_spec = "quadratic",
                                   trial_period_spec = "none"))
  b_em <- obj$model$b_treat

  # TrialEmulation
  d_te <- d
  names(d_te)[names(d_te) == "seqn"] <- "id"
  te <- suppressWarnings(TrialEmulation::initiators(
    data.frame(d_te), id = "id", period = "period", treatment = "treatment",
    outcome = "outcome", eligible = "eligible", estimand_type = "ITT",
    outcome_cov = as.formula(paste("~", paste(covs, collapse = " + "))),
    include_followup_time = ~ followup_time + I(followup_time^2),
    include_trial_period = ~ 1,
    model_var = "assigned_treatment", use_censor_weights = FALSE,
    data_dir = tempdir(), quiet = TRUE
  ))
  te_coefs <- te$robust$summary
  te_row <- te_coefs[te_coefs$names == "assigned_treatment", ]
  b_te <- te_row$estimate

  # Absolute agreement within 0.15
  expect_true(abs(b_em - b_te) < 0.15,
              info = sprintf("NHEFS emulate=%.4f, TE=%.4f, diff=%.4f",
                             b_em, b_te, b_em - b_te))
})
