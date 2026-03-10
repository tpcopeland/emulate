# V1: Cross-validation against R TrialEmulation reference values
# Mirrors: Stata-Tools/tte/qa/validate_trialemulation.do
#
# Dataset: trial_example.csv (503 patients, 48,400 person-periods)
#
# R TrialEmulation reference results (Maringe et al. 2024, arXiv:2402.12083):
#   ITT: Coefficient = -0.2829, Robust SE = 0.3138
#   PP:  Coefficient = -0.4143, Robust SE = 0.4152
#
# Tolerance rationale:
#   - ITT coefficient: 10% relative (no weighting, only SE formula differs)
#   - PP coefficient: 20% relative (weight model strata differ: R uses 4, emulate uses 2)
#   - SE: wider tolerance due to HC1 + G/(G-1) correction differences

# ---------------------------------------------------------------------------
# Test 1: ITT coefficient matches R TrialEmulation reference
# ---------------------------------------------------------------------------
test_that("V01-1: ITT coefficient matches R TrialEmulation (-0.2829 +/- 10%)", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                       "nvarA", "nvarB", "nvarC"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))

  # R TrialEmulation benchmark: -0.282945
  r_coef <- -0.2829
  r_se   <- 0.3138

  b  <- obj$model$b_treat
  se <- obj$model$se_treat

  # Coefficient within 10% relative tolerance
  rel_diff <- abs(b - r_coef) / abs(r_coef)
  expect_true(rel_diff < 0.10,
              info = sprintf("ITT coef %.4f not within 10%% of %.4f (rel diff: %.1f%%)",
                             b, r_coef, rel_diff * 100))

  # SE within 15% relative tolerance (HC1 + G/(G-1) differences)
  se_diff <- abs(se - r_se) / r_se
  expect_true(se_diff < 0.15,
              info = sprintf("ITT SE %.4f not within 15%% of %.4f (rel diff: %.1f%%)",
                             se, r_se, se_diff * 100))
})

# ---------------------------------------------------------------------------
# Test 2: ITT expansion produces >100k person-periods
# ---------------------------------------------------------------------------
test_that("V01-2: ITT expansion produces >100k person-periods", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                       "nvarA", "nvarB", "nvarC"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))

  n_expanded <- nrow(obj$data)

  expect_true(obj$state$expanded,
              info = "Object should be in expanded state")
  expect_true(n_expanded > 100000,
              info = sprintf("Expanded obs = %d, expected > 100,000", n_expanded))
  expect_true(obj$expansion$n_trials > 1,
              info = "Should produce multiple emulated trials")
})

# ---------------------------------------------------------------------------
# Test 3: ITT cumulative incidence is monotonically non-decreasing and in [0,1]
# ---------------------------------------------------------------------------
test_that("V01-3: ITT cumulative incidence is monotonic and bounded [0,1]", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                       "nvarA", "nvarB", "nvarC"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))
  obj <- suppressMessages(emulate_predict(obj, times = 0:8, type = "cum_inc",
                                       samples = 100, seed = 12345,
                                       difference = TRUE))

  pred <- obj$predictions

  # Control arm (est_0) should be monotonically non-decreasing
  diffs_0 <- diff(pred$est_0)
  expect_true(all(diffs_0 >= -0.001),
              info = "Control arm cumulative incidence not monotonically non-decreasing")

  # Treated arm (est_1) should be monotonically non-decreasing
  diffs_1 <- diff(pred$est_1)
  expect_true(all(diffs_1 >= -0.001),
              info = "Treated arm cumulative incidence not monotonically non-decreasing")

  # All values in [0, 1]
  expect_true(all(pred$est_0 >= 0 & pred$est_0 <= 1),
              info = "Control arm values outside [0,1]")
  expect_true(all(pred$est_1 >= 0 & pred$est_1 <= 1),
              info = "Treated arm values outside [0,1]")

  # CIs should bracket estimates (allowing small MC noise)
  expect_true(all(pred$ci_lo_0 <= pred$est_0 + 0.001),
              info = "Control arm CI lower bound exceeds estimate")
  expect_true(all(pred$ci_hi_0 >= pred$est_0 - 0.001),
              info = "Control arm CI upper bound below estimate")
})

# ---------------------------------------------------------------------------
# Test 4: PP coefficient matches R TrialEmulation reference (20% tolerance)
# ---------------------------------------------------------------------------
test_that("V01-4: PP coefficient matches R TrialEmulation (-0.4143 +/- 20%)", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                       "nvarA", "nvarB", "nvarC"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      switch_n_cov = c("nvarA", "nvarB"),
                                      quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))

  # R TrialEmulation PP benchmark: -0.414264
  r_pp_coef <- -0.4143
  r_pp_se   <- 0.4152

  b <- obj$model$b_treat

  # PP has more variance due to weighting; sign + magnitude check
  # plus 20% relative tolerance on coefficient
  expect_true(b < 0,
              info = sprintf("PP coefficient should be negative, got %.4f", b))

  pp_rdiff <- abs(b - r_pp_coef) / abs(r_pp_coef)
  expect_true(pp_rdiff < 0.20,
              info = sprintf("PP coef %.4f not within 20%% of %.4f (rel diff: %.1f%%)",
                             b, r_pp_coef, pp_rdiff * 100))

  # Magnitude in plausible range
  expect_true(abs(b) > 0.1 & abs(b) < 2.0,
              info = sprintf("PP coef magnitude %.4f outside plausible range [0.1, 2.0]",
                             abs(b)))
})

# ---------------------------------------------------------------------------
# Test 5: PP weights are non-degenerate (ESS > 100, mean weight 0.5-2.0)
# ---------------------------------------------------------------------------
test_that("V01-5: PP weights non-degenerate (ESS > 100, mean 0.5-2.0)", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                       "nvarA", "nvarB", "nvarC"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      switch_n_cov = c("nvarA", "nvarB"),
                                      quiet = TRUE))

  ess     <- obj$weights$ess
  mean_wt <- obj$weights$mean

  # ESS should be a reasonable fraction of total sample
  expect_true(ess > 100,
              info = sprintf("ESS = %.1f, expected > 100", ess))

  # Mean weight should be near 1 (stabilized weights)
  expect_true(mean_wt > 0.5 & mean_wt < 2.0,
              info = sprintf("Mean weight = %.4f, expected in [0.5, 2.0]", mean_wt))

  # Weights should be strictly positive
  expect_true(obj$weights$min > 0,
              info = "Minimum weight should be positive")
})

# ---------------------------------------------------------------------------
# Test 6: PP predictions valid (risk difference exists, values in [0,1])
# ---------------------------------------------------------------------------
test_that("V01-6: PP predictions valid (risk diff exists, values in [0,1])", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                       "nvarA", "nvarB", "nvarC"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      switch_n_cov = c("nvarA", "nvarB"),
                                      quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))
  obj <- suppressMessages(emulate_predict(obj, times = 0:8, type = "cum_inc",
                                       samples = 100, seed = 12345,
                                       difference = TRUE))

  pred <- obj$predictions

  # Risk difference should exist and be non-zero at the last time point
  rd_max <- pred$diff[nrow(pred)]
  expect_true(abs(rd_max) > 0.0001,
              info = sprintf("Risk difference at max follow-up = %.6f, expected non-zero",
                             rd_max))

  # All values in [0, 1] (allowing small MC noise)
  expect_true(all(pred$est_0 >= -0.01 & pred$est_0 <= 1.01),
              info = "Control arm values outside [0,1]")
  expect_true(all(pred$est_1 >= -0.01 & pred$est_1 <= 1.01),
              info = "Treated arm values outside [0,1]")
})

# ---------------------------------------------------------------------------
# Test 7: ITT vs PP directional consistency (both negative, PP magnitude >= ITT)
# ---------------------------------------------------------------------------
test_that("V01-7: ITT vs PP directional consistency", {
  skip_on_cran()

  d <- load_trial_example()

  # ITT pipeline
  obj_itt <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = c("catvarA", "catvarB",
                                                           "nvarA", "nvarB", "nvarC"),
                                           estimand = "ITT"))
  obj_itt <- suppressMessages(emulate_expand(obj_itt))
  obj_itt <- suppressMessages(emulate_weight(obj_itt))
  obj_itt <- suppressMessages(emulate_fit(obj_itt,
                                       outcome_cov = c("catvarA", "catvarB",
                                                       "nvarA", "nvarB", "nvarC")))

  # PP pipeline
  obj_pp <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                          treatment = "treatment",
                                          outcome = "outcome",
                                          eligible = "eligible",
                                          covariates = c("catvarA", "catvarB",
                                                          "nvarA", "nvarB", "nvarC"),
                                          estimand = "PP"))
  obj_pp <- suppressMessages(emulate_expand(obj_pp))
  obj_pp <- suppressMessages(emulate_weight(obj_pp,
                                         switch_d_cov = c("nvarA", "nvarB"),
                                         switch_n_cov = c("nvarA", "nvarB"),
                                         quiet = TRUE))
  obj_pp <- suppressMessages(emulate_fit(obj_pp,
                                      outcome_cov = c("catvarA", "catvarB",
                                                      "nvarA", "nvarB", "nvarC")))

  b_itt <- obj_itt$model$b_treat
  b_pp  <- obj_pp$model$b_treat

  # Both should be negative (treatment is protective in this dataset)
  expect_true(b_itt < 0,
              info = sprintf("ITT coefficient should be negative, got %.4f", b_itt))
  expect_true(b_pp < 0,
              info = sprintf("PP coefficient should be negative, got %.4f", b_pp))

  # PP should be at least as strong as ITT (treatment switching dilutes ITT)
  # Use 80% threshold to account for sampling variability
  expect_true(abs(b_pp) >= abs(b_itt) * 0.8,
              info = sprintf("PP magnitude (%.4f) should be >= 80%% of ITT magnitude (%.4f)",
                             abs(b_pp), abs(b_itt)))
})
