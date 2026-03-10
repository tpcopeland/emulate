# Cross-validation: emulate vs Stata tte reference values
# Mirrors: Stata-Tools/tte/qa/crossval_emulate_vs_r.do
#
# This file reverses the cross-validation direction: R validates against
# hardcoded Stata reference values from known results on trial_example.csv
# (503 patients, 48,400 person-periods).
#
# Stata reference values come from crossval_emulate_vs_r.do run with:
#   - tte version aligned with emulate v0.1.x
#   - Harrell RCS splines (same formula as emulate)
#   - vce(cluster id) with HC1 + cadjust
#
# Known algorithmic differences between R (emulate) and Stata (tte):
#   1. SE computation: emulate uses sandwich::vcovCL(type="HC1", cadjust=TRUE)
#      which applies HC1 + G/(G-1) correction. Stata uses vce(cluster) which
#      applies G/(G-1) only. This produces slightly different SEs.
#   2. Weight truncation percentiles: emulate uses quantile(type=2) matching
#      Stata's _pctile, so truncation cutoffs should be identical.
#   3. Spline basis: Both use Harrell RCS with identical formulas.
#      emulate's .emulate_rcs_basis matches Stata's _emulate_natural_spline exactly.
#
# R TrialEmulation reference (from 01_r_analysis.R output):
#   Config 1 (ITT): coef = -0.282945, SE = 0.313778
#   Config 2 (PP):  coef = -0.414264, SE = 0.415158
#   Config 3 (PP truncated): coef = -0.414264, SE = 0.415158
#     (Config 3 coefficients match Config 2 because R TrialEmulation's p99
#      truncation on this dataset has minimal effect.)
#
# Stata tte reference (from crossval_emulate_vs_r.do):
#   Config 1 (ITT): coef ~ -0.282, SE ~ 0.312
#   Config 2 (PP, stabilized): coef ~ -0.420
#   Config 3 (PP, truncated 1/99): coef ~ -0.420

# Hardcoded Stata reference values
# These are from known Stata tte runs on trial_example dataset
stata_c1_coef <- -0.282   # ITT quadratic, no weights
stata_c1_se   <- 0.312
stata_c2_coef <- -0.420   # PP quadratic, stabilized IPTW
stata_c3_coef <- -0.420   # PP quadratic, stabilized + truncated 1/99

# ---------------------------------------------------------------------------
# Test 1: Config 1 ITT coefficient within tolerance
# ---------------------------------------------------------------------------
test_that("Crossval-1: Config 1 ITT coefficient matches Stata (+/- 0.02)", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))

  b <- obj$model$b_treat

  # Tight tolerance: ITT has no weighting, so differences come only from
  # SE formula and numeric precision
  expect_true(abs(b - stata_c1_coef) < 0.02,
              info = sprintf("emulate ITT coef %.4f vs Stata %.4f, diff = %.4f (tol = 0.02)",
                             b, stata_c1_coef, abs(b - stata_c1_coef)))
})

# ---------------------------------------------------------------------------
# Test 2: Config 1 ITT SE within tolerance (wider due to HC1 vs G/(G-1))
# ---------------------------------------------------------------------------
test_that("Crossval-2: Config 1 ITT SE matches Stata (+/- 0.05)", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))

  se <- obj$model$se_treat

  # Wider tolerance for SE: HC1 + G/(G-1) in emulate vs G/(G-1) only in Stata
  # With ~500 clusters, the HC1 correction N/(N-p) is close to 1, so
  # differences should be small. Allow +/- 0.05 absolute.
  expect_true(abs(se - stata_c1_se) < 0.05,
              info = sprintf("emulate ITT SE %.4f vs Stata %.4f, diff = %.4f (tol = 0.05)",
                             se, stata_c1_se, abs(se - stata_c1_se)))
})

# ---------------------------------------------------------------------------
# Test 3: Config 2 PP coefficient within tolerance
# ---------------------------------------------------------------------------
test_that("Crossval-3: Config 2 PP coefficient matches Stata (+/- 0.15)", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      switch_n_cov = c("nvarA", "nvarB"),
                                      quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))

  b <- obj$model$b_treat

  # PP has wider tolerance: weight model stratification differs
  # (emulate uses 2 strata by arm, R TrialEmulation uses 4 strata by arm x lag_treat)
  expect_true(abs(b - stata_c2_coef) < 0.15,
              info = sprintf("emulate PP coef %.4f vs Stata %.4f, diff = %.4f (tol = 0.15)",
                             b, stata_c2_coef, abs(b - stata_c2_coef)))
})

# ---------------------------------------------------------------------------
# Test 4: Config 3 PP truncated coefficient within tolerance
# ---------------------------------------------------------------------------
test_that("Crossval-4: Config 3 PP truncated coefficient matches Stata (+/- 0.15)", {
  skip_on_cran()

  d <- load_trial_example()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      switch_n_cov = c("nvarA", "nvarB"),
                                      truncate = c(1, 99),
                                      quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))

  b <- obj$model$b_treat

  # Truncation uses quantile(type=2) matching Stata's _pctile, so
  # truncation cutoffs should be identical. Coefficient differences
  # come from weight model strata (same as Config 2).
  expect_true(abs(b - stata_c3_coef) < 0.15,
              info = sprintf("emulate PP-trunc coef %.4f vs Stata %.4f, diff = %.4f (tol = 0.15)",
                             b, stata_c3_coef, abs(b - stata_c3_coef)))
})

# ---------------------------------------------------------------------------
# Test 5: All configs — risk difference at t=10 in correct direction (negative)
# ---------------------------------------------------------------------------
test_that("Crossval-5: All configs have negative risk difference at t=10", {
  skip_on_cran()

  d <- load_trial_example()

  rds <- numeric(3)

  # Config 1: ITT
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))
  # Predict at times 0-10 to include t=10
  max_fu <- max(obj$data[["_emulate_followup"]], na.rm = TRUE)
  pred_times <- 0:min(10, max_fu)
  obj <- suppressMessages(emulate_predict(obj, times = pred_times,
                                       type = "cum_inc",
                                       samples = 50, seed = 12345,
                                       difference = TRUE))
  # Risk difference at the last available time point
  rds[1] <- obj$predictions$diff[nrow(obj$predictions)]

  # Config 2: PP, stabilized
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      switch_n_cov = c("nvarA", "nvarB"),
                                      quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))
  max_fu <- max(obj$data[["_emulate_followup"]], na.rm = TRUE)
  pred_times <- 0:min(10, max_fu)
  obj <- suppressMessages(emulate_predict(obj, times = pred_times,
                                       type = "cum_inc",
                                       samples = 50, seed = 12345,
                                       difference = TRUE))
  rds[2] <- obj$predictions$diff[nrow(obj$predictions)]

  # Config 3: PP, stabilized + truncated
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      switch_n_cov = c("nvarA", "nvarB"),
                                      truncate = c(1, 99),
                                      quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))
  max_fu <- max(obj$data[["_emulate_followup"]], na.rm = TRUE)
  pred_times <- 0:min(10, max_fu)
  obj <- suppressMessages(emulate_predict(obj, times = pred_times,
                                       type = "cum_inc",
                                       samples = 50, seed = 12345,
                                       difference = TRUE))
  rds[3] <- obj$predictions$diff[nrow(obj$predictions)]

  # Treatment is protective in trial_example: treated arm has LOWER
  # cumulative incidence, so risk_diff = CI_treated - CI_control should
  # be negative.
  for (i in seq_along(rds)) {
    config_label <- c("Config 1 (ITT)", "Config 2 (PP)", "Config 3 (PP-trunc)")[i]
    expect_true(rds[i] < 0,
                info = sprintf("%s: risk diff at last time = %.6f, expected < 0",
                               config_label, rds[i]))
  }
})
