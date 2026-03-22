# V26: Risk Ratio and Risk Difference Hand-Computed Validation
#
# Validates emulate_predict() risk differences and risk ratios against
# manually computed values from the marginal survival curves.

# ---------------------------------------------------------------------------
# Test 1: Risk difference = treated CI minus control CI
# ---------------------------------------------------------------------------
test_that("V26-1: risk difference equals est_1 - est_0 at each time", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2601)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
  obj <- suppressMessages(emulate_predict(obj, times = c(0, 2, 4, 6),
    difference = TRUE, samples = 50, seed = 2601))

  preds <- obj$predictions

  # diff should be est_1 - est_0 at every row
  hand_diff <- preds$est_1 - preds$est_0
  expect_equal(preds$diff, hand_diff, tolerance = 1e-10,
    label = "diff column = est_1 - est_0")
})

# ---------------------------------------------------------------------------
# Test 2: Cumulative incidence is monotonically non-decreasing
# ---------------------------------------------------------------------------
test_that("V26-2: cumulative incidence is monotonically non-decreasing", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 2602)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
  obj <- suppressMessages(emulate_predict(obj, times = 0:8, type = "cum_inc",
    samples = 50, seed = 2602))

  preds <- obj$predictions

  # Cumulative incidence must be non-decreasing for both arms
  for (arm in c("est_0", "est_1")) {
    vals <- preds[[arm]]
    diffs <- diff(vals)
    expect_true(all(diffs >= -1e-10),
      label = paste("cum_inc non-decreasing for", arm))
  }
})

# ---------------------------------------------------------------------------
# Test 3: Survival = 1 - cumulative incidence
# ---------------------------------------------------------------------------
test_that("V26-3: survival + cumulative incidence = 1 at each time", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2603)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  obj_ci <- suppressMessages(emulate_predict(obj, times = c(0, 2, 4, 6),
    type = "cum_inc", samples = 50, seed = 2603))
  obj_surv <- suppressMessages(emulate_predict(obj, times = c(0, 2, 4, 6),
    type = "survival", samples = 50, seed = 2603))

  for (arm in c("est_0", "est_1")) {
    total <- obj_ci$predictions[[arm]] + obj_surv$predictions[[arm]]
    expect_equal(total, rep(1, nrow(obj_ci$predictions)), tolerance = 1e-10,
      label = paste("CI + Survival = 1 for", arm))
  }
})

# ---------------------------------------------------------------------------
# Test 4: CI at time 0 is near zero (very few events at baseline)
# ---------------------------------------------------------------------------
test_that("V26-4: cumulative incidence at time 0 is near zero", {
  skip_on_cran()

  d <- dgp_simple(n = 5000, periods = 10, effect = -0.50, seed = 2604,
                   outcome_intercept = -4.0)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
  obj <- suppressMessages(emulate_predict(obj, times = c(0, 4, 8),
    type = "cum_inc", samples = 50, seed = 2604))

  preds <- obj$predictions

  # At time 0, CI should be very small (just P(Y=1|t=0))
  expect_lt(preds$est_0[1], 0.10, label = "CI_0(t=0) small")
  expect_lt(preds$est_1[1], 0.10, label = "CI_1(t=0) small")
})

# ---------------------------------------------------------------------------
# Test 5: Protective treatment produces negative risk difference
# ---------------------------------------------------------------------------
test_that("V26-5: protective treatment yields negative risk difference", {
  skip_on_cran()

  # Strong protective effect
  d <- dgp_simple(n = 5000, periods = 10, effect = -1.0, seed = 2605)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
  obj <- suppressMessages(emulate_predict(obj, times = c(0, 4, 8),
    difference = TRUE, samples = 80, seed = 2605))

  preds <- obj$predictions

  # At later times, treated should have lower CI (protective)
  # Risk diff = treated - control should be negative
  expect_lt(preds$diff[preds$time == 8], 0,
    label = "risk difference negative for protective treatment at t=8")
})

# ---------------------------------------------------------------------------
# Test 6: MC CI bounds bracket point estimate
# ---------------------------------------------------------------------------
test_that("V26-6: MC confidence intervals bracket point estimates", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2606)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
  obj <- suppressMessages(emulate_predict(obj, times = c(0, 2, 4, 6),
    difference = TRUE, samples = 200, seed = 2606))

  preds <- obj$predictions

  # Point estimates should be within CIs (or very close)
  for (i in seq_len(nrow(preds))) {
    expect_true(preds$est_0[i] >= preds$ci_lo_0[i] - 1e-6 &&
                preds$est_0[i] <= preds$ci_hi_0[i] + 1e-6,
      label = paste("est_0 within CI at time", preds$time[i]))
    expect_true(preds$est_1[i] >= preds$ci_lo_1[i] - 1e-6 &&
                preds$est_1[i] <= preds$ci_hi_1[i] + 1e-6,
      label = paste("est_1 within CI at time", preds$time[i]))
  }
})
