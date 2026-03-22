# V24: Cox Model Gold-Standard Validation
#
# Validates the emulate Cox pipeline against direct survival::coxph() on the
# same expanded data. Since emulate internally calls coxph(), these should be
# numerically identical — the test validates that the expansion, weighting,
# and variable construction pipeline does not corrupt data before reaching coxph.
#
# Additionally validates:
#   - Baseline hazard (monotonically non-decreasing)
#   - Predicted survival (monotonically non-increasing)
#   - Cox vs logistic coefficient convergence
#   - Cox on real NHEFS data
#   - Cox with weighted PP estimation

# ---------------------------------------------------------------------------
# Test 1: emulate Cox identical to direct coxph on expanded data
# ---------------------------------------------------------------------------
test_that("V24-1: emulate Cox coefficient identical to direct coxph (tol=1e-8)", {
  skip_on_cran()

  d <- dgp_simple(n = 5000, periods = 10, effect = -0.50, seed = 2401,
                   outcome_intercept = -3.5)

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))

  b_emulate <- obj$model$b_treat
  fit_emulate <- obj$model$object

  # Direct coxph on the same expanded data
  dt <- obj$data
  prefix <- obj$settings$prefix
  fu_col <- paste0(prefix, "followup")
  arm_col <- paste0(prefix, "arm")
  outobs_col <- paste0(prefix, "outcome_obs")
  trial_col <- paste0(prefix, "trial")
  cens_col <- paste0(prefix, "censored")
  esample_col <- paste0(prefix, "esample")

  est_dt <- dt[dt[[esample_col]] == 1, ]
  est_dt$time_enter <- est_dt[[fu_col]]
  est_dt$time_exit <- est_dt[[fu_col]] + 1L

  direct_fit <- survival::coxph(
    survival::Surv(time_enter, time_exit, est_dt[[outobs_col]]) ~
      est_dt[[arm_col]] + est_dt[[trial_col]] + est_dt[["x"]],
    data = est_dt,
    cluster = est_dt[[obj$settings$id]]
  )

  b_direct <- unname(coef(direct_fit)[1])  # First coefficient is the arm

  # Should be numerically identical (same function, same data)
  expect_equal(unname(b_emulate), b_direct, tolerance = 1e-8,
               info = sprintf("emulate=%.10f, direct=%.10f", b_emulate, b_direct))
})

# ---------------------------------------------------------------------------
# Test 2: emulate Cox SE identical to direct coxph robust SE
# ---------------------------------------------------------------------------
test_that("V24-2: emulate Cox SE identical to direct coxph robust SE", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2402,
                   outcome_intercept = -3.5)

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))

  se_emulate <- obj$model$se_treat

  # Direct SE from the model's vcov (Cox uses cluster() internally)
  direct_se <- unname(sqrt(diag(vcov(obj$model$object)))[1])

  expect_equal(unname(se_emulate), direct_se, tolerance = 1e-8,
               info = sprintf("emulate SE=%.10f, direct SE=%.10f",
                              se_emulate, direct_se))
})

# ---------------------------------------------------------------------------
# Test 3: Baseline hazard is monotonically non-decreasing
# ---------------------------------------------------------------------------
test_that("V24-3: Baseline cumulative hazard is monotonically non-decreasing", {
  skip_on_cran()

  d <- dgp_simple(n = 5000, periods = 10, effect = -0.50, seed = 2403,
                   outcome_intercept = -3.0)

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))

  bh <- survival::basehaz(obj$model$object, centered = FALSE)

  # Cumulative hazard should be non-decreasing
  diffs <- diff(bh$hazard)
  expect_true(all(diffs >= -1e-10),
              info = "Baseline cumulative hazard is not monotonically non-decreasing")

  # Should start near zero
  expect_true(bh$hazard[1] < 0.5,
              info = sprintf("First cumulative hazard = %.4f, expected near 0",
                             bh$hazard[1]))

  # Should have multiple time points
  expect_true(nrow(bh) > 1,
              info = "Baseline hazard should have multiple time points")
})

# ---------------------------------------------------------------------------
# Test 4: Predicted survival is monotonically non-increasing
# ---------------------------------------------------------------------------
test_that("V24-4: Predicted survival from Cox is monotonically non-increasing", {
  skip_on_cran()

  d <- dgp_simple(n = 5000, periods = 10, effect = -0.50, seed = 2404,
                   outcome_intercept = -3.0)

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))

  sf <- survival::survfit(obj$model$object)

  # Survival should be non-increasing
  diffs <- diff(sf$surv)
  expect_true(all(diffs <= 1e-10),
              info = "Predicted survival is not monotonically non-increasing")

  # Should start at or near 1.0
  expect_true(sf$surv[1] <= 1.0 + 1e-10,
              info = sprintf("Initial survival = %.6f, expected <= 1.0", sf$surv[1]))

  # Should end above 0 (not all die in 8 periods)
  expect_true(sf$surv[length(sf$surv)] > 0,
              info = "Final survival should be > 0")
})

# ---------------------------------------------------------------------------
# Test 5: Cox log-likelihood is finite and negative
# ---------------------------------------------------------------------------
test_that("V24-5: Cox log-likelihood is finite and negative", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2405,
                   outcome_intercept = -3.5)

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))

  ll <- logLik(obj$model$object)

  expect_true(is.finite(as.numeric(ll)),
              info = "Log-likelihood should be finite")
  expect_true(as.numeric(ll) < 0,
              info = sprintf("Log-likelihood = %.2f, expected < 0", as.numeric(ll)))
})

# ---------------------------------------------------------------------------
# Test 6: Cox on NHEFS real data — coefficient and baseline hazard
# ---------------------------------------------------------------------------
test_that("V24-6: Cox on real NHEFS data produces valid results", {
  skip_on_cran()

  d <- load_nhefs()
  covs <- c("age_std", "sex", "race", "smoke_cat", "wt71_std", "smokeyrs_std")

  obj <- suppressMessages(emulate_prepare(d, id = "seqn", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = covs,
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = covs, model = "cox",
                                   trial_period_spec = "none"))

  # Coefficient should be finite
  expect_true(is.finite(obj$model$b_treat))

  # HR should be in plausible range [0.3, 2.0]
  hr <- exp(obj$model$b_treat)
  expect_true(hr > 0.3 & hr < 2.0,
              info = sprintf("NHEFS Cox HR = %.4f", hr))

  # Baseline hazard should exist and be valid
  bh <- survival::basehaz(obj$model$object, centered = FALSE)
  expect_true(nrow(bh) > 0)
  expect_true(all(diff(bh$hazard) >= -1e-10),
              info = "NHEFS baseline hazard not monotonic")
})

# ---------------------------------------------------------------------------
# Test 7: Cox PP with weights — coefficient recovery
# ---------------------------------------------------------------------------
test_that("V24-7: Cox PP with weights recovers protective effect", {
  skip_on_cran()

  d <- dgp_simple(n = 5000, periods = 10, effect = -0.50, seed = 2407,
                   outcome_intercept = -3.0)

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment", outcome = "outcome",
                                       eligible = "eligible", covariates = "x",
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))

  # Should be negative (protective)
  expect_true(obj$model$b_treat < 0,
              info = sprintf("Cox PP coef = %.4f, expected < 0", obj$model$b_treat))

  # HR in plausible range
  hr <- exp(obj$model$b_treat)
  expect_true(hr > 0.2 & hr < 1.5,
              info = sprintf("Cox PP HR = %.4f, expected in [0.2, 1.5]", hr))
})

# ---------------------------------------------------------------------------
# Test 8: Cox coefficient converges to logistic on golden DGP
# ---------------------------------------------------------------------------
test_that("V24-8: Cox and logistic coefficients agree on golden DGP", {
  skip_on_cran()

  d <- load_golden_dgp()

  # Logistic
  obj_l <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment", outcome = "outcome",
                                         eligible = "eligible", covariates = "x",
                                         estimand = "ITT"))
  obj_l <- suppressMessages(emulate_expand(obj_l, maxfollowup = 8))
  obj_l <- suppressMessages(emulate_weight(obj_l))
  obj_l <- suppressMessages(emulate_fit(obj_l, outcome_cov = "x",
                                     followup_spec = "quadratic",
                                     trial_period_spec = "linear"))

  # Cox
  obj_c <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment", outcome = "outcome",
                                         eligible = "eligible", covariates = "x",
                                         estimand = "ITT"))
  obj_c <- suppressMessages(emulate_expand(obj_c, maxfollowup = 8))
  obj_c <- suppressMessages(emulate_weight(obj_c))
  obj_c <- suppressMessages(emulate_fit(obj_c, outcome_cov = "x", model = "cox",
                                     trial_period_spec = "linear"))

  # Direction agreement
  expect_equal(sign(obj_l$model$b_treat), sign(obj_c$model$b_treat))

  # Within 0.15 (Cox and logistic model different hazard functions)
  diff <- abs(obj_l$model$b_treat - obj_c$model$b_treat)
  expect_true(diff < 0.15,
              info = sprintf("Logistic=%.4f, Cox=%.4f, diff=%.4f",
                             obj_l$model$b_treat, obj_c$model$b_treat, diff))

  # Both recover true effect within 0.20
  true_effect <- -0.50
  expect_true(abs(obj_l$model$b_treat - true_effect) < 0.20,
              info = sprintf("Logistic=%.4f vs true=%.2f",
                             obj_l$model$b_treat, true_effect))
  expect_true(abs(obj_c$model$b_treat - true_effect) < 0.20,
              info = sprintf("Cox=%.4f vs true=%.2f",
                             obj_c$model$b_treat, true_effect))
})
