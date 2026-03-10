# V12: Sensitivity Sweeps and Stress Tests (5 tests)

test_that("V12.1: truncation sweep PP all negative", {
  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 1201)

  for (trunc in list(c(1, 99), c(5, 95), c(10, 90))) {
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment",
                                         outcome = "outcome",
                                         eligible = "eligible",
                                         covariates = "x", estimand = "PP"))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
    obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                        truncate = trunc, quiet = TRUE))
    obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                     followup_spec = "quadratic",
                                     trial_period_spec = "linear"))
    expect_true(obj$model$b_treat < 0,
                info = sprintf("Trunc(%d,%d) coef %.4f should be < 0",
                               trunc[1], trunc[2], obj$model$b_treat))
  }
})

test_that("V12.2: time spec sweep ITT all negative", {
  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 1202)

  for (spec in c("linear", "quadratic", "cubic", "ns(3)")) {
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment",
                                         outcome = "outcome",
                                         eligible = "eligible",
                                         covariates = "x", estimand = "ITT"))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
    obj <- suppressMessages(emulate_weight(obj))
    obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                     followup_spec = spec,
                                     trial_period_spec = "linear"))
    expect_true(obj$model$b_treat < 0,
                info = sprintf("followup_spec=%s coef %.4f should be < 0",
                               spec, obj$model$b_treat))
  }
})

test_that("V12.3: follow-up length sweep ITT all negative", {
  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 1203)

  for (mfu in c(4, 6, 8)) {
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment",
                                         outcome = "outcome",
                                         eligible = "eligible",
                                         covariates = "x", estimand = "ITT"))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = mfu))
    obj <- suppressMessages(emulate_weight(obj))
    obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                     followup_spec = "quadratic",
                                     trial_period_spec = "linear"))
    expect_true(obj$model$b_treat < 0,
                info = sprintf("maxfollowup=%d coef %.4f should be < 0",
                               mfu, obj$model$b_treat))
  }
})

test_that("V12.4: large N=50000 ITT stress test completes", {
  skip_on_cran()

  d <- dgp_simple(n = 50000, periods = 8, effect = -0.50, seed = 1204)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "linear",
                                   trial_period_spec = "linear"))
  expect_true(obj$state$fitted)
  expect_true(obj$model$b_treat < 0)
})

test_that("V12.5: large N=50000 PP stress test completes", {
  skip_on_cran()

  d <- dgp_simple(n = 50000, periods = 8, effect = -0.50, seed = 1205)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "linear",
                                   trial_period_spec = "linear"))
  expect_true(obj$state$fitted)
  expect_true(obj$model$b_treat < 0)
})
