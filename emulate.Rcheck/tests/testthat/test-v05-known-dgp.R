# V5: Known DGP Recovery (6 tests)

test_that("V5.1: large-sample ITT recovers negative effect", {
  d <- dgp_simple(n = 10000, periods = 10, effect = -0.50, seed = 501)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("ITT coef %.4f should be < 0", obj$model$b_treat))
})

test_that("V5.2: large-sample PP recovers negative effect", {
  d <- dgp_simple(n = 10000, periods = 10, effect = -0.50, seed = 502)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("PP coef %.4f should be < 0", obj$model$b_treat))
})

test_that("V5.3: both ITT and PP are negative", {
  d <- dgp_simple(n = 5000, periods = 10, effect = -0.50, seed = 503)

  # ITT
  obj_itt <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = "x", estimand = "ITT"))
  obj_itt <- suppressMessages(emulate_expand(obj_itt, maxfollowup = 8))
  obj_itt <- suppressMessages(emulate_weight(obj_itt))
  obj_itt <- suppressMessages(emulate_fit(obj_itt, outcome_cov = "x",
                                       followup_spec = "quadratic",
                                       trial_period_spec = "linear"))

  # PP
  obj_pp <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                          treatment = "treatment",
                                          outcome = "outcome",
                                          eligible = "eligible",
                                          covariates = "x", estimand = "PP"))
  obj_pp <- suppressMessages(emulate_expand(obj_pp, maxfollowup = 8))
  obj_pp <- suppressMessages(emulate_weight(obj_pp, switch_d_cov = "x",
                                         truncate = c(1, 99), quiet = TRUE))
  obj_pp <- suppressMessages(emulate_fit(obj_pp, outcome_cov = "x",
                                      followup_spec = "quadratic",
                                      trial_period_spec = "linear"))

  expect_true(obj_itt$model$b_treat < 0)
  expect_true(obj_pp$model$b_treat < 0)
})

test_that("V5.4: Monte Carlo PP mean is negative", {
  skip_on_cran()

  n_reps <- 50
  coefs <- numeric(n_reps)
  for (r in seq_len(n_reps)) {
    d <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 5040 + r)
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
    coefs[r] <- obj$model$b_treat
  }
  expect_true(mean(coefs) < 0,
              info = sprintf("MC mean PP coef %.4f should be < 0", mean(coefs)))
})

test_that("V5.5: natural spline spec produces non-zero coefficient", {
  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 505)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "ns(3)",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat != 0)
  expect_true(obj$model$se_treat > 0)
})

test_that("V5.6: cubic spec completes with negative coefficient", {
  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 506)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "cubic",
                                   trial_period_spec = "linear"))
  expect_true(obj$state$fitted)
  expect_true(obj$model$b_treat < 0)
})
