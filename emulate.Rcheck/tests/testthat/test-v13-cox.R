# V13: Cox Model Ground Truth (6 tests)

test_that("V13.1: Cox ITT pipeline completes", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 1301,
                   outcome_intercept = -3)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))
  expect_true(obj$state$fitted)
  expect_equal(obj$model$type, "cox")
})

test_that("V13.2: Cox ITT coefficient is negative", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 1302,
                   outcome_intercept = -3)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("Cox ITT coef %.4f should be < 0",
                             obj$model$b_treat))
})

test_that("V13.3: Cox ITT close to logistic ITT", {
  skip_on_cran()

  d <- dgp_simple(n = 5000, periods = 10, effect = -0.50, seed = 1303,
                   outcome_intercept = -3.5)

  # Logistic
  obj_log <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = "x", estimand = "ITT"))
  obj_log <- suppressMessages(emulate_expand(obj_log, maxfollowup = 8))
  obj_log <- suppressMessages(emulate_weight(obj_log))
  obj_log <- suppressMessages(emulate_fit(obj_log, outcome_cov = "x",
                                       followup_spec = "quadratic",
                                       trial_period_spec = "linear"))

  # Cox
  obj_cox <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = "x", estimand = "ITT"))
  obj_cox <- suppressMessages(emulate_expand(obj_cox, maxfollowup = 8))
  obj_cox <- suppressMessages(emulate_weight(obj_cox))
  obj_cox <- suppressMessages(emulate_fit(obj_cox, outcome_cov = "x", model = "cox",
                                       trial_period_spec = "linear"))

  expect_true(abs(obj_cox$model$b_treat - obj_log$model$b_treat) < 0.3,
              info = sprintf("Cox(%.4f) and logistic(%.4f) should be within 0.3",
                             obj_cox$model$b_treat, obj_log$model$b_treat))
})

test_that("V13.4: Cox PP pipeline completes", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 1304,
                   outcome_intercept = -3)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))
  expect_true(obj$state$fitted)
})

test_that("V13.5: Cox PP coefficient is negative", {
  skip_on_cran()

  d <- dgp_simple(n = 3000, periods = 10, effect = -0.50, seed = 1305,
                   outcome_intercept = -3)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("Cox PP coef %.4f should be < 0",
                             obj$model$b_treat))
})

test_that("V13.6: emulate_predict after Cox errors", {
  skip_on_cran()

  d <- dgp_simple(n = 500, periods = 8, seed = 1306, outcome_intercept = -3)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x", model = "cox",
                                   trial_period_spec = "linear"))
  expect_error(suppressMessages(emulate_predict(obj, times = 0:4)), "logistic")
})
