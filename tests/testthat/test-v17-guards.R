# V17: Pipeline Guards (6 tests)

test_that("V17.1: emulate_fit before expand errors", {
  d <- dgp_simple(n = 100, periods = 6, seed = 1701)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  expect_error(suppressMessages(emulate_fit(obj)), "expanded")
})

test_that("V17.2: emulate_predict before fit errors", {
  d <- dgp_simple(n = 100, periods = 6, seed = 1702)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  expect_error(suppressMessages(emulate_predict(obj, times = 0:3)),
               "fitted|expanded")
})

test_that("V17.3: emulate_diagnose before expand errors", {
  d <- dgp_simple(n = 100, periods = 6, seed = 1703)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  expect_error(suppressMessages(emulate_diagnose(obj)), "expanded")
})

test_that("V17.4: emulate_weight PP without switch_d_cov errors", {
  d <- dgp_simple(n = 100, periods = 6, seed = 1704)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 4))
  expect_error(suppressMessages(emulate_weight(obj)), "switch_d_cov")
})

test_that("V17.5: emulate_predict after Cox errors", {
  d <- dgp_simple(n = 500, periods = 8, seed = 1705, outcome_intercept = -3)
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

test_that("V17.6: emulate_expand with negative grace errors", {
  d <- dgp_simple(n = 100, periods = 6, seed = 1706)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  expect_error(suppressMessages(emulate_expand(obj, grace = -1)), "non-negative")
})
