# V9: Edge Cases and Strict Validation (8 tests)

test_that("V9.1: small N ITT completes", {
  d <- dgp_simple(n = 50, periods = 8, effect = -0.50, seed = 901,
                   outcome_intercept = -2.5)
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
  expect_true(is.finite(obj$model$b_treat))
})

test_that("V9.2: sparse events ITT completes", {
  d <- dgp_simple(n = 500, periods = 10, effect = -0.50, seed = 902,
                   outcome_intercept = -6)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  # Use linear spec to reduce parameters for sparse data
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "linear",
                                   trial_period_spec = "none"))
  expect_true(obj$state$fitted)
  expect_true(is.finite(obj$model$b_treat))
})

test_that("V9.3: single eligible period produces exactly 1 trial", {
  d <- dgp_simple(n = 200, periods = 8, effect = -0.50, seed = 903)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, trials = c(0)))
  expect_equal(obj$expansion$n_trials, 1L)
})

test_that("V9.4: all binary covariates PP has valid weights", {
  d <- dgp_simple(n = 300, periods = 8, effect = -0.50, seed = 904)
  # Convert x to binary
  d$x <- as.integer(d$x > 0)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  expect_true(obj$weights$ess > 10)
  expect_true(obj$weights$mean > 0)
})

test_that("V9.5: strict validation catches period gaps", {
  d <- dgp_simple(n = 100, periods = 6, seed = 905)
  # Introduce a gap: remove period 2 for first 10 individuals
  d <- d[!(d$id <= 10 & d$period == 2), ]
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  expect_error(suppressMessages(emulate_validate(obj, strict = TRUE)),
               "validation failed")
})

test_that("V9.6: strict validation catches post-outcome rows", {
  d <- dgp_simple(n = 100, periods = 6, seed = 906, outcome_intercept = -2)
  # Find an individual with an event and add rows after it
  events <- d[d$outcome == 1, ]
  if (nrow(events) > 0) {
    victim <- events$id[1]
    event_period <- events$period[events$id == victim][1]
    extra <- data.frame(id = victim, period = event_period + 1L,
                        treatment = 1L, outcome = 0L, eligible = 1L,
                        x = d$x[d$id == victim][1])
    d <- rbind(d, extra)
  }
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  expect_error(suppressMessages(emulate_validate(obj, strict = TRUE)),
               "validation failed")
})

test_that("V9.7: strict validation catches missing treatment", {
  d <- dgp_simple(n = 100, periods = 6, seed = 907)
  # Introduce NA in treatment
  d$treatment[5] <- NA
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  expect_error(suppressMessages(emulate_validate(obj, strict = TRUE)),
               "validation failed")
})

test_that("V9.8: non-strict validation with issues gives warnings not errors", {
  d <- dgp_simple(n = 100, periods = 6, seed = 908)
  # Introduce a gap
  d <- d[!(d$id <= 5 & d$period == 2), ]
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  # Non-strict should not error
  expect_no_error(suppressMessages(emulate_validate(obj, strict = FALSE)))
  expect_true(obj$state$prepared)
})
