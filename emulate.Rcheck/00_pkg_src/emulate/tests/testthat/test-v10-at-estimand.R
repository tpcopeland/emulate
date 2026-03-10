# V10: As-Treated Estimand (6 tests)

test_that("V10.1: AT pipeline completes", {
  skip_on_cran()

  d <- dgp_at(n = 5000, periods = 10, effect = -0.50, seed = 1001)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "AT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$state$fitted)
})

test_that("V10.2: AT coefficient is negative and bounded", {
  skip_on_cran()

  d <- dgp_at(n = 5000, periods = 10, effect = -0.50, seed = 1002)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "AT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("AT coef %.4f should be < 0", obj$model$b_treat))
  expect_true(abs(obj$model$b_treat) < 3,
              info = sprintf("|AT coef| %.4f should be < 3",
                             abs(obj$model$b_treat)))
})

test_that("V10.3: AT weights are valid", {
  skip_on_cran()

  d <- dgp_at(n = 5000, periods = 10, effect = -0.50, seed = 1003)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "AT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  expect_true(obj$weights$mean >= 0.1 && obj$weights$mean <= 10,
              info = sprintf("AT weight mean %.4f should be in [0.1, 10]",
                             obj$weights$mean))
  # No NAs in weights
  w <- obj$data[[paste0(obj$settings$prefix, "weight")]]
  expect_true(all(!is.na(w)))
})

test_that("V10.4: AT approximates PP", {
  skip_on_cran()

  d <- dgp_simple(n = 5000, periods = 10, effect = -0.50, seed = 1004)

  run_estimand <- function(est) {
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment",
                                         outcome = "outcome",
                                         eligible = "eligible",
                                         covariates = "x", estimand = est))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
    obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                        truncate = c(1, 99), quiet = TRUE))
    obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                     followup_spec = "quadratic",
                                     trial_period_spec = "linear"))
    obj$model$b_treat
  }

  b_at <- run_estimand("AT")
  b_pp <- run_estimand("PP")

  expect_true(abs(b_at - b_pp) < 0.5,
              info = sprintf("|AT(%.4f) - PP(%.4f)| should be < 0.5",
                             b_at, b_pp))
})

test_that("V10.5: AT with pool_switch runs with negative coefficient", {
  skip_on_cran()

  d <- dgp_at(n = 5000, periods = 10, effect = -0.50, seed = 1005)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "AT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      pool_switch = TRUE,
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0)
})

test_that("V10.6: AT predictions have cumulative incidence in [0,1]", {
  skip_on_cran()

  d <- dgp_at(n = 5000, periods = 10, effect = -0.50, seed = 1006)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "AT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "linear",
                                   trial_period_spec = "linear"))
  obj <- suppressMessages(emulate_predict(obj, times = c(2, 4, 6),
                                       samples = 20, seed = 1006))
  preds <- obj$predictions
  expect_true(all(preds$est_0 >= 0 & preds$est_0 <= 1))
  expect_true(all(preds$est_1 >= 0 & preds$est_1 <= 1))
})
