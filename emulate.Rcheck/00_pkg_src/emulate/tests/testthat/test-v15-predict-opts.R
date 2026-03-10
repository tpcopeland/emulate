# V15: emulate_predict Options (7 tests)

# Helper: build a fitted ITT object for predict tests
.v15_fit <- function(seed = 1500) {
  d <- dgp_simple(n = 500, periods = 8, effect = -0.50, seed = seed,
                   outcome_intercept = -3)
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
  obj
}

test_that("V15.1: survival predictions in [0,1]", {
  obj <- .v15_fit(seed = 1501)
  obj <- suppressMessages(emulate_predict(obj, times = 0:5, type = "survival",
                                       samples = 20, seed = 1501))
  surv <- obj$predictions
  expect_true(all(surv$est_0 >= 0 & surv$est_0 <= 1))
  expect_true(all(surv$est_1 >= 0 & surv$est_1 <= 1))
})

test_that("V15.2: survival + cum_inc complementary", {
  obj <- .v15_fit(seed = 1502)
  obj_surv <- suppressMessages(emulate_predict(obj, times = 0:5, type = "survival",
                                            samples = 20, seed = 99))
  obj_ci <- suppressMessages(emulate_predict(obj, times = 0:5, type = "cum_inc",
                                          samples = 20, seed = 99))

  s0 <- obj_surv$predictions$est_0
  c0 <- obj_ci$predictions$est_0
  expect_true(all(abs(s0 + c0 - 1.0) < 0.01))

  s1 <- obj_surv$predictions$est_1
  c1 <- obj_ci$predictions$est_1
  expect_true(all(abs(s1 + c1 - 1.0) < 0.01))
})

test_that("V15.3: difference=TRUE stores diff column", {
  obj <- .v15_fit(seed = 1503)
  obj <- suppressMessages(emulate_predict(obj, times = 0:5, type = "cum_inc",
                                       samples = 20, seed = 1503,
                                       difference = TRUE))
  expect_true("diff" %in% names(obj$predictions))
  expect_true("diff_lo" %in% names(obj$predictions))
  expect_true("diff_hi" %in% names(obj$predictions))
})

test_that("V15.4: risk difference sign is correct", {
  # Use larger N for stable risk difference estimate
  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 1504,
                   outcome_intercept = -3)
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
  obj <- suppressMessages(emulate_predict(obj, times = c(3, 6), type = "cum_inc",
                                       samples = 30, seed = 1504,
                                       difference = TRUE))
  # With negative treatment effect, treated should have lower cum_inc
  # diff = treated - control, so diff should be < 0 at later times
  last_diff <- obj$predictions$diff[nrow(obj$predictions)]
  expect_true(last_diff < 0,
              info = sprintf("Risk difference %.4f should be negative", last_diff))
})

test_that("V15.5: same seed gives identical predictions", {
  obj <- .v15_fit(seed = 1505)
  obj1 <- suppressMessages(emulate_predict(obj, times = 0:4, samples = 20,
                                        seed = 777))
  obj2 <- suppressMessages(emulate_predict(obj, times = 0:4, samples = 20,
                                        seed = 777))
  expect_equal(obj1$predictions$est_0, obj2$predictions$est_0)
  expect_equal(obj1$predictions$est_1, obj2$predictions$est_1)
})

test_that("V15.6: level=90 narrower CIs than level=99", {
  obj <- .v15_fit(seed = 1506)
  obj_90 <- suppressMessages(emulate_predict(obj, times = c(3, 5), samples = 30,
                                          seed = 42, level = 90))
  obj_99 <- suppressMessages(emulate_predict(obj, times = c(3, 5), samples = 30,
                                          seed = 42, level = 99))
  # Width of CI for arm 0 at last time
  w90 <- obj_90$predictions$ci_hi_0 - obj_90$predictions$ci_lo_0
  w99 <- obj_99$predictions$ci_hi_0 - obj_99$predictions$ci_lo_0
  # 90% CI should be narrower than 99% CI
  expect_true(all(w90 <= w99 + 1e-10))
})

test_that("V15.7: samples=10 minimum runs", {
  obj <- .v15_fit(seed = 1507)
  obj <- suppressMessages(emulate_predict(obj, times = 0:3, samples = 10,
                                       seed = 1507))
  expect_true(!is.null(obj$predictions))
  expect_equal(nrow(obj$predictions), 4L)
})
