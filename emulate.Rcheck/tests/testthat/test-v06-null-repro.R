# V6: Null Effect and Reproducibility (5 tests)

test_that("V6.1: PP 95% CI covers 0 under null effect", {
  d <- dgp_null(n = 5000, periods = 10, seed = 601)
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

  b <- obj$model$b_treat
  se <- obj$model$se_treat
  ci_lo <- b - 1.96 * se
  ci_hi <- b + 1.96 * se
  expect_true(ci_lo <= 0 && ci_hi >= 0,
              info = sprintf("95%% CI [%.4f, %.4f] should cover 0", ci_lo, ci_hi))
})

test_that("V6.2: ITT 95% CI covers 0 under null effect", {
  d <- dgp_null(n = 5000, periods = 10, seed = 602)
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

  b <- obj$model$b_treat
  se <- obj$model$se_treat
  ci_lo <- b - 1.96 * se
  ci_hi <- b + 1.96 * se
  expect_true(ci_lo <= 0 && ci_hi >= 0,
              info = sprintf("95%% CI [%.4f, %.4f] should cover 0", ci_lo, ci_hi))
})

test_that("V6.3: type-I error rate under null is controlled", {
  skip_on_cran()

  n_reps <- 100
  rejections <- 0L
  for (r in seq_len(n_reps)) {
    d <- dgp_null(n = 1000, periods = 8, seed = 6030 + r)
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
    z <- obj$model$b_treat / obj$model$se_treat
    if (abs(z) > 1.96) rejections <- rejections + 1L
  }
  expect_true(rejections <= 15,
              info = sprintf("Rejection rate %d/100 should be <= 15", rejections))
})

test_that("V6.4: same seed produces identical coefficients", {
  run_one <- function(seed) {
    d <- dgp_simple(n = 1000, periods = 8, effect = -0.50, seed = seed)
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
    obj$model$b_treat
  }

  b1 <- run_one(604)
  b2 <- run_one(604)
  expect_equal(b1, b2)
})

test_that("V6.5: different seed produces different coefficients", {
  run_one <- function(seed) {
    d <- dgp_simple(n = 1000, periods = 8, effect = -0.50, seed = seed)
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
    obj$model$b_treat
  }

  b1 <- run_one(605)
  b2 <- run_one(606)
  expect_false(b1 == b2)
})
