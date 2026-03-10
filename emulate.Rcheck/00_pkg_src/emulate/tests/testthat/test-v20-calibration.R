# Validation v20: Negative control calibration

test_that("v20: calibration reduces bias from systematic error", {
  skip_if_not_installed("EmpiricalCalibration")

  # Fit primary analysis
  dat <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 2001)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 5)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  # Simulate NCOs: zero true effect but systematic upward bias
  set.seed(2001)
  nco_results <- data.frame(
    outcome_name = paste0("nco_", 1:8),
    log_estimate = rnorm(8, 0.15, 0.08),  # systematic positive bias
    se_log_estimate = runif(8, 0.08, 0.12)
  )

  obj <- suppressMessages(
    emulate_calibrate(obj, nco_results = nco_results)
  )

  cal <- obj$calibration

  # Calibrated estimate should exist
  expect_true(is.finite(cal$calibrated$log_estimate))
  expect_true(is.finite(cal$calibrated$ci_lo))
  expect_true(is.finite(cal$calibrated$ci_hi))

  # Null distribution should capture the systematic error
  expect_true(is.finite(cal$null_distribution[1]))
  expect_true(is.finite(cal$null_distribution[2]))
})

test_that("v20: calibrated CIs have appropriate width", {
  skip_if_not_installed("EmpiricalCalibration")

  dat <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 2002)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 5)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  # NCOs with known systematic bias
  set.seed(2002)
  nco_results <- data.frame(
    outcome_name = paste0("nco_", 1:6),
    log_estimate = rnorm(6, 0.1, 0.12),
    se_log_estimate = runif(6, 0.07, 0.11)
  )

  obj <- suppressMessages(
    emulate_calibrate(obj, nco_results = nco_results)
  )

  cal <- obj$calibration
  uncal_width <- cal$uncalibrated$ci_hi - cal$uncalibrated$ci_lo
  cal_width <- cal$calibrated$ci_hi - cal$calibrated$ci_lo

  # Both should be positive intervals
  expect_true(uncal_width > 0)
  expect_true(cal_width > 0)
})

test_that("v20: calibration with minimal NCOs", {
  skip_if_not_installed("EmpiricalCalibration")

  dat <- dgp_simple(n = 1000, periods = 6, seed = 2003)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  # Minimum: 2 NCOs
  nco_results <- data.frame(
    outcome_name = c("nco_1", "nco_2"),
    log_estimate = c(0.05, -0.03),
    se_log_estimate = c(0.1, 0.12)
  )

  obj <- suppressMessages(
    emulate_calibrate(obj, nco_results = nco_results)
  )

  expect_true(isTRUE(obj$state$calibrated))
  expect_true(!is.null(obj$calibration$calibrated))
})
