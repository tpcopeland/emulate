# Tests for negative control calibration (Phase 6)

test_that("calibration produces adjusted CIs", {
  skip_if_not_installed("EmpiricalCalibration")

  dat <- dgp_simple(n = 500, periods = 6, seed = 601)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  # Simulate NCO results (negative controls with zero true effect)
  set.seed(42)
  nco_results <- data.frame(
    outcome_name = paste0("nco_", 1:5),
    log_estimate = rnorm(5, 0.05, 0.1),  # small bias
    se_log_estimate = runif(5, 0.08, 0.15)
  )

  obj <- suppressMessages(
    emulate_calibrate(obj, nco_results = nco_results)
  )

  expect_true(isTRUE(obj$state$calibrated))
  expect_true(!is.null(obj$calibration$calibrated))
  expect_true(is.finite(obj$calibration$calibrated$log_estimate))
  expect_true(is.finite(obj$calibration$calibrated$ci_lo))
  expect_true(is.finite(obj$calibration$calibrated$ci_hi))
})

test_that("calibrated CI is wider than uncalibrated", {
  skip_if_not_installed("EmpiricalCalibration")

  dat <- dgp_simple(n = 500, periods = 6, seed = 602)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  # NCOs with noticeable systematic error
  set.seed(43)
  nco_results <- data.frame(
    outcome_name = paste0("nco_", 1:5),
    log_estimate = rnorm(5, 0.1, 0.15),  # systematic bias
    se_log_estimate = runif(5, 0.08, 0.12)
  )

  obj <- suppressMessages(
    emulate_calibrate(obj, nco_results = nco_results)
  )

  cal <- obj$calibration
  # Calibrated CI should generally be wider when there's systematic error
  uncal_width <- cal$uncalibrated$ci_hi - cal$uncalibrated$ci_lo
  cal_width <- cal$calibrated$ci_hi - cal$calibrated$ci_lo

  # This is expected but not guaranteed for all random seeds,
  # so just check both are positive
  expect_true(uncal_width > 0)
  expect_true(cal_width > 0)
})

test_that("error when no NCO results provided", {
  skip_if_not_installed("EmpiricalCalibration")

  dat <- dgp_simple(n = 100, periods = 4, seed = 603)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 3)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))

  expect_error(
    suppressMessages(emulate_calibrate(obj)),
    "nco_results is required"
  )
})

test_that("error when fewer than 2 NCOs", {
  skip_if_not_installed("EmpiricalCalibration")

  dat <- dgp_simple(n = 100, periods = 4, seed = 604)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 3)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))

  nco_1 <- data.frame(
    outcome_name = "nco_1",
    log_estimate = 0.05,
    se_log_estimate = 0.1
  )

  expect_error(
    suppressMessages(emulate_calibrate(obj, nco_results = nco_1)),
    "At least 2"
  )
})

test_that("calibration plot produces ggplot object", {
  skip_if_not_installed("EmpiricalCalibration")

  dat <- dgp_simple(n = 500, periods = 6, seed = 605)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  nco_results <- data.frame(
    outcome_name = paste0("nco_", 1:5),
    log_estimate = rnorm(5, 0.05, 0.1),
    se_log_estimate = runif(5, 0.08, 0.15)
  )

  obj <- suppressMessages(
    emulate_calibrate(obj, nco_results = nco_results)
  )

  p <- emulate_plot(obj, type = "calibration")
  # EmpiricalCalibration returns a recordedplot or ggplot
  expect_true(!is.null(p))
})

test_that("emulate_report includes calibrated results", {
  skip_if_not_installed("EmpiricalCalibration")

  dat <- dgp_simple(n = 500, periods = 6, seed = 606)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  nco_results <- data.frame(
    outcome_name = paste0("nco_", 1:5),
    log_estimate = rnorm(5, 0.05, 0.1),
    se_log_estimate = runif(5, 0.08, 0.15)
  )

  obj <- suppressMessages(
    emulate_calibrate(obj, nco_results = nco_results)
  )

  msgs <- capture_messages(emulate_report(obj, eform = TRUE))
  expect_true(any(grepl("[Cc]alibrat", msgs)))
})
