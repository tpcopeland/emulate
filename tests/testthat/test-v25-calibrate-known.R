# V25: Calibrate Known-Answer Correctness
#
# Validates emulate_calibrate() against hand-computed values and known
# statistical properties of the empirical calibration procedure.

# ---------------------------------------------------------------------------
# Test 1: Zero-bias NCOs produce near-zero bias estimate
# ---------------------------------------------------------------------------
test_that("V25-1: zero-bias NCOs yield bias near zero", {
  skip_on_cran()
  skip_if_not_installed("EmpiricalCalibration")
  # Install with: install.packages("EmpiricalCalibration")

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2501)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  # NCOs with true effect = 0 (centered around zero)
  set.seed(2501)
  nco <- data.frame(
    outcome_name = paste0("nco_", 1:10),
    log_estimate = rnorm(10, mean = 0, sd = 0.05),
    se_log_estimate = runif(10, 0.03, 0.08)
  )

  obj <- suppressMessages(emulate_calibrate(obj, nco_results = nco))

  # Bias should be near zero (within 0.15 of zero)
  expect_true(abs(obj$calibration$null_distribution[1]) < 0.15,
    label = "bias estimate near zero for unbiased NCOs")
})

# ---------------------------------------------------------------------------
# Test 2: Calibrated SE is always >= uncalibrated SE
# ---------------------------------------------------------------------------
test_that("V25-2: calibrated SE >= uncalibrated SE", {
  skip_on_cran()
  skip_if_not_installed("EmpiricalCalibration")
  # Install with: install.packages("EmpiricalCalibration")

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2502)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  nco <- data.frame(
    outcome_name = paste0("nco_", 1:5),
    log_estimate = c(0.10, -0.05, 0.15, -0.08, 0.03),
    se_log_estimate = c(0.05, 0.06, 0.04, 0.07, 0.05)
  )

  obj <- suppressMessages(emulate_calibrate(obj, nco_results = nco))

  cal <- obj$calibration
  expect_gte(cal$calibrated$se_log_estimate, cal$uncalibrated$se_log_estimate)
})

# ---------------------------------------------------------------------------
# Test 3: Calibrated CI is wider than uncalibrated CI
# ---------------------------------------------------------------------------
test_that("V25-3: calibrated CI wider than uncalibrated CI", {
  skip_on_cran()
  skip_if_not_installed("EmpiricalCalibration")
  # Install with: install.packages("EmpiricalCalibration")

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2503)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  nco <- data.frame(
    outcome_name = paste0("nco_", 1:6),
    log_estimate = c(0.20, -0.15, 0.10, -0.25, 0.05, 0.18),
    se_log_estimate = c(0.05, 0.06, 0.04, 0.07, 0.05, 0.06)
  )

  obj <- suppressMessages(emulate_calibrate(obj, nco_results = nco))

  cal <- obj$calibration
  uncal_width <- cal$uncalibrated$ci_hi - cal$uncalibrated$ci_lo
  cal_width <- cal$calibrated$ci_hi - cal$calibrated$ci_lo
  expect_gt(cal_width, uncal_width)
})

# ---------------------------------------------------------------------------
# Test 4: Known systematic bias is corrected
# ---------------------------------------------------------------------------
test_that("V25-4: known systematic bias is corrected in calibrated estimate", {
  skip_on_cran()
  skip_if_not_installed("EmpiricalCalibration")
  # Install with: install.packages("EmpiricalCalibration")

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2504)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  # NCOs with systematic positive bias (all shifted +0.3)
  nco <- data.frame(
    outcome_name = paste0("nco_", 1:8),
    log_estimate = rnorm(8, mean = 0.30, sd = 0.05),
    se_log_estimate = rep(0.05, 8)
  )
  set.seed(2504)

  obj <- suppressMessages(emulate_calibrate(obj, nco_results = nco))

  cal <- obj$calibration
  # Calibrated estimate should be shifted toward null relative to uncalibrated
  # (bias ~ +0.3, so calibrated = uncalibrated - 0.3ish)
  expect_true(cal$calibrated$log_estimate < cal$uncalibrated$log_estimate,
    label = "calibration corrects positive bias")
})

# ---------------------------------------------------------------------------
# Test 5: Level parameter affects CI width
# ---------------------------------------------------------------------------
test_that("V25-5: level(90) produces narrower CI than level(95)", {
  skip_on_cran()
  skip_if_not_installed("EmpiricalCalibration")
  # Install with: install.packages("EmpiricalCalibration")

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2505)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  nco <- data.frame(
    outcome_name = paste0("nco_", 1:5),
    log_estimate = c(0.05, -0.03, 0.08, -0.02, 0.04),
    se_log_estimate = rep(0.05, 5)
  )

  obj95 <- suppressMessages(emulate_calibrate(obj, nco_results = nco, level = 95))
  obj90 <- suppressMessages(emulate_calibrate(obj, nco_results = nco, level = 90))

  w95 <- obj95$calibration$calibrated$ci_hi - obj95$calibration$calibrated$ci_lo
  w90 <- obj90$calibration$calibrated$ci_hi - obj90$calibration$calibrated$ci_lo

  expect_gt(w95, w90)
})

# ---------------------------------------------------------------------------
# Test 6: Error on invalid NCO inputs
# ---------------------------------------------------------------------------
test_that("V25-6: error on fewer than 2 NCOs", {
  skip_on_cran()
  skip_if_not_installed("EmpiricalCalibration")
  # Install with: install.packages("EmpiricalCalibration")

  d <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 2506)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  nco_bad <- data.frame(
    outcome_name = "nco_1",
    log_estimate = 0.1,
    se_log_estimate = 0.05
  )

  expect_error(emulate_calibrate(obj, nco_results = nco_bad))
})

# ---------------------------------------------------------------------------
# Test 7: All return structure fields present
# ---------------------------------------------------------------------------
test_that("V25-7: calibration return structure complete", {
  skip_on_cran()
  skip_if_not_installed("EmpiricalCalibration")
  # Install with: install.packages("EmpiricalCalibration")

  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 2507)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  nco <- data.frame(
    outcome_name = paste0("nco_", 1:5),
    log_estimate = c(0.05, -0.03, 0.08, -0.02, 0.04),
    se_log_estimate = rep(0.05, 5)
  )

  obj <- suppressMessages(emulate_calibrate(obj, nco_results = nco))

  expect_true(obj$state$calibrated)
  expect_true(is.list(obj$calibration))
  expect_true(all(c("null_distribution", "uncalibrated", "calibrated", "level")
    %in% names(obj$calibration)))
  expect_true(all(c("log_estimate", "se_log_estimate", "ci_lo", "ci_hi")
    %in% names(obj$calibration$uncalibrated)))
  expect_true(all(c("log_estimate", "se_log_estimate", "ci_lo", "ci_hi")
    %in% names(obj$calibration$calibrated)))
})
