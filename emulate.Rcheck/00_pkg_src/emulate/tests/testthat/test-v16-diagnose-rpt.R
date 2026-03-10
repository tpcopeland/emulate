# V16: emulate_diagnose and emulate_report (8 tests)

# Helper: build a fitted PP object for diagnose/report tests
.v16_pp_fit <- function(seed = 1600) {
  d <- dgp_simple(n = 500, periods = 8, effect = -0.50, seed = seed,
                   outcome_intercept = -3)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressWarnings(suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                                    followup_spec = "linear",
                                                    trial_period_spec = "linear")))
  obj
}

test_that("V16.1: emulate_diagnose weight distribution is valid", {
  obj <- .v16_pp_fit(seed = 1601)
  obj <- suppressMessages(emulate_diagnose(obj))
  expect_true(obj$diagnostics$ess > 0)
  expect_true(obj$diagnostics$weight_mean >= 0.5 &&
                obj$diagnostics$weight_mean <= 2.0)
  expect_true(obj$diagnostics$weight_sd > 0)
})

test_that("V16.2: emulate_diagnose balance covariates produces SMD", {
  obj <- .v16_pp_fit(seed = 1602)
  obj <- suppressMessages(emulate_diagnose(obj, balance_covariates = "x"))
  expect_true(!is.null(obj$diagnostics$balance))
  expect_true(all(!is.na(obj$diagnostics$balance$smd_wt)))
  expect_true(all(!is.na(obj$diagnostics$balance$smd_unwt)))
})

test_that("V16.3: balance data.frame has correct structure", {
  obj <- .v16_pp_fit(seed = 1603)
  obj <- suppressMessages(emulate_diagnose(obj, balance_covariates = "x"))
  bal <- obj$diagnostics$balance
  expect_true(nrow(bal) >= 1)
  expect_true(ncol(bal) >= 2)
  expect_true("covariate" %in% names(bal))
  expect_true("smd_wt" %in% names(bal))
})

test_that("V16.4: emulate_diagnose by_trial completes", {
  obj <- .v16_pp_fit(seed = 1604)
  expect_no_error(suppressMessages(emulate_diagnose(obj, by_trial = TRUE)))
})

test_that("V16.5: emulate_diagnose on ITT completes", {
  d <- dgp_simple(n = 500, periods = 8, seed = 1605, outcome_intercept = -3)
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
  expect_no_error(suppressMessages(emulate_diagnose(obj)))
})

test_that("V16.6: emulate_report after fit runs without error", {
  obj <- .v16_pp_fit(seed = 1606)
  expect_no_error(suppressMessages(emulate_report(obj)))
})

test_that("V16.7: emulate_report with eform completes", {
  obj <- .v16_pp_fit(seed = 1607)
  expect_no_error(suppressMessages(emulate_report(obj, eform = TRUE)))
})

test_that("V16.8: emulate_report CSV export creates file", {
  obj <- .v16_pp_fit(seed = 1608)
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  suppressMessages(emulate_report(obj, format = "csv", export = tmp))
  expect_true(file.exists(tmp))
  contents <- read.csv(tmp)
  expect_true(nrow(contents) > 0)
})
