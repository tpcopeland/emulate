# Tests for propensity score stratification (Phase 4)

test_that("stratification creates correct number of strata", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 401)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x", n_strata = 5L)
  )

  expect_true(isTRUE(obj$state$stratified))
  expect_true(obj$state$fitted)
  expect_equal(obj$stratification$n_strata, 5)
})

test_that("stratum ID added to data", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 402)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x")
  )

  expect_true("ps_stratum" %in% names(obj$data))
  expect_s3_class(obj$data$ps_stratum, "factor")
})

test_that("fitted model includes stratum terms for logistic", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 403)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x", model = "logistic")
  )

  coef_names <- names(coef(obj$model$object))
  # Should have ps_stratum terms
  expect_true(any(grepl("ps_stratum", coef_names)))
})

test_that("Cox model uses strata() correctly", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 404)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x", model = "cox", n_strata = 3L)
  )

  expect_true(obj$state$fitted)
  expect_equal(obj$model$type, "cox")
  # Cox strata() doesn't show as coefficients
  expect_true(is.finite(obj$model$b_treat))
})

test_that("emulate_report works with stratification", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 405)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x")
  )

  msgs <- capture_messages(emulate_report(obj))
  expect_true(any(grepl("logistic", msgs, ignore.case = TRUE)))
})

test_that("adjustment_method is set to stratification", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 406)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x")
  )

  expect_equal(obj$model$adjustment_method, "stratification")
})

test_that("error when no strat_cov provided", {
  dat <- dgp_simple(n = 100, periods = 4, seed = 407)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 3)

  expect_error(
    suppressMessages(emulate_stratify(obj)),
    "strat_cov is required"
  )
})
