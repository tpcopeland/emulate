# Tests for propensity score matching (Phase 3)

test_that("matching produces fewer rows than original", {
  skip_if_not_installed("MatchIt")

  dat <- dgp_simple(n = 500, periods = 6, seed = 301)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)
  n_expanded <- nrow(obj$data)

  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", match_ratio = 3L, caliper = 0.5)
  )

  expect_true(isTRUE(obj$state$matched))
  expect_true(obj$state$fitted)
  expect_true(nrow(obj$data) <= n_expanded)
  expect_true(obj$matching$n_matched <= n_expanded)
})

test_that("MatchIt object is stored", {
  skip_if_not_installed("MatchIt")

  dat <- dgp_simple(n = 300, periods = 5, seed = 302)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 3)

  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", match_ratio = 2L, caliper = 0.5)
  )

  expect_true(!is.null(obj$matching$matchit_object))
  expect_s3_class(obj$matching$matchit_object, "matchit")
  expect_true(!is.null(obj$matching$ps_values))
})

test_that("emulate_predict works on matched object", {
  skip_if_not_installed("MatchIt")

  dat <- dgp_simple(n = 500, periods = 6, seed = 303)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", outcome_cov = "x",
                  match_ratio = 3L, caliper = 0.5)
  )

  obj <- suppressMessages(
    emulate_predict(obj, times = 1:3, samples = 20, seed = 42)
  )

  expect_true(!is.null(obj$predictions))
  expect_equal(nrow(obj$predictions), 3)
})

test_that("emulate_report works on matched object", {
  skip_if_not_installed("MatchIt")

  dat <- dgp_simple(n = 500, periods = 6, seed = 304)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", outcome_cov = "x",
                  match_ratio = 3L, caliper = 0.5)
  )

  msgs <- capture_messages(emulate_report(obj))
  expect_true(any(grepl("logistic", msgs, ignore.case = TRUE)))
})

test_that("matching with lasso PS works", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("glmnet")

  dat <- dgp_simple(n = 500, periods = 6, seed = 305)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", ps_method = "lasso",
                  match_ratio = 3L, caliper = 0.5)
  )

  expect_true(isTRUE(obj$state$matched))
  expect_equal(obj$matching$ps_method, "lasso")
})

test_that("error when no match covariates provided", {
  dat <- dgp_simple(n = 100, periods = 4, seed = 306)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 3)

  expect_error(
    suppressMessages(emulate_match(obj)),
    "match_cov is required"
  )
})

test_that("adjustment_method is set to matching", {
  skip_if_not_installed("MatchIt")

  dat <- dgp_simple(n = 300, periods = 5, seed = 307)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 3)

  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", match_ratio = 2L, caliper = 0.5)
  )

  expect_equal(obj$model$adjustment_method, "matching")
})
