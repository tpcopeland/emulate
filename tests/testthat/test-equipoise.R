# Tests for equipoise diagnostics (Phase 5)

test_that("preference scores are in [0, 1]", {
  # Test the internal function directly
  ps <- runif(100, 0.1, 0.9)
  pref <- emulate:::.compute_preference_score(ps, 0.3)

  expect_true(all(pref >= 0 & pref <= 1))
})

test_that("equipoise percentage is sensible", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 501)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(
    emulate_weight(obj, switch_d_cov = "x", ps_trim = c(0.01, 0.99))
  )

  obj <- suppressMessages(
    emulate_diagnose(obj, equipoise = TRUE, equipoise_bounds = c(0.3, 0.7))
  )

  eq <- obj$diagnostics$equipoise
  expect_true(!is.null(eq))
  expect_true(eq$pct_equipoise >= 0 && eq$pct_equipoise <= 100)
  expect_true(all(eq$preference_scores >= 0 & eq$preference_scores <= 1))
})

test_that("PS plot type produces ggplot object", {
  skip_if_not_installed("MatchIt")

  dat <- dgp_simple(n = 300, periods = 5, seed = 502)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 3)
  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", match_ratio = 2L, caliper = 0.5)
  )

  p <- emulate_plot(obj, type = "ps")
  expect_s3_class(p, "ggplot")
})

test_that("equipoise plot type produces ggplot object", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 503)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(
    emulate_weight(obj, switch_d_cov = "x", ps_trim = c(0.01, 0.99))
  )
  obj <- suppressMessages(
    emulate_diagnose(obj, equipoise = TRUE)
  )

  p <- emulate_plot(obj, type = "equipoise")
  expect_s3_class(p, "ggplot")
})

test_that("equipoise works with stratified objects", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 504)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x")
  )

  # Need weight var for diagnose, add dummy
  prefix <- obj$settings$prefix
  weight_col <- paste0(prefix, "weight")
  if (!weight_col %in% names(obj$data)) {
    obj$data[, (weight_col) := 1.0]
    obj$state$weighted <- TRUE
    obj$weights$weight_var <- weight_col
  }

  obj <- suppressMessages(
    emulate_diagnose(obj, equipoise = TRUE)
  )

  expect_true(!is.null(obj$diagnostics$equipoise))
})
