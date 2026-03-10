# Validation v18: Propensity score matching recovers treatment effect

test_that("v18: matching recovers ATT within tolerance", {
  skip_if_not_installed("MatchIt")

  # DGP with known effect
  dat <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 1801)

  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 5)

  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", outcome_cov = "x",
                  match_ratio = 5L, caliper = 0.5)
  )

  # Treatment effect should be negative (protective)
  expect_true(obj$model$b_treat < 0)

  # Should be in the right ballpark (log-OR around -0.50)
  expect_true(obj$model$b_treat > -2.0)
  expect_true(obj$model$b_treat < 0.5)
})

test_that("v18: matching vs IPTW on same data", {
  skip_if_not_installed("MatchIt")

  dat <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 1802)

  # IPTW path
  obj_iptw <- emulate_prepare(dat, id = "id", period = "period",
                               treatment = "treatment", outcome = "outcome",
                               eligible = "eligible", covariates = "x",
                               estimand = "PP")
  obj_iptw <- emulate_expand(obj_iptw, maxfollowup = 5)
  obj_iptw <- suppressMessages(
    emulate_weight(obj_iptw, switch_d_cov = "x")
  )
  obj_iptw <- suppressMessages(
    emulate_fit(obj_iptw, outcome_cov = "x")
  )

  # Matching path
  obj_match <- emulate_prepare(dat, id = "id", period = "period",
                                treatment = "treatment", outcome = "outcome",
                                eligible = "eligible", covariates = "x",
                                estimand = "PP")
  obj_match <- emulate_expand(obj_match, maxfollowup = 5)
  obj_match <- suppressMessages(
    emulate_match(obj_match, match_cov = "x", outcome_cov = "x",
                  match_ratio = 5L, caliper = 0.5)
  )

  # Both should estimate negative treatment effects
  expect_true(obj_iptw$model$b_treat < 0)
  expect_true(obj_match$model$b_treat < 0)

  # They should be in the same direction but may differ in magnitude
  # (different estimands: ATE vs ATT, different data)
  expect_true(sign(obj_iptw$model$b_treat) == sign(obj_match$model$b_treat))
})

test_that("v18: matching preserves pipeline integrity", {
  skip_if_not_installed("MatchIt")

  dat <- dgp_simple(n = 500, periods = 6, seed = 1803)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_match(obj, match_cov = "x", outcome_cov = "x",
                  match_ratio = 3L, caliper = 0.5)
  )

  # State flags

  expect_true(obj$state$prepared)
  expect_true(obj$state$expanded)
  expect_true(isTRUE(obj$state$matched))
  expect_true(obj$state$fitted)

  # Model metadata
  expect_equal(obj$model$adjustment_method, "matching")

  # Report should work
  msgs <- capture_messages(emulate_report(obj))
  expect_true(length(msgs) > 0)
})
