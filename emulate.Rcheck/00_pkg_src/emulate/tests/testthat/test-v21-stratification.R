# Validation v21: Propensity score stratification

test_that("v21: stratification recovers treatment effect", {
  dat <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 2101)

  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 5)

  obj <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x", n_strata = 5L, outcome_cov = "x")
  )

  # Treatment effect should be negative (protective)
  expect_true(obj$model$b_treat < 0)

  # Should be in the right ballpark
  expect_true(obj$model$b_treat > -2.0)
  expect_true(obj$model$b_treat < 0.5)
})

test_that("v21: stratification vs IPTW on same data", {
  dat <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 2102)

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

  # Stratification path
  obj_strat <- emulate_prepare(dat, id = "id", period = "period",
                                treatment = "treatment", outcome = "outcome",
                                eligible = "eligible", covariates = "x",
                                estimand = "PP")
  obj_strat <- emulate_expand(obj_strat, maxfollowup = 5)
  obj_strat <- suppressMessages(
    emulate_stratify(obj_strat, strat_cov = "x", n_strata = 5L,
                     outcome_cov = "x")
  )

  # Both should estimate negative treatment effects
  expect_true(obj_iptw$model$b_treat < 0)
  expect_true(obj_strat$model$b_treat < 0)

  # Same direction
  expect_true(sign(obj_iptw$model$b_treat) == sign(obj_strat$model$b_treat))
})

test_that("v21: stratification with different strata counts", {
  dat <- dgp_simple(n = 1000, periods = 6, seed = 2103)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  # 3 strata
  obj3 <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x", n_strata = 3L, outcome_cov = "x")
  )

  # 10 strata
  obj10 <- suppressMessages(
    emulate_stratify(obj, strat_cov = "x", n_strata = 10L, outcome_cov = "x")
  )

  # Both should produce valid estimates
  expect_true(is.finite(obj3$model$b_treat))
  expect_true(is.finite(obj10$model$b_treat))

  # More strata = finer adjustment, but both should be in same direction
  expect_true(sign(obj3$model$b_treat) == sign(obj10$model$b_treat))
})
