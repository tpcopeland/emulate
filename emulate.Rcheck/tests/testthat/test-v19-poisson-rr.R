# Validation v19: Poisson regression risk ratios

test_that("v19: Poisson RR converges toward truth", {
  # DGP with known effect (protective treatment, log-OR = -0.50)
  dat <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 1901)

  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 5)
  obj <- suppressMessages(emulate_weight(obj))

  obj <- suppressMessages(
    emulate_fit(obj, outcome_cov = "x", model = "poisson")
  )

  # Should be negative (protective)
  expect_true(obj$model$b_treat < 0)

  # RR should be < 1
  rr <- exp(obj$model$b_treat)
  expect_true(rr < 1)
  expect_true(rr > 0.2)  # not implausibly small
})

test_that("v19: Poisson RR vs logistic OR - RR attenuated for common outcomes", {
  # Use higher intercept for more common outcomes where OR != RR
  dat <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = 1902,
                     outcome_intercept = -2)

  obj_pois <- emulate_prepare(dat, id = "id", period = "period",
                               treatment = "treatment", outcome = "outcome",
                               eligible = "eligible", covariates = "x",
                               estimand = "ITT")
  obj_pois <- emulate_expand(obj_pois, maxfollowup = 5)
  obj_pois <- suppressMessages(emulate_weight(obj_pois))
  obj_pois <- suppressMessages(
    emulate_fit(obj_pois, outcome_cov = "x", model = "poisson")
  )

  obj_log <- emulate_prepare(dat, id = "id", period = "period",
                              treatment = "treatment", outcome = "outcome",
                              eligible = "eligible", covariates = "x",
                              estimand = "ITT")
  obj_log <- emulate_expand(obj_log, maxfollowup = 5)
  obj_log <- suppressMessages(emulate_weight(obj_log))
  obj_log <- suppressMessages(
    emulate_fit(obj_log, outcome_cov = "x", model = "logistic")
  )

  rr <- exp(obj_pois$model$b_treat)
  or <- exp(obj_log$model$b_treat)

  # For common outcomes with protective effect:
  # |log(RR)| < |log(OR)| (RR is attenuated toward 1 compared to OR)
  expect_true(abs(obj_pois$model$b_treat) < abs(obj_log$model$b_treat))
})

test_that("v19: Poisson with clustered SEs", {
  dat <- dgp_simple(n = 1000, periods = 6, seed = 1903)

  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))

  obj <- suppressMessages(
    emulate_fit(obj, outcome_cov = "x", model = "poisson")
  )

  # Cluster-robust variance should exist
  expect_true(!is.null(obj$model$vcov))
  # SE should be positive
  expect_true(obj$model$se_treat > 0)
})
