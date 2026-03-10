# V11: RCT vs Observational Benchmark (3 tests)

test_that("V11.1: RCT ITT shows negative effect", {
  skip_on_cran()

  d <- dgp_rct(n = 5000, periods = 10, effect = -0.50, seed = 1101)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("RCT ITT coef %.4f should be < 0",
                             obj$model$b_treat))
})

test_that("V11.2: observational PP approximates RCT ITT", {
  skip_on_cran()

  # RCT data
  d_rct <- dgp_rct(n = 5000, periods = 10, effect = -0.50, seed = 1102)
  obj_rct <- suppressMessages(emulate_prepare(d_rct, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = "x", estimand = "ITT"))
  obj_rct <- suppressMessages(emulate_expand(obj_rct, maxfollowup = 8))
  obj_rct <- suppressMessages(emulate_weight(obj_rct))
  obj_rct <- suppressMessages(emulate_fit(obj_rct, outcome_cov = "x",
                                       followup_spec = "quadratic",
                                       trial_period_spec = "linear"))

  # Observational data with confounding
  d_obs <- dgp_obs(n = 5000, periods = 10, effect = -0.50, seed = 1102)
  obj_obs <- suppressMessages(emulate_prepare(d_obs, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = "x", estimand = "PP"))
  obj_obs <- suppressMessages(emulate_expand(obj_obs, maxfollowup = 8))
  obj_obs <- suppressMessages(emulate_weight(obj_obs, switch_d_cov = "x",
                                          truncate = c(1, 99), quiet = TRUE))
  obj_obs <- suppressMessages(emulate_fit(obj_obs, outcome_cov = "x",
                                       followup_spec = "quadratic",
                                       trial_period_spec = "linear"))

  # Same direction
  expect_true(sign(obj_rct$model$b_treat) == sign(obj_obs$model$b_treat),
              info = "RCT and observational PP should have same direction")
  # Within 0.5
  expect_true(abs(obj_rct$model$b_treat - obj_obs$model$b_treat) < 0.5,
              info = sprintf("RCT(%.4f) and Obs PP(%.4f) should be within 0.5",
                             obj_rct$model$b_treat, obj_obs$model$b_treat))
})

test_that("V11.3: observational ITT attenuated relative to PP", {
  skip_on_cran()

  d <- dgp_obs(n = 5000, periods = 10, effect = -0.50, seed = 1103)

  # ITT
  obj_itt <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = "x", estimand = "ITT"))
  obj_itt <- suppressMessages(emulate_expand(obj_itt, maxfollowup = 8))
  obj_itt <- suppressMessages(emulate_weight(obj_itt))
  obj_itt <- suppressMessages(emulate_fit(obj_itt, outcome_cov = "x",
                                       followup_spec = "quadratic",
                                       trial_period_spec = "linear"))

  # PP
  obj_pp <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                          treatment = "treatment",
                                          outcome = "outcome",
                                          eligible = "eligible",
                                          covariates = "x", estimand = "PP"))
  obj_pp <- suppressMessages(emulate_expand(obj_pp, maxfollowup = 8))
  obj_pp <- suppressMessages(emulate_weight(obj_pp, switch_d_cov = "x",
                                         truncate = c(1, 99), quiet = TRUE))
  obj_pp <- suppressMessages(emulate_fit(obj_pp, outcome_cov = "x",
                                      followup_spec = "quadratic",
                                      trial_period_spec = "linear"))

  # ITT should be attenuated: |ITT| <= |PP| + 0.2
  expect_true(abs(obj_itt$model$b_treat) <= abs(obj_pp$model$b_treat) + 0.2,
              info = sprintf("|ITT|(%.4f) should be <= |PP|(%.4f) + 0.2",
                             abs(obj_itt$model$b_treat),
                             abs(obj_pp$model$b_treat)))
})
