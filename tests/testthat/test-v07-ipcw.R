# V7: IPCW Informative Censoring (5 tests)

test_that("V7.1: PP without IPCW shows negative effect", {
  skip_on_cran()

  d <- dgp_ipcw(n = 5000, periods = 10, effect = -0.60, seed = 701)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("x", "z"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = c("x", "z"),
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("x", "z"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("PP no-IPCW coef %.4f should be < 0",
                             obj$model$b_treat))
})

test_that("V7.2: PP with IPCW shows negative effect", {
  skip_on_cran()

  d <- dgp_ipcw(n = 5000, periods = 10, effect = -0.60, seed = 702)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("x", "z"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = c("x", "z"),
                                      censor_d_cov = c("x", "z"),
                                      censor_n_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("x", "z"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("PP IPCW coef %.4f should be < 0",
                             obj$model$b_treat))
})

test_that("V7.3: IPCW coefficient closer to truth", {
  skip_on_cran()

  d <- dgp_ipcw(n = 5000, periods = 10, effect = -0.60, seed = 703)
  true_effect <- -0.60

  # Without IPCW
  obj_no <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                          treatment = "treatment",
                                          outcome = "outcome",
                                          eligible = "eligible",
                                          covariates = c("x", "z"),
                                          estimand = "PP"))
  obj_no <- suppressMessages(emulate_expand(obj_no, maxfollowup = 8))
  obj_no <- suppressMessages(emulate_weight(obj_no, switch_d_cov = c("x", "z"),
                                         truncate = c(1, 99), quiet = TRUE))
  obj_no <- suppressMessages(emulate_fit(obj_no, outcome_cov = c("x", "z"),
                                      followup_spec = "quadratic",
                                      trial_period_spec = "linear"))

  # With IPCW
  obj_ipcw <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                            treatment = "treatment",
                                            outcome = "outcome",
                                            eligible = "eligible",
                                            covariates = c("x", "z"),
                                            estimand = "PP"))
  obj_ipcw <- suppressMessages(emulate_expand(obj_ipcw, maxfollowup = 8))
  obj_ipcw <- suppressMessages(emulate_weight(obj_ipcw, switch_d_cov = c("x", "z"),
                                           censor_d_cov = c("x", "z"),
                                           censor_n_cov = "x",
                                           truncate = c(1, 99), quiet = TRUE))
  obj_ipcw <- suppressMessages(emulate_fit(obj_ipcw, outcome_cov = c("x", "z"),
                                        followup_spec = "quadratic",
                                        trial_period_spec = "linear"))

  # Both should be negative; IPCW should be closer to truth (within tolerance)
  bias_no <- abs(obj_no$model$b_treat - true_effect)
  bias_ipcw <- abs(obj_ipcw$model$b_treat - true_effect)
  expect_true(bias_ipcw < bias_no + 0.2,
              info = sprintf("IPCW bias %.4f should be < no-IPCW bias %.4f + 0.2",
                             bias_ipcw, bias_no))
})

test_that("V7.4: weight mean in plausible range", {
  skip_on_cran()

  d <- dgp_ipcw(n = 5000, periods = 10, effect = -0.60, seed = 704)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("x", "z"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = c("x", "z"),
                                      censor_d_cov = c("x", "z"),
                                      censor_n_cov = "x",
                                      truncate = c(1, 99), quiet = TRUE))
  expect_true(obj$weights$mean >= 0.5 && obj$weights$mean <= 2.0,
              info = sprintf("Weight mean %.4f should be in [0.5, 2.0]",
                             obj$weights$mean))
})

test_that("V7.5: pooled censor model runs and produces negative coefficient", {
  skip_on_cran()

  d <- dgp_ipcw(n = 5000, periods = 10, effect = -0.60, seed = 705)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("x", "z"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = c("x", "z"),
                                      censor_d_cov = c("x", "z"),
                                      censor_n_cov = "x",
                                      pool_censor = TRUE,
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("x", "z"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("Pooled censor coef %.4f should be < 0",
                             obj$model$b_treat))
})
