# V4: G-Formula HIV/ART (5 tests)

test_that("V4.1: ITT with time-varying confounding shows negative effect", {
  skip_on_cran()

  d <- dgp_gformula(n = 5000, periods = 15, seed = 401)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("cd4_std", "age_cat", "male"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("cd4_std", "age_cat", "male"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("ITT coef %.4f should be < 0", obj$model$b_treat))
})

test_that("V4.2: PP with IPTW shows negative effect", {
  skip_on_cran()

  d <- dgp_gformula(n = 5000, periods = 15, seed = 402)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("cd4_std", "age_cat", "male"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("cd4_std", "age_cat", "male"),
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("cd4_std", "age_cat", "male"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("PP coef %.4f should be < 0", obj$model$b_treat))
})

test_that("V4.3: PP shows stronger or comparable effect to naive", {
  skip_on_cran()

  d <- dgp_gformula(n = 5000, periods = 15, seed = 403)

  # Naive unadjusted analysis
  naive_fit <- glm(outcome ~ treatment, data = d, family = binomial())
  naive_b <- coef(naive_fit)["treatment"]

  # PP pipeline
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("cd4_std", "age_cat", "male"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("cd4_std", "age_cat", "male"),
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("cd4_std", "age_cat", "male"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))

  # PP should show a stronger or comparable effect (more negative or similar)
  expect_true(obj$model$b_treat <= naive_b + 0.2,
              info = sprintf("PP coef %.4f should be <= naive %.4f + 0.2",
                             obj$model$b_treat, naive_b))
})

test_that("V4.4: weight diagnostics show ESS > 500", {
  skip_on_cran()

  d <- dgp_gformula(n = 5000, periods = 15, seed = 404)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("cd4_std", "age_cat", "male"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("cd4_std", "age_cat", "male"),
                                      truncate = c(1, 99), quiet = TRUE))
  expect_true(obj$weights$ess > 500,
              info = sprintf("ESS %.1f should be > 500", obj$weights$ess))
})

test_that("V4.5: predictions show treated < control", {
  skip_on_cran()

  d <- dgp_gformula(n = 5000, periods = 15, seed = 405)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("cd4_std", "age_cat", "male"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("cd4_std", "age_cat", "male"),
                                   followup_spec = "linear",
                                   trial_period_spec = "linear"))
  obj <- suppressMessages(emulate_predict(obj, times = c(4, 8),
                                       samples = 30, seed = 405))
  last <- nrow(obj$predictions)
  expect_true(obj$predictions$est_1[last] < obj$predictions$est_0[last])
})
