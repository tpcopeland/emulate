# V3: CCW/Immortal-Time Bias (5 tests)

test_that("V3.1: naive logistic shows exaggerated protective effect", {
  skip_on_cran()

  d <- dgp_ccw(n = 2000, seed = 301)

  # Naive analysis: simple logistic on raw data (ignores immortal-time bias)
  naive_fit <- glm(outcome ~ treatment + age_std + ps + stage,
                   data = d, family = binomial())
  naive_b <- coef(naive_fit)["treatment"]

  # Should show an exaggerated protective effect (more negative than truth)
  # True HR ~ 0.60, so log(0.60) ~ -0.51
  expect_true(naive_b < 0,
              info = sprintf("Naive coef %.4f should be < 0", naive_b))
})

test_that("V3.2: PP/CCW pipeline produces negative coefficient", {
  skip_on_cran()

  d <- dgp_ccw(n = 2000, seed = 302)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("age_std", "ps", "stage"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 12))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("age_std", "ps", "stage"),
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("age_std", "ps", "stage"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("PP coef %.4f should be < 0", obj$model$b_treat))
})

test_that("V3.3: ITT pipeline produces negative coefficient", {
  skip_on_cran()

  d <- dgp_ccw(n = 2000, seed = 303)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("age_std", "ps", "stage"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 12))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("age_std", "ps", "stage"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  expect_true(obj$model$b_treat < 0,
              info = sprintf("ITT coef %.4f should be < 0", obj$model$b_treat))
})

test_that("V3.4: diagnose shows ESS > 100 and max SMD < 0.5", {
  skip_on_cran()

  d <- dgp_ccw(n = 2000, seed = 304)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("age_std", "ps", "stage"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 12))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("age_std", "ps", "stage"),
                                      truncate = c(1, 99), quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("age_std", "ps", "stage"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "linear"))
  obj <- suppressMessages(emulate_diagnose(obj,
                                        balance_covariates = c("age_std", "ps", "stage")))
  expect_true(obj$diagnostics$ess > 100)
  max_smd <- max(abs(obj$diagnostics$balance$smd_wt))
  expect_true(max_smd < 0.5,
              info = sprintf("Max weighted SMD %.4f should be < 0.5", max_smd))
})

test_that("V3.5: predictions show treated CI < control CI", {
  skip_on_cran()

  d <- dgp_ccw(n = 2000, seed = 305)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("age_std", "ps", "stage"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("age_std", "ps", "stage"),
                                   followup_spec = "linear",
                                   trial_period_spec = "linear"))
  obj <- suppressMessages(emulate_predict(obj, times = c(5, 10),
                                       samples = 30, seed = 305))
  # Treated (arm 1) should have lower cumulative incidence
  last <- nrow(obj$predictions)
  expect_true(obj$predictions$est_1[last] < obj$predictions$est_0[last],
              info = "Treated cum_inc should be < control at end of follow-up")
})
