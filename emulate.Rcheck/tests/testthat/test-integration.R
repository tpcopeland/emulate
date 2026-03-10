test_that("ITT pipeline on trial_example matches R TrialEmulation benchmark", {
  skip_on_cran()

  # Load the trial_example dataset
  f <- system.file("extdata", "trial_example.csv", package = "emulate")
  if (f == "" || !file.exists(f)) {
    f <- file.path(normalizePath("../../inst/extdata", mustWork = FALSE),
                   "trial_example.csv")
  }
  skip_if(!file.exists(f), "trial_example.csv not found")

  d <- read.csv(f)

  # Config 1: ITT, quadratic time, no weights
  # R TrialEmulation: outcome_cov = ~ catvarA + catvarB + nvarA + nvarB + nvarC
  # All eligible trial periods, unlimited follow-up
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                      "nvarA", "nvarB", "nvarC"),
                                       estimand = "ITT"))

  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))

  # R TrialEmulation benchmark: treatment coef = -0.2829, SE = 0.3138
  b <- obj$model$b_treat
  se <- obj$model$se_treat

  # Within 5% relative tolerance of benchmark
  expect_true(abs(b - (-0.2829)) / abs(-0.2829) < 0.05,
              info = sprintf("ITT coef %.4f not within 5%% of -0.2829", b))

  # SE check (wider tolerance due to sandwich estimator differences)
  expect_true(abs(se - 0.3138) / 0.3138 < 0.15,
              info = sprintf("ITT SE %.4f not within 15%% of 0.3138", se))
})

test_that("PP pipeline on trial_example matches benchmark", {
  skip_on_cran()

  f <- system.file("extdata", "trial_example.csv", package = "emulate")
  if (f == "" || !file.exists(f)) {
    f <- file.path(normalizePath("../../inst/extdata", mustWork = FALSE),
                   "trial_example.csv")
  }
  skip_if(!file.exists(f), "trial_example.csv not found")

  d <- read.csv(f)

  # Config 2: PP, quadratic time, stabilized IPTW
  # R TrialEmulation: switch_n_cov = ~ nvarA + nvarB, switch_d_cov = ~ nvarA + nvarB
  # outcome_cov = ~ catvarA + catvarB + nvarA + nvarB + nvarC
  # All eligible trial periods, unlimited follow-up, no censor weights
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("catvarA", "catvarB",
                                                      "nvarA", "nvarB", "nvarC"),
                                       estimand = "PP"))

  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      switch_n_cov = c("nvarA", "nvarB"),
                                      quiet = TRUE))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("catvarA", "catvarB",
                                                   "nvarA", "nvarB", "nvarC")))

  # R TrialEmulation benchmark: treatment coef = -0.4143, SE = 0.4152
  b <- obj$model$b_treat
  se <- obj$model$se_treat

  expect_true(abs(b - (-0.4143)) / abs(-0.4143) < 0.10,
              info = sprintf("PP coef %.4f not within 10%% of -0.4143", b))
})
