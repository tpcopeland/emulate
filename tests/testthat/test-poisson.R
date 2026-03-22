# Tests for Poisson regression (Phase 2)

test_that("Poisson model fits and produces log-RR coefficient", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 201)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))

  obj <- suppressMessages(
    emulate_fit(obj, outcome_cov = "x", model = "poisson")
  )

  expect_true(obj$state$fitted)
  expect_equal(obj$model$type, "poisson")
  expect_true(is.finite(obj$model$b_treat))
  expect_true(is.finite(obj$model$se_treat))

  # Check it's a glm with poisson family
  expect_equal(family(obj$model$object)$family, "poisson")
})

test_that("Poisson RR differs from logistic OR", {
  dat <- dgp_simple(n = 1000, periods = 6, seed = 202,
                     outcome_intercept = -2)  # higher event rate

  obj_pois <- emulate_prepare(dat, id = "id", period = "period",
                               treatment = "treatment", outcome = "outcome",
                               eligible = "eligible", covariates = "x",
                               estimand = "ITT")
  obj_pois <- emulate_expand(obj_pois, maxfollowup = 4)
  obj_pois <- suppressMessages(emulate_weight(obj_pois))
  obj_pois <- suppressMessages(
    emulate_fit(obj_pois, outcome_cov = "x", model = "poisson")
  )

  obj_log <- emulate_prepare(dat, id = "id", period = "period",
                              treatment = "treatment", outcome = "outcome",
                              eligible = "eligible", covariates = "x",
                              estimand = "ITT")
  obj_log <- emulate_expand(obj_log, maxfollowup = 4)
  obj_log <- suppressMessages(emulate_weight(obj_log))
  obj_log <- suppressMessages(
    emulate_fit(obj_log, outcome_cov = "x", model = "logistic")
  )

  # RR and OR should differ (especially for common outcomes)
  rr <- exp(obj_pois$model$b_treat)
  or <- exp(obj_log$model$b_treat)
  expect_false(isTRUE(all.equal(rr, or, tolerance = 0.01)))
})

test_that("Poisson works with weights", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 203)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(
    emulate_weight(obj, switch_d_cov = "x")
  )

  obj <- suppressMessages(
    emulate_fit(obj, outcome_cov = "x", model = "poisson")
  )

  expect_true(obj$state$fitted)
  expect_equal(obj$model$type, "poisson")
})

test_that("emulate_report shows RR with eform for Poisson", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 204)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(
    emulate_fit(obj, outcome_cov = "x", model = "poisson")
  )

  # Should show "RR" column when eform = TRUE
  msgs <- capture_messages(emulate_report(obj, eform = TRUE))
  expect_true(any(grepl("poisson", msgs, ignore.case = TRUE)))
})

test_that("emulate_predict rejects Poisson with informative error", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 205)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 4)
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(
    emulate_fit(obj, outcome_cov = "x", model = "poisson")
  )

  expect_error(
    emulate_predict(obj, times = 1:3),
    "Poisson"
  )
})

test_that("invalid model type errors", {
  dat <- dgp_simple(n = 100, periods = 4, seed = 206)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", estimand = "ITT")
  obj <- emulate_expand(obj, maxfollowup = 3)
  obj <- suppressMessages(emulate_weight(obj))

  expect_error(
    suppressMessages(emulate_fit(obj, model = "bad")),
    "model must be"
  )
})
