# V14: emulate_expand Options (4 tests)

test_that("V14.1: selective trials produces correct count", {
  d <- dgp_simple(n = 200, periods = 10, seed = 1401)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, trials = c(0, 2, 4, 6, 8)))
  expect_equal(obj$expansion$n_trials, 5L)
})

test_that("V14.2: single trial produces n_trials == 1", {
  d <- dgp_simple(n = 200, periods = 10, seed = 1402)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, trials = c(0)))
  expect_equal(obj$expansion$n_trials, 1L)
})

test_that("V14.3: selective trials coefficient same direction as full", {
  d <- dgp_simple(n = 500, periods = 10, effect = -0.50, seed = 1403)
  obj_full <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                            treatment = "treatment",
                                            outcome = "outcome",
                                            eligible = "eligible",
                                            covariates = "x", estimand = "ITT"))
  obj_full <- suppressMessages(emulate_expand(obj_full, maxfollowup = 8))
  obj_full <- suppressMessages(emulate_weight(obj_full))
  obj_full <- suppressMessages(emulate_fit(obj_full, outcome_cov = "x"))

  obj_sel <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = "x", estimand = "ITT"))
  obj_sel <- suppressMessages(emulate_expand(obj_sel, trials = c(0, 2, 4, 6, 8),
                                          maxfollowup = 8))
  obj_sel <- suppressMessages(emulate_weight(obj_sel))
  obj_sel <- suppressMessages(emulate_fit(obj_sel, outcome_cov = "x"))

  # Same direction
  expect_true(sign(obj_full$model$b_treat) == sign(obj_sel$model$b_treat))
})

test_that("V14.4: shorter maxfollowup produces fewer rows", {
  d <- dgp_simple(n = 300, periods = 10, seed = 1404)
  obj_long <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                            treatment = "treatment",
                                            outcome = "outcome",
                                            eligible = "eligible",
                                            covariates = "x", estimand = "ITT"))
  obj_long <- suppressMessages(emulate_expand(obj_long, maxfollowup = 0))

  obj_short <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                             treatment = "treatment",
                                             outcome = "outcome",
                                             eligible = "eligible",
                                             covariates = "x", estimand = "ITT"))
  obj_short <- suppressMessages(emulate_expand(obj_short, maxfollowup = 3))

  expect_true(obj_short$expansion$n_expanded < obj_long$expansion$n_expanded)
})
