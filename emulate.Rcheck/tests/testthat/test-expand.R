test_that("emulate_expand creates expanded dataset for ITT", {
  d <- make_test_data(n_ids = 30, n_periods = 8)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))

  expect_true(obj$state$expanded)
  expect_true(nrow(obj$data) > nrow(d))
  expect_true("_emulate_followup" %in% names(obj$data))
  expect_true("_emulate_trial" %in% names(obj$data))
  expect_true("_emulate_arm" %in% names(obj$data))
  expect_true("_emulate_censored" %in% names(obj$data))
  expect_true("_emulate_outcome_obs" %in% names(obj$data))

  # ITT: no artificial censoring
  expect_equal(sum(obj$data[["_emulate_censored"]]), 0)
})

test_that("emulate_expand clones for PP estimand", {
  d <- make_test_data(n_ids = 30, n_periods = 8)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))

  expect_true(obj$state$expanded)
  # Both arms should exist
  expect_true(all(c(0L, 1L) %in% obj$data[["_emulate_arm"]]))
})

test_that("emulate_expand respects maxfollowup", {
  d <- make_test_data(n_ids = 30, n_periods = 10)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 3))

  expect_true(max(obj$data[["_emulate_followup"]]) <= 3)
})

test_that("emulate_expand respects trials argument", {
  d <- make_test_data(n_ids = 30, n_periods = 8)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, trials = c(0, 1, 2)))

  trials_used <- unique(obj$data[["_emulate_trial"]])
  expect_true(all(trials_used %in% c(0, 1, 2)))
})

test_that("emulate_expand requires prepared object", {
  expect_error(
    emulate_expand(list(state = list(prepared = FALSE))),
    "emulate"
  )
})
