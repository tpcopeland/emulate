test_that("pipeline functions enforce correct order", {
  d <- make_test_data()
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible"))

  # Can't expand before prepare (already prepared, but test weight/fit guards)
  expect_error(emulate_fit(obj), "expanded")
  expect_error(emulate_predict(obj, times = c(1)), "expanded")
  expect_error(emulate_diagnose(obj), "expanded")
})

test_that("full ITT pipeline runs end to end", {
  d <- make_test_data(n_ids = 50, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_validate(obj))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))
  obj <- suppressMessages(emulate_predict(obj, times = c(2, 4),
                                       samples = 20, seed = 42))

  expect_true(obj$state$prepared)
  expect_true(obj$state$expanded)
  expect_true(obj$state$weighted)
  expect_true(obj$state$fitted)
  expect_false(is.null(obj$predictions))
})

test_that("print.emulate works at each pipeline stage", {
  d <- make_test_data(n_ids = 30, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  expect_output(print(obj), "Target Trial")

  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  expect_output(print(obj), "expanded")
})
