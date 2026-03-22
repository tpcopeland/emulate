test_that("emulate_fit fits logistic model for ITT", {
  d <- make_test_data(n_ids = 50, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))

  expect_true(obj$state$fitted)
  expect_equal(obj$model$type, "logistic")
  expect_true(!is.null(obj$model$object))
  expect_true(!is.null(obj$model$vcov))
  expect_true("_emulate_esample" %in% names(obj$data))
})

test_that("emulate_fit supports followup_spec options", {
  d <- make_test_data(n_ids = 50, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj))

  # Linear
  obj_l <- suppressMessages(emulate_fit(obj, followup_spec = "linear"))
  expect_true(obj_l$state$fitted)

  # Quadratic (default)
  obj_q <- suppressMessages(emulate_fit(obj, followup_spec = "quadratic"))
  expect_true(obj_q$state$fitted)
})

test_that("emulate_fit requires expanded data", {
  d <- make_test_data()
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible"))
  expect_error(emulate_fit(obj), "expanded")
})
