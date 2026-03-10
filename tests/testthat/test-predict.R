test_that("emulate_predict produces predictions", {
  d <- make_test_data(n_ids = 50, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))
  obj <- suppressMessages(emulate_predict(obj, times = c(2, 4),
                                       samples = 20, seed = 42))

  expect_false(is.null(obj$predictions))
  expect_equal(nrow(obj$predictions), 2)
  expect_true(all(c("time", "est_0", "est_1") %in% names(obj$predictions)))
  expect_true(all(obj$predictions$est_0 >= 0 & obj$predictions$est_0 <= 1))
  expect_true(all(obj$predictions$est_1 >= 0 & obj$predictions$est_1 <= 1))
})

test_that("emulate_predict supports difference option", {
  d <- make_test_data(n_ids = 50, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))
  obj <- suppressMessages(emulate_predict(obj, times = c(3),
                                       samples = 20, seed = 42,
                                       difference = TRUE))

  expect_true("diff" %in% names(obj$predictions))
})

test_that("emulate_predict requires fitted model", {
  d <- make_test_data()
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible"))
  expect_error(emulate_predict(obj, times = c(1)), "expanded")
})
