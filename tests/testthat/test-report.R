test_that("emulate_report displays results", {
  d <- make_test_data(n_ids = 50, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))

  expect_no_error(suppressMessages(emulate_report(obj)))
})

test_that("emulate_report exports CSV", {
  d <- make_test_data(n_ids = 50, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))

  tmp <- tempfile(fileext = ".csv")
  suppressMessages(emulate_report(obj, format = "csv", export = tmp))
  expect_true(file.exists(tmp))
  result <- read.csv(tmp)
  expect_true(nrow(result) > 0)
  unlink(tmp)
})
