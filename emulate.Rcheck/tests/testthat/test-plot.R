test_that("emulate_plot creates cumhaz plot", {
  d <- make_test_data(n_ids = 50, n_periods = 8, event_prob = 0.03)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj))
  obj <- suppressMessages(emulate_predict(obj, times = c(1, 3, 5),
                                       samples = 20, seed = 42))

  p <- emulate_plot(obj, type = "cumhaz")
  expect_s3_class(p, "ggplot")
})

test_that("emulate_plot creates weights plot", {
  d <- make_test_data(n_ids = 50, n_periods = 8, treat_prob = 0.3)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA"),
                                      quiet = TRUE))
  p <- emulate_plot(obj, type = "weights")
  expect_s3_class(p, "ggplot")
})
