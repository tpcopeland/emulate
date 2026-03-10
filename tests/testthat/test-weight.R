test_that("emulate_weight sets all weights to 1 for ITT", {
  d <- make_test_data(n_ids = 30, n_periods = 8)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))

  expect_true(obj$state$weighted)
  expect_true(all(obj$data[["_emulate_weight"]] == 1))
})

test_that("emulate_weight computes weights for PP", {
  d <- make_test_data(n_ids = 50, n_periods = 8, treat_prob = 0.3)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB"),
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 5))
  obj <- suppressMessages(emulate_weight(obj,
                                      switch_d_cov = c("nvarA", "nvarB"),
                                      quiet = TRUE))

  expect_true(obj$state$weighted)
  w <- obj$data[["_emulate_weight"]]
  expect_true(all(!is.na(w)))
  expect_true(all(w > 0))
  # Weights should not all be 1 for PP
  expect_false(all(w == 1))
})

test_that("emulate_weight truncation works", {
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
                                      truncate = c(1, 99),
                                      quiet = TRUE))

  expect_true(obj$state$weighted)
})

test_that("emulate_weight requires switch_d_cov for PP", {
  d <- make_test_data(n_ids = 30, n_periods = 8)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj))
  expect_error(
    suppressMessages(emulate_weight(obj)),
    "switch_d_cov required"
  )
})
