test_that("emulate_diagnose reports weight distribution", {
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
  obj <- suppressMessages(emulate_diagnose(obj))

  expect_true(!is.null(obj$diagnostics$ess))
  expect_true(obj$diagnostics$ess > 0)
})

test_that("emulate_diagnose computes covariate balance", {
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
  obj <- suppressMessages(emulate_diagnose(obj,
                                        balance_covariates = c("nvarA", "nvarB")))

  bal <- obj$diagnostics$balance
  expect_equal(nrow(bal), 2)
  expect_true("smd_unwt" %in% names(bal))
  expect_true("smd_wt" %in% names(bal))
})
