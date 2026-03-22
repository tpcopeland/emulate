test_that("emulate_prepare creates valid emulate object", {
  d <- make_test_data()
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible"))
  expect_s3_class(obj, "emulate")
  expect_true(obj$state$prepared)
  expect_false(obj$state$expanded)
  expect_equal(obj$settings$estimand, "PP")
})

test_that("emulate_prepare validates required columns", {
  d <- make_test_data()
  expect_error(
    suppressMessages(emulate_prepare(d, id = "id", period = "MISSING",
                                  treatment = "treatment",
                                  outcome = "outcome",
                                  eligible = "eligible")),
    "Missing columns"
  )
})

test_that("emulate_prepare rejects non-binary treatment", {
  d <- make_test_data()
  d$treatment[1] <- 2
  expect_error(
    suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                  treatment = "treatment",
                                  outcome = "outcome",
                                  eligible = "eligible")),
    "binary"
  )
})

test_that("emulate_prepare rejects duplicate id-period", {
  d <- make_test_data()
  d <- rbind(d, d[1, ])
  expect_error(
    suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                  treatment = "treatment",
                                  outcome = "outcome",
                                  eligible = "eligible")),
    "person-period"
  )
})

test_that("emulate_prepare accepts covariates", {
  d <- make_test_data()
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("nvarA", "nvarB")))
  expect_equal(obj$settings$covariates, c("nvarA", "nvarB"))
})

test_that("emulate_prepare validates estimand", {
  d <- make_test_data()
  expect_error(
    suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                  treatment = "treatment",
                                  outcome = "outcome",
                                  eligible = "eligible",
                                  estimand = "WRONG")),
    "estimand"
  )
})
