test_that("emulate_validate runs 10 checks on clean data", {
  d <- make_test_data()
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible"))
  obj <- suppressMessages(emulate_validate(obj))
  expect_equal(obj$validation$n_checks, 10)
  expect_equal(obj$validation$n_errors, 0)
})

test_that("emulate_validate detects missing data in strict mode", {
  d <- make_test_data()
  d$treatment[5] <- NA
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible"))
  # Strict mode should error
  expect_error(
    suppressMessages(emulate_validate(obj, strict = TRUE)),
    "validation failed"
  )
})

test_that("emulate_validate returns object invisibly", {
  d <- make_test_data()
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible"))
  result <- suppressMessages(emulate_validate(obj))
  expect_s3_class(result, "emulate")
})
