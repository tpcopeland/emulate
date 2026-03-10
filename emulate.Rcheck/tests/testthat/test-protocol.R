test_that("emulate_protocol creates protocol table", {
  result <- suppressMessages(emulate_protocol(
    eligibility = "Adults aged 18+",
    treatment = "Drug A vs placebo",
    assignment = "Random assignment at baseline",
    followup_start = "Start of treatment",
    outcome = "All-cause mortality",
    causal_contrast = "ITT effect of assignment to Drug A",
    analysis = "Pooled logistic regression"
  ))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 7)
  expect_true("Component" %in% names(result))
  expect_true("Specification" %in% names(result))
})

test_that("emulate_protocol exports CSV", {
  tmp <- tempfile(fileext = ".csv")
  suppressMessages(emulate_protocol(
    eligibility = "Adults aged 18+",
    treatment = "Drug A vs placebo",
    assignment = "Random assignment",
    followup_start = "Baseline",
    outcome = "Mortality",
    causal_contrast = "ITT",
    analysis = "Logistic",
    format = "csv",
    export = tmp
  ))
  expect_true(file.exists(tmp))
  unlink(tmp)
})
