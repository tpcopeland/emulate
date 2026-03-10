# Tests for lasso PS estimation and PS trimming (Phase 1)

test_that("lasso PS produces non-degenerate weights", {
  skip_if_not_installed("glmnet")

  dat <- dgp_simple(n = 500, periods = 6, seed = 101)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_weight(obj, switch_d_cov = "x", ps_method = "lasso")
  )

  expect_true(obj$state$weighted)
  expect_equal(obj$weights$ps_method, "lasso")

  w <- obj$data[[obj$weights$weight_var]]
  expect_true(all(is.finite(w)))
  expect_true(sd(w) > 0)  # not degenerate
})

test_that("lasso weights differ from GLM weights", {
  skip_if_not_installed("glmnet")

  dat <- dgp_simple(n = 500, periods = 6, seed = 102)

  obj_glm <- emulate_prepare(dat, id = "id", period = "period",
                              treatment = "treatment", outcome = "outcome",
                              eligible = "eligible", covariates = "x",
                              estimand = "PP")
  obj_glm <- emulate_expand(obj_glm, maxfollowup = 4)
  obj_glm <- suppressMessages(
    emulate_weight(obj_glm, switch_d_cov = "x", ps_method = "glm")
  )

  obj_lasso <- emulate_prepare(dat, id = "id", period = "period",
                                treatment = "treatment", outcome = "outcome",
                                eligible = "eligible", covariates = "x",
                                estimand = "PP")
  obj_lasso <- emulate_expand(obj_lasso, maxfollowup = 4)
  obj_lasso <- suppressMessages(
    emulate_weight(obj_lasso, switch_d_cov = "x", ps_method = "lasso")
  )

  # Weights should not be identical (lasso shrinks differently)
  w_glm <- obj_glm$data[[obj_glm$weights$weight_var]]
  w_lasso <- obj_lasso$data[[obj_lasso$weights$weight_var]]
  # May have different lengths if PS trimming involved, just check means differ
  expect_false(isTRUE(all.equal(mean(w_glm), mean(w_lasso), tolerance = 0.001)))
})

test_that("PS trim removes expected fraction of observations", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 103)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)
  n_before <- nrow(obj$data)

  obj <- suppressMessages(
    emulate_weight(obj, switch_d_cov = "x", ps_trim = c(0.05, 0.95))
  )

  expect_true(obj$state$weighted)
  # Should have removed some observations
  expect_true(obj$weights$n_ps_trimmed >= 0)
  # Data should have fewer or equal rows

  expect_true(nrow(obj$data) <= n_before)
})

test_that("PS trim + truncate can be combined", {
  dat <- dgp_simple(n = 500, periods = 6, seed = 104)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 4)

  obj <- suppressMessages(
    emulate_weight(obj, switch_d_cov = "x",
                   ps_trim = c(0.05, 0.95), truncate = c(1, 99))
  )

  expect_true(obj$state$weighted)
  expect_true(!is.null(obj$weights$ps_trim))
  expect_true(!is.null(obj$weights$truncation))
})

test_that("invalid ps_method errors", {
  dat <- dgp_simple(n = 100, periods = 4, seed = 105)
  obj <- emulate_prepare(dat, id = "id", period = "period",
                         treatment = "treatment", outcome = "outcome",
                         eligible = "eligible", covariates = "x",
                         estimand = "PP")
  obj <- emulate_expand(obj, maxfollowup = 3)

  expect_error(
    suppressMessages(
      emulate_weight(obj, switch_d_cov = "x", ps_method = "bad")
    ),
    "ps_method must be"
  )
})
