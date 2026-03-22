# V27: ATT vs ATE Predictions
#
# Validates that ATT (average treatment effect on the treated) and ATE
# (average treatment effect) produce different reference populations and
# that ATT correctly restricts to the treated subgroup.
#
# NOTE: emulate_predict() does not yet support the `att` parameter.
# Tests that require `att` skip gracefully. Test 2 validates the
# reference population structure without calling predict(att=TRUE).

# ---------------------------------------------------------------------------
# Helper: check if att parameter is supported
# ---------------------------------------------------------------------------
.has_att_param <- function() {

  "att" %in% names(formals(emulate_predict))
}

# ---------------------------------------------------------------------------
# Helper: run full pipeline
# ---------------------------------------------------------------------------
.v27_pipeline <- function(seed = 2700) {
  d <- dgp_simple(n = 5000, periods = 8, effect = -0.50, seed = seed)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
  obj
}

# ---------------------------------------------------------------------------
# Test 1: ATE and ATT produce different predictions (when att is supported)
# ---------------------------------------------------------------------------
test_that("V27-1: ATE and ATT predictions differ", {
  skip_on_cran()
  skip_if(!.has_att_param(), "emulate_predict() does not have att parameter yet")

  obj <- .v27_pipeline(seed = 2701)

  obj_ate <- suppressMessages(emulate_predict(obj, times = c(0, 3, 6),
    type = "cum_inc", samples = 50, seed = 2701))
  obj_att <- suppressMessages(emulate_predict(obj, times = c(0, 3, 6),
    type = "cum_inc", samples = 50, seed = 2701, att = TRUE))

  ate_est0 <- obj_ate$predictions$est_0
  att_est0 <- obj_att$predictions$est_0

  max_diff <- max(abs(ate_est0 - att_est0))
  expect_gt(max_diff, 1e-6,
    label = "ATE and ATT produce different control-arm predictions")
})

# ---------------------------------------------------------------------------
# Test 2: ATT reference population is subset of ATE population
# ---------------------------------------------------------------------------
test_that("V27-2: treated baseline is a strict subset of all baseline", {
  skip_on_cran()

  obj <- .v27_pipeline(seed = 2702)

  prefix <- obj$settings$prefix
  dt <- obj$data
  fu_col <- paste0(prefix, "followup")
  treat_col <- obj$settings$treatment

  baseline <- dt[dt[[fu_col]] == 0, ]
  n_ate_ref <- nrow(baseline)
  n_att_ref <- sum(baseline[[treat_col]] == 1)

  expect_lt(n_att_ref, n_ate_ref,
    label = "ATT reference pop is smaller than ATE reference pop")
  expect_gt(n_att_ref, 0,
    label = "ATT reference pop is non-empty")
})

# ---------------------------------------------------------------------------
# Test 3: Predictions have correct dimensions (when att is supported)
# ---------------------------------------------------------------------------
test_that("V27-3: ATE and ATT return same-shaped predictions", {
  skip_on_cran()
  skip_if(!.has_att_param(), "emulate_predict() does not have att parameter yet")

  obj <- .v27_pipeline(seed = 2703)

  times <- c(0, 2, 4, 6)
  obj_ate <- suppressMessages(emulate_predict(obj, times = times,
    difference = TRUE, samples = 50, seed = 2703))
  obj_att <- suppressMessages(emulate_predict(obj, times = times,
    difference = TRUE, samples = 50, seed = 2703, att = TRUE))

  expect_equal(nrow(obj_ate$predictions), length(times))
  expect_equal(nrow(obj_att$predictions), length(times))
  expect_equal(ncol(obj_ate$predictions), ncol(obj_att$predictions))
})

# ---------------------------------------------------------------------------
# Test 4: ATE ~ ATT under randomization (when att is supported)
# ---------------------------------------------------------------------------
test_that("V27-4: ATE ~ ATT under randomization (RCT data)", {
  skip_on_cran()
  skip_if(!.has_att_param(), "emulate_predict() does not have att parameter yet")

  d <- dgp_rct(n = 8000, periods = 8, effect = -0.50, seed = 2704)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))

  obj_ate <- suppressMessages(emulate_predict(obj, times = c(0, 3, 6),
    difference = TRUE, samples = 80, seed = 2704))
  obj_att <- suppressMessages(emulate_predict(obj, times = c(0, 3, 6),
    difference = TRUE, samples = 80, seed = 2704, att = TRUE))

  rd_ate <- obj_ate$predictions$diff[obj_ate$predictions$time == 6]
  rd_att <- obj_att$predictions$diff[obj_att$predictions$time == 6]

  expect_lt(abs(rd_ate - rd_att), 0.02,
    label = "ATE ~ ATT risk difference under randomization")
})

# ---------------------------------------------------------------------------
# Test 5: ATE prediction seed reproducibility
# ---------------------------------------------------------------------------
test_that("V27-5: ATE predictions are reproducible with same seed", {
  skip_on_cran()

  obj <- .v27_pipeline(seed = 2705)

  obj_a <- suppressMessages(emulate_predict(obj, times = c(0, 3, 6),
    samples = 50, seed = 42))
  obj_b <- suppressMessages(emulate_predict(obj, times = c(0, 3, 6),
    samples = 50, seed = 42))

  expect_equal(obj_a$predictions$est_0, obj_b$predictions$est_0,
    tolerance = 1e-12, label = "same seed produces identical est_0")
  expect_equal(obj_a$predictions$est_1, obj_b$predictions$est_1,
    tolerance = 1e-12, label = "same seed produces identical est_1")
})
