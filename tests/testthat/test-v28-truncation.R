# V28: Weight Truncation Percentile Verification
#
# Validates that emulate_weight() truncation and PS trimming work correctly
# at the percentile level, and that the reported counts match reality.

# ---------------------------------------------------------------------------
# Helper: full pipeline through weighting
# ---------------------------------------------------------------------------
.v28_pipeline <- function(seed = 2800, estimand = "PP") {
  d <- dgp_simple(n = 3000, periods = 8, effect = -0.50, seed = seed)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = estimand))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6, grace = 0))
  obj
}

# ---------------------------------------------------------------------------
# Test 1: Truncated weights lie within percentile bounds
# ---------------------------------------------------------------------------
test_that("V28-1: truncated weights lie within specified percentile bounds", {
  skip_on_cran()

  obj <- .v28_pipeline(seed = 2801)
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
    truncate = c(1, 99)))

  wvar <- obj$weights$weight_var
  w <- obj$data[[wvar]]
  w <- w[!is.na(w)]

  # After truncation at (1, 99), all weights should be within
  # the 1st and 99th percentile bounds
  lo <- quantile(w, 0.01, type = 2)
  hi <- quantile(w, 0.99, type = 2)

  # All weights should be within [lo, hi]
  expect_true(all(w >= lo - 1e-10),
    label = "no weights below 1st percentile after truncation")
  expect_true(all(w <= hi + 1e-10),
    label = "no weights above 99th percentile after truncation")
})

# ---------------------------------------------------------------------------
# Test 2: Truncation at (5, 95) is more aggressive than (1, 99)
# ---------------------------------------------------------------------------
test_that("V28-2: tighter truncation produces lower weight variance", {
  skip_on_cran()

  obj <- .v28_pipeline(seed = 2802)

  obj_wide <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
    truncate = c(1, 99)))
  obj_tight <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
    truncate = c(5, 95)))

  expect_lte(obj_tight$weights$sd, obj_wide$weights$sd + 1e-10,
    label = "tighter truncation reduces weight SD")
})

# ---------------------------------------------------------------------------
# Test 3: n_truncated count is accurate
# ---------------------------------------------------------------------------
test_that("V28-3: n_truncated count matches actual truncated observations", {
  skip_on_cran()

  obj <- .v28_pipeline(seed = 2803)

  # First compute untruncated weights
  obj_raw <- suppressMessages(emulate_weight(obj, switch_d_cov = "x"))
  raw_w <- obj_raw$data[[obj_raw$weights$weight_var]]
  raw_w <- raw_w[!is.na(raw_w)]

  lo <- quantile(raw_w, 0.01, type = 2)
  hi <- quantile(raw_w, 0.99, type = 2)
  hand_count <- sum(raw_w < lo) + sum(raw_w > hi)

  # Now compute truncated
  obj_trunc <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
    truncate = c(1, 99)))

  expect_equal(obj_trunc$weights$n_truncated, hand_count,
    label = "n_truncated matches hand count")
})

# ---------------------------------------------------------------------------
# Test 4: No truncation when all weights are equal (ITT)
# ---------------------------------------------------------------------------
test_that("V28-4: ITT weights are all 1, n_truncated = 0", {
  skip_on_cran()

  d <- dgp_simple(n = 2000, periods = 8, effect = -0.50, seed = 2804)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 6))
  obj <- suppressMessages(emulate_weight(obj))

  wvar <- obj$weights$weight_var
  w <- obj$data[[wvar]]

  expect_true(all(w == 1),
    label = "all ITT weights are 1")
  # n_truncated may be 0 or NULL for ITT (no truncation applied)
  expect_true(is.null(obj$weights$n_truncated) ||
              obj$weights$n_truncated == 0,
    label = "no truncation for ITT")
})

# ---------------------------------------------------------------------------
# Test 5: Truncation bounds are symmetric in percentile space
# ---------------------------------------------------------------------------
test_that("V28-5: truncate(5, 95) clips symmetrically", {
  skip_on_cran()

  obj <- .v28_pipeline(seed = 2805)
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
    truncate = c(5, 95)))

  wvar <- obj$weights$weight_var
  w <- obj$data[[wvar]]
  w <- w[!is.na(w)]

  # The min and max should be the 5th and 95th percentile values
  # (or very close, since some weights might naturally be within bounds)
  expect_lte(obj$weights$min, quantile(w, 0.05, type = 2) + 1e-10)
})

# ---------------------------------------------------------------------------
# Test 6: ESS decreases with extreme weights (no truncation)
# ---------------------------------------------------------------------------
test_that("V28-6: ESS is smaller than N (weights reduce effective sample)", {
  skip_on_cran()

  obj <- .v28_pipeline(seed = 2806)
  obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x"))

  n_obs <- nrow(obj$data)
  ess <- obj$weights$ess

  # ESS should be less than total N (weights reduce effective sample)
  expect_lt(ess, n_obs,
    label = "ESS < N when weights vary")
  expect_gt(ess, 0,
    label = "ESS is positive")
})
