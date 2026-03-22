# V30: Grace Period Monotonicity and Edge Cases
#
# Validates grace period handling in emulate_expand(): censoring timing,
# monotonicity properties, and boundary conditions.

# ---------------------------------------------------------------------------
# Test 1: grace=0 censors immediately on deviation
# ---------------------------------------------------------------------------
test_that("V30-1: grace=0 censors at first deviation", {
  skip_on_cran()

  d <- dgp_grace(n = 2000, periods = 10, seed = 3001)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8, grace = 0))

  prefix <- obj$settings$prefix
  dt <- obj$data
  arm_col <- paste0(prefix, "arm")
  cens_col <- paste0(prefix, "censored")
  fu_col <- paste0(prefix, "followup")

  # In treatment arm, censored observations should have treatment == 0
  # (they deviated at or before the censoring period)
  treat_cens <- dt[dt[[arm_col]] == 1 & dt[[cens_col]] == 1, ]
  if (nrow(treat_cens) > 0) {
    expect_true(all(treat_cens[[obj$settings$treatment]] == 0),
      label = "censored treatment-arm obs have treatment=0")
  }
})

# ---------------------------------------------------------------------------
# Test 2: Larger grace period produces more uncensored observations
# ---------------------------------------------------------------------------
test_that("V30-2: larger grace period produces more uncensored observations", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 10, seed = 3002)

  run_expand <- function(grace_val) {
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
      treatment = "treatment", outcome = "outcome", eligible = "eligible",
      covariates = "x", estimand = "PP"))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8,
      grace = grace_val))
    obj$expansion
  }

  exp_g0 <- run_expand(0)
  exp_g1 <- run_expand(1)
  exp_g2 <- run_expand(2)

  # More grace -> fewer censored -> more total observations
  expect_lte(exp_g0$n_expanded, exp_g1$n_expanded + 1,
    label = "grace=1 has >= obs than grace=0")
  expect_lte(exp_g1$n_expanded, exp_g2$n_expanded + 1,
    label = "grace=2 has >= obs than grace=1")

  # More grace -> fewer censored
  expect_gte(exp_g0$n_censored, exp_g1$n_censored,
    label = "grace=1 has <= censored than grace=0")
})

# ---------------------------------------------------------------------------
# Test 3: grace period does not affect ITT
# ---------------------------------------------------------------------------
test_that("V30-3: grace parameter has no effect on ITT expansion", {
  skip_on_cran()

  d <- dgp_grace(n = 2000, periods = 10, seed = 3003)

  run_itt <- function(grace_val) {
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
      treatment = "treatment", outcome = "outcome", eligible = "eligible",
      covariates = "x", estimand = "ITT"))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8,
      grace = grace_val))
    obj$expansion
  }

  exp_g0 <- run_itt(0)
  exp_g3 <- run_itt(3)

  expect_equal(exp_g0$n_expanded, exp_g3$n_expanded,
    label = "ITT expansion identical regardless of grace")
  expect_equal(exp_g0$n_censored, 0,
    label = "ITT has zero censored observations")
  expect_equal(exp_g3$n_censored, 0,
    label = "ITT with grace=3 still has zero censored")
})

# ---------------------------------------------------------------------------
# Test 4: No rows exist AFTER censoring point
# ---------------------------------------------------------------------------
test_that("V30-4: no follow-up rows exist after censoring point", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 10, seed = 3004)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8, grace = 1))

  prefix <- obj$settings$prefix
  dt <- obj$data
  id_col <- obj$settings$id
  trial_col <- paste0(prefix, "trial")
  arm_col <- paste0(prefix, "arm")
  cens_col <- paste0(prefix, "censored")
  fu_col <- paste0(prefix, "followup")

  # For each censored observation, no rows should exist at later follow-up
  censored <- dt[dt[[cens_col]] == 1, ]
  if (nrow(censored) > 0) {
    # Check a sample of censored observations
    sample_idx <- head(seq_len(nrow(censored)), 100)
    for (i in sample_idx) {
      row <- censored[i, ]
      post <- dt[dt[[id_col]] == row[[id_col]] &
                 dt[[trial_col]] == row[[trial_col]] &
                 dt[[arm_col]] == row[[arm_col]] &
                 dt[[fu_col]] > row[[fu_col]], ]
      expect_equal(nrow(post), 0,
        label = paste("no rows after censoring for id",
          row[[id_col]], "trial", row[[trial_col]]))
    }
  }
})

# ---------------------------------------------------------------------------
# Test 5: Outcome is 0 at censoring point
# ---------------------------------------------------------------------------
test_that("V30-5: outcome_obs is 0 at censoring point", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 10, seed = 3005)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 8, grace = 2))

  prefix <- obj$settings$prefix
  dt <- obj$data
  cens_col <- paste0(prefix, "censored")
  outobs_col <- paste0(prefix, "outcome_obs")

  censored <- dt[dt[[cens_col]] == 1, ]
  if (nrow(censored) > 0) {
    expect_true(all(censored[[outobs_col]] == 0),
      label = "outcome_obs is 0 at all censoring points")
  }
})

# ---------------------------------------------------------------------------
# Test 6: Very large grace period effectively disables censoring
# ---------------------------------------------------------------------------
test_that("V30-6: grace >= maxfollowup effectively disables censoring", {
  skip_on_cran()

  d <- dgp_grace(n = 2000, periods = 10, seed = 3006)

  obj_pp <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "PP"))
  obj_pp <- suppressMessages(emulate_expand(obj_pp, maxfollowup = 6,
    grace = 100))

  obj_itt <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
    treatment = "treatment", outcome = "outcome", eligible = "eligible",
    covariates = "x", estimand = "ITT"))
  obj_itt <- suppressMessages(emulate_expand(obj_itt, maxfollowup = 6))

  # With grace >= maxfollowup, PP should have very few censored observations
  # (only those who deviate at exactly the last allowed period)
  # Total obs should be close to ITT (which has no censoring)
  # PP with large grace is close to ITT in terms of observations
  ratio <- obj_pp$expansion$n_expanded / obj_itt$expansion$n_expanded
  expect_gt(ratio, 0.8,
    label = "PP with huge grace has similar obs count to ITT")
})
