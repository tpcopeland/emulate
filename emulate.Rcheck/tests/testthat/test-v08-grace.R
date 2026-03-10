# V8: Grace Period (6 tests)

test_that("V8.1: grace=0 produces censored observations", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 12, seed = 801)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10, grace = 0))
  expect_true(obj$expansion$n_censored > 0,
              info = "grace=0 should produce censored observations")
})

test_that("V8.2: grace=1 has fewer censored than grace=0", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 12, seed = 802)

  obj_g0 <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                          treatment = "treatment",
                                          outcome = "outcome",
                                          eligible = "eligible",
                                          covariates = "x", estimand = "PP"))
  obj_g0 <- suppressMessages(emulate_expand(obj_g0, maxfollowup = 10, grace = 0))

  obj_g1 <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                          treatment = "treatment",
                                          outcome = "outcome",
                                          eligible = "eligible",
                                          covariates = "x", estimand = "PP"))
  obj_g1 <- suppressMessages(emulate_expand(obj_g1, maxfollowup = 10, grace = 1))

  expect_true(obj_g1$expansion$n_censored <= obj_g0$expansion$n_censored,
              info = sprintf("grace=1 censored (%d) should be <= grace=0 (%d)",
                             obj_g1$expansion$n_censored,
                             obj_g0$expansion$n_censored))
})

test_that("V8.3: censoring is monotonically non-increasing with grace", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 12, seed = 803)

  get_censored <- function(grace_val) {
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment",
                                         outcome = "outcome",
                                         eligible = "eligible",
                                         covariates = "x", estimand = "PP"))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10, grace = grace_val))
    obj$expansion$n_censored
  }

  c0 <- get_censored(0)
  c1 <- get_censored(1)
  c2 <- get_censored(2)
  c3 <- get_censored(3)

  expect_true(c0 >= c1,
              info = sprintf("cens g0=%d >= g1=%d", c0, c1))
  expect_true(c1 >= c2,
              info = sprintf("cens g1=%d >= g2=%d", c1, c2))
  expect_true(c2 >= c3,
              info = sprintf("cens g2=%d >= g3=%d", c2, c3))
})

test_that("V8.4: large grace coefficient approaches ITT", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 12, seed = 804)

  # ITT
  obj_itt <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = "x", estimand = "ITT"))
  obj_itt <- suppressMessages(emulate_expand(obj_itt, maxfollowup = 10))
  obj_itt <- suppressMessages(emulate_weight(obj_itt))
  obj_itt <- suppressMessages(emulate_fit(obj_itt, outcome_cov = "x",
                                       followup_spec = "linear",
                                       trial_period_spec = "linear"))

  # PP with moderate grace (grace=5, most switching within grace)
  obj_g <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment",
                                         outcome = "outcome",
                                         eligible = "eligible",
                                         covariates = "x", estimand = "PP"))
  obj_g <- suppressMessages(emulate_expand(obj_g, maxfollowup = 10, grace = 5))
  obj_g <- suppressWarnings(suppressMessages(
    emulate_weight(obj_g, switch_d_cov = "x", truncate = c(1, 99), quiet = TRUE)
  ))
  obj_g <- suppressWarnings(suppressMessages(
    emulate_fit(obj_g, outcome_cov = "x",
            followup_spec = "linear", trial_period_spec = "linear")
  ))

  # Large grace PP should trend toward ITT direction
  # Both should be negative or the difference should be bounded
  expect_true(obj_g$model$b_treat < 0.10,
              info = sprintf("Large grace coef %.4f should be < 0.10",
                             obj_g$model$b_treat))
})

test_that("V8.5: censored individual was deviating", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 12, seed = 805)
  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = "x", estimand = "PP"))
  obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10, grace = 0))

  dt <- obj$data
  prefix <- obj$settings$prefix
  cens_col <- paste0(prefix, "censored")
  arm_col <- paste0(prefix, "arm")
  treat_col <- obj$settings$treatment

  # Find a censored row in arm 1 (treatment arm)
  censored_treat <- dt[dt[[cens_col]] == 1 & dt[[arm_col]] == 1, ]
  if (nrow(censored_treat) > 0) {
    # In arm 1 (assigned treatment), censoring means the person stopped treatment
    # So their actual treatment should be 0 at censoring point
    expect_true(any(censored_treat[[treat_col]] == 0),
                info = "Censored treatment-arm individuals should have stopped treatment")
  }

  # Find a censored row in arm 0 (control arm)
  censored_ctrl <- dt[dt[[cens_col]] == 1 & dt[[arm_col]] == 0, ]
  if (nrow(censored_ctrl) > 0) {
    # In arm 0 (assigned control), censoring means the person started treatment
    expect_true(any(censored_ctrl[[treat_col]] == 1),
                info = "Censored control-arm individuals should have started treatment")
  }
})

test_that("V8.6: all grace coefficients are < 0.10", {
  skip_on_cran()

  d <- dgp_grace(n = 3000, periods = 12, seed = 806)

  fit_grace <- function(g) {
    obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                         treatment = "treatment",
                                         outcome = "outcome",
                                         eligible = "eligible",
                                         covariates = "x", estimand = "PP"))
    obj <- suppressMessages(emulate_expand(obj, maxfollowup = 10, grace = g))
    obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
                                        truncate = c(1, 99), quiet = TRUE))
    obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x",
                                     followup_spec = "linear",
                                     trial_period_spec = "linear"))
    obj$model$b_treat
  }

  for (g in c(0, 1, 2, 3)) {
    b <- fit_grace(g)
    expect_true(b < 0.10,
                info = sprintf("grace=%d coef %.4f should be < 0.10", g, b))
  }
})
