# V2: Synthetic NHEFS-style smoking cessation validation
# Mirrors: Stata-Tools/tte/qa/validate_nhefs.do
#
# Since we may not have NHEFS data shipped, this uses a synthetic dataset
# mimicking the NHEFS structure: smoking cessation, mortality over 10 years,
# ~1600 subjects.
#
# Reference: Hernan MA, Robins JM. Causal Inference: What If. 2020.
#   IP-weighted HR for smoking cessation on mortality: ~0.80-0.90
#   Chapter 12: IP weighting; Chapter 17: Causal survival analysis

# ---------------------------------------------------------------------------
# DGP: Synthetic NHEFS (smoking cessation & 10-year mortality)
# ---------------------------------------------------------------------------
dgp_nhefs <- function(n = 1600, periods = 10, seed = 20260309) {
  set.seed(seed)

  rows <- vector("list", n * periods)
  idx <- 0L

  for (i in seq_len(n)) {
    # Baseline confounders
    age <- rnorm(1, 45, 12)       # mean age ~45
    age_std <- (age - 45) / 12
    sex <- sample(0:1, 1)         # binary sex
    wt_base <- rnorm(1, 75, 15)   # baseline weight (kg)
    wt_std <- (wt_base - 75) / 15

    # Treatment: quit smoking (~25% prevalence, confounded by age and sex)
    # Older and female subjects slightly more likely to quit
    p_quit <- plogis(-1.1 + 0.15 * age_std + 0.3 * sex)
    treat <- as.integer(runif(1) < p_quit)

    alive <- TRUE

    for (t in 0:(periods - 1L)) {
      if (!alive) break

      # Outcome: death (true protective effect of quitting, log-OR ~ -0.4)
      # Hazard increases with age, higher baseline weight, and time
      log_haz <- -4.0 + 0.4 * age_std + 0.15 * wt_std +
                 0.05 * t - 0.40 * treat
      p_death <- plogis(log_haz)
      event <- as.integer(runif(1) < p_death)

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        id = i,
        period = t,
        treatment = treat,    # time-invariant (quit at baseline)
        outcome = event,
        eligible = as.integer(t == 0),  # single enrollment at baseline
        age_std = age_std,
        sex = sex,
        wt_std = wt_std,
        stringsAsFactors = FALSE
      )

      if (event == 1L) alive <- FALSE
    }
  }

  do.call(rbind, rows[seq_len(idx)])
}

# ---------------------------------------------------------------------------
# Test 1: ITT shows protective effect (negative coefficient, OR in [0.3, 1.5])
# ---------------------------------------------------------------------------
test_that("V02-1: NHEFS-style ITT shows protective effect of quitting", {
  skip_on_cran()

  d <- dgp_nhefs()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("age_std", "sex", "wt_std"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("age_std", "sex", "wt_std"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "none"))

  b  <- obj$model$b_treat
  or <- exp(b)

  # Direction: should be negative (quitting is protective)
  expect_true(b < 0,
              info = sprintf("ITT coefficient should be negative (protective), got %.4f", b))

  # Magnitude: OR should be in plausible range [0.3, 1.5]
  expect_true(or > 0.3 & or < 1.5,
              info = sprintf("ITT odds ratio = %.4f, expected in [0.3, 1.5]", or))
})

# ---------------------------------------------------------------------------
# Test 2: Survival curves are valid (CI increases over 10 years, both in [0,1])
# ---------------------------------------------------------------------------
test_that("V02-2: NHEFS-style survival curves are valid over 10 years", {
  skip_on_cran()

  d <- dgp_nhefs()

  obj <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                       treatment = "treatment",
                                       outcome = "outcome",
                                       eligible = "eligible",
                                       covariates = c("age_std", "sex", "wt_std"),
                                       estimand = "ITT"))
  obj <- suppressMessages(emulate_expand(obj))
  obj <- suppressMessages(emulate_weight(obj))
  obj <- suppressMessages(emulate_fit(obj,
                                   outcome_cov = c("age_std", "sex", "wt_std"),
                                   followup_spec = "quadratic",
                                   trial_period_spec = "none"))
  obj <- suppressMessages(emulate_predict(obj, times = 0:9, type = "cum_inc",
                                       samples = 50, seed = 42,
                                       difference = TRUE))

  pred <- obj$predictions

  # Cumulative incidence should increase over time
  ci_start_0 <- pred$est_0[1]
  ci_end_0   <- pred$est_0[nrow(pred)]
  ci_start_1 <- pred$est_1[1]
  ci_end_1   <- pred$est_1[nrow(pred)]

  expect_true(ci_end_0 > ci_start_0,
              info = "Control arm CI should increase over 10 years")
  expect_true(ci_end_1 > ci_start_1,
              info = "Treated arm CI should increase over 10 years")

  # Both arms in [0, 1]
  expect_true(all(pred$est_0 >= 0 & pred$est_0 <= 1),
              info = "Control arm predictions outside [0,1]")
  expect_true(all(pred$est_1 >= 0 & pred$est_1 <= 1),
              info = "Treated arm predictions outside [0,1]")

  # Final CI should be reasonable for 10-year mortality (~5-30%)
  expect_true(ci_end_0 > 0.01 & ci_end_0 < 0.90,
              info = sprintf("Control CI at t=9 = %.4f, expected in [0.01, 0.90]",
                             ci_end_0))
})

# ---------------------------------------------------------------------------
# Test 3: Pooled logistic vs Cox — both in plausible range (OR 0.5-2.0)
# ---------------------------------------------------------------------------
test_that("V02-3: Logistic and Cox both produce plausible estimates", {
  skip_on_cran()

  d <- dgp_nhefs()

  # Pooled logistic
  obj_logit <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                             treatment = "treatment",
                                             outcome = "outcome",
                                             eligible = "eligible",
                                             covariates = c("age_std", "sex", "wt_std"),
                                             estimand = "ITT"))
  obj_logit <- suppressMessages(emulate_expand(obj_logit))
  obj_logit <- suppressMessages(emulate_weight(obj_logit))
  obj_logit <- suppressMessages(emulate_fit(obj_logit,
                                         outcome_cov = c("age_std", "sex", "wt_std"),
                                         model = "logistic",
                                         followup_spec = "quadratic",
                                         trial_period_spec = "none"))

  # Cox model
  obj_cox <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = c("age_std", "sex", "wt_std"),
                                           estimand = "ITT"))
  obj_cox <- suppressMessages(emulate_expand(obj_cox))
  obj_cox <- suppressMessages(emulate_weight(obj_cox))
  obj_cox <- suppressMessages(emulate_fit(obj_cox,
                                       outcome_cov = c("age_std", "sex", "wt_std"),
                                       model = "cox",
                                       trial_period_spec = "none"))

  or_logit <- exp(obj_logit$model$b_treat)
  hr_cox   <- exp(obj_cox$model$b_treat)

  # Both in plausible range [0.5, 2.0]
  expect_true(or_logit > 0.5 & or_logit < 2.0,
              info = sprintf("Logistic OR = %.4f, expected in [0.5, 2.0]", or_logit))
  expect_true(hr_cox > 0.5 & hr_cox < 2.0,
              info = sprintf("Cox HR = %.4f, expected in [0.5, 2.0]", hr_cox))
})

# ---------------------------------------------------------------------------
# Test 4: Cox and logistic models agree on direction
# ---------------------------------------------------------------------------
test_that("V02-4: Cox and logistic agree on direction of effect", {
  skip_on_cran()

  d <- dgp_nhefs()

  # Pooled logistic
  obj_logit <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                             treatment = "treatment",
                                             outcome = "outcome",
                                             eligible = "eligible",
                                             covariates = c("age_std", "sex", "wt_std"),
                                             estimand = "ITT"))
  obj_logit <- suppressMessages(emulate_expand(obj_logit))
  obj_logit <- suppressMessages(emulate_weight(obj_logit))
  obj_logit <- suppressMessages(emulate_fit(obj_logit,
                                         outcome_cov = c("age_std", "sex", "wt_std"),
                                         model = "logistic",
                                         followup_spec = "quadratic",
                                         trial_period_spec = "none"))

  # Cox model
  obj_cox <- suppressMessages(emulate_prepare(d, id = "id", period = "period",
                                           treatment = "treatment",
                                           outcome = "outcome",
                                           eligible = "eligible",
                                           covariates = c("age_std", "sex", "wt_std"),
                                           estimand = "ITT"))
  obj_cox <- suppressMessages(emulate_expand(obj_cox))
  obj_cox <- suppressMessages(emulate_weight(obj_cox))
  obj_cox <- suppressMessages(emulate_fit(obj_cox,
                                       outcome_cov = c("age_std", "sex", "wt_std"),
                                       model = "cox",
                                       trial_period_spec = "none"))

  b_logit <- obj_logit$model$b_treat
  b_cox   <- obj_cox$model$b_treat

  # Both should point in the same direction
  expect_equal(sign(b_logit), sign(b_cox),
               info = sprintf("Logistic (%.4f) and Cox (%.4f) disagree on direction",
                              b_logit, b_cox))
})
