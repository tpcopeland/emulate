# Shared data-generating process (DGP) functions for emulate validation suite
#
# All generators return a data.frame in person-period format ready for
# emulate_prepare(). Treatment is absorbing, outcome is terminal.

# ---------------------------------------------------------------------------
# dgp_simple: Core confounded DGP
# ---------------------------------------------------------------------------
dgp_simple <- function(n = 2000, periods = 10, effect = -0.50,
                        seed = 42, outcome_intercept = -3.5) {
  set.seed(seed)

  rows <- vector("list", n * periods)
  idx <- 0L

  for (i in seq_len(n)) {
    x <- rnorm(1)
    treat <- 0L
    alive <- TRUE

    for (t in 0:(periods - 1L)) {
      if (!alive) break

      # Treatment initiation (absorbing): P(T=1) = invlogit(-2 + 0.3*x)
      if (treat == 0L) {
        p_treat <- plogis(-2 + 0.3 * x)
        if (runif(1) < p_treat) treat <- 1L
      }

      # Outcome: P(Y=1) = invlogit(intercept + 0.3*x + effect*treat)
      p_out <- plogis(outcome_intercept + 0.3 * x + effect * treat)
      event <- as.integer(runif(1) < p_out)

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        id = i, period = t, treatment = treat, outcome = event,
        eligible = 1L, x = x, stringsAsFactors = FALSE
      )

      if (event == 1L) alive <- FALSE
    }
  }

  do.call(rbind, rows[seq_len(idx)])
}

# ---------------------------------------------------------------------------
# dgp_null: Null effect (effect = 0)
# ---------------------------------------------------------------------------
dgp_null <- function(n = 5000, periods = 10, seed = 42) {
  dgp_simple(n = n, periods = periods, effect = 0, seed = seed,
             outcome_intercept = -3.5)
}

# ---------------------------------------------------------------------------
# dgp_ccw: Surgery timing DGP for CCW/immortal-time bias tests
# ---------------------------------------------------------------------------
dgp_ccw <- function(n = 2000, seed = 20260303) {
  set.seed(seed)
  periods <- 24L

  rows <- vector("list", n * periods)
  idx <- 0L

  for (i in seq_len(n)) {
    age <- rnorm(1, 60, 10)
    age_std <- (age - 60) / 10
    ps <- rnorm(1)            # propensity score confounder
    stage <- sample(0:2, 1)   # disease stage

    # Surgery timing: some get surgery early, some late, some never
    surgery_time <- if (runif(1) < 0.4) {
      sample(0:12, 1, prob = pmax(0.01, dnorm(0:12, 4, 3)))
    } else {
      NA_integer_
    }

    treat <- 0L
    alive <- TRUE

    for (t in 0:(periods - 1L)) {
      if (!alive) break

      # Surgery is absorbing
      if (treat == 0L && !is.na(surgery_time) && t >= surgery_time) {
        treat <- 1L
      }

      # Exponential survival: hazard depends on age, stage, and treatment
      # True HR for treatment ~ 0.60
      log_haz <- -4 + 0.2 * age_std + 0.3 * stage + log(0.60) * treat
      p_out <- 1 - exp(-exp(log_haz))
      event <- as.integer(runif(1) < p_out)

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        id = i, period = t, treatment = treat, outcome = event,
        eligible = 1L, age_std = age_std, ps = ps,
        stage = as.numeric(stage), stringsAsFactors = FALSE
      )

      if (event == 1L) alive <- FALSE
    }
  }

  do.call(rbind, rows[seq_len(idx)])
}

# ---------------------------------------------------------------------------
# dgp_gformula: Time-varying CD4, ART, confounding by indication
# ---------------------------------------------------------------------------
dgp_gformula <- function(n = 5000, periods = 15, seed = 20260304) {
  set.seed(seed)

  rows <- vector("list", n * periods)
  idx <- 0L

  for (i in seq_len(n)) {
    age_cat <- sample(0:2, 1)
    male <- sample(0:1, 1)
    cd4 <- rnorm(1, 500, 150)
    treat <- 0L
    alive <- TRUE

    for (t in 0:(periods - 1L)) {
      if (!alive) break

      # Time-varying CD4 (declines, treatment helps)
      if (t > 0) {
        cd4 <- cd4 - 10 + 20 * treat + rnorm(1, 0, 30)
        cd4 <- max(cd4, 50)
      }
      cd4_std <- (cd4 - 500) / 150

      # ART initiation (absorbing): more likely with low CD4
      if (treat == 0L) {
        p_treat <- plogis(-1.5 - 0.8 * cd4_std + 0.2 * age_cat)
        if (runif(1) < p_treat) treat <- 1L
      }

      # Outcome: treatment is protective (effect = -0.80)
      p_out <- plogis(-4.5 + 0.3 * age_cat + 0.1 * male -
                        0.4 * cd4_std - 0.80 * treat)
      event <- as.integer(runif(1) < p_out)

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        id = i, period = t, treatment = treat, outcome = event,
        eligible = 1L, cd4_std = cd4_std,
        age_cat = as.numeric(age_cat), male = male,
        stringsAsFactors = FALSE
      )

      if (event == 1L) alive <- FALSE
    }
  }

  do.call(rbind, rows[seq_len(idx)])
}

# ---------------------------------------------------------------------------
# dgp_ipcw: Adds informative censoring
# ---------------------------------------------------------------------------
dgp_ipcw <- function(n = 5000, periods = 10, effect = -0.60,
                      seed = 20260305) {
  set.seed(seed)

  rows <- vector("list", n * periods)
  idx <- 0L

  for (i in seq_len(n)) {
    x <- rnorm(1)
    z <- rnorm(1)
    treat <- 0L
    alive <- TRUE
    censored <- FALSE

    for (t in 0:(periods - 1L)) {
      if (!alive || censored) break

      # Informative censoring: P(C=1) = invlogit(-3 + 0.5*x + 0.4*z)
      if (t > 0) {
        p_cens <- plogis(-3 + 0.5 * x + 0.4 * z)
        if (runif(1) < p_cens) {
          censored <- TRUE
          break
        }
      }

      # Treatment (absorbing)
      if (treat == 0L) {
        p_treat <- plogis(-2 + 0.3 * x)
        if (runif(1) < p_treat) treat <- 1L
      }

      # Outcome
      p_out <- plogis(-3.5 + 0.3 * x + 0.2 * z + effect * treat)
      event <- as.integer(runif(1) < p_out)

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        id = i, period = t, treatment = treat, outcome = event,
        eligible = 1L, x = x, z = z, stringsAsFactors = FALSE
      )

      if (event == 1L) alive <- FALSE
    }
  }

  do.call(rbind, rows[seq_len(idx)])
}

# ---------------------------------------------------------------------------
# dgp_grace: Deterministic switching groups for grace period tests
# ---------------------------------------------------------------------------
dgp_grace <- function(n = 3000, periods = 12, seed = 20260306) {
  set.seed(seed)

  rows <- vector("list", n * periods)
  idx <- 0L

  for (i in seq_len(n)) {
    x <- rnorm(1)
    # Assign switching behavior group:
    #   15% switch at period 1, 10% at period 2, 5% at period 3, 70% never
    u <- runif(1)
    switch_time <- if (u < 0.15) 1L
                   else if (u < 0.25) 2L
                   else if (u < 0.30) 3L
                   else NA_integer_

    treat <- 0L
    alive <- TRUE

    for (t in 0:(periods - 1L)) {
      if (!alive) break

      # Treatment initiation (absorbing)
      if (treat == 0L) {
        p_treat <- plogis(-2 + 0.3 * x)
        if (runif(1) < p_treat) treat <- 1L
      }

      # After treatment starts, some switch OFF (non-absorbing for treated)
      # This is the switching behavior that grace period handles
      if (treat == 1L && !is.na(switch_time)) {
        # Calculate periods since treatment started
        # Reset to 0 at switch_time periods after starting
        # For simplicity: if treated and t >= switch_time, stop treatment
        if (t >= switch_time) treat <- 0L
      }

      # Outcome: effect = -0.50
      p_out <- plogis(-3.5 + 0.3 * x - 0.50 * treat)
      event <- as.integer(runif(1) < p_out)

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        id = i, period = t, treatment = treat, outcome = event,
        eligible = 1L, x = x, stringsAsFactors = FALSE
      )

      if (event == 1L) alive <- FALSE
    }
  }

  do.call(rbind, rows[seq_len(idx)])
}

# ---------------------------------------------------------------------------
# dgp_rct: Randomized trial (no confounding)
# ---------------------------------------------------------------------------
dgp_rct <- function(n = 5000, periods = 10, effect = -0.50,
                     seed = 20260307) {
  set.seed(seed)

  rows <- vector("list", n * periods)
  idx <- 0L

  for (i in seq_len(n)) {
    x <- rnorm(1)
    # Random assignment at baseline: P(T=1) = 0.3 (no confounding)
    treat <- as.integer(runif(1) < 0.3)
    alive <- TRUE

    for (t in 0:(periods - 1L)) {
      if (!alive) break

      # Outcome: affected by x and treatment, but treatment is randomized
      p_out <- plogis(-3.5 + 0.3 * x + effect * treat)
      event <- as.integer(runif(1) < p_out)

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        id = i, period = t, treatment = treat, outcome = event,
        eligible = 1L, x = x, stringsAsFactors = FALSE
      )

      if (event == 1L) alive <- FALSE
    }
  }

  do.call(rbind, rows[seq_len(idx)])
}

# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------
dgp_at <- function(...) dgp_simple(...)
dgp_obs <- function(...) dgp_simple(...)
