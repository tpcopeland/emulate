# TOST (Two One-Sided Tests) equivalence testing helpers
# No dependency on TOSTER package — implements z-based TOST directly

#' Z-based TOST for equivalence of two estimates
#'
#' Tests H0: |theta1 - theta2| >= delta vs H1: |theta1 - theta2| < delta.
#' Rejection supports the conclusion that the two estimates are equivalent
#' within the specified margin.
#'
#' @param est1 First estimate
#' @param est2 Second estimate
#' @param se1 Standard error of first estimate
#' @param se2 Standard error of second estimate
#' @param delta Equivalence margin (symmetric)
#' @param alpha Significance level (default 0.05)
#' @return List with: equivalent (logical), p_tost, diff, se_diff, z_lo, z_hi
tost_z <- function(est1, est2, se1, se2, delta, alpha = 0.05) {
  diff <- est1 - est2
  se_diff <- sqrt(se1^2 + se2^2)

  # Two one-sided tests
  z_lo <- (diff + delta) / se_diff   # H0: diff <= -delta, reject if large
  z_hi <- (diff - delta) / se_diff   # H0: diff >= +delta, reject if small
  p_lo <- pnorm(z_lo, lower.tail = FALSE)  # upper-tail p-value
  p_hi <- pnorm(z_hi, lower.tail = TRUE)   # lower-tail p-value
  p_tost <- max(p_lo, p_hi)

  list(
    equivalent = (p_tost < alpha),
    p_tost = p_tost,
    diff = diff,
    se_diff = se_diff,
    delta = delta,
    z_lo = z_lo,
    z_hi = z_hi
  )
}

#' MC-based TOST for equivalence using replicated differences
#'
#' Given a vector of paired differences (e.g., from MC replications),
#' tests equivalence using a one-sample TOST on the mean difference.
#'
#' @param diffs Numeric vector of paired differences
#' @param delta Equivalence margin
#' @param alpha Significance level (default 0.05)
#' @return List with: equivalent (logical), p_tost, mean_diff, se_diff, n
tost_mc <- function(diffs, delta, alpha = 0.05) {
  diffs <- diffs[!is.na(diffs)]
  n <- length(diffs)
  if (n < 3) return(list(equivalent = FALSE, p_tost = 1, mean_diff = NA,
                          se_diff = NA, n = n))
  mean_diff <- mean(diffs)
  se_diff <- sd(diffs) / sqrt(n)

  # TOST: test |mean_diff| < delta
  t_lo <- (mean_diff + delta) / se_diff
  t_hi <- (mean_diff - delta) / se_diff
  p_lo <- pt(t_lo, df = n - 1, lower.tail = FALSE)  # upper-tail
  p_hi <- pt(t_hi, df = n - 1, lower.tail = TRUE)    # lower-tail
  p_tost <- max(p_lo, p_hi)

  list(
    equivalent = (p_tost < alpha),
    p_tost = p_tost,
    mean_diff = mean_diff,
    se_diff = se_diff,
    n = n,
    delta = delta
  )
}

#' Load the NHEFS person-period dataset
#' @return data.frame
load_nhefs <- function() {
  f <- system.file("extdata", "nhefs_personperiod.csv", package = "emulate")
  if (f == "") {
    f <- normalizePath("../../inst/extdata/nhefs_personperiod.csv",
                       mustWork = FALSE)
  }
  if (!file.exists(f)) {
    stop("nhefs_personperiod.csv not found. Is the package installed?")
  }
  read.csv(f)
}

#' Load the golden reference DGP dataset
#' @return data.frame
load_golden_dgp <- function() {
  f <- system.file("extdata", "known_dgp_golden.csv", package = "emulate")
  if (f == "") {
    f <- normalizePath("../../inst/extdata/known_dgp_golden.csv",
                       mustWork = FALSE)
  }
  if (!file.exists(f)) {
    stop("known_dgp_golden.csv not found. Is the package installed?")
  }
  read.csv(f)
}
