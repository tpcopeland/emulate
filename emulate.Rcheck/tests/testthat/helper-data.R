# Test data generators for emulate tests

#' Create a minimal valid person-period dataset for testing
#' @param n_ids Number of individuals
#' @param n_periods Number of periods per individual
#' @param treat_prob Probability of treatment at each period
#' @param event_prob Probability of outcome event at each period
#' @param seed Random seed
#' @return data.frame in person-period format
make_test_data <- function(n_ids = 50, n_periods = 10, treat_prob = 0.3,
                           event_prob = 0.02, seed = 12345) {
  set.seed(seed)

  rows <- list()
  for (i in seq_len(n_ids)) {
    treat_status <- 0
    event_happened <- FALSE
    for (t in seq_len(n_periods)) {
      if (event_happened) break

      # Treatment can start (absorbing for simplicity)
      if (treat_status == 0 && runif(1) < treat_prob) {
        treat_status <- 1
      }

      # Outcome event
      event <- as.integer(runif(1) < event_prob)
      if (event == 1) event_happened <- TRUE

      rows[[length(rows) + 1]] <- data.frame(
        id = i,
        period = t - 1L,
        treatment = treat_status,
        outcome = event,
        eligible = 1L,
        nvarA = rnorm(1, 50, 10),
        nvarB = rnorm(1, 30, 5),
        catvarA = sample(0:2, 1)
      )
    }
  }

  do.call(rbind, rows)
}

#' Load the trial_example dataset shipped with the package
#' @return data.frame
load_trial_example <- function() {
  f <- system.file("extdata", "trial_example.csv", package = "emulate")
  if (f == "") {
    # Fallback for when package is not installed (development)
    f <- file.path(find.package("emulate", quiet = TRUE), "inst", "extdata",
                   "trial_example.csv")
    if (!file.exists(f)) {
      f <- normalizePath("../../inst/extdata/trial_example.csv",
                         mustWork = FALSE)
    }
  }
  if (!file.exists(f)) {
    stop("trial_example.csv not found. Is the package installed?")
  }
  read.csv(f)
}
