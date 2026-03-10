#!/usr/bin/env Rscript
# ==============================================================================
# emulate Validation Suite Runner
# Mirrors: Stata-Tools/tte/qa/run_all_validations.do
# ==============================================================================
#
# Master runner for all emulate validation test modules. Runs devtools::test()
# with filter patterns for each V01-V17 module, times each one, and produces
# a summary report.
#
# Validations:
#   V01. R TrialEmulation cross-validation (trial_example.csv)
#   V02. NHEFS-style smoking cessation (synthetic DGP)
#   V03. Clone-censor-weight / immortal-time bias
#   V04. G-formula / time-varying confounding
#   V05. Known DGP Monte Carlo
#   V06. Null effect & reproducibility
#   V07. IPCW / informative censoring
#   V08. Grace period correctness
#   V09. Edge cases & strict validation
#   V10. As-treated (AT) estimand
#   V11. RCT benchmarks
#   V12. Sensitivity sweep & stress tests
#   V13. Cox model ground truth
#   V14. emulate_expand options
#   V15. emulate_predict options
#   V16. emulate_diagnose and emulate_report
#   V17. Pipeline guards
#
# Additional cross-validation:
#   crossval: emulate vs Stata tte reference values
#
# Usage:
#   Rscript inst/scripts/run_all_validations.R             # run all
#   Rscript inst/scripts/run_all_validations.R 1 5 13      # run V01, V05, V13
#   Rscript inst/scripts/run_all_validations.R crossval    # run cross-validation only
#
# Output: VALIDATION_REPORT.md in the current working directory

# ==============================================================================
# Parse command-line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Validation module definitions
all_modules <- list(
  list(num = 1,  filter = "v01",           name = "R TrialEmulation"),
  list(num = 2,  filter = "v02",           name = "NHEFS"),
  list(num = 3,  filter = "v03",           name = "CCW / Immortal-Time Bias"),
  list(num = 4,  filter = "v04",           name = "G-Formula / Time-Varying"),
  list(num = 5,  filter = "v05",           name = "Known DGP Monte Carlo"),
  list(num = 6,  filter = "v06",           name = "Null Effect & Reproducibility"),
  list(num = 7,  filter = "v07",           name = "IPCW / Informative Censoring"),
  list(num = 8,  filter = "v08",           name = "Grace Period"),
  list(num = 9,  filter = "v09",           name = "Edge Cases"),
  list(num = 10, filter = "v10",           name = "As-Treated (AT)"),
  list(num = 11, filter = "v11",           name = "RCT Benchmarks"),
  list(num = 12, filter = "v12",           name = "Sensitivity & Stress"),
  list(num = 13, filter = "v13",           name = "Cox Model"),
  list(num = 14, filter = "v14",           name = "Expand Options"),
  list(num = 15, filter = "v15",           name = "Predict Options"),
  list(num = 16, filter = "v16",           name = "Diagnose & Report"),
  list(num = 17, filter = "v17",           name = "Pipeline Guards")
)

# Extra module for cross-validation
crossval_module <- list(num = 99, filter = "crossval", name = "Cross-Val (emulate vs Stata)")

# Determine which modules to run
if (length(args) == 0) {
  # Run all
  run_modules <- all_modules
  cat("Running ALL validation modules (V01-V17)\n")
} else {
  # Parse selective arguments
  run_modules <- list()
  for (arg in args) {
    if (tolower(arg) == "crossval") {
      run_modules <- c(run_modules, list(crossval_module))
    } else if (tolower(arg) == "all") {
      run_modules <- c(all_modules, list(crossval_module))
      break
    } else {
      num <- suppressWarnings(as.integer(arg))
      if (!is.na(num) && num >= 1 && num <= length(all_modules)) {
        run_modules <- c(run_modules, list(all_modules[[num]]))
      } else {
        cat(sprintf("WARNING: Unknown module '%s', skipping\n", arg))
      }
    }
  }
  nums <- vapply(run_modules, function(m) m$num, integer(1))
  cat(sprintf("Running selective modules: %s\n",
              paste(ifelse(nums == 99, "crossval",
                           sprintf("V%02d", nums)), collapse = ", ")))
}

if (length(run_modules) == 0) {
  cat("No valid modules specified. Exiting.\n")
  quit(status = 1)
}

# ==============================================================================
# Setup
# ==============================================================================

cat("\n")
cat(strrep("=", 72), "\n")
cat("emulate VALIDATION SUITE\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("R version: %s\n", R.version.string))
cat(strrep("=", 72), "\n\n")

# Check that emulate is loadable
tryCatch({
  library(emulate)
  cat(sprintf("emulate version: %s\n\n",
              as.character(packageVersion("emulate"))))
}, error = function(e) {
  cat("ERROR: emulate package not found. Install or load with devtools::load_all().\n")
  quit(status = 1)
})

# Check devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
  cat("ERROR: devtools package required. Install with install.packages('devtools').\n")
  quit(status = 1)
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  cat("ERROR: testthat package required. Install with install.packages('testthat').\n")
  quit(status = 1)
}

# Detect package root (look for DESCRIPTION file)
pkg_root <- NULL
candidates <- c(
  getwd(),
  normalizePath(file.path(getwd(), ".."), mustWork = FALSE),
  normalizePath(file.path(getwd(), "../.."), mustWork = FALSE),
  Sys.getenv("EMULATE_PKG_DIR", unset = "")
)
for (cand in candidates) {
  if (nchar(cand) > 0 && file.exists(file.path(cand, "DESCRIPTION"))) {
    desc <- readLines(file.path(cand, "DESCRIPTION"), n = 3)
    if (any(grepl("^Package:\\s*emulate", desc))) {
      pkg_root <- cand
      break
    }
  }
}
if (is.null(pkg_root)) {
  cat("WARNING: Could not find emulate package root. Using installed tests.\n")
  cat("         Set EMULATE_PKG_DIR or run from the package directory.\n\n")
}

# ==============================================================================
# Run validation modules
# ==============================================================================

results <- data.frame(
  Module   = character(),
  Name     = character(),
  Tests    = integer(),
  Passed   = integer(),
  Failed   = integer(),
  Skipped  = integer(),
  Warnings = integer(),
  Time_sec = numeric(),
  Status   = character(),
  stringsAsFactors = FALSE
)

total_start <- proc.time()

for (mod in run_modules) {
  label <- if (mod$num == 99) "crossval" else sprintf("V%02d", mod$num)
  cat(strrep("-", 72), "\n")
  cat(sprintf("Running %s: %s\n", label, mod$name))
  cat(strrep("-", 72), "\n")

  t0 <- proc.time()

  # Run tests using devtools::test with filter
  test_result <- tryCatch({
    if (!is.null(pkg_root)) {
      devtools::test(pkg_root, filter = mod$filter, reporter = "summary",
                     stop_on_failure = FALSE)
    } else {
      # Fallback: use testthat directly on installed package
      testthat::test_package("emulate", filter = mod$filter,
                              reporter = "summary",
                              stop_on_failure = FALSE)
    }
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    NULL
  })

  elapsed <- (proc.time() - t0)["elapsed"]

  # Parse results
  n_tests   <- 0L
  n_passed  <- 0L
  n_failed  <- 0L
  n_skipped <- 0L
  n_warns   <- 0L
  status    <- "ERROR"

  if (!is.null(test_result)) {
    # testthat returns a list of test results
    test_df <- as.data.frame(test_result)
    if (nrow(test_df) > 0) {
      n_tests   <- nrow(test_df)
      n_passed  <- sum(test_df$failed == 0 & test_df$skipped == FALSE, na.rm = TRUE)
      n_failed  <- sum(test_df$failed > 0, na.rm = TRUE)
      n_skipped <- sum(test_df$skipped == TRUE, na.rm = TRUE)
      n_warns   <- sum(test_df$warning > 0, na.rm = TRUE)

      if (n_failed == 0 && n_tests > 0) {
        status <- "PASS"
      } else if (n_failed > 0) {
        status <- "FAIL"
      } else {
        status <- "SKIP"
      }
    } else {
      status <- "NO TESTS"
    }
  }

  cat(sprintf("  Result: %s (%d tests, %d passed, %d failed, %d skipped) [%.1fs]\n\n",
              status, n_tests, n_passed, n_failed, n_skipped, elapsed))

  results <- rbind(results, data.frame(
    Module   = label,
    Name     = mod$name,
    Tests    = n_tests,
    Passed   = n_passed,
    Failed   = n_failed,
    Skipped  = n_skipped,
    Warnings = n_warns,
    Time_sec = round(elapsed, 1),
    Status   = status,
    stringsAsFactors = FALSE
  ))
}

total_elapsed <- (proc.time() - total_start)["elapsed"]

# ==============================================================================
# Summary table
# ==============================================================================

cat("\n")
cat(strrep("=", 72), "\n")
cat("VALIDATION SUMMARY\n")
cat(strrep("=", 72), "\n\n")

# Print summary table
cat(sprintf("%-10s %-30s %5s %6s %6s %7s %8s %7s\n",
            "Module", "Name", "Tests", "Pass", "Fail", "Skip", "Time(s)", "Status"))
cat(strrep("-", 90), "\n")

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-10s %-30s %5d %6d %6d %7d %8.1f %7s\n",
              r$Module, r$Name, r$Tests, r$Passed, r$Failed,
              r$Skipped, r$Time_sec, r$Status))
}

cat(strrep("-", 90), "\n")
cat(sprintf("%-10s %-30s %5d %6d %6d %7d %8.1f\n",
            "TOTAL", "",
            sum(results$Tests), sum(results$Passed), sum(results$Failed),
            sum(results$Skipped), total_elapsed))

n_modules_pass <- sum(results$Status == "PASS")
n_modules_fail <- sum(results$Status == "FAIL")
n_modules_err  <- sum(results$Status == "ERROR")

cat("\n")
if (n_modules_fail == 0 && n_modules_err == 0) {
  cat("OVERALL: ALL VALIDATIONS PASSED\n")
} else {
  if (n_modules_fail > 0)
    cat(sprintf("OVERALL: %d MODULE(S) FAILED\n", n_modules_fail))
  if (n_modules_err > 0)
    cat(sprintf("         %d MODULE(S) HAD ERRORS\n", n_modules_err))
}

# ==============================================================================
# Generate VALIDATION_REPORT.md
# ==============================================================================

report_file <- file.path(getwd(), "VALIDATION_REPORT.md")

lines <- c(
  "# emulate Validation Report",
  "",
  sprintf("**Date:** %s", Sys.time()),
  sprintf("**R version:** %s", R.version.string),
  sprintf("**emulate version:** %s", as.character(packageVersion("emulate"))),
  sprintf("**Platform:** %s", R.version$platform),
  "",
  "## Summary",
  "",
  sprintf("- **Modules run:** %d", nrow(results)),
  sprintf("- **Total tests:** %d", sum(results$Tests)),
  sprintf("- **Passed:** %d", sum(results$Passed)),
  sprintf("- **Failed:** %d", sum(results$Failed)),
  sprintf("- **Skipped:** %d", sum(results$Skipped)),
  sprintf("- **Total time:** %.1f seconds", total_elapsed),
  sprintf("- **Overall status:** %s",
          if (n_modules_fail == 0 && n_modules_err == 0) "PASS" else "FAIL"),
  "",
  "## Results by Module",
  "",
  "| Module | Name | Tests | Pass | Fail | Skip | Time (s) | Status |",
  "|--------|------|------:|-----:|-----:|-----:|---------:|--------|"
)

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  lines <- c(lines, sprintf("| %s | %s | %d | %d | %d | %d | %.1f | %s |",
                             r$Module, r$Name, r$Tests, r$Passed, r$Failed,
                             r$Skipped, r$Time_sec, r$Status))
}

lines <- c(lines, "",
  "## Known Algorithmic Differences",
  "",
  "When comparing emulate to Stata tte and R TrialEmulation:",
  "",
  "1. **Standard errors:** emulate uses `sandwich::vcovCL(type='HC1', cadjust=TRUE)`",
  "   which applies both HC1 and G/(G-1) corrections. Stata uses `vce(cluster)`",
  "   which applies G/(G-1) only. R TrialEmulation uses HC1 + G/(G-1).",
  "2. **Weight model strata:** emulate and Stata use 2 strata (by arm).",
  "   R TrialEmulation uses 4 strata (by arm x lagged treatment).",
  "3. **Truncation percentiles:** All three use equivalent quantile methods",
  "   (emulate: `quantile(type=2)`, Stata: `_pctile`).",
  "4. **Spline basis:** emulate and Stata use identical Harrell RCS.",
  "   R TrialEmulation uses `splines::ns()` (different basis, different coefficients).",
  "",
  sprintf("---\n*Generated by `run_all_validations.R` on %s*", Sys.time())
)

writeLines(lines, report_file)
cat(sprintf("\nValidation report written to: %s\n", report_file))

# Exit with appropriate code
if (n_modules_fail > 0 || n_modules_err > 0) {
  quit(status = 1)
} else {
  quit(status = 0)
}
