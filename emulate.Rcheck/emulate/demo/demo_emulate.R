#!/usr/bin/env Rscript
# ==============================================================================
# emulate Package Demo: Target Trial Emulation
# Mirrors: Stata-Tools/tte/demo/demo_emulate.do
# ==============================================================================
#
# This script demonstrates the full emulate pipeline on the trial_example dataset
# (503 patients, 48,400 person-periods, 11 variables). It runs two analyses
# (ITT and PP) and benchmarks results against known R TrialEmulation and
# Stata tte reference values.
#
# Produces:
#   - Console: protocol overview, ITT + PP results, benchmark comparison
#   - PNG:     cumulative incidence plots (ITT and PP)
#   - PNG:     weight distribution plot (PP)
#   - Excel:   protocol table, ITT report, PP report
#   - CSV:     coefficient comparison table
#
# Usage:
#   Rscript demo/demo_emulate.R               # from emulate package root
#   Rscript inst/scripts/demo_emulate.R       # or via inst path
#
# Requirements: emulate, ggplot2, openxlsx (optional, for Excel export)

# ==============================================================================
# STEP 0: Setup and library loading
# ==============================================================================

cat("==============================================================================\n")
cat("emulate Package Demo: Target Trial Emulation\n")
cat("==============================================================================\n\n")

# Load the package
library(emulate)
library(ggplot2)

# Create output directory for demo artifacts
output_dir <- file.path(getwd(), "demo_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Output directory: %s\n\n", output_dir))

# ==============================================================================
# STEP 1: Data loading and exploration
# ==============================================================================

cat("STEP 1: Loading trial_example dataset\n")
cat(strrep("-", 70), "\n")

# Load bundled dataset
f <- system.file("extdata", "trial_example.csv", package = "emulate")
if (f == "" || !file.exists(f)) {
  # Development fallback
  f <- normalizePath("inst/extdata/trial_example.csv", mustWork = FALSE)
}
if (!file.exists(f)) {
  stop("trial_example.csv not found. Is emulate installed?")
}

d <- read.csv(f)

cat(sprintf("  File: %s\n", f))
cat(sprintf("  Dimensions: %d rows x %d columns\n", nrow(d), ncol(d)))
cat(sprintf("  Unique patients: %d\n", length(unique(d$id))))
cat(sprintf("  Period range: %d to %d\n", min(d$period), max(d$period)))
cat(sprintf("  Events: %d (%.2f%%)\n",
            sum(d$outcome == 1), 100 * mean(d$outcome)))
cat(sprintf("  Treatment prevalence: %.1f%%\n", 100 * mean(d$treatment)))

cat("\nData structure:\n")
str(d, give.attr = FALSE)

cat("\nFirst 6 rows:\n")
print(head(d))
cat("\n")

# ==============================================================================
# STEP 2: Protocol specification (Hernan 7-component framework)
# ==============================================================================
# The target trial protocol should be defined BEFORE any analysis,
# following Hernan & Robins (2016). This formalizes the causal question.

cat("STEP 2: Target Trial Protocol\n")
cat(strrep("-", 70), "\n")

protocol <- emulate_protocol(
  eligibility    = "Eligible at period start (eligible == 1); no prior outcome",
  treatment      = "Initiate treatment vs. do not initiate treatment",
  assignment     = "At each eligible period, based on observed treatment decision",
  followup_start = "Start of the period when eligibility criteria are met",
  outcome        = "Binary outcome event (outcome == 1)",
  causal_contrast = "Intention-to-treat (ITT) and per-protocol (PP)",
  analysis       = "Pooled logistic regression with robust SE, clustered by id",
  format         = "display"
)

# Export protocol to Excel (if openxlsx available)
if (requireNamespace("openxlsx", quietly = TRUE)) {
  protocol_file <- file.path(output_dir, "protocol.xlsx")
  emulate_protocol(
    eligibility    = "Eligible at period start (eligible == 1); no prior outcome",
    treatment      = "Initiate treatment vs. do not initiate treatment",
    assignment     = "At each eligible period, based on observed treatment decision",
    followup_start = "Start of the period when eligibility criteria are met",
    outcome        = "Binary outcome event (outcome == 1)",
    causal_contrast = "Intention-to-treat (ITT) and per-protocol (PP)",
    analysis       = "Pooled logistic regression with robust SE, clustered by id",
    format         = "excel",
    export         = protocol_file,
    title          = "Target Trial Protocol: trial_example"
  )
  cat(sprintf("\nProtocol exported to: %s\n", protocol_file))
}
cat("\n")

# ==============================================================================
# ANALYSIS 1: Intention-to-Treat (ITT)
# ==============================================================================
# ITT ignores treatment switching -- everyone is analyzed as initially
# assigned. No IP weights needed. This is the simplest and most robust
# analysis, but diluted toward the null when there is treatment switching.

cat("==============================================================================\n")
cat("ANALYSIS 1: Intention-to-Treat (ITT)\n")
cat("==============================================================================\n\n")

# Step 1: Prepare -- map variables, validate structure, set estimand
obj_itt <- emulate_prepare(d, id = "id", period = "period",
                        treatment = "treatment",
                        outcome = "outcome",
                        eligible = "eligible",
                        covariates = c("catvarA", "catvarB",
                                       "nvarA", "nvarB", "nvarC"),
                        estimand = "ITT")

# Step 2: Validate -- run 10 data quality checks
cat("\nRunning data validation...\n")
obj_itt <- emulate_validate(obj_itt)

# Step 3: Expand -- create sequential emulated trials
# Each eligible time point becomes a new trial. Patients are enrolled
# into the trial at their eligible period with their baseline treatment.
cat("\nExpanding into sequential trials...\n")
obj_itt <- emulate_expand(obj_itt)

# Step 4: Weight -- for ITT, all weights are 1 (no reweighting needed)
cat("\nComputing weights (ITT: all = 1)...\n")
obj_itt <- emulate_weight(obj_itt)

# Step 5: Fit -- pooled logistic regression with clustered SEs
# outcome_cov: covariates in the outcome model (confounders)
# followup_spec: functional form for follow-up time (quadratic = time + time^2)
# trial_period_spec: functional form for trial period (quadratic)
cat("\nFitting outcome model...\n")
obj_itt <- emulate_fit(obj_itt,
                    outcome_cov = c("catvarA", "catvarB",
                                    "nvarA", "nvarB", "nvarC"),
                    followup_spec = "quadratic",
                    trial_period_spec = "quadratic")

# Step 6: Predict -- G-formula standardized marginal predictions
# Computes cumulative incidence by averaging individual-level survival
# curves across the reference population. Monte Carlo CIs from MVN draws.
cat("\nComputing marginal predictions with MC confidence intervals...\n")
obj_itt <- emulate_predict(obj_itt, times = 0:8, type = "cum_inc",
                        samples = 200, seed = 12345, difference = TRUE)

# Step 7: Plot -- cumulative incidence curves
cat("\nGenerating cumulative incidence plot...\n")
p_itt <- emulate_plot(obj_itt, type = "cumhaz", ci = TRUE,
                   title = "Cumulative Incidence (ITT)")
p_itt <- p_itt + theme_minimal(base_size = 12)

# Save plot
itt_plot_file <- file.path(output_dir, "cumulative_incidence_itt.png")
ggsave(itt_plot_file, p_itt, width = 8, height = 5, dpi = 150)
cat(sprintf("  Plot saved: %s\n", itt_plot_file))

# Step 8: Report -- export to Excel
if (requireNamespace("openxlsx", quietly = TRUE)) {
  itt_report_file <- file.path(output_dir, "emulate_report_itt.xlsx")
  emulate_report(obj_itt, format = "excel", export = itt_report_file,
             predictions = TRUE)
  cat(sprintf("  Report saved: %s\n", itt_report_file))
}

# Also export to CSV
itt_csv_file <- file.path(output_dir, "emulate_report_itt.csv")
emulate_report(obj_itt, format = "csv", export = itt_csv_file)
cat(sprintf("  CSV saved: %s\n\n", itt_csv_file))

# Store ITT results for comparison
itt_coef <- obj_itt$model$b_treat
itt_se   <- obj_itt$model$se_treat

# ==============================================================================
# ANALYSIS 2: Per-Protocol (PP) with IPTW
# ==============================================================================
# PP censors treatment switchers and reweights via inverse probability of
# treatment weights (IPTW) to adjust for the informative censoring induced
# by the clone-censor-weight approach. This gives an unbiased estimate of
# the per-protocol effect.

cat("==============================================================================\n")
cat("ANALYSIS 2: Per-Protocol (PP) with IPTW\n")
cat("==============================================================================\n\n")

# Prepare with PP estimand
obj_pp <- emulate_prepare(d, id = "id", period = "period",
                       treatment = "treatment",
                       outcome = "outcome",
                       eligible = "eligible",
                       covariates = c("catvarA", "catvarB",
                                      "nvarA", "nvarB", "nvarC"),
                       estimand = "PP")

# Validate
obj_pp <- emulate_validate(obj_pp)

# Expand -- for PP, this clones each person into two arms (treat/control)
# and applies artificial censoring when they deviate from their assigned arm
cat("\nExpanding with clone-censor-weight...\n")
obj_pp <- emulate_expand(obj_pp)

# Weight -- compute stabilized IPTW with truncation at 1st/99th percentiles
# switch_d_cov: covariates for switch denominator (full model)
# switch_n_cov: covariates for switch numerator (reduced model)
# truncate: winsorize extreme weights at these percentiles
cat("\nComputing stabilized IPTW (truncated at 1/99 percentiles)...\n")
obj_pp <- emulate_weight(obj_pp,
                      switch_d_cov = c("catvarA", "catvarB",
                                       "nvarA", "nvarB", "nvarC"),
                      switch_n_cov = c("catvarA", "nvarA"),
                      stabilized = TRUE,
                      truncate = c(1, 99),
                      quiet = TRUE)

# Diagnose -- weight diagnostics (ESS, distribution, balance)
cat("\nWeight diagnostics:\n")
obj_pp <- emulate_diagnose(obj_pp,
                        balance_covariates = c("nvarA", "nvarB", "nvarC"))

# Weight distribution plot
cat("\nGenerating weight distribution plot...\n")
p_wt <- emulate_plot(obj_pp, type = "weights",
                  title = "IPTW Distribution by Arm")
p_wt <- p_wt + theme_minimal(base_size = 12)

wt_plot_file <- file.path(output_dir, "weight_distribution.png")
ggsave(wt_plot_file, p_wt, width = 8, height = 5, dpi = 150)
cat(sprintf("  Plot saved: %s\n", wt_plot_file))

# Fit PP outcome model
cat("\nFitting PP outcome model...\n")
obj_pp <- emulate_fit(obj_pp,
                   outcome_cov = c("catvarA", "catvarB",
                                   "nvarA", "nvarB", "nvarC"),
                   followup_spec = "quadratic",
                   trial_period_spec = "quadratic")

# Predict
cat("\nComputing PP marginal predictions...\n")
obj_pp <- emulate_predict(obj_pp, times = 0:8, type = "cum_inc",
                       samples = 200, seed = 12345, difference = TRUE)

# PP cumulative incidence plot
p_pp <- emulate_plot(obj_pp, type = "cumhaz", ci = TRUE,
                  title = "Cumulative Incidence (Per-Protocol)")
p_pp <- p_pp + theme_minimal(base_size = 12)

pp_plot_file <- file.path(output_dir, "cumulative_incidence_pp.png")
ggsave(pp_plot_file, p_pp, width = 8, height = 5, dpi = 150)
cat(sprintf("  Plot saved: %s\n", pp_plot_file))

# Report
if (requireNamespace("openxlsx", quietly = TRUE)) {
  pp_report_file <- file.path(output_dir, "emulate_report_pp.xlsx")
  emulate_report(obj_pp, format = "excel", export = pp_report_file,
             predictions = TRUE)
  cat(sprintf("  Report saved: %s\n", pp_report_file))
}

pp_csv_file <- file.path(output_dir, "emulate_report_pp.csv")
emulate_report(obj_pp, format = "csv", export = pp_csv_file)
cat(sprintf("  CSV saved: %s\n\n", pp_csv_file))

# Store PP results
pp_coef <- obj_pp$model$b_treat
pp_se   <- obj_pp$model$se_treat

# ==============================================================================
# STEP 9: Summary comparison of ITT vs PP
# ==============================================================================

cat("==============================================================================\n")
cat("SUMMARY: ITT vs PP Comparison\n")
cat("==============================================================================\n\n")

cat(sprintf("%-25s %12s %12s\n", "", "ITT", "PP"))
cat(strrep("-", 50), "\n")
cat(sprintf("%-25s %12.4f %12.4f\n", "Treatment coefficient", itt_coef, pp_coef))
cat(sprintf("%-25s %12.4f %12.4f\n", "Robust SE", itt_se, pp_se))
cat(sprintf("%-25s %12.4f %12.4f\n", "Odds ratio", exp(itt_coef), exp(pp_coef)))
cat(sprintf("%-25s %12.4f %12.4f\n", "95% CI lower", exp(itt_coef - 1.96 * itt_se),
            exp(pp_coef - 1.96 * pp_se)))
cat(sprintf("%-25s %12.4f %12.4f\n", "95% CI upper", exp(itt_coef + 1.96 * itt_se),
            exp(pp_coef + 1.96 * pp_se)))
cat(sprintf("%-25s %12.4f %12.4f\n", "p-value",
            2 * pnorm(-abs(itt_coef / itt_se)),
            2 * pnorm(-abs(pp_coef / pp_se))))
cat(strrep("-", 50), "\n\n")

# Theory: PP effect should be further from null than ITT because treatment
# switching dilutes the ITT estimate toward the null.
cat("Note: PP coefficient magnitude should be >= ITT (treatment switching dilutes ITT).\n")
cat(sprintf("  |PP|/|ITT| ratio: %.2f\n\n", abs(pp_coef) / abs(itt_coef)))

# ==============================================================================
# STEP 10: Benchmark comparison to known reference values
# ==============================================================================

cat("==============================================================================\n")
cat("BENCHMARK COMPARISON\n")
cat("==============================================================================\n\n")

# R TrialEmulation reference (Maringe et al. 2024)
r_itt_coef <- -0.2829
r_itt_se   <- 0.3138
r_pp_coef  <- -0.4143
r_pp_se    <- 0.4152

# Stata tte reference (from crossval_emulate_vs_r.do)
stata_itt_coef <- -0.282
stata_itt_se   <- 0.312
stata_pp_coef  <- -0.420

cat(sprintf("%-25s %12s %12s %12s %12s\n",
            "", "emulate", "R TrialEmul", "Stata tte", "emulate vs R (%)"))
cat(strrep("-", 75), "\n")

# ITT comparison
itt_rdiff <- 100 * abs(itt_coef - r_itt_coef) / abs(r_itt_coef)
cat(sprintf("%-25s %12.4f %12.4f %12.4f %11.1f%%\n",
            "ITT coefficient", itt_coef, r_itt_coef, stata_itt_coef, itt_rdiff))

itt_se_rdiff <- 100 * abs(itt_se - r_itt_se) / r_itt_se
cat(sprintf("%-25s %12.4f %12.4f %12.4f %11.1f%%\n",
            "ITT SE", itt_se, r_itt_se, stata_itt_se, itt_se_rdiff))

# PP comparison
pp_rdiff <- 100 * abs(pp_coef - r_pp_coef) / abs(r_pp_coef)
cat(sprintf("%-25s %12.4f %12.4f %12.4f %11.1f%%\n",
            "PP coefficient", pp_coef, r_pp_coef, stata_pp_coef, pp_rdiff))

pp_se_rdiff <- 100 * abs(pp_se - r_pp_se) / r_pp_se
cat(sprintf("%-25s %12.4f %12.4f %12s %11.1f%%\n",
            "PP SE", pp_se, r_pp_se, "--", pp_se_rdiff))

cat(strrep("-", 75), "\n")

# Status
itt_pass <- itt_rdiff < 10
pp_pass  <- pp_rdiff < 20
cat(sprintf("\nITT benchmark (10%% tolerance): %s\n",
            if (itt_pass) "PASS" else "FAIL"))
cat(sprintf("PP  benchmark (20%% tolerance): %s\n",
            if (pp_pass) "PASS" else "FAIL"))

# Save comparison to CSV
comparison_df <- data.frame(
  Estimand = c("ITT", "ITT", "PP", "PP"),
  Metric = c("Coefficient", "SE", "Coefficient", "SE"),
  emulate = c(itt_coef, itt_se, pp_coef, pp_se),
  R_TrialEmulation = c(r_itt_coef, r_itt_se, r_pp_coef, r_pp_se),
  Stata_tte = c(stata_itt_coef, stata_itt_se, stata_pp_coef, NA),
  Rel_diff_pct = c(itt_rdiff, itt_se_rdiff, pp_rdiff, pp_se_rdiff)
)
comp_file <- file.path(output_dir, "benchmark_comparison.csv")
write.csv(comparison_df, comp_file, row.names = FALSE)
cat(sprintf("\nComparison table saved: %s\n", comp_file))

# ==============================================================================
# Cleanup and summary
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("DEMO COMPLETE\n")
cat("==============================================================================\n")
cat(sprintf("\nOutput files in: %s\n", output_dir))
cat(sprintf("  %s\n", paste(list.files(output_dir), collapse = "\n  ")))
cat("\nKnown algorithmic differences vs R TrialEmulation:\n")
cat("  1. SE: emulate uses sandwich::vcovCL(HC1, cadjust=TRUE) vs R's HC1+G/(G-1)\n")
cat("  2. Weight strata: emulate uses 2 (by arm) vs R's 4 (by arm x lag_treat)\n")
cat("  3. Splines: both use Harrell RCS (identical formula)\n")
cat("  4. Truncation: both use quantile type=2 (matches Stata _pctile)\n")
