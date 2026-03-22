#!/usr/bin/env Rscript
# =============================================================================
# Cross-Validation: emulate vs TrialEmulation vs Known DGP Ground Truth
# =============================================================================
#
# Compares two R implementations of sequential target trial emulation:
#   - emulate (v0.2.0+): ~/R-Packages/emulate
#   - TrialEmulation (v0.0.4.9, CRAN): Maringe et al. 2024, arXiv:2402.12083
#
# Three datasets with known data-generating processes:
#   A) known_dgp      -- 10K patients, 10 periods, true log-OR = -0.50
#   B) gformula_sim    -- 5K patients, 15 periods, true ART log-OR = -0.80
#   C) ipcw_dgp        -- 5K patients, 10 periods, true log-OR = -0.60
#
# Each dataset is analysed under ITT and/or PP estimands. Coefficients,
# robust SEs, and risk differences are compared across implementations
# and against the true DGP parameter.
#
# Known algorithmic differences:
#   1. Weight model stratification: emulate uses 2 strata (by arm),
#      TrialEmulation uses 4 strata (by arm x lagged treatment)
#   2. SE computation: emulate uses sandwich::vcovCL(HC1 + cadjust),
#      TrialEmulation uses sandwich::vcovCL with HC1
#   3. IPCW: emulate stratifies censoring models by arm, producing
#      non-trivial weights even when numerator = denominator covariates.
#      TrialEmulation produces unit weights in that case.
#
# Output:
#   qa/crossval_results.csv   -- Machine-readable results table
#   Console summary           -- Human-readable comparison + assessment
#
# Usage:
#   cd ~/R-Packages/emulate
#   Rscript qa/crossval_emulate_vs_trialemulation.R
#
# Prerequisites:
#   install.packages("TrialEmulation")
#   devtools::install(".")  # or install emulate from local source
#   # Datasets from Stata-Tools/tte/qa/data/ (read via haven)
# =============================================================================

suppressPackageStartupMessages({
  library(emulate)
  library(TrialEmulation)
  library(haven)
  library(data.table)
})

# --- Paths ---
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_arg) > 0) dirname(normalizePath(script_arg)) else getwd()
}, error = function(e) getwd())

# Dataset location: Stata tte QA data
data_dir <- normalizePath("~/Stata-Tools/tte/qa/data", mustWork = TRUE)
output_csv <- file.path(script_dir, "crossval_results.csv")
te_tempdir <- tempdir()

# --- Helpers ---
results <- data.frame(
  dataset    = character(),
  config     = character(),
  metric     = character(),
  emulate    = numeric(),
  trialemul  = numeric(),
  true_value = numeric(),
  em_true_diff  = numeric(),
  te_true_diff  = numeric(),
  em_te_diff    = numeric(),
  stringsAsFactors = FALSE
)

add_result <- function(dataset, config, metric, em_val, te_val, true_val = NA) {
  results <<- rbind(results, data.frame(
    dataset       = dataset,
    config        = config,
    metric        = metric,
    emulate       = round(em_val, 6),
    trialemul     = round(te_val, 6),
    true_value    = if (is.na(true_val)) NA_real_ else round(true_val, 6),
    em_true_diff  = if (is.na(true_val)) NA_real_ else round(abs(em_val - true_val), 6),
    te_true_diff  = if (is.na(true_val)) NA_real_ else round(abs(te_val - true_val), 6),
    em_te_diff    = round(abs(em_val - te_val), 6),
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# Header
# =============================================================================
cat("\n")
cat(strrep("=", 78), "\n")
cat("CROSS-VALIDATION: emulate vs TrialEmulation vs DGP Ground Truth\n")
cat(strrep("=", 78), "\n")
cat(sprintf("Date:               %s\n", Sys.time()))
cat(sprintf("emulate version:    %s\n", as.character(packageVersion("emulate"))))
cat(sprintf("TrialEmulation:     %s\n", as.character(packageVersion("TrialEmulation"))))
cat(sprintf("R version:          %s\n", R.version.string))
cat(sprintf("Data source:        %s\n", data_dir))
cat("\n")

# =============================================================================
# DATASET A: known_dgp
# True treatment log-OR = -0.50
# 10,000 patients, 10 periods, binary covariate x
# =============================================================================
cat(strrep("-", 78), "\n")
cat("DATASET A: known_dgp (N=10,000, true log-OR = -0.50)\n")
cat(strrep("-", 78), "\n\n")

d_known <- as.data.frame(read_dta(file.path(data_dir, "known_dgp.dta")))
TRUE_KNOWN <- -0.50
cat(sprintf("  Loaded: %d rows, %d IDs, %d periods\n",
            nrow(d_known), length(unique(d_known$id)),
            length(unique(d_known$period))))

# A1: ITT -- emulate
cat("  A1: ITT (emulate)...\n")
obj <- suppressMessages(emulate_prepare(d_known, id = "id", period = "period",
  treatment = "treatment", outcome = "outcome", eligible = "eligible",
  covariates = "x", estimand = "ITT"))
obj <- suppressMessages(emulate_expand(obj))
obj <- suppressMessages(emulate_weight(obj))
obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
em_a1_coef <- obj$model$b_treat
em_a1_se   <- obj$model$se_treat
cat(sprintf("    Coef: %.6f  SE: %.6f\n", em_a1_coef, em_a1_se))

# A1: ITT -- TrialEmulation
cat("  A1: ITT (TrialEmulation)...\n")
te <- suppressWarnings(initiators(
  data = d_known, id = "id", period = "period", treatment = "treatment",
  outcome = "outcome", eligible = "eligible", estimand_type = "ITT",
  outcome_cov = ~ x,
  include_followup_time = ~ followup_time + I(followup_time^2),
  include_trial_period = ~ trial_period + I(trial_period^2),
  data_dir = te_tempdir, quiet = TRUE
))
te_a1_coef <- te$robust$summary$estimate[2]
te_a1_se   <- te$robust$summary$robust_se[2]
cat(sprintf("    Coef: %.6f  SE: %.6f\n", te_a1_coef, te_a1_se))

add_result("known_dgp", "ITT", "Coef", em_a1_coef, te_a1_coef, TRUE_KNOWN)
add_result("known_dgp", "ITT", "SE",   em_a1_se,   te_a1_se)

# A1: ITT predictions
cat("  A1: Predictions...\n")
max_fu <- max(obj$data[["_emulate_followup"]], na.rm = TRUE)
obj <- suppressMessages(emulate_predict(obj, times = 0:min(5, max_fu),
  type = "cum_inc", samples = 100, seed = 12345, difference = TRUE))
em_a_rd <- obj$predictions$diff[nrow(obj$predictions)]

set.seed(12345)
pred <- predict(te, predict_times = 0:5, conf_int_type = "normal")
te_a_rd <- pred$difference$cum_inc_diff[nrow(pred$difference)]
cat(sprintf("    RD(t=5): emulate=%.6f  TrialEmulation=%.6f\n", em_a_rd, te_a_rd))
add_result("known_dgp", "ITT", "RD(t=5)", em_a_rd, te_a_rd)

# A2: PP -- emulate
cat("  A2: PP (emulate)...\n")
obj <- suppressMessages(emulate_prepare(d_known, id = "id", period = "period",
  treatment = "treatment", outcome = "outcome", eligible = "eligible",
  covariates = "x", estimand = "PP"))
obj <- suppressMessages(emulate_expand(obj))
obj <- suppressMessages(emulate_weight(obj, switch_d_cov = "x",
  switch_n_cov = "x", quiet = TRUE))
obj <- suppressMessages(emulate_fit(obj, outcome_cov = "x"))
em_a2_coef <- obj$model$b_treat
em_a2_se   <- obj$model$se_treat
cat(sprintf("    Coef: %.6f  SE: %.6f\n", em_a2_coef, em_a2_se))

# A2: PP -- TrialEmulation
cat("  A2: PP (TrialEmulation)...\n")
te <- suppressWarnings(initiators(
  data = d_known, id = "id", period = "period", treatment = "treatment",
  outcome = "outcome", eligible = "eligible", estimand_type = "PP",
  outcome_cov = ~ x, switch_d_cov = ~ x, switch_n_cov = ~ x,
  include_followup_time = ~ followup_time + I(followup_time^2),
  include_trial_period = ~ trial_period + I(trial_period^2),
  data_dir = te_tempdir, quiet = TRUE
))
te_a2_coef <- te$robust$summary$estimate[2]
te_a2_se   <- te$robust$summary$robust_se[2]
cat(sprintf("    Coef: %.6f  SE: %.6f\n", te_a2_coef, te_a2_se))

add_result("known_dgp", "PP", "Coef", em_a2_coef, te_a2_coef, TRUE_KNOWN)
add_result("known_dgp", "PP", "SE",   em_a2_se,   te_a2_se)


# =============================================================================
# DATASET B: gformula_simulated
# True ART log-OR = -0.80
# 5,000 patients, 15 periods, time-varying CD4 confounding
# =============================================================================
cat("\n")
cat(strrep("-", 78), "\n")
cat("DATASET B: gformula_simulated (N=5,000, true ART log-OR = -0.80)\n")
cat(strrep("-", 78), "\n\n")

d_gform <- as.data.frame(read_dta(file.path(data_dir, "gformula_simulated.dta")))
TRUE_GFORM <- -0.80
cat(sprintf("  Loaded: %d rows, %d IDs, %d periods\n",
            nrow(d_gform), length(unique(d_gform$id)),
            length(unique(d_gform$period))))

outcome_covs <- c("age_cat", "male", "cd4_std")
switch_covs  <- c("age_cat", "cd4_std")

# B1: ITT -- emulate
cat("  B1: ITT (emulate)...\n")
obj <- suppressMessages(emulate_prepare(d_gform, id = "id", period = "period",
  treatment = "treatment", outcome = "outcome", eligible = "eligible",
  covariates = switch_covs, estimand = "ITT"))
obj <- suppressMessages(emulate_expand(obj))
obj <- suppressMessages(emulate_weight(obj))
obj <- suppressMessages(emulate_fit(obj, outcome_cov = outcome_covs))
em_b1_coef <- obj$model$b_treat
em_b1_se   <- obj$model$se_treat
cat(sprintf("    Coef: %.6f  SE: %.6f\n", em_b1_coef, em_b1_se))

# B1: ITT -- TrialEmulation
cat("  B1: ITT (TrialEmulation)...\n")
te <- suppressWarnings(initiators(
  data = d_gform, id = "id", period = "period", treatment = "treatment",
  outcome = "outcome", eligible = "eligible", estimand_type = "ITT",
  outcome_cov = ~ age_cat + male + cd4_std,
  include_followup_time = ~ followup_time + I(followup_time^2),
  include_trial_period = ~ trial_period + I(trial_period^2),
  data_dir = te_tempdir, quiet = TRUE
))
te_b1_coef <- te$robust$summary$estimate[2]
te_b1_se   <- te$robust$summary$robust_se[2]
cat(sprintf("    Coef: %.6f  SE: %.6f\n", te_b1_coef, te_b1_se))

add_result("gformula", "ITT", "Coef", em_b1_coef, te_b1_coef, TRUE_GFORM)
add_result("gformula", "ITT", "SE",   em_b1_se,   te_b1_se)

# B1: ITT predictions
cat("  B1: Predictions...\n")
max_fu <- max(obj$data[["_emulate_followup"]], na.rm = TRUE)
obj <- suppressMessages(emulate_predict(obj, times = 0:min(8, max_fu),
  type = "cum_inc", samples = 100, seed = 12345, difference = TRUE))
em_b_rd <- obj$predictions$diff[nrow(obj$predictions)]

set.seed(12345)
pred <- predict(te, predict_times = 0:8, conf_int_type = "normal")
te_b_rd <- pred$difference$cum_inc_diff[nrow(pred$difference)]
cat(sprintf("    RD(t=8): emulate=%.6f  TrialEmulation=%.6f\n", em_b_rd, te_b_rd))
add_result("gformula", "ITT", "RD(t=8)", em_b_rd, te_b_rd)

# B2: PP -- emulate
cat("  B2: PP (emulate)...\n")
obj <- suppressMessages(emulate_prepare(d_gform, id = "id", period = "period",
  treatment = "treatment", outcome = "outcome", eligible = "eligible",
  covariates = switch_covs, estimand = "PP"))
obj <- suppressMessages(emulate_expand(obj))
obj <- suppressWarnings(suppressMessages(emulate_weight(obj,
  switch_d_cov = switch_covs, switch_n_cov = switch_covs, quiet = TRUE)))
obj <- suppressMessages(emulate_fit(obj, outcome_cov = outcome_covs))
em_b2_coef <- obj$model$b_treat
em_b2_se   <- obj$model$se_treat
cat(sprintf("    Coef: %.6f  SE: %.6f\n", em_b2_coef, em_b2_se))

# B2: PP -- TrialEmulation
cat("  B2: PP (TrialEmulation)...\n")
te <- suppressWarnings(initiators(
  data = d_gform, id = "id", period = "period", treatment = "treatment",
  outcome = "outcome", eligible = "eligible", estimand_type = "PP",
  outcome_cov = ~ age_cat + male + cd4_std,
  switch_d_cov = ~ age_cat + cd4_std, switch_n_cov = ~ age_cat + cd4_std,
  include_followup_time = ~ followup_time + I(followup_time^2),
  include_trial_period = ~ trial_period + I(trial_period^2),
  data_dir = te_tempdir, quiet = TRUE
))
te_b2_coef <- te$robust$summary$estimate[2]
te_b2_se   <- te$robust$summary$robust_se[2]
cat(sprintf("    Coef: %.6f  SE: %.6f\n", te_b2_coef, te_b2_se))

add_result("gformula", "PP", "Coef", em_b2_coef, te_b2_coef, TRUE_GFORM)
add_result("gformula", "PP", "SE",   em_b2_se,   te_b2_se)


# =============================================================================
# DATASET C: ipcw_dgp
# True treatment log-OR = -0.60, informative censoring
# 5,000 patients, 10 periods, covariates x (binary), z (continuous)
# =============================================================================
cat("\n")
cat(strrep("-", 78), "\n")
cat("DATASET C: ipcw_dgp (N=5,000, true log-OR = -0.60, informative censoring)\n")
cat(strrep("-", 78), "\n\n")

d_ipcw <- as.data.frame(read_dta(file.path(data_dir, "ipcw_dgp.dta")))
TRUE_IPCW <- -0.60
cat(sprintf("  Loaded: %d rows, %d IDs, %d periods\n",
            nrow(d_ipcw), length(unique(d_ipcw$id)),
            length(unique(d_ipcw$period))))

# C1: PP naive (no IPCW) -- emulate
cat("  C1: PP naive (emulate)...\n")
obj <- suppressMessages(emulate_prepare(d_ipcw, id = "id", period = "period",
  treatment = "treatment", outcome = "outcome", eligible = "eligible",
  covariates = c("x", "z"), estimand = "PP"))
obj <- suppressMessages(emulate_expand(obj))
obj <- suppressWarnings(suppressMessages(emulate_weight(obj,
  switch_d_cov = c("x", "z"), switch_n_cov = c("x", "z"), quiet = TRUE)))
obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("x", "z")))
em_c1_coef <- obj$model$b_treat
em_c1_se   <- obj$model$se_treat
cat(sprintf("    Coef: %.6f  SE: %.6f\n", em_c1_coef, em_c1_se))

# C1: PP naive -- TrialEmulation
cat("  C1: PP naive (TrialEmulation)...\n")
te <- suppressWarnings(initiators(
  data = d_ipcw, id = "id", period = "period", treatment = "treatment",
  outcome = "outcome", eligible = "eligible", estimand_type = "PP",
  outcome_cov = ~ x + z, switch_d_cov = ~ x + z, switch_n_cov = ~ x + z,
  include_followup_time = ~ followup_time + I(followup_time^2),
  include_trial_period = ~ trial_period + I(trial_period^2),
  data_dir = te_tempdir, quiet = TRUE
))
te_c1_coef <- te$robust$summary$estimate[2]
te_c1_se   <- te$robust$summary$robust_se[2]
cat(sprintf("    Coef: %.6f  SE: %.6f\n", te_c1_coef, te_c1_se))

add_result("ipcw_dgp", "PP-naive", "Coef", em_c1_coef, te_c1_coef, TRUE_IPCW)
add_result("ipcw_dgp", "PP-naive", "SE",   em_c1_se,   te_c1_se)

# C2: PP with IPCW -- emulate
# Proper stabilization: denominator includes x + z, numerator includes x only
cat("  C2: PP + IPCW (emulate)...\n")
obj <- suppressMessages(emulate_prepare(d_ipcw, id = "id", period = "period",
  treatment = "treatment", outcome = "outcome", eligible = "eligible",
  censor = "censored", covariates = c("x", "z"), estimand = "PP"))
obj <- suppressMessages(emulate_expand(obj))
obj <- suppressWarnings(suppressMessages(emulate_weight(obj,
  switch_d_cov = c("x", "z"), switch_n_cov = c("x", "z"),
  censor_d_cov = c("x", "z"), censor_n_cov = "x", quiet = TRUE)))
obj <- suppressMessages(emulate_fit(obj, outcome_cov = c("x", "z")))
em_c2_coef <- obj$model$b_treat
em_c2_se   <- obj$model$se_treat
cat(sprintf("    Coef: %.6f  SE: %.6f\n", em_c2_coef, em_c2_se))

# C2: PP with IPCW -- TrialEmulation
cat("  C2: PP + IPCW (TrialEmulation)...\n")
te <- suppressWarnings(initiators(
  data = d_ipcw, id = "id", period = "period", treatment = "treatment",
  outcome = "outcome", eligible = "eligible", estimand_type = "PP",
  outcome_cov = ~ x + z, switch_d_cov = ~ x + z, switch_n_cov = ~ x + z,
  use_censor_weights = TRUE, cense = "censored",
  cense_d_cov = ~ x + z, cense_n_cov = ~ x,
  include_followup_time = ~ followup_time + I(followup_time^2),
  include_trial_period = ~ trial_period + I(trial_period^2),
  data_dir = te_tempdir, quiet = TRUE
))
te_c2_coef <- te$robust$summary$estimate[2]
te_c2_se   <- te$robust$summary$robust_se[2]
cat(sprintf("    Coef: %.6f  SE: %.6f\n", te_c2_coef, te_c2_se))

add_result("ipcw_dgp", "PP-IPCW", "Coef", em_c2_coef, te_c2_coef, TRUE_IPCW)
add_result("ipcw_dgp", "PP-IPCW", "SE",   em_c2_se,   te_c2_se)

# C2: IPCW predictions
cat("  C2: Predictions...\n")
max_fu <- max(obj$data[["_emulate_followup"]], na.rm = TRUE)
obj <- suppressMessages(emulate_predict(obj, times = 0:min(5, max_fu),
  type = "cum_inc", samples = 100, seed = 12345, difference = TRUE))
em_c_rd <- obj$predictions$diff[nrow(obj$predictions)]

set.seed(12345)
pred <- predict(te, predict_times = 0:5, conf_int_type = "normal")
te_c_rd <- pred$difference$cum_inc_diff[nrow(pred$difference)]
cat(sprintf("    RD(t=5): emulate=%.6f  TrialEmulation=%.6f\n", em_c_rd, te_c_rd))
add_result("ipcw_dgp", "PP-IPCW", "RD(t=5)", em_c_rd, te_c_rd)


# =============================================================================
# RESULTS TABLE
# =============================================================================
cat("\n")
cat(strrep("=", 94), "\n")
cat("CROSS-VALIDATION RESULTS\n")
cat(strrep("=", 94), "\n\n")

cat(sprintf("%-12s %-10s %-8s %10s %10s %10s %9s %9s %9s\n",
            "Dataset", "Config", "Metric", "emulate", "TrialEmul", "True",
            "|em-true|", "|te-true|", "|em-te|"))
cat(strrep("-", 94), "\n")

prev_dataset <- ""
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  true_str  <- if (is.na(r$true_value))   sprintf("%10s", "---") else sprintf("%10.4f", r$true_value)
  em_d_str  <- if (is.na(r$em_true_diff)) sprintf("%9s", "---")  else sprintf("%9.4f", r$em_true_diff)
  te_d_str  <- if (is.na(r$te_true_diff)) sprintf("%9s", "---")  else sprintf("%9.4f", r$te_true_diff)

  # Separator between datasets
  if (r$dataset != prev_dataset && prev_dataset != "") {
    cat(strrep("-", 94), "\n")
  }
  prev_dataset <- r$dataset

  cat(sprintf("%-12s %-10s %-8s %10.4f %10.4f %s %s %s %9.4f\n",
              r$dataset, r$config, r$metric,
              r$emulate, r$trialemul, true_str,
              em_d_str, te_d_str, r$em_te_diff))
}

# =============================================================================
# ASSESSMENT
# =============================================================================
cat("\n")
cat(strrep("=", 94), "\n")
cat("ASSESSMENT\n")
cat(strrep("=", 94), "\n\n")

# Tolerances for inter-implementation agreement
tol_coef_itt <- 0.10   # ITT coefficients (no weighting differences)
tol_coef_pp  <- 0.20   # PP coefficients (weight stratification differs)
tol_se       <- 0.15   # Standard errors
tol_rd       <- 0.05   # Risk differences
tol_true_itt <- 0.25   # Recovery of true effect (ITT, finite-sample)
tol_true_pp  <- 0.35   # Recovery of true effect (PP, more variable)

n_pass <- 0L; n_note <- 0L; n_fail <- 0L

for (i in seq_len(nrow(results))) {
  r <- results[i, ]

  # Inter-implementation agreement
  if (r$metric == "Coef") {
    tol <- if (grepl("ITT", r$config)) tol_coef_itt else tol_coef_pp
  } else if (r$metric == "SE") {
    tol <- tol_se
  } else {
    tol <- tol_rd
  }

  status <- if (r$em_te_diff <= tol) "PASS" else if (r$em_te_diff <= tol * 1.5) "NOTE" else "FAIL"

  # True value recovery
  true_info <- ""
  if (!is.na(r$true_value) && r$metric == "Coef") {
    t <- if (grepl("ITT", r$config)) tol_true_itt else tol_true_pp
    em_ok <- if (r$em_true_diff <= t) "OK" else "WIDE"
    te_ok <- if (r$te_true_diff <= t) "OK" else "WIDE"
    true_info <- sprintf(" | DGP: em=%s te=%s", em_ok, te_ok)
  }

  cat(sprintf("  %-12s %-10s %-8s  |em-te|=%.4f (tol=%.2f)  %s%s\n",
              r$dataset, r$config, r$metric, r$em_te_diff, tol, status, true_info))

  if (status == "PASS") n_pass <- n_pass + 1L
  else if (status == "NOTE") n_note <- n_note + 1L
  else n_fail <- n_fail + 1L
}

cat(sprintf("\nTOTALS: PASS=%d  NOTE=%d  FAIL=%d\n", n_pass, n_note, n_fail))
if (n_fail == 0) {
  cat("OVERALL: ALL COMPARISONS WITHIN TOLERANCE\n")
} else {
  cat(sprintf("OVERALL: %d COMPARISON(S) EXCEEDED TOLERANCE\n", n_fail))
}

cat("\n")
cat("Algorithmic differences (expected sources of divergence):\n")
cat("  1. Weight stratification: emulate 2 strata (arm), TrialEmulation 4 (arm x lag_treat)\n")
cat("  2. Robust SE: emulate sandwich::vcovCL(HC1+cadjust), TrialEmulation vcovCL(HC1)\n")
cat("  3. IPCW: emulate stratifies by arm -> non-trivial weights; TrialEmulation\n")
cat("     produces unit weights when numerator = denominator covariates\n")
cat("  4. Data manipulation: TrialEmulation removes post-outcome observations,\n")
cat("     which can affect the expanded dataset size\n")
cat("\n")

# =============================================================================
# EXPORT CSV
# =============================================================================
write.csv(results, output_csv, row.names = FALSE)
cat(sprintf("Results saved to %s\n", output_csv))
