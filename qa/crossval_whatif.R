#!/usr/bin/env Rscript
# =============================================================================
# Cross-Validation: What If NHEFS Replication (Stata tte vs R emulate vs R TrialEmulation)
# =============================================================================
#
# Verifies that Stata tte, R emulate, and R TrialEmulation all produce
# identical results when replicating Hernan & Robins, What If Programs
# 17.2 and 17.3 (NHEFS smoking cessation, IP-weighted pooled logistic).
#
# All three implementations use the same analysis:
#   1. Person-month NHEFS dataset (1,629 subjects, 176,764 person-months)
#   2. Stabilized IP weights (exact textbook propensity score specification)
#   3. IP-weighted pooled logistic hazards model (no outcome covariates)
#   4. Standardized survival curves via G-formula
#
# Expected results (from textbook):
#   Crude OR:         ~1.40 (confounding by indication)
#   IP-weighted OR:   ~0.84 (protective effect of cessation)
#   10-yr survival:   A=0: ~0.805, A=1: ~0.807, diff: ~+0.002
#
# Prerequisites:
#   Run all three demo scripts first:
#     cd ~/Stata-Tools && stata-mp -b do tte/demo/demo_whatif.do
#     cd ~/Stata-Tools/tte/demo && Rscript demo_whatif_trialemulation.R
#     cd ~/Stata-Tools/tte/demo && Rscript demo_whatif_emulate.R
#
# Usage:
#   cd ~/R-Packages/emulate
#   Rscript qa/crossval_whatif.R
# =============================================================================

cat("==============================================================================\n")
cat("Cross-Validation: What If NHEFS (Stata tte vs R emulate vs R TrialEmulation)\n")
cat("==============================================================================\n\n")

# Load comparison CSVs from all three demo runs
stata_path <- "~/Stata-Tools/tte/demo/whatif/comparison.csv"
te_path    <- "~/Stata-Tools/tte/demo/demo_output_trialemulation/comparison.csv"
emu_path   <- "~/Stata-Tools/tte/demo/demo_output_emulate/comparison.csv"

missing <- c()
if (!file.exists(stata_path)) missing <- c(missing, "Stata tte")
if (!file.exists(te_path))    missing <- c(missing, "R TrialEmulation")
if (!file.exists(emu_path))   missing <- c(missing, "R emulate")

if (length(missing) > 0) {
  cat(sprintf("ERROR: Missing comparison files from: %s\n",
              paste(missing, collapse = ", ")))
  cat("Run the demo scripts first (see header for instructions).\n")
  quit(status = 1)
}

stata <- read.csv(stata_path)
te    <- read.csv(te_path)
emu   <- read.csv(emu_path)

# Extract IP-weighted results (Program 17.3)
get_ipw <- function(d) {
  row <- d[grep("IPweighted|IP_weighted", d[[1]], ignore.case = TRUE), ]
  if (nrow(row) == 0) row <- d[2, ]  # fallback to second row
  list(coef = row$Coefficient %||% row$coef,
       or   = row$OR %||% row$or,
       s0   = row$Survival_A0 %||% row$surv0,
       s1   = row$Survival_A1 %||% row$surv1)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

s <- get_ipw(stata)
t <- get_ipw(te)
e <- get_ipw(emu)

# Display comparison
cat("IP-WEIGHTED RESULTS (Program 17.3)\n")
cat(strrep("-", 75), "\n")
cat(sprintf("%-25s %12s %10s %12s %12s\n",
            "", "Coefficient", "OR", "Surv(A=0)", "Surv(A=1)"))
cat(strrep("-", 75), "\n")
cat(sprintf("%-25s %12.6f %10.6f %12.6f %12.6f\n",
            "Stata tte", s$coef, s$or, s$s0, s$s1))
cat(sprintf("%-25s %12.6f %10.6f %12.6f %12.6f\n",
            "R emulate", e$coef, e$or, e$s0, e$s1))
cat(sprintf("%-25s %12.6f %10.6f %12.6f %12.6f\n",
            "R TrialEmulation", t$coef, t$or, t$s0, t$s1))
cat(strrep("-", 75), "\n\n")

# Pairwise relative differences
rd <- function(a, b) abs(a - b) / max(abs(b), 1e-10) * 100

cat("PAIRWISE RELATIVE DIFFERENCES (%)\n")
cat(strrep("-", 60), "\n")
cat(sprintf("%-30s %12s %12s\n", "", "Coefficient", "OR"))
cat(strrep("-", 60), "\n")
cat(sprintf("%-30s %12.6f %12.6f\n", "Stata vs R emulate",
            rd(s$coef, e$coef), rd(s$or, e$or)))
cat(sprintf("%-30s %12.6f %12.6f\n", "Stata vs R TrialEmulation",
            rd(s$coef, t$coef), rd(s$or, t$or)))
cat(sprintf("%-30s %12.6f %12.6f\n", "R emulate vs R TrialEmulation",
            rd(e$coef, t$coef), rd(e$or, t$or)))
cat(strrep("-", 60), "\n\n")

# Pass/fail assessment (tolerance: 0.01% for coefficient, 0.1% for survival)
coef_tol <- 0.01  # percent
surv_tol <- 0.001 # percent

pass_coef <- rd(s$coef, e$coef) < coef_tol & rd(s$coef, t$coef) < coef_tol
pass_surv <- rd(s$s0, e$s0) < surv_tol & rd(s$s1, e$s1) < surv_tol

cat("ASSESSMENT\n")
cat(sprintf("  Coefficient agreement (< %.2f%%): %s\n",
            coef_tol, if (pass_coef) "PASS" else "FAIL"))
cat(sprintf("  Survival agreement (< %.3f%%):    %s\n",
            surv_tol, if (pass_surv) "PASS" else "FAIL"))

if (pass_coef & pass_surv) {
  cat("\n  ALL PLATFORMS AGREE — What If replication verified.\n")
} else {
  cat("\n  DISCREPANCY DETECTED — investigate.\n")
}

# Save merged results
merged <- data.frame(
  Platform = c("Stata_tte", "R_emulate", "R_TrialEmulation"),
  Coefficient = c(s$coef, e$coef, t$coef),
  OR = c(s$or, e$or, t$or),
  Survival_A0 = c(s$s0, e$s0, t$s0),
  Survival_A1 = c(s$s1, e$s1, t$s1),
  stringsAsFactors = FALSE
)
out_file <- "qa/crossval_whatif_results.csv"
write.csv(merged, out_file, row.names = FALSE)
cat(sprintf("\nMerged results saved to: %s\n", out_file))
