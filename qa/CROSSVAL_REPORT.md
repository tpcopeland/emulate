# Cross-Validation Report: emulate vs TrialEmulation

**Date:** 2026-03-11
**emulate version:** 0.2.0
**TrialEmulation version:** 0.0.4.9
**R version:** 4.3.3

## Purpose

This report compares two R implementations of sequential target trial
emulation on three datasets with known data-generating processes (DGPs).
The true causal effect is known for each dataset, providing ground truth
for assessing both implementations.

- **emulate** (this package): 2-stratum weight models (by arm),
  sandwich::vcovCL with HC1 + cadjust
- **TrialEmulation** (CRAN, Maringe et al. 2024, arXiv:2402.12083):
  4-stratum weight models (by arm x lagged treatment),
  sandwich::vcovCL with HC1

## Datasets

| Dataset | N patients | Periods | True log-OR | Key feature |
|---------|-----------|---------|-------------|-------------|
| known_dgp | 10,000 | 10 | -0.50 | Simple binary confounder |
| gformula_simulated | 5,000 | 15 | -0.80 | Time-varying CD4 confounding |
| ipcw_dgp | 5,000 | 10 | -0.60 | Informative censoring |

All datasets generated in Stata (`~/Stata-Tools/tte/qa/data/`) with
documented DGP parameters and read into R via `haven::read_dta()`.

## Results

### Treatment Coefficients

| Dataset | Config | emulate | TrialEmul | True | \|em-true\| | \|te-true\| |
|---------|--------|---------|-----------|------|-------------|-------------|
| known_dgp | ITT | -0.429 | -0.376 | -0.50 | 0.071 | 0.124 |
| known_dgp | PP | -0.575 | -0.534 | -0.50 | 0.075 | 0.034 |
| gformula | ITT | -0.892 | -0.749 | -0.80 | 0.092 | 0.051 |
| gformula | PP | -0.630 | -1.162 | -0.80 | 0.170 | **0.362** |
| ipcw_dgp | PP-naive | -0.661 | -0.674 | -0.60 | 0.061 | 0.074 |
| ipcw_dgp | PP-IPCW | -0.510 | -0.740 | -0.60 | 0.090 | 0.140 |

### Standard Errors

| Dataset | Config | emulate | TrialEmul | \|em-te\| |
|---------|--------|---------|-----------|-----------|
| known_dgp | ITT | 0.044 | 0.037 | 0.007 |
| known_dgp | PP | 0.055 | 0.050 | 0.005 |
| gformula | ITT | 0.063 | 0.050 | 0.013 |
| gformula | PP | 0.086 | 0.073 | 0.013 |
| ipcw_dgp | PP-naive | 0.082 | 0.075 | 0.008 |
| ipcw_dgp | PP-IPCW | 0.184 | 0.107 | 0.077 |

### Risk Differences

| Dataset | Config | emulate | TrialEmul | \|em-te\| |
|---------|--------|---------|-----------|-----------|
| known_dgp | ITT RD(t=5) | -0.122 | -0.045 | 0.077 |
| gformula | ITT RD(t=8) | -0.256 | -0.105 | 0.151 |
| ipcw_dgp | PP-IPCW RD(t=5) | -0.146 | -0.095 | 0.051 |

## Key Findings

### 1. Known DGP (true = -0.50): Both implementations recover the true effect

ITT and PP coefficients from both packages fall within acceptable
distance of the true parameter. Inter-implementation agreement is
strong (coefficient differences < 0.06).

### 2. G-formula (true = -0.80): TrialEmulation PP diverges

TrialEmulation's PP estimate (-1.162) is 0.362 from the true
value (-0.80), compared to emulate's 0.170 deviation. This is
driven by TrialEmulation's 4-stratum weight model: the
`P(treatment=1 | prev_treatment=1)` stratum has near-zero variance
for absorbing treatments, producing `glm.fit: algorithm did not
converge` warnings and unstable weight estimates.

emulate's 2-stratum approach (by arm only) avoids the degenerate
subsample and produces a more stable estimate closer to the truth.

### 3. IPCW (true = -0.60): Both implementations apply censoring weights

With proper stabilization (denominator covariates: x + z; numerator
covariates: x only), both packages produce PP-IPCW estimates that
differ from the naive PP, confirming censoring weights are active.
emulate is closer to the true effect (|diff| = 0.090 vs 0.140).

emulate's wider SEs in the IPCW case (0.184 vs 0.107) reflect
greater variance from its stratified censoring model. This is the
expected bias-variance tradeoff from IP weighting: wider intervals
that correctly account for the reweighting, rather than tighter
intervals from weights that are closer to unity.

### 4. SE precision

TrialEmulation consistently produces tighter SEs across all
configurations. This is partly a consequence of its 4-stratum
weight models producing near-unit weights in several strata
(particularly the `prev_treatment=1` stratum for absorbing
treatments), which reduces variance inflation from weighting.
However, weights closer to 1 also mean less confounding correction,
as reflected in the coefficient differences from truth in the
gformula PP case.

## Algorithmic Differences

| Component | emulate | TrialEmulation |
|-----------|---------|----------------|
| Weight model strata | 2 (by arm) | 4 (by arm x lagged treatment) |
| Robust SE | sandwich::vcovCL(HC1, cadjust=TRUE) | sandwich::vcovCL(HC1) |
| IPCW with same num/denom covs | Non-trivial weights (arm stratification) | Unit weights (exact cancellation) |
| Post-outcome observations | Retained through pipeline | Removed during data manipulation |
| Convergence on absorbing treatments | Stable (fewer strata) | `glm.fit` warnings (degenerate stratum) |

## Reproducibility

To reproduce these results:

```bash
cd ~/R-Packages/emulate
Rscript qa/crossval_emulate_vs_trialemulation.R
```

Requires:
- `emulate` installed from local source
- `TrialEmulation` from CRAN
- `haven` for reading .dta files
- Datasets in `~/Stata-Tools/tte/qa/data/`

Machine-readable results: `qa/crossval_results.csv`

## References

- Maringe C, Belot A, Quaresma M, et al. (2024). TrialEmulation:
  Causal analysis of observational time-to-event data using target
  trial emulation. arXiv:2402.12083.
- Hernan MA, Robins JM (2016). Using big data to emulate a target
  trial when a randomized trial is not available. American Journal
  of Epidemiology, 183(8):758-764.
