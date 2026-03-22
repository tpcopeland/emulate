# emulate 0.2.0

## New features

* **Lasso propensity scores**: `emulate_weight()` gains `ps_method = "lasso"` for
  cross-validated L1 logistic regression via `glmnet`. Falls back to GLM for
  single-predictor models.

* **Propensity score trimming**: `emulate_weight()` gains `ps_trim` parameter
  to remove observations with extreme propensity scores before weight computation.
  Distinct from `truncate`, which caps weight values after computation.

* **Poisson regression**: `emulate_fit()` now supports `model = "poisson"` for
  log-linear risk models with sandwich standard errors. Reports risk ratios (RR)
  instead of odds ratios. G-formula predictions (`emulate_predict()`) are not
  supported for Poisson models.

* **Propensity score matching**: New `emulate_match()` function provides an
  alternative to IPTW using nearest-neighbor matching via `MatchIt`. Estimates
  ATT by default. Supports both GLM and lasso PS estimation.

* **Propensity score stratification**: New `emulate_stratify()` function creates
  PS quantile strata and fits the outcome model with stratum as a fixed effect
  (logistic/Poisson) or `strata()` term (Cox).

* **Equipoise diagnostics**: `emulate_diagnose()` gains `equipoise` parameter
  to compute preference scores and the percentage of observations in clinical
  equipoise. New plot types `"ps"` and `"equipoise"` in `emulate_plot()`.

* **Negative control calibration**: New `emulate_calibrate()` function adjusts
  treatment effect estimates using negative control outcomes via
  `EmpiricalCalibration`. New `"calibration"` plot type.

## Pipeline changes

* The pipeline state machine now supports three paths after expansion:
  weighting (`emulate_weight` -> `emulate_fit`), matching (`emulate_match`),
  or stratification (`emulate_stratify`). All paths lead to a fitted model
  compatible with downstream functions.

* New state flags: `matched`, `stratified`, `calibrated`.

* `emulate_report()` automatically includes calibrated results when available.

## Dependencies

* Added `glmnet` and `MatchIt` to Imports.
* Added `EmpiricalCalibration` to Suggests.

# emulate 0.1.1

* Initial release with IPTW/IPCW, pooled logistic regression, Cox PH,
  G-formula predictions, weight diagnostics, and publication-quality output.
