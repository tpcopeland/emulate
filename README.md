# emulate: Target Trial Emulation via Sequential Trials

<!-- badges: start -->
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-0.2.0-orange.svg)]()
<!-- badges: end -->

An R implementation of the sequential trials framework for target trial
emulation from observational data. Produces numerically equivalent results
to the Stata `tte` package.

---

## What is Target Trial Emulation?

Randomized controlled trials (RCTs) are the gold standard for estimating
causal treatment effects, but they are not always feasible. Ethical
constraints, cost, time, and generalizability concerns often make it
impossible to randomize patients to treatment strategies of interest. Yet
clinical decisions still need to be made, and observational data from
electronic health records, registries, and claims databases are abundant.

Target trial emulation is a framework that bridges this gap. The idea,
formalized by Hernan and Robins, is straightforward: before analyzing
observational data, explicitly specify the hypothetical randomized trial you
*wish* you could run -- the "target trial" -- and then design the
observational analysis to emulate that trial as closely as possible. This
forces the analyst to articulate eligibility criteria, treatment strategies,
follow-up start, outcomes, and the causal contrast of interest *before*
touching the data, reducing the risk of ad hoc analytic decisions that can
introduce bias.

The key insight is that many biases in observational studies (immortal-time
bias, selection bias, confounding by indication) arise because the
observational analysis does not correspond to any well-defined trial. By
mapping each component of the analysis to a component of the target trial,
these biases become visible and correctable.

Target trial emulation is relevant to epidemiologists, biostatisticians,
pharmacoepidemiologists, health economists, and clinical researchers who work
with longitudinal observational data and need to estimate causal effects of
treatments, exposures, or interventions.

## The Sequential Trials Framework

The approach implemented in emulate follows Hernan and Robins (2016). In an
observational dataset with longitudinal follow-up, patients become eligible
for the target trial at different calendar times. Rather than picking a
single enrollment date, the sequential trials method treats *each* eligible
time period as the start of a new emulated trial.

At each trial period, every eligible individual is "cloned" into two copies:
one assigned to the treatment strategy, one assigned to the control strategy.
Each clone is then followed forward in time and artificially censored if
their observed behavior deviates from the strategy they were assigned to.
For example, a clone assigned to "initiate treatment" who never actually
starts treatment is censored at the time of deviation.

This artificial censoring introduces informative censoring bias, which is
corrected using inverse probability of censoring weights (IPCW). The final
analysis pools all the emulated trials together, fitting a weighted outcome
model to the combined data. This "clone-censor-weight" (CCW) approach
handles time-varying confounding, immortal-time bias, and treatment
switching in a principled way.

The result: a marginal causal effect estimate that approximates what you
would have obtained from the target trial, with valid confidence intervals.

## Installation

### From GitHub (recommended)

```r
# Install devtools if needed
install.packages("devtools")

# Install emulate from GitHub
devtools::install_github("tpcopeland/emulate")
```

### From a local source

```r
# If you have the source directory
install.packages("/path/to/emulate", repos = NULL, type = "source")
```

### Dependencies

emulate depends on:

- **data.table** -- fast data manipulation
- **survival** -- Cox proportional hazards models
- **sandwich** -- clustered standard errors
- **MASS** -- multivariate normal sampling for Monte Carlo CIs
- **ggplot2** -- publication-quality plots

Optional:

- **openxlsx** -- Excel export in `emulate_report()` and `emulate_protocol()`
- **knitr** / **rmarkdown** -- building vignettes

## Quick Start

A complete analysis in about 20 lines:

```r
library(emulate)

# Load the bundled example dataset (503 patients, 48,400 person-periods)
trial_data <- read.csv(system.file("extdata", "trial_example.csv", package = "emulate"))

# --- Intention-to-Treat Analysis ---

# Step 1: Prepare data
obj <- emulate_prepare(trial_data,
                   id = "id", period = "period",
                   treatment = "treatment", outcome = "outcome",
                   eligible = "eligible",
                   covariates = c("catvarA", "catvarB", "catvarC",
                                  "nvarA", "nvarB", "nvarC"),
                   estimand = "ITT")

# Step 2: Validate
obj <- emulate_validate(obj)

# Step 3: Expand into sequential trials
obj <- emulate_expand(obj, maxfollowup = 8)

# Step 4: Fit pooled logistic regression
obj <- emulate_fit(obj,
               outcome_cov = c("catvarA", "catvarB", "catvarC",
                               "nvarA", "nvarB", "nvarC"))

# Step 5: Predict cumulative incidence
obj <- emulate_predict(obj, times = 0:8, type = "cum_inc",
                   difference = TRUE, samples = 100, seed = 12345)

# Step 6: Plot
emulate_plot(obj, type = "cumhaz", title = "ITT Cumulative Incidence")

# Step 7: Report
emulate_report(obj, format = "display", predictions = TRUE)
```

## The Pipeline

emulate implements a step-by-step pipeline. Each function takes a `emulate` object
and returns an updated `emulate` object, so calls are chained naturally:

```
emulate_prepare --> emulate_validate --> emulate_expand --> emulate_weight --> emulate_fit --> emulate_predict
                                                                       \-> emulate_diagnose
                                                                       \-> emulate_report
                                                                       \-> emulate_plot
                                                                       \-> emulate_protocol
```

### Pipeline Steps

| Step | Function | Purpose | Required? |
|------|----------|---------|-----------|
| 1 | `emulate_prepare()` | Map variables, validate structure, create emulate object | Yes |
| 2 | `emulate_validate()` | Run 10 data quality checks | Recommended |
| 3 | `emulate_expand()` | Clone-censor-weight expansion into sequential trials | Yes |
| 4 | `emulate_weight()` | Compute inverse probability weights (IPTW/IPCW) | PP/AT only |
| 5 | `emulate_fit()` | Fit pooled logistic or Cox PH model | Yes |
| 6 | `emulate_predict()` | G-formula marginal predictions with Monte Carlo CIs | Optional |
| 7 | `emulate_diagnose()` | Weight diagnostics and covariate balance | Optional |
| 8 | `emulate_report()` | Publication-quality results table | Optional |
| 9 | `emulate_plot()` | Visualization (cumulative incidence, KM, weights, balance) | Optional |
| 0 | `emulate_protocol()` | Standalone: 7-component target trial protocol table | Optional |

Guards enforce the correct order: calling `emulate_fit()` before `emulate_expand()`
will produce a clear error message telling you what to run first.

## Estimands Explained

The `estimand` argument in `emulate_prepare()` controls the causal contrast.
emulate supports three estimands:

| Estimand | Full Name | Question Answered | Weighting? |
|----------|-----------|-------------------|------------|
| **ITT** | Intention-to-Treat | What happens to outcomes if patients are *assigned* to treatment, regardless of whether they actually take it? | No (weights = 1) |
| **PP** | Per-Protocol | What happens to outcomes if everyone *adheres* to their assigned treatment strategy? | Yes (IPTW) |
| **AT** | As-Treated | What happens to outcomes given the treatment patients *actually received* over time? | Yes (IPTW) |

### When to Use Each

- **ITT** is the most robust estimand. It does not require modeling treatment
  switching and is less susceptible to model misspecification. Use it as your
  primary analysis when treatment switching is rare, or as a complement to PP.

- **PP** answers the clinically relevant question "what if everyone followed
  the protocol?" but requires correct modeling of treatment switching via
  inverse probability weights. Use it when treatment adherence varies and you
  want to estimate the effect under full compliance.

- **AT** is appropriate when treatment is absorbing (once started, never
  stopped) or when you are interested in the effect of treatment as actually
  received. For absorbing treatment, AT and PP produce identical results.

## Data Format Requirements

emulate expects data in **person-period** (long) format: one row per individual
per time period. Here is an example:

| id | period | treatment | outcome | eligible | age | biomarker |
|----|--------|-----------|---------|----------|-----|-----------|
| 1  | 0      | 0         | 0       | 1        | 55  | 3.2       |
| 1  | 1      | 1         | 0       | 0        | 55  | 3.5       |
| 1  | 2      | 1         | 0       | 0        | 55  | 3.8       |
| 1  | 3      | 1         | 1       | 0        | 55  | 4.1       |
| 2  | 0      | 0         | 0       | 1        | 62  | 2.8       |
| 2  | 1      | 0         | 0       | 1        | 62  | 2.9       |
| 2  | 2      | 1         | 0       | 0        | 62  | 3.0       |
| 2  | 3      | 1         | 0       | 0        | 62  | 3.1       |

### Required Variables

- **id**: Unique patient identifier (any type)
- **period**: Integer time period (0, 1, 2, ...). Must be consecutive within each individual.
- **treatment**: Binary (0/1). Current treatment status at each period.
- **outcome**: Binary (0/1). Terminal event indicator. Once outcome = 1, the individual should have no further rows.
- **eligible**: Binary (0/1). Whether the individual is eligible to enter a new trial at this period.

### Optional Variables

- **censor**: Binary (0/1). External censoring indicator (e.g., loss to follow-up, administrative end).
- **covariates**: Time-varying covariates (e.g., biomarkers, lab values). Frozen at baseline of each trial.
- **baseline_covariates**: Time-fixed covariates (e.g., sex, birth year). Also frozen at baseline.

## Detailed Function Reference

### `emulate_prepare()`

Maps variable names and creates the `emulate` object that flows through the pipeline.

```r
obj <- emulate_prepare(
  data,                           # data.frame or data.table
  id = "id",                      # patient identifier
  period = "period",              # time period variable
  treatment = "treatment",        # binary treatment indicator
  outcome = "outcome",            # binary outcome indicator
  eligible = "eligible",          # binary eligibility indicator
  censor = NULL,                  # optional: external censoring
  covariates = c("x1", "x2"),    # optional: time-varying covariates
  baseline_covariates = c("z1"),  # optional: time-fixed covariates
  estimand = "PP",                # "ITT", "PP", or "AT"
  prefix = "_emulate_"                # prefix for generated variables
)
```

The function validates that all named columns exist, checks that binary
variables are truly 0/1, confirms person-period structure (no duplicate
id-period rows), and prints a summary of the data.

### `emulate_validate()`

Runs 10 data quality checks that mirror the Stata `emulate_validate` command:

```r
obj <- emulate_validate(obj, strict = FALSE, verbose = FALSE)
```

| Check | What it tests |
|-------|---------------|
| 1 | Person-period format (no duplicates) |
| 2 | No gaps in period sequences |
| 3 | Outcome is terminal (no rows after event) |
| 4 | Treatment consistency (switching patterns) |
| 5 | Missing data in core variables |
| 6 | Eligibility consistency (not eligible after outcome) |
| 7 | Sufficient eligible observations per period (>= 10) |
| 8 | Positivity (treatment variation among eligible) |
| 9 | Period numbering (starts at 0 or 1) |
| 10 | Event rate (>= 5 events) |

With `strict = TRUE`, warnings become errors and the function stops. With
`verbose = TRUE`, per-variable missing data counts are shown.

### `emulate_expand()`

Performs the clone-censor-weight expansion into sequential emulated trials:

```r
obj <- emulate_expand(obj, trials = NULL, maxfollowup = 0, grace = 0)
```

- **trials**: Integer vector of specific trial periods to use. Default `NULL`
  uses all periods where at least one individual is eligible.
- **maxfollowup**: Maximum follow-up periods after trial entry. `0` means
  unlimited.
- **grace**: Number of periods of deviation allowed before censoring
  (PP/AT only). A grace period of 2 means a patient can deviate for up to
  2 periods without being censored.

For ITT, individuals are *not* cloned; they keep their observed treatment
and are never artificially censored. For PP/AT, each eligible individual is
cloned into a treatment arm and a control arm, then censored when their
observed behavior deviates from the assigned strategy.

The expansion can dramatically increase dataset size. A dataset with 503
individuals and 4 eligible periods produces over 100,000 expanded rows.

### `emulate_weight()`

Computes inverse probability weights to correct for the informative
censoring introduced by the clone-censor-weight expansion:

```r
obj <- emulate_weight(
  obj,
  switch_d_cov = c("x1", "x2"),   # denominator covariates for switch model
  switch_n_cov = c("x1"),          # numerator covariates (stabilized)
  censor_d_cov = NULL,             # denominator covariates for IPCW
  censor_n_cov = NULL,             # numerator covariates for IPCW
  pool_switch = FALSE,             # pool switch models across arms?
  pool_censor = FALSE,             # pool censor models across arms?
  truncate = c(1, 99),             # percentile truncation
  stabilized = TRUE                # use stabilized weights?
)
```

**Switch weights (IPTW)**: Model the probability of treatment switching to
correct for artificial censoring. The denominator model includes time-varying
covariates; the numerator model uses only baseline covariates (or none) for
stabilization.

**Censoring weights (IPCW)**: Optionally model external censoring. Useful
when loss to follow-up is informative (correlated with treatment or outcome).

**Truncation**: Extreme weights can make estimates unstable. The `truncate`
argument clips weights at the specified percentiles. Common choices are
`c(1, 99)` or `c(5, 95)`.

For ITT analyses, `emulate_weight()` simply sets all weights to 1.

### `emulate_fit()`

Fits the outcome model with cluster-robust standard errors:

```r
obj <- emulate_fit(
  obj,
  outcome_cov = c("x1", "x2"),        # covariates in outcome model
  model = "logistic",                   # "logistic" or "cox"
  followup_spec = "quadratic",         # time specification for follow-up
  trial_period_spec = "quadratic",     # time specification for trial period
  cluster = NULL,                       # cluster variable (default: id)
  level = 95                            # confidence level
)
```

The model includes the treatment arm indicator, follow-up time terms, trial
period terms, and any specified covariates. Standard errors are clustered by
patient ID using the HC1 sandwich estimator from the `sandwich` package.

See the [Time Specifications](#time-specifications) section below for
guidance on choosing `followup_spec` and `trial_period_spec`.

### `emulate_predict()`

Computes G-formula marginal predictions with Monte Carlo confidence
intervals:

```r
obj <- emulate_predict(
  obj,
  times = 0:8,              # follow-up times for prediction
  type = "cum_inc",          # "cum_inc" or "survival"
  samples = 100,             # Monte Carlo samples for CIs
  seed = 12345,              # reproducibility
  level = 95,                # confidence level
  difference = FALSE         # compute risk differences?
)
```

Predictions are computed by averaging individual-level survival curves across
a reference population (all individuals at follow-up time 0). Confidence
intervals are obtained by drawing coefficients from their multivariate
normal sampling distribution and recomputing predictions for each draw.

With `difference = TRUE`, the function also computes risk differences
(treatment minus control) at each time point with CIs.

Currently only supports the logistic model. Calling after a Cox fit produces
a clear error.

### `emulate_diagnose()`

Reports weight diagnostics and optionally assesses covariate balance:

```r
obj <- emulate_diagnose(
  obj,
  balance_covariates = c("x1", "x2"),  # covariates for SMD calculation
  by_trial = FALSE                      # show per-trial weight stats?
)
```

Outputs:

- Weight distribution (mean, SD, percentiles, min, max)
- Effective sample size (ESS) overall and by arm
- Standardized mean differences (SMD) before and after weighting
- Per-trial weight statistics (if `by_trial = TRUE`)

A maximum weighted SMD below 0.1 indicates adequate balance.

### `emulate_plot()`

Creates ggplot2 visualizations:

```r
p <- emulate_plot(obj, type = "cumhaz", ci = TRUE, title = "My Plot")
```

| Type | What it shows | Requires |
|------|---------------|----------|
| `"cumhaz"` | Marginal cumulative incidence curves with CIs | `emulate_predict()` |
| `"km"` | Kaplan-Meier survival curves | `emulate_expand()` |
| `"weights"` | Density plot of IP weight distribution | `emulate_weight()` |
| `"balance"` | Love plot of standardized mean differences | `emulate_diagnose()` with balance |

### `emulate_report()`

Generates a publication-quality results table:

```r
emulate_report(
  obj,
  format = "display",          # "display", "csv", or "excel"
  export = NULL,               # file path for csv/excel
  decimals = 3,                # decimal places
  eform = FALSE,               # exponentiate (OR/HR)?
  predictions = FALSE,         # include prediction table?
  title = NULL                 # table title
)
```

The display format prints to console. CSV and Excel formats export to a
file. Excel output includes separate worksheets for summary, coefficients,
and predictions.

### `emulate_protocol()`

Generates the 7-component target trial protocol table (standalone, does not
require data):

```r
protocol <- emulate_protocol(
  eligibility    = "Adults aged 18+ with condition X, no prior treatment",
  treatment      = "Strategy A: initiate drug at diagnosis; Strategy B: no drug",
  assignment     = "Patients will be assigned at the start of each eligible period",
  followup_start = "Start of the period when eligibility criteria are met",
  outcome        = "First occurrence of event Y within 2 years",
  causal_contrast = "Intention-to-treat effect and per-protocol effect",
  analysis       = "Pooled logistic regression with IPTW; G-formula predictions",
  format = "display"           # also "csv", "excel", or "latex"
)
```

## Time Specifications

The `followup_spec` and `trial_period_spec` arguments in `emulate_fit()` control
how time is modeled in the outcome regression. The goal is to flexibly
capture the baseline hazard shape without overfitting.

| Specification | Terms Added | When to Use |
|---------------|-------------|-------------|
| `"linear"` | t | Short follow-up, approximately constant hazard |
| `"quadratic"` | t, t^2 | **Default.** Good for most applications. |
| `"cubic"` | t, t^2, t^3 | Longer follow-up with complex hazard shapes |
| `"ns(3)"` | 3 restricted cubic spline basis variables | Long follow-up, highly non-linear hazard. Change `3` to desired df. |
| `"none"` | (nothing) | When you do not want time adjustment |

**Guidance**: Start with `"quadratic"`. If results are sensitive to the time
specification, try `"ns(3)"` or `"cubic"`. The validation suite confirms
that quadratic, cubic, and ns(3) produce nearly identical estimates on the
test datasets.

**Important**: emulate uses Harrell restricted cubic splines (RCS), not B-splines.
This matches the Stata `tte` package and differs from `splines::ns()` in R.
Do not use `splines::ns()` as a substitute -- the basis functions, knot
placement, and resulting coefficients will differ.

## Weight Models

### Inverse Probability of Treatment Weights (IPTW)

IPTW corrects for the artificial censoring introduced by the clone-censor-weight
expansion. The idea is to up-weight observations that were likely to be
censored (because they were "swimming against the tide" of their covariates)
and down-weight observations that were unlikely to be censored.

**Denominator model**: P(treatment_t | treatment_{t-1}, covariates, follow-up time).
This is the probability of the observed treatment given the full covariate
history. It captures what "actually happened."

**Numerator model** (stabilized): P(treatment_t | treatment_{t-1}).
This is a simpler model that only conditions on prior treatment. Using a
numerator model produces *stabilized* weights that have mean closer to 1
and lower variance.

The stabilized weight at each time point is:

```
sw_t = P(A_t | A_{t-1}) / P(A_t | A_{t-1}, L_t)
```

The cumulative weight is the running product of sw_t over follow-up time.

### Inverse Probability of Censoring Weights (IPCW)

IPCW corrects for external censoring (loss to follow-up, administrative end
of study) when it is informative -- that is, when the probability of being
censored depends on covariates related to the outcome.

The IPCW weight at each time point is:

```
cw_t = P(uncensored_t) / P(uncensored_t | L_t)
```

When both IPTW and IPCW are used, the final weight is their product.

## Diagnostics

After computing weights, use `emulate_diagnose()` to check:

1. **Effective Sample Size (ESS)**: ESS = (sum(w))^2 / sum(w^2). If ESS is
   much smaller than the actual sample size, extreme weights are driving the
   estimates. Consider tighter truncation.

2. **Standardized Mean Differences (SMD)**: After weighting, the weighted SMD
   for each covariate should be below 0.1 (the conventional threshold). If
   any covariate has SMD > 0.1 after weighting, the weight model may be
   misspecified.

3. **Weight Distribution**: Look for extreme weights. Mean should be near 1
   for stabilized weights. Large max/min ratios indicate potential instability.
   The 1st and 99th percentiles help identify outliers.

4. **Per-Trial Statistics**: Use `by_trial = TRUE` to check whether weight
   distributions are consistent across trial periods. Large differences
   suggest time-varying model misspecification.

## Output Formats

`emulate_report()` and `emulate_protocol()` support multiple output formats:

| Format | Description | Requires |
|--------|-------------|----------|
| `"display"` | Print to R console | Nothing |
| `"csv"` | Comma-separated values file | `export` path |
| `"excel"` | Multi-sheet Excel workbook | `export` path + `openxlsx` |
| `"latex"` | LaTeX table (protocol only) | `export` path |

Example:

```r
# Console display
emulate_report(obj, format = "display", eform = TRUE)

# Excel with predictions
emulate_report(obj, format = "excel", export = "results.xlsx", predictions = TRUE)

# Protocol as LaTeX
emulate_protocol(..., format = "latex", export = "protocol.tex")
```

## Comparison with Stata tte

emulate is designed to produce numerically equivalent results to the Stata `tte`
package (v1.0.4). The two implementations share the same algorithmic design.

### Equivalent Commands

| emulate (R) | tte (Stata) |
|----------|-------------|
| `emulate_prepare()` | `emulate_prepare` |
| `emulate_validate()` | `emulate_validate` |
| `emulate_expand()` | `emulate_expand` |
| `emulate_weight()` | `emulate_weight` |
| `emulate_fit()` | `emulate_fit` |
| `emulate_predict()` | `emulate_predict` |
| `emulate_diagnose()` | `emulate_diagnose` |
| `emulate_plot()` | `emulate_plot` |
| `emulate_report()` | `emulate_report` |
| `emulate_protocol()` | `emulate_protocol` |

### Known Algorithmic Differences

While results are designed to match, small numerical differences arise from:

| Difference | Stata | R | Impact |
|-----------|-------|---|--------|
| **Clustered SE small-sample correction** | G/(G-1) correction | HC1 via `sandwich::vcovCL` | ~3% SE difference on small datasets |
| **Weighted variance** | `_emulate_diagnose` formula | `sum(w*(x-wm)^2)/sum(w)` | Minor differences in SMD |
| **Quantile type** | `_pctile` (type 2 equivalent) | `quantile(..., type = 2)` | Matched by design |
| **GLM solver** | IRLS (Stata `glm`) | IRLS (R `glm`) | Machine-precision differences |

In cross-validation on the trial_example dataset (503 patients), the ITT
coefficient difference between Stata and R TrialEmulation is 3.2%.

## Comparison with R TrialEmulation

The CRAN `TrialEmulation` package (Maringe et al., 2024) is another R
implementation of sequential trial emulation. Key differences:

| Feature | emulate | TrialEmulation |
|---------|------|----------------|
| **Spline basis** | Harrell RCS (matches Stata) | B-splines via `splines::ns()` |
| **Weight strata** | Per-arm by default, pooled optional | Pooled by default |
| **Data expansion** | In-memory data.table | Data sampling for large datasets |
| **Pipeline style** | Object-based (emulate object flows through) | Function-based (separate data objects) |
| **Cox model** | Supported (`model = "cox"`) | Not supported |
| **Protocol table** | Built-in `emulate_protocol()` | Not included |
| **Diagnostics** | Built-in `emulate_diagnose()` with ESS, SMD, Love plots | Separate tools needed |
| **Output** | Display, CSV, Excel, LaTeX | Display |
| **Stata parity** | Designed for numerical equivalence | Independent implementation |

Because emulate uses Harrell RCS while TrialEmulation uses B-splines,
coefficients on spline terms will differ when using natural spline time
specifications. Point estimates and predictions should be similar but not
identical.

## Validation

emulate includes a comprehensive test suite:

- **127 tests** in `tests/testthat/` covering all exported functions
- Unit tests for each pipeline step (prepare, validate, expand, weight, fit, predict, diagnose, plot, report, protocol)
- Integration tests for full ITT and PP pipelines
- Validation tests against known data-generating processes
- Cross-validation against the R TrialEmulation package
- Edge case tests (small N, few events, single trial period)
- Spline basis tests (Harrell RCS correctness)
- Pipeline guard tests (correct error on out-of-order calls)

The companion Stata `tte` package has a separate 103-test validation suite
(see `VALIDATION_REPORT.md` in the Stata tte package) covering 17
independent validation exercises including cross-validation with R
TrialEmulation, NHEFS real data, known DGP Monte Carlo, null effect
calibration, and stress tests.

## References

- Hernan MA, Robins JM. Using big data to emulate a target trial when a
  randomized trial is not available. *American Journal of Epidemiology*.
  2016;183(8):758-764. doi:10.1093/aje/kwv254

- Hernan MA, Robins JM. *Causal Inference: What If*. Boca Raton: Chapman &
  Hall/CRC, 2020.

- Maringe C, Benitez Majano S, et al. TrialEmulation: An R package for
  target trial emulation. *arXiv*. 2024;2402.12083.

- Hernan MA, Sauer BC, Hernandez-Diaz S, Platt R, Shrier I, Robins JM.
  Specifying a target trial prevents immortal time bias and other
  self-inflicted injuries in observational analyses. *Journal of Clinical
  Epidemiology*. 2016;79:70-75.

- Danaei G, Rodriguez LAG, Cantero OF, Logan R, Hernan MA. Observational
  data for comparative effectiveness research: An emulation of randomised
  trials of statins and primary prevention of coronary heart disease.
  *Statistical Methods in Medical Research*. 2013;22(1):70-96.

- Maringe C, et al. Reflection on modern methods: trial emulation in the
  presence of immortal-time bias. *International Journal of Epidemiology*.
  2020;49(5):1719-1729.

## Authors and License

**Author**: Timothy P. Copeland (<timothy.copeland@ki.se>)

**License**: MIT

Copyright 2026 Timothy P. Copeland
