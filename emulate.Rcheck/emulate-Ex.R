pkgname <- "emulate"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('emulate')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("emulate_diagnose")
### * emulate_diagnose

flush(stderr()); flush(stdout())

### Name: emulate_diagnose
### Title: Weight diagnostics and covariate balance assessment
### Aliases: emulate_diagnose

### ** Examples





cleanEx()
nameEx("emulate_expand")
### * emulate_expand

flush(stderr()); flush(stdout())

### Name: emulate_expand
### Title: Expand person-period data into sequential emulated trials
### Aliases: emulate_expand

### ** Examples





cleanEx()
nameEx("emulate_fit")
### * emulate_fit

flush(stderr()); flush(stdout())

### Name: emulate_fit
### Title: Fit outcome model for target trial emulation
### Aliases: emulate_fit

### ** Examples





cleanEx()
nameEx("emulate_plot")
### * emulate_plot

flush(stderr()); flush(stdout())

### Name: emulate_plot
### Title: Plot results from target trial emulation
### Aliases: emulate_plot

### ** Examples





cleanEx()
nameEx("emulate_predict")
### * emulate_predict

flush(stderr()); flush(stdout())

### Name: emulate_predict
### Title: Marginal predictions with Monte Carlo confidence intervals
### Aliases: emulate_predict

### ** Examples





cleanEx()
nameEx("emulate_prepare")
### * emulate_prepare

flush(stderr()); flush(stdout())

### Name: emulate_prepare
### Title: Prepare data for target trial emulation
### Aliases: emulate_prepare

### ** Examples

# Create a small person-period dataset
set.seed(42)
n_ids <- 20
n_periods <- 8
dat <- data.frame(
  id = rep(seq_len(n_ids), each = n_periods),
  period = rep(0:(n_periods - 1), times = n_ids),
  treatment = 0L,
  outcome = 0L,
  eligible = 1L,
  age = rep(sample(40:70, n_ids, replace = TRUE), each = n_periods),
  sex = rep(sample(0:1, n_ids, replace = TRUE), each = n_periods)
)
# Simulate some treatment initiation
dat$treatment <- ifelse(dat$period >= sample(2:6, nrow(dat), replace = TRUE), 1L, 0L)
# Simulate a few events
dat$outcome[sample(which(dat$period > 3), 5)] <- 1L

# Prepare the data for a per-protocol analysis
obj <- emulate_prepare(
  data = dat,
  id = "id",
  period = "period",
  treatment = "treatment",
  outcome = "outcome",
  eligible = "eligible",
  covariates = "age",
  baseline_covariates = "sex",
  estimand = "PP"
)
print(obj)




cleanEx()
nameEx("emulate_protocol")
### * emulate_protocol

flush(stderr()); flush(stdout())

### Name: emulate_protocol
### Title: Generate a target trial protocol table
### Aliases: emulate_protocol

### ** Examples

# Define and display a protocol for a statin trial emulation
protocol <- emulate_protocol(
  eligibility     = "Adults 40-80, newly diagnosed hyperlipidemia, no prior CVD",
  treatment       = "Initiate statin within 6 months vs. no statin",
  assignment      = "Observational: assigned by physician prescribing decision",
  followup_start  = "Date of first eligible visit",
  outcome         = "Composite of MI, stroke, or CV death",
  causal_contrast = "Per-protocol effect on 5-year risk",
  analysis        = "Sequential trials, IPTW, pooled logistic regression"
)

# The result is a data.frame
str(protocol)





cleanEx()
nameEx("emulate_report")
### * emulate_report

flush(stderr()); flush(stdout())

### Name: emulate_report
### Title: Publication-quality results table
### Aliases: emulate_report

### ** Examples





cleanEx()
nameEx("emulate_validate")
### * emulate_validate

flush(stderr()); flush(stdout())

### Name: emulate_validate
### Title: Validate prepared data for target trial emulation
### Aliases: emulate_validate

### ** Examples





cleanEx()
nameEx("emulate_weight")
### * emulate_weight

flush(stderr()); flush(stdout())

### Name: emulate_weight
### Title: Compute inverse probability weights for target trial emulation
### Aliases: emulate_weight

### ** Examples





cleanEx()
nameEx("print.emulate")
### * print.emulate

flush(stderr()); flush(stdout())

### Name: print.emulate
### Title: Print method for emulate objects
### Aliases: print.emulate

### ** Examples





cleanEx()
nameEx("summary.emulate")
### * summary.emulate

flush(stderr()); flush(stdout())

### Name: summary.emulate
### Title: Summary method for emulate objects
### Aliases: summary.emulate

### ** Examples





### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
