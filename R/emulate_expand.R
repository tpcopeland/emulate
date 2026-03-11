#' Expand person-period data into sequential emulated trials
#'
#' This is the core step that implements the \strong{clone-censor-weight}
#' approach (Danaei et al., 2013) for sequential trial emulation. The
#' function takes your prepared person-period data and creates a set of
#' emulated trials -- one starting at each eligible time period -- then
#' stacks them into a single expanded dataset.
#'
#' @section How it works (clone-censor-weight):
#' The expansion proceeds in three conceptual steps:
#' \enumerate{
#'   \item \strong{Clone}: At each trial start time, every eligible individual
#'     is duplicated ("cloned") into two copies -- one assigned to the
#'     treatment arm (\code{arm = 1}) and one to the control arm
#'     (\code{arm = 0}). This is necessary because we cannot observe what
#'     would have happened if the individual had followed the opposite
#'     strategy.
#'   \item \strong{Censor}: Each clone is followed forward in time. If a
#'     clone deviates from its assigned strategy (e.g., the treatment-arm
#'     clone stops taking treatment, or the control-arm clone starts
#'     treatment), it is \strong{artificially censored} at the time of
#'     deviation. This ensures that each arm only contains person-time that
#'     is consistent with the assigned strategy.
#'   \item \strong{Weight}: The artificial censoring introduces selection bias
#'     (deviators may differ from adherers). This is corrected in the next
#'     step (\code{\link{emulate_weight}}) using inverse probability weights.
#' }
#'
#' For the \strong{ITT estimand}, no cloning or artificial censoring is
#' performed. Each eligible individual enters one arm based on their observed
#' treatment at the trial start time.
#'
#' @section Generated columns:
#' The function adds the following columns to the expanded dataset (prefixed
#' by the \code{prefix} set in \code{\link{emulate_prepare}}, default
#' \code{"_emulate_"}):
#' \describe{
#'   \item{\code{_emulate_trial}}{Integer. The period at which this emulated trial
#'     started. For example, if \code{_emulate_trial = 3}, this row belongs to the
#'     trial that started at period 3.}
#'   \item{\code{_emulate_arm}}{Integer (0 or 1). The assigned treatment arm.
#'     \code{1} = treatment strategy (expected to remain on treatment),
#'     \code{0} = control strategy (expected to remain off treatment).}
#'   \item{\code{_emulate_followup}}{Integer. The follow-up time within this trial,
#'     starting at 0. Computed as \code{period - trial_start_period}.}
#'   \item{\code{_emulate_censored}}{Integer (0 or 1). Whether this observation
#'     was artificially censored due to deviation from the assigned strategy.
#'     \code{1} = censored at this time point, \code{0} = not censored.
#'     Always \code{0} for ITT.}
#'   \item{\code{_emulate_outcome_obs}}{Integer (0 or 1). The observed outcome,
#'     but set to \code{0} at the censoring time point (since the outcome
#'     is not "observed" when an individual deviates). For uncensored
#'     observations, this equals the original outcome variable.}
#' }
#'
#' @section Grace period:
#' The \code{grace} parameter allows a specified number of follow-up periods
#' before treatment deviation triggers censoring. For example, with
#' \code{grace = 2}, an individual in the treatment arm who stops treatment at
#' follow-up time 1 would \emph{not} be censored (deviation within the grace
#' period), but stopping at follow-up time 3 \emph{would} trigger censoring.
#' This is useful when treatment adherence cannot be instantaneous. The grace
#' period is ignored for the ITT estimand.
#'
#' @param obj A \code{emulate} object returned by \code{\link{emulate_prepare}} (or
#'   \code{\link{emulate_validate}}). The \code{prepared} state flag must be
#'   \code{TRUE}.
#' @param trials An optional integer vector specifying which time periods to
#'   use as trial start times. If \code{NULL} (the default), all periods
#'   where at least one individual is eligible (\code{eligible = 1}) are
#'   used. Specifying a subset can be useful for computational reasons or
#'   sensitivity analyses.
#' @param maxfollowup An integer specifying the maximum number of follow-up
#'   periods to include per trial. If \code{0} (the default), follow-up is
#'   unlimited (runs until end of data, outcome, or censoring). For example,
#'   \code{maxfollowup = 10} restricts each trial to at most 10 follow-up
#'   periods.
#' @param grace A non-negative integer specifying the grace period (in number
#'   of periods) allowed before treatment deviation triggers artificial
#'   censoring. Default is \code{0} (no grace period -- any deviation is
#'   immediately censored). Only applies to PP and AT estimands. See the
#'   "Grace period" section for details.
#'
#' @return The input \code{emulate} object with updated data (the expanded
#'   dataset) and \code{state$expanded} set to \code{TRUE}. The
#'   \code{expansion} slot is populated with metadata:
#'   \describe{
#'     \item{\code{n_trials}}{Number of emulated trials created.}
#'     \item{\code{trial_periods}}{Integer vector of trial start periods.}
#'     \item{\code{n_expanded}}{Total number of rows in the expanded dataset.}
#'     \item{\code{n_treat}}{Number of rows in the treatment arm.}
#'     \item{\code{n_control}}{Number of rows in the control arm.}
#'     \item{\code{n_censored}}{Number of artificially censored observations.}
#'     \item{\code{n_events}}{Number of observed outcome events.}
#'     \item{\code{expansion_ratio}}{Ratio of expanded rows to original rows.}
#'   }
#'   The object is returned invisibly.
#'
#' @seealso \code{\link{emulate_prepare}} for the previous step,
#'   \code{\link{emulate_weight}} for the next step (computing IP weights).
#'
#' @references
#' Danaei G, Rodriguez LAG, Cantero OF, Logan R, Hernan MA (2013).
#' "Observational data for comparative effectiveness research: An emulation
#' of randomised trials of statins and primary prevention of coronary heart
#' disease." \emph{Statistical Methods in Medical Research}, 22(1), 70-96.
#'
#' Hernan MA, Robins JM (2016). "Using Big Data to Emulate a Target Trial
#' When a Randomized Trial Is Not Available." \emph{American Journal of
#' Epidemiology}, 183(8), 758-764.
#'
#' @examples
#' \donttest{
#' # Create a small test dataset
#' set.seed(42)
#' n_ids <- 20
#' n_per <- 8
#' dat <- data.frame(
#'   id = rep(seq_len(n_ids), each = n_per),
#'   period = rep(0:(n_per - 1), times = n_ids),
#'   treatment = 0L, outcome = 0L, eligible = 1L
#' )
#' dat$treatment <- ifelse(dat$period >= sample(3:7, nrow(dat), replace = TRUE), 1L, 0L)
#' dat$outcome[sample(which(dat$period > 4), 3)] <- 1L
#'
#' obj <- emulate_prepare(dat, id = "id", period = "period",
#'                    treatment = "treatment", outcome = "outcome",
#'                    eligible = "eligible", estimand = "PP")
#'
#' # Expand with all eligible periods, max 5 follow-up periods
#' obj <- emulate_expand(obj, maxfollowup = 5)
#'
#' # Expand with specific trial periods
#' obj2 <- emulate_prepare(dat, id = "id", period = "period",
#'                     treatment = "treatment", outcome = "outcome",
#'                     eligible = "eligible", estimand = "PP")
#' obj2 <- emulate_expand(obj2, trials = c(0, 2, 4), maxfollowup = 5)
#'
#' # ITT estimand: no cloning
#' obj3 <- emulate_prepare(dat, id = "id", period = "period",
#'                     treatment = "treatment", outcome = "outcome",
#'                     eligible = "eligible", estimand = "ITT")
#' obj3 <- emulate_expand(obj3, maxfollowup = 5)
#' }
#'
#' @export
emulate_expand <- function(obj, trials = NULL, maxfollowup = 0L, grace = 0L) {
  .check_prepared(obj)

  dt <- data.table::copy(obj$data)
  s <- obj$settings
  id_col <- s$id
  period_col <- s$period
  treat_col <- s$treatment
  out_col <- s$outcome
  elig_col <- s$eligible
  covariates <- s$covariates
  bl_covs <- s$baseline_covariates
  estimand <- s$estimand
  prefix <- s$prefix

  if (grace < 0) stop("grace must be non-negative", call. = FALSE)
  if (maxfollowup < 0) stop("maxfollowup must be non-negative", call. = FALSE)

  # Determine trial periods
  if (is.null(trials)) {
    trial_periods <- sort(unique(dt[get(elig_col) == 1][[period_col]]))
  } else {
    trial_periods <- sort(as.integer(trials))
  }

  n_trial_periods <- length(trial_periods)

  message("emulate_expand - Sequential Trial Expansion")
  message(strrep("-", 50))
  message("Estimand:        ", estimand)
  message("Trial periods:   ", n_trial_periods)
  message("Max follow-up:   ", if (maxfollowup > 0) paste(maxfollowup, "periods") else "unlimited")
  if (estimand != "ITT") message("Grace period:    ", grace, " periods")

  # Sort data
  data.table::setorderv(dt, c(id_col, period_col))

  # Process each trial period
  trial_list <- list()
  trial_count <- 0L
  used_trial_periods <- integer(0)

  fu_col <- paste0(prefix, "followup")
  trial_col <- paste0(prefix, "trial")
  arm_col <- paste0(prefix, "arm")
  cens_col <- paste0(prefix, "censored")
  outobs_col <- paste0(prefix, "outcome_obs")

  for (t in trial_periods) {
    # Find eligible individuals at period t
    elig_at_t <- dt[get(period_col) == t & get(elig_col) == 1][[id_col]]
    if (length(elig_at_t) == 0) next

    trial_count <- trial_count + 1L
    used_trial_periods <- c(used_trial_periods, t)

    # Keep obs for eligible IDs from period t onward
    trial_dt <- dt[get(id_col) %in% elig_at_t & get(period_col) >= t]

    # Baseline treatment at period t
    bl_treat <- trial_dt[get(period_col) == t,
                         .(bl_treat = get(treat_col)),
                         by = c(id_col)]
    trial_dt <- merge(trial_dt, bl_treat, by = id_col, all.x = TRUE)

    # Freeze covariates at baseline values
    all_covs <- c(covariates, bl_covs)
    if (length(all_covs) > 0) {
      data.table::setorderv(trial_dt, c(id_col, period_col))
      for (v in all_covs) {
        trial_dt[, (v) := get(v)[1L], by = c(id_col)]
      }
    }

    # Follow-up time
    trial_dt[, (fu_col) := get(period_col) - t]

    # Apply max follow-up
    if (maxfollowup > 0) {
      trial_dt <- trial_dt[get(fu_col) <= maxfollowup]
    }

    # Trial identifier
    trial_dt[, (trial_col) := t]

    if (estimand == "ITT") {
      # ITT: no cloning, use actual baseline treatment, no censoring
      trial_dt[, (arm_col) := bl_treat]
      trial_dt[, (cens_col) := 0L]
      trial_dt[, (outobs_col) := get(out_col)]
      trial_dt[, bl_treat := NULL]
      trial_list[[trial_count]] <- trial_dt
    } else {
      # PP/AT: clone into two arms
      arm1 <- data.table::copy(trial_dt)
      arm0 <- data.table::copy(trial_dt)

      arm1[, (arm_col) := 1L]
      arm0[, (arm_col) := 0L]

      # Censor arm 1 (treatment arm): censored when stops treatment
      arm1 <- .expand_censor(arm1, id_col, treat_col, out_col,
                             arm_val = 1L, grace = grace,
                             fu_col = fu_col, cens_col = cens_col,
                             outobs_col = outobs_col)

      # Censor arm 0 (control arm): censored when starts treatment
      arm0 <- .expand_censor(arm0, id_col, treat_col, out_col,
                             arm_val = 0L, grace = grace,
                             fu_col = fu_col, cens_col = cens_col,
                             outobs_col = outobs_col)

      arm1[, bl_treat := NULL]
      arm0[, bl_treat := NULL]

      trial_list[[trial_count]] <- data.table::rbindlist(list(arm1, arm0),
                                                         use.names = TRUE)
    }

    if (trial_count %% 10 == 0) {
      message("  ... processed ", trial_count, " of ", n_trial_periods, " trial periods")
    }
  }

  if (trial_count == 0) {
    stop("No trials produced; check eligibility criteria", call. = FALSE)
  }

  # Combine all trials
  expanded <- data.table::rbindlist(trial_list, use.names = TRUE, fill = TRUE)
  data.table::setorderv(expanded, c(trial_col, arm_col, id_col, fu_col))

  # Statistics
  n_orig <- nrow(dt)
  n_expanded <- nrow(expanded)
  n_treat <- sum(expanded[[arm_col]] == 1)
  n_control <- sum(expanded[[arm_col]] == 0)
  n_cens <- sum(expanded[[cens_col]] == 1)
  n_events <- sum(expanded[[outobs_col]] == 1)

  message("\nExpansion complete:")
  message("  Original obs:   ", format(n_orig, big.mark = ","))
  message("  Expanded obs:   ", format(n_expanded, big.mark = ","))
  message("  Expansion ratio:", sprintf(" %.1fx", n_expanded / n_orig))
  message("  Trials created: ", trial_count)
  message("  Treatment arm:  ", format(n_treat, big.mark = ","))
  message("  Control arm:    ", format(n_control, big.mark = ","))
  message("  Censored:       ", format(n_cens, big.mark = ","))
  message("  Events:         ", format(n_events, big.mark = ","))
  message(strrep("-", 50))

  # Update object
  obj$data <- expanded
  obj$state$expanded <- TRUE
  obj$expansion <- list(
    n_trials = trial_count,
    trial_periods = used_trial_periods,
    n_expanded = n_expanded,
    n_treat = n_treat,
    n_control = n_control,
    n_censored = n_cens,
    n_events = n_events,
    expansion_ratio = n_expanded / n_orig
  )

  invisible(obj)
}

#' Apply artificial censoring for one arm (internal helper)
#'
#' For a single treatment arm, identifies the first time point where the
#' individual deviates from the assigned strategy (accounting for the grace
#' period), marks that observation as censored, drops all subsequent rows,
#' and sets the outcome to 0 at the censoring time point.
#'
#' @param dt A \code{data.table} of trial data for one arm.
#' @param id_col Character: name of the patient identifier column.
#' @param treat_col Character: name of the treatment variable.
#' @param out_col Character: name of the outcome variable.
#' @param arm_val Integer (0 or 1): the assigned arm value.
#' @param grace Integer: the grace period in number of follow-up periods.
#' @param fu_col Character: name of the follow-up time column.
#' @param cens_col Character: name of the censoring indicator column to create.
#' @param outobs_col Character: name of the observed outcome column to create.
#'
#' @return The input \code{data.table}, modified in place, with rows after
#'   the censoring point removed and the \code{cens_col} and \code{outobs_col}
#'   columns populated.
#'
#' @keywords internal
.expand_censor <- function(dt, id_col, treat_col, out_col,
                           arm_val, grace, fu_col, cens_col, outobs_col) {
  # Initialize
  dt[, (cens_col) := 0L]
  dt[, (outobs_col) := 0L]

  data.table::setorderv(dt, c(id_col, fu_col))

  if (arm_val == 1L) {
    # Treatment arm: censored when STOPS treatment (treat==0 after grace)
    dt[, deviated := as.integer(get(treat_col) == 0 & get(fu_col) >= grace)]
  } else {
    # Control arm: censored when STARTS treatment (treat==1 after grace)
    dt[, deviated := as.integer(get(treat_col) == 1 & get(fu_col) >= grace)]
  }

  # Find first deviation time per individual
  dt[, first_dev := as.integer(deviated == 1 &
                                 (shift(deviated, 1L, fill = 0L, type = "lag") == 0 |
                                    seq_len(.N) == 1L)),
     by = c(id_col)]

  dt[first_dev == 1, cens_time := get(fu_col)]
  dt[, min_cens := min(cens_time, na.rm = TRUE), by = c(id_col)]
  # If no deviation, min will be Inf
  dt[is.infinite(min_cens), min_cens := NA_real_]

  # Mark censoring point
  dt[get(fu_col) == min_cens & !is.na(min_cens), (cens_col) := 1L]

  # Drop rows after censoring
  dt <- dt[is.na(min_cens) | get(fu_col) <= min_cens]

  # Outcome observed only if not censored
  dt[get(cens_col) == 0L, (outobs_col) := get(out_col)]

  # Clean up temp columns
  dt[, c("deviated", "first_dev", "cens_time", "min_cens") := NULL]

  dt
}
