# Harrell restricted cubic spline basis
# Matches Stata's _emulate_natural_spline.ado exactly.
# Do NOT use splines::ns() — different basis, different coefficients.

#' Compute Harrell RCS basis at given x values with given knots
#'
#' @param x Numeric vector of values
#' @param knots Numeric vector of knot positions (length = df + 1).
#'   Includes boundary knots at min/max.
#' @param df Degrees of freedom
#' @return Matrix with df columns (basis variables)
#' @keywords internal
.emulate_rcs_basis <- function(x, knots, df) {
  # Basis 1 is always linear (x itself)
  basis <- matrix(NA_real_, nrow = length(x), ncol = df)
  basis[, 1] <- x

  if (df == 1) return(basis)

  n_knots <- df + 1L
  t_last <- knots[n_knots]
  n_internal <- df - 1L
  t_pen <- knots[n_internal + 1L]  # second-to-last knot = knot[n_internal]
  # Note: knots[1] = boundary min, knots[2..n_internal+1] = internal,
  #       knots[n_knots] = boundary max

  if (n_internal == 1L) {
    # df=2: single nonlinear basis
    # d_1(x) = ((x - k1)_+^3 - (x - k_last)_+^3) / (k_last - k1)
    k1 <- knots[2]  # the one internal knot
    basis[, 2] <- (pmax(0, x - k1)^3 - pmax(0, x - t_last)^3) /
      (t_last - k1)
  } else {
    # df >= 3: Harrell RCS formula
    # h_j(x) = d_j(x) - d_pen(x)  for j = 1, ..., n_internal
    # (k-2 nonlinear bases for k = df+1 knots, matching Harrell)
    n_nonlinear <- n_internal
    # d_pen(x) uses penultimate knot
    d_pen <- (pmax(0, x - t_pen)^3 - pmax(0, x - t_last)^3) /
      (t_last - t_pen)

    for (j in seq_len(n_nonlinear)) {
      k_j <- knots[j + 1L]  # internal knots at positions 2..n_internal+1
      d_j <- (pmax(0, x - k_j)^3 - pmax(0, x - t_last)^3) /
        (t_last - k_j)
      basis[, j + 1L] <- d_j - d_pen
    }
  }

  basis
}

#' Compute knot positions for natural splines (matching Stata)
#'
#' @param x Numeric vector
#' @param df Degrees of freedom
#' @return Numeric vector of knot positions (length df + 1)
#' @keywords internal
.emulate_compute_knots <- function(x, df) {
  x <- x[!is.na(x)]
  x_min <- min(x)
  x_max <- max(x)

  if (x_max == x_min) stop("Variable has no variation", call. = FALSE)

  if (df == 1) return(c(x_min, x_max))

  n_internal <- df - 1L
  # Internal knots at equally spaced quantiles
  pcts <- seq_len(n_internal) / (n_internal + 1L)
  internal <- as.numeric(quantile(x, probs = pcts, type = 2))
  # type=2 matches Stata's _pctile (uses (lower+upper)/2 at discontinuities)

  c(x_min, internal, x_max)
}

#' Generate natural spline basis variables and add to data.table
#'
#' @param dt data.table (modified in place)
#' @param varname Column name of the variable to spline
#' @param df Degrees of freedom
#' @param prefix Prefix for generated column names
#' @param touse Optional logical vector for subsetting when computing knots
#' @return List with knots, df, and variable names
#' @keywords internal
.emulate_natural_spline <- function(dt, varname, df, prefix, touse = NULL) {
  x <- dt[[varname]]
  x_for_knots <- if (!is.null(touse)) x[touse] else x

  knots <- .emulate_compute_knots(x_for_knots, df)
  basis <- .emulate_rcs_basis(x, knots, df)

  varnames <- paste0(prefix, seq_len(df))
  for (j in seq_len(df)) {
    data.table::set(dt, j = varnames[j], value = basis[, j])
  }

  list(knots = knots, df = df, varnames = varnames)
}
