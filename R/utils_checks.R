# Internal pipeline state validators
# These mirror Stata's _emulate_check_prepared, _emulate_check_expanded, etc.

.check_emulate <- function(obj) {
  if (!inherits(obj, "emulate")) {
    stop("Expected a 'emulate' object. Run emulate_prepare() first.", call. = FALSE)
  }
}

.check_prepared <- function(obj) {
  .check_emulate(obj)
  if (!isTRUE(obj$state$prepared)) {
    stop("Data has not been prepared. Run emulate_prepare() first.", call. = FALSE)
  }
}

.check_expanded <- function(obj) {
  .check_prepared(obj)
  if (!isTRUE(obj$state$expanded)) {
    stop("Data has not been expanded. Run emulate_expand() first.", call. = FALSE)
  }
}

.check_weighted <- function(obj) {
  .check_expanded(obj)
  if (!isTRUE(obj$state$weighted)) {
    stop("Data has not been weighted. Run emulate_weight() first.", call. = FALSE)
  }
}

.check_fitted <- function(obj) {
  .check_expanded(obj)
  if (!isTRUE(obj$state$fitted)) {
    stop("Model has not been fitted. Run emulate_fit() first.", call. = FALSE)
  }
}
