test_that(".emulate_compute_knots places knots correctly", {
  x <- 0:100
  knots <- emulate:::.emulate_compute_knots(x, df = 3)
  expect_length(knots, 4)  # df + 1
  expect_equal(knots[1], 0)
  expect_equal(knots[4], 100)
  expect_true(knots[2] > 0 && knots[2] < knots[3])
  expect_true(knots[3] < 100)
})

test_that(".emulate_rcs_basis returns correct dimensions", {
  x <- 0:20
  knots <- emulate:::.emulate_compute_knots(x, df = 3)
  basis <- emulate:::.emulate_rcs_basis(x, knots, df = 3)
  expect_equal(nrow(basis), 21)
  expect_equal(ncol(basis), 3)

  # First column should be x itself
  expect_equal(basis[, 1], as.numeric(x))
})

test_that(".emulate_rcs_basis df=1 returns linear only", {
  x <- 0:10
  knots <- emulate:::.emulate_compute_knots(x, df = 1)
  basis <- emulate:::.emulate_rcs_basis(x, knots, df = 1)
  expect_equal(ncol(basis), 1)
  expect_equal(basis[, 1], as.numeric(x))
})

test_that(".emulate_rcs_basis df=2 returns linear + 1 nonlinear", {
  x <- 0:20
  knots <- emulate:::.emulate_compute_knots(x, df = 2)
  basis <- emulate:::.emulate_rcs_basis(x, knots, df = 2)
  expect_equal(ncol(basis), 2)
  # Nonlinear column should be 0 for x < first internal knot
  expect_true(basis[1, 2] == 0)
})

test_that("Harrell RCS is continuous", {
  x <- seq(0, 20, by = 0.1)
  knots <- emulate:::.emulate_compute_knots(x, df = 3)
  basis <- emulate:::.emulate_rcs_basis(x, knots, df = 3)

  # Check no NAs and columns are finite
  expect_true(all(!is.na(basis)))
  expect_true(all(is.finite(basis)))
  # First column is linear
  expect_equal(basis[, 1], x)
})
