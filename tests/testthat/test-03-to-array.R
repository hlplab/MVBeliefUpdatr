test_that("to_array converts vectors into arrays with the requested inner dimensions", {
  x <- to_array(c(1, 2), inner_dims = 2, outer_dims = 2, simplify = FALSE)

  expect_equal(dim(x), c(2, 2))
  expect_equal(x, matrix(c(1, 2, 1, 2), nrow = 2, byrow = TRUE))
})

test_that("to_array stacks list inputs into the expected higher-dimensional array", {
  x <- to_array(
    list(matrix(1:4, nrow = 2, byrow = TRUE), matrix(5:8, nrow = 2, byrow = TRUE)),
    inner_dims = c(2, 2),
    outer_dims = 2,
    simplify = FALSE
  )

  expect_equal(dim(x), c(2, 2, 2))
  expect_equal(x[, , 1], matrix(1:4, nrow = 2, byrow = TRUE))
  expect_equal(x[, , 2], matrix(5:8, nrow = 2, byrow = TRUE))
})
