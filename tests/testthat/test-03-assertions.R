test_that(".assert_that handles simple and multi-condition checks", {
  expect_true(MVBeliefUpdatr:::.assert_that(TRUE))
  expect_true(MVBeliefUpdatr:::.assert_that(1 > 0, 2 > 1))
  expect_error(
    MVBeliefUpdatr:::.assert_that(FALSE, msg = "custom message"),
    "custom message"
  )
})

test_that(".assert_true and .assert_false provide clear wrappers", {
  expect_true(MVBeliefUpdatr:::.assert_true(1 == 1, "must hold"))
  expect_true(MVBeliefUpdatr:::.assert_false(FALSE, "must be false"))
  expect_error(
    MVBeliefUpdatr:::.assert_true(FALSE, "must fail"),
    "must fail"
  )
  expect_error(
    MVBeliefUpdatr:::.assert_false(TRUE, "must fail"),
    "must fail"
  )
})
