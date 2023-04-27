context("cov, css, uss")

test_that("uss2css, css2cov - does sum-of-square to cov conversion work?", {
  expect_equal(
    .data %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      pull(x_cov),
    .data %>%
      mutate(x_css_from_uss = pmap(list(x_uss, x_N, x_mean), ~ uss2css(..1, ..2, ..3)),
             x_cov_from_uss = map2(x_css_from_uss, x_N, ~ css2cov(.x, n = .y))) %>%
      pull(x_cov_from_uss))
})
