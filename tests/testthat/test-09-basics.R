.io <- suppressMessages(suppressWarnings(example_mvg_ideal_observer(n_cues = 2)))
.cues <- get_cue_labels(.io)
.data <- sample_observations(.io, Ns = 50)

test_that("uss2css, css2cov - does sum-of-square to cov conversion work?", {
  expect_equal(
    .data %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      dplyr::pull(x_cov),
    .data %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      dplyr::mutate(
        x_css_from_uss = purrr::pmap(list(x_uss, x_N, x_mean), ~ uss2css(..1, ..2, ..3)),
        x_cov_from_uss = purrr::map2(x_css_from_uss, x_N, ~ css2cov(.x, n = .y))
      ) %>%
      dplyr::pull(x_cov_from_uss),
    ignore_attr = TRUE
  )
})

