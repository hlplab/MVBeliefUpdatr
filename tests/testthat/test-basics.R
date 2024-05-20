#' @import tidyverse
require(tidyverse)

context("cov, css, uss")

.cues <- c("cue1", "cue2")
.io <- example_MVG_ideal_observer(5)
.data <- sample_MVG_data_from_model(model = .io, Ns = 50, keep.input_parameters = F)

test_that("uss2css, css2cov - does sum-of-square to cov conversion work?", {
  expect_equivalent(
    .data %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      pull(x_cov),
    .data %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      mutate(x_css_from_uss = pmap(list(x_uss, x_N, x_mean), ~ uss2css(..1, ..2, ..3)),
             x_cov_from_uss = map2(x_css_from_uss, x_N, ~ css2cov(.x, n = .y))) %>%
      pull(x_cov_from_uss))
})


