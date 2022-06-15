context("Transforming and untransforming cues and sufficient statistics")

.cues <- c("cue1", "cue2")
.io <- example_MVG_ideal_observer(1)
.data <- make_MVG_data_from_model(model = .io, Ns = 50, keep.input_parameters = F)
test_that("transform_cues() - type check", {
  expect_true(is_tibble(transform_cues(data = .data, cues = .cues, center = T, scale = T)))
  expect_true(is.list(transform_cues(data = .data, cues = .cues, center = T, scale = T, return.transformed.data = F, return.transform.parameters = T)))
})


.data.transformed <- transform_cues(data = .data, cues = .cues, center = T, scale = T)
.transform <- transform_cues(data = .data, cues = .cues, center = T, scale = T, return.transformed.data = F, return.transform.parameters = T)
test_that("transform_cues() - do cue values have changed?", {
  expect_false(all(.data$cue1 == .data.transformed$cue1))
  expect_false(all(.data$cue2 == .data.transformed$cue2))
})


.data.untransformed <- untransform_cues(cues = .cues, transform.parameters = .transform)
test_that("(un)transform_cues() - does reverting to original cues work?", {
  expect_equal(.data$cue1, .data.untransformed$cue1)
  expect_equal(.data$cue2 == .data.untransformed$cue2)
})


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


test_that("untransform_category_mean - does back-transformation of category means work?", {
  expect_equal(
    .data %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      pull(x_mean),
    .data.transformed %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      pull(x_mean) %>%
      map(~ untransform_category_mean(.x, transform = .transform)))
})


test_that("untransform_category_cov - does back-transformation of category covariance matrix work?", {
  expect_equal(
    .data %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      pull(x_cov),
    .data.transformed %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      pull(x_cov) %>%
      map(~ untransform_category_cov(.x, transform = .transform)))
})
