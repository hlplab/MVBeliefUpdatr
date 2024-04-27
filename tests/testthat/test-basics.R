#' @import tidyverse
require(tidyverse)

context("cov, css, uss")

.cues <- c("cue1", "cue2")
.io <- example_MVG_ideal_observer(1)
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

test_that("transform_cues output", {
  expect_message(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_true(
    is.data.frame(
      transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F)))
  expect_true(
    is.list(
      transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F)))
  expect_true(
    is.function(
      transform_cues(
        data = .data,
        cues = c("cue1"),
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = F,
        return.transform.parameters = F,
        return.transform.function = T,
        return.untransform.function = F)))
  expect_true(
    is.function(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = T)))
  expect_true(
    is.list(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F)))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = c("cue1"),
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = F)))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = c("cue1"),
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = T)))
  expect_true(
    is.data.frame(
      transform_cues(
        data = .data,
        cues = c("cue1"),
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = T)[["data"]]))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = c("cue1"),
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = T)[["transform.parameters"]]))
  expect_true(
    is.function(
      transform_cues(
        data = .data,
        cues = c("cue1"),
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = T)[["transform.function"]]))
  expect_true(
    is.function(
      transform_cues(
        data = .data,
        cues = c("cue1"),
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = T)[["untransform.function"]]))
  tc <-
    transform_cues(
    data = .data,
    cues = c("cue1"),
    center = T,
    scale = T,
    pca = F,
    return.transformed.data = T,
    return.transform.parameters = T,
    return.transform.function = T,
    return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue1, .data$cue1)
  tc <-
    transform_cues(
      data = .data,
      cues = c("cue1", "cue2"),
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue1, .data$cue1)
})
