context("Transforming and untransforming cues and sufficient statistics")

.cues <- c("cue1", "cue2")
.io <- example_MVG_ideal_observer(1)
.data <- sample_MVG_data_from_model(model = .io, Ns = 50, keep.input_parameters = F)

test_that("transform_cues - input (data)", {
  expect_error(
    transform_cues(
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = NA,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = NULL,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = "data",
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
})

test_that("transform_cues - input (cues)", {
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue2"),
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1", "cue2"),
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = cue1,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = NA,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = NULL,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = c("cue3"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
})

test_that("transform_cues - input (center, scale, pca, attach)", {
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = F,
      scale = F,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = F,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = F,
      scale = T,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = F,
      scale = F,
      pca = T,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = F,
      scale = F,
      pca = F,
      attach = T,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = T,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = T,
      attach = T,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = "yes",
      scale = F,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = F,
      scale = "yes",
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = F,
      scale = F,
      pca = "yes",
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = F,
      scale = F,
      pca = F,
      attach = "yes",
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
})

test_that("transform_cues - input (return arguments)", {
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = c("cue1"),
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T))
})

test_that("transform_cues - basic output", {
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
})

.data.transformed <- transform_cues(data = .data, cues = .cues, center = T, scale = T)
.transform <- transform_cues(data = .data, cues = .cues, center = T, scale = T, return.transformed.data = F, return.transform.parameters = T)
test_that("transform_cues - do transformed cue values have changed?", {
  expect_false(all(.data$cue1 == .data.transformed$cue1))
  expect_false(all(.data$cue2 == .data.transformed$cue2))
})

test_that("transform_cues - output equality of input and un-transformation", {
  # one cue, without PCA
  tc <-
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue1, .data$cue1)
  tc <-
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue1, .data$cue1)
  tc <-
    transform_cues(
      data = .data,
      cues = c("cue1"),
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue1, .data$cue1)
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
  # two cues, with PCA
  tc <-
    transform_cues(
      data = .data,
      cues = c("cue1", "cue2"),
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue1, .data$cue1)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue2, .data$cue2)
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
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue2, .data$cue2)
  tc <-
    transform_cues(
      data = .data,
      cues = c("cue1", "cue2"),
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue1, .data$cue1)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue2, .data$cue2)
  tc <-
    transform_cues(
      data = .data,
      cues = c("cue1", "cue2"),
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue1, .data$cue1)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))$cue2, .data$cue2)
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
  expect_equivalent(
    .data %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      pull(x_cov),
    .data.transformed %>%
      get_sufficient_category_statistics(cues = .cues) %>%
      pull(x_cov) %>%
      map(~ untransform_category_cov(.x, transform = .transform)))
})
