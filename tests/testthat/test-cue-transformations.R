context("Transforming and untransforming cues and sufficient statistics")

.io <- example_MVG_ideal_observer(3)
.cues <- get_cue_labels_from_model(.io)
.data <- sample_data_from_model(model = .io, Ns = 50)

test_that("transform_cues - input (data)", {
  expect_error(
    transform_cues(
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = NA,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = NULL,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = "data",
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
})

test_that("transform_cues - input (cues)", {
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1:2],
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
      cues = c("not a cue"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
})

test_that("transform_cues - input (center, scale, pca, attach)", {
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T))
})

test_that("transform_cues - basic output", {
  expect_message(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_true(
    is.data.frame(
      transform_cues(
        data = .data,
        cues = .cues[1],
        return.transformed.data = T,
        return.transform.parameters = F,
        return.transform.function = F,
        return.untransform.function = F)))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
        cues = .cues[1],
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
        cues = .cues[1],
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
        cues = .cues[1],
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
        cues = .cues[1],
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
        cues = .cues[1],
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
      cues = .cues[1],
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
      cues = .cues[1],
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
  expect_false(all(.data[[.cues[1]]] == .data.transformed[[.cues[1]]]))
  expect_false(all(.data[[.cues[2]]] == .data.transformed[[.cues[2]]]))
  expect_false(all(.data[[.cues[3]]] == .data.transformed[[.cues[3]]]))
})

test_that("transform_cues - output equality of input and un-transformation", {
  # one cue, without PCA
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  # two cues, with PCA
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[2]]], .data[[.cues[2]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[2]]], .data[[.cues[2]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[2]]], .data[[.cues[2]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[2]]], .data[[.cues[2]]])
})


test_that("untransform_category_mean - does back-transformation of category means work?", {
  # one cue
  expect_equal(
    .data %>%
      get_sufficient_category_statistics(cues = .cues[1]) %>%
      pull(x_mean),
    transform_cues(data = .data, cues = .cues[1], center = T, scale = T, return.transformed.data = T, return.transform.parameters = F) %>%
      get_sufficient_category_statistics(cues = .cues[1]) %>%
      pull(x_mean) %>%
      map(~ untransform_category_mean(.x, transform = transform_cues(data = .data, cues = .cues[1], center = T, scale = T, return.transformed.data = F, return.transform.parameters = T))))
  # two cues
  expect_equal(
    .data %>%
      get_sufficient_category_statistics(cues = .cues[1:2]) %>%
      pull(x_mean),
    transform_cues(data = .data, cues = .cues[1:2], center = T, scale = T, return.transformed.data = T, return.transform.parameters = F) %>%
      get_sufficient_category_statistics(cues = .cues[1:2]) %>%
      pull(x_mean) %>%
      map(~ untransform_category_mean(.x, transform = transform_cues(data = .data, cues = .cues[1:2], center = T, scale = T, return.transformed.data = F, return.transform.parameters = T))))
  # three cues
  expect_equal(
    .data %>%
      get_sufficient_category_statistics(cues = .cues[1:3]) %>%
      pull(x_mean),
    transform_cues(data = .data, cues = .cues[1:3], center = T, scale = T, return.transformed.data = T, return.transform.parameters = F) %>%
      get_sufficient_category_statistics(cues = .cues[1:3]) %>%
      pull(x_mean) %>%
      map(~ untransform_category_mean(.x, transform = transform_cues(data = .data, cues = .cues[1:3], center = T, scale = T, return.transformed.data = F, return.transform.parameters = T))))
})

test_that("untransform_category_cov - does back-transformation of category covariance matrix work?", {
  # one cue
  expect_equivalent(
    .data %>%
      get_sufficient_category_statistics(cues = .cues[1]) %>%
      pull(x_cov),
    transform_cues(data = .data, cues = .cues[1], center = T, scale = T, return.transformed.data = T, return.transform.parameters = F) %>%
      get_sufficient_category_statistics(cues = .cues[1]) %>%
      pull(x_cov) %>%
      map(~ untransform_category_cov(.x, transform = transform_cues(data = .data, cues = .cues[1], center = T, scale = T, return.transformed.data = F, return.transform.parameters = T))))
  # two cues
  expect_equivalent(
    .data %>%
      get_sufficient_category_statistics(cues = .cues[1:2]) %>%
      pull(x_cov),
    transform_cues(data = .data, cues = .cues[1:2], center = T, scale = T, return.transformed.data = T, return.transform.parameters = F) %>%
      get_sufficient_category_statistics(cues = .cues[1:2]) %>%
      pull(x_cov) %>%
      map(~ untransform_category_cov(.x, transform = transform_cues(data = .data, cues = .cues[1:2], center = T, scale = T, return.transformed.data = F, return.transform.parameters = T))))
  # three cues
  expect_equivalent(
    .data %>%
      get_sufficient_category_statistics(cues = .cues[1:3]) %>%
      pull(x_cov),
    transform_cues(data = .data, cues = .cues[1:3], center = T, scale = T, return.transformed.data = T, return.transform.parameters = F) %>%
      get_sufficient_category_statistics(cues = .cues[1:3]) %>%
      pull(x_cov) %>%
      map(~ untransform_category_cov(.x, transform = transform_cues(data = .data, cues = .cues[1:3], center = T, scale = T, return.transformed.data = F, return.transform.parameters = T))))
})

