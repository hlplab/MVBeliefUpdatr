context("fit stanfit")

source("../functions-to-make-or-load-models.R")
.data <- make_data_for_stanfit(5)
.exposure <- .data %>% filter(Phase == "exposure")
.test <- .data %>% filter(Phase == "test")

test_that("Test compose stanfit", {
  expect_no_error(
    infer_NIW_ideal_adaptor(
      exposure = .exposure,
      test = .test,
      cues = c("cue1"),
      category = "category",
      response = "Response",
      group = "Subject",
      group.unique = "Condition",
      sample = F))
  expect_no_error(
    infer_NIW_ideal_adaptor(
      exposure = .exposure,
      test = .test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "Response",
      group = "Subject",
      group.unique = "Condition",
      sample = F))
})


test_that("Test fitting stanfit", {
  # check that default works
  expect_no_error(
    suppressMessages(suppressWarnings(
      infer_NIW_ideal_adaptor(
        exposure = .exposure,
        test = .test,
        cues = c("cue1"),
        category = "category",
        response = "Response",
        group = "Subject",
        group.unique = "Condition",
        sample = T,
        cores = 4,
        refresh = -1,
        iter = 100))))
  expect_no_error(
    suppressMessages(suppressWarnings(
      infer_NIW_ideal_adaptor(
        exposure = .exposure,
        test = .test,
        cues = c("cue1", "cue2"),
        category = "category",
        response = "Response",
        group = "Subject",
        group.unique = "Condition",
        sample = T,
        cores = 4,
        refresh = -1,
        iter = 100))))
  # check that forcing multivariate updating works even if there is only one cue
  expect_no_error(
    suppressMessages(suppressWarnings(
      infer_NIW_ideal_adaptor(
        exposure = .exposure,
        test = .test,
        cues = c("cue1"),
        category = "category",
        response = "Response",
        group = "Subject",
        group.unique = "Condition",
        sample = T,
        use_univariate_updating = F,
        cores = 4,
        refresh = -1,
        iter = 100))))
  # forcing univariate updating should throw an error if and only if there is more than 1 cue
  expect_no_error(
    suppressMessages(suppressWarnings(
      infer_NIW_ideal_adaptor(
        exposure = .exposure,
        test = .test,
        cues = c("cue1"),
        category = "category",
        response = "Response",
        group = "Subject",
        group.unique = "Condition",
        sample = T,
        use_univariate_updating = T,
        cores = 4,
        refresh = -1,
        iter = 100))))
  expect_error(
    suppressWarnings(
      infer_NIW_ideal_adaptor(
        exposure = .exposure,
        test = .test,
        cues = c("cue1", "cue2"),
        category = "category",
        response = "Response",
        group = "Subject",
        group.unique = "Condition",
        sample = T,
        use_univariate_updating = T,
        cores = 4,
        refresh = -1,
        iter = 100)))
  })

.data <- make_data_for_stanfit(2)
.exposure <- .data %>% filter(Phase == "exposure")
.test <- .data %>% filter(Phase == "test")
test_that("stanfit output", {
  expect_no_error(
    suppressMessages(suppressWarnings(
      infer_NIW_ideal_adaptor(
        exposure = .exposure,
        test = .test,
        cues = c("VOT", "f0_semitones"),
        category = "category",
        response = "Response",
        group = "Subject",
        group.unique = "Condition",
        sample = T,
        cores = 4,
        refresh = -1,
        iter = 100))))
})
