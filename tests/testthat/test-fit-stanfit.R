source("../functions-to-make-or-load-models.R")

.silent <- 0
.verbose <- F
.chains <- 1
.warmup <- 100
.iter <- 200
.file_refit <- "always"

# Note: Many of the tests below expect warnings because it would be too time-consuming to fit the models with adequate warmup, etc.
# (so I'm fitting them with very low warmup, and swallow the warnings)

context("to_array")

x <- NULL
test_that("NULL inputs with simplify = TRUE", {
  expect_identical(
    to_array(x, simplify = TRUE), array(numeric()))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = TRUE),
    array(numeric()))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = 1, simplify = TRUE),
    array(numeric(), dim = c(0)))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = 1, simplify = TRUE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = 2, simplify = TRUE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = 1, simplify = TRUE),
    array(numeric(), dim = c(0)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(1, 1), simplify = TRUE),
    array(numeric(), dim = c(0, 0)))
})

test_that("NULL inputs with simplify = FALSE", {
  expect_identical(
    to_array(x, simplify = FALSE),
    array(numeric()))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = FALSE),
    array(numeric()))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = 1, simplify = FALSE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = 1, simplify = FALSE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = 2, simplify = FALSE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = 1, simplify = FALSE),
    array(numeric(), dim = c(0, 0, 0)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(1, 1), simplify = FALSE),
    array(numeric(), dim = c(0, 0, 0, 0)))
})

x <- 5
test_that("scalar inputs with simplify = TRUE", {
  expect_identical(
    to_array(x, simplify = TRUE),
    x)
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = TRUE),
    x)
  expect_error(
    to_array(x, inner_dims = 3, outer_dims = NULL, simplify = TRUE))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = 3, simplify = TRUE),
    array(x, dim = c(3)))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = TRUE))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = TRUE))
  # No simplification for atomic values
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = TRUE),
    array(x, dim = c(3)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = TRUE),
    array(x, dim = c(3, 4)))
})

test_that("scalar inputs with simplify = FALSE", {
  expect_identical(
    to_array(x, simplify = FALSE),
    array(x))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = FALSE),
    array(x))
  expect_error(
    to_array(x, inner_dims = 3, outer_dims = NULL, simplify = FALSE))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = 3, simplify = FALSE),
    array(x, dim = c(3, 1)))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = FALSE))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = FALSE))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = FALSE),
    array(x, dim = c(3, 1, 1)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = FALSE),
    array(x, dim = c(3, 4, 1, 1)))
})

x <- 1:2
test_that("vector inputs with simplify = TRUE", {
  expect_identical(
    to_array(x, simplify = TRUE), array(x))
  expect_error(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = TRUE))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = NULL, simplify = TRUE),
    array(x))
  expect_equal(
    dim(to_array(x, inner_dims = 2, outer_dims = 3, simplify = TRUE)),
    c(3, 2))
  # No simplification for atomic values
  expect_error(
    to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = TRUE))
  expect_error(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = TRUE))
})

test_that("vector inputs with simplify = FALSE", {
  expect_identical(
    to_array(x, simplify = FALSE), array(x))
  expect_error(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = FALSE))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = NULL, simplify = FALSE),
    array(x))
  expect_equal(
    dim(to_array(x, inner_dims = 2, outer_dims = 3, simplify = FALSE)),
    c(3, 2))
  expect_error(
    to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = FALSE))
  expect_error(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = FALSE))
})

x <- matrix(1, nrow = 1)
test_that("matrix inputs with simplify = TRUE", {
  # No simplification for atomic values
  expect_equal(
    to_array(x, simplify = TRUE),
    1)
  # Coercion
  expect_equal(
   to_array(x, inner_dims = 1, outer_dims = NULL, simplify = TRUE),
    1)
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = NULL, simplify = TRUE))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = TRUE))
  # No simplification for atomic values
  expect_equal(
    dim(to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = TRUE)),
    c(3))
  expect_equal(
    dim(to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = TRUE)),
    c(3, 4))
})

test_that("matrix inputs with simplify = FALSE", {
  # No simplification for atomic values
  expect_equal(
    dim(to_array(x, simplify = FALSE)),
    c(1, 1))
  # Coercion
  expect_equal(
    dim(to_array(x, inner_dims = 1, outer_dims = NULL, simplify = FALSE)),
    c(1))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = NULL, simplify = FALSE))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = FALSE))
  expect_equal(
    dim(to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = FALSE)),
    c(3, 1, 1))
  expect_equal(
    dim(to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = FALSE)),
    c(3, 4, 1, 1))
})

# ADD LIST INPUT TEST

context("make_staninput_for_ideal_adaptor")

# Test whether input formatting works before testing fitting
test_that("test make_staninput_for_NIX_ideal_adaptor (one cue)", {
  expect_no_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "center"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "PCA whiten"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "other"))
})

test_that("test make_staninput_for_NIW_ideal_adaptor (one cue)", {
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "center"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "other"))
})

test_that("test make_staninput_for_MNIX_ideal_adaptor (one cue)", {
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "identity"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "center"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "PCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "other"))
})

test_that("test make_staninput_for_NIX_ideal_adaptor (two cues)", {
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "identity"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "center"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "PCA whiten"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "other"))
})

test_that("test make_staninput_for_NIW_ideal_adaptor (two cues)", {
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "center"))
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize"))
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten"))
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "other"))
})

test_that("test make_staninput_for_MNIX_ideal_adaptor (two cues)", {
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "center"))
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize"))
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "PCA whiten"))
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "other"))
})

# Test whether conditions with empty exposure information also work
test_that("test make_staninput_for_NIX_ideal_adaptor (two cues)", {
  # expect_no_error(get_example_staninput(4, stanmodel = "NIX_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(4, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize"))
})

test_that("test make_staninput_for_NIW_ideal_adaptor (two cues)", {
  # expect_no_error(get_example_staninput(5, stanmodel = "NIW_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(5, stanmodel = "NIW_ideal_adaptor", transform_type = "ZCA whiten"))
})

test_that("test make_staninput_for_MNIX_ideal_adaptor (two cues)", {
  # expect_no_error(get_example_staninput(5, stanmodel = "MNIX_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(5, stanmodel = "MNIX_ideal_adaptor", transform_type = "ZCA whiten"))
})


context("get_category_statistics_as_list_of_arrays")

# FOR NOW, ONLY SOME VALUE CHECKS ARE INCLUDED HERE, AND ONLY FOR NIW_IDEAL_ADAPTORS
# COULD ADD MORE TESTS OF THIS TYPE BELOW, E.G., FOR NIX AND MNIX
.staninput <- get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "identity")
test_that("test return values of get_category_statistics_as_list_of_arrays (NIW_ideal_adaptor)", {
  expect_true(
    {
      .x_mean_exposure <-
        .staninput$data %>%
        filter(Phase == "exposure") %>%
        group_by(Condition, category) %>%
        summarise(x_mean_exposure = list(c(mean(VOT), mean(f0_semitones)))) %>%
        ungroup() %>%
        pull(x_mean_exposure)
      .match <- c()
      for (j in 1:.staninput$staninput$transformed$L)
        for (i in 1:.staninput$staninput$transformed$M) {
          .match %<>%
            append(
              all.equal(
                .staninput$staninput$untransformed$x_mean_exposure[i,j,],
                .x_mean_exposure[[(j - 1) * 2 + i]]))
        }
      all(.match)
    })
  expect_true(
    {
      .x_ss_exposure <-
        .staninput$data %>%
        filter(Phase == "exposure") %>%
        group_by(Condition, category) %>%
        group_map(.f = ~ get_sum_of_uncentered_squares_from_df(.x %>% select(VOT, f0_semitones) %>% as.matrix()))
      .match <- c()
      for (j in 1:.staninput$staninput$transformed$L)
        for (i in 1:.staninput$staninput$transformed$M) {
          .match %<>%
            append(
              all.equal(
                .staninput$staninput$untransformed$x_ss_exposure[i,j,,],
                .x_ss_exposure[[(j - 1) * 2 + i]]))
        }
      all(.match)
    })
})

# Test fitting
context("fit_ideal_adaptor (unknown parameters)")

# Running with low sampling sizes, which should elicit warnings (insufficient samples but no errors)

test_that("test fitting NIX (one cue)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "identity",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "center",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "PCA whiten",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "ZCA whiten",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit)))
})

test_that("test fitting NIW (one cue)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "identity",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "center",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "ZCA whiten",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting NIX (two+ cues)", {
  expect_error(fit <- get_example_stanfit(2, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_error(fit <- get_example_stanfit(3, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
})

test_that("test fitting NIW (two+ cues)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(3, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(3, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting MNIX (two+ cues)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "PCA whiten",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(3, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
  expect_warning(expect_no_error(fit <- get_example_stanfit(3, stanmodel = "MNIX_ideal_adaptor", transform_type = "PCA whiten",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})



context("fit_ideal_adaptor (known lapse rate)")

# Get MVG ideal observers that the example data is generated from to get the ground truth for the parameters
m_prior1 <- example_MVG_ideal_observer(1)
m_prior2 <- example_MVG_ideal_observer(2)

test_that("test fitting NIX (one cue)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                                            lapse_rate = unique(m_prior1$lapse_rate),
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting NIW (one cue)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            lapse_rate = unique(m_prior1$lapse_rate),
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting NIW (two cues)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            lapse_rate = unique(m_prior2$lapse_rate),
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting MNIX (two cues)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize",
                                                            lapse_rate = unique(m_prior2$lapse_rate),
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})


context("fit_ideal_adaptor (known mu_0)")

test_that("test fitting NIX (one cue)", {
  expect_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                          mu_0 = m_prior2$mu,
                                          filename = "temp",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                                            mu_0 = m_prior1$mu,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting NIW (one cue)", {
  expect_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                          mu_0 = m_prior2$mu,
                                          filename = "temp",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_warning(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                            mu_0 = m_prior1$mu,
                                            filename = "temp",
                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                            file_refit = .file_refit))
})

test_that("test fitting NIW (two cues)", {
  expect_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                          mu_0 = m_prior1$mu,
                                          filename = "temp",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            mu_0 = m_prior2$mu,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting MNIX (two cues)", {
  expect_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize",
                                          mu_0 = m_prior1$mu,
                                          filename = "temp",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize",
                                                            mu_0 = m_prior2$mu,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})


context("fit_ideal_adaptor (known Sigma_0)")

test_that("test fitting NIX (one cue)", {
  expect_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                          Sigma_0 = m_prior2$Sigma,
                                          filename = "temp",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                                            Sigma_0 = m_prior1$Sigma,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting NIW (one cue)", {
  expect_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                          Sigma_0 = m_prior2$Sigma,
                                          filename = "temp",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            Sigma_0 = m_prior1$Sigma,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting NIW (two cues)", {
  expect_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                          Sigma_0 = m_prior1$Sigma,
                                          filename = "temp",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            Sigma_0 = m_prior2$Sigma,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting MNIX (two cues)", {
  expect_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize",
                                          Sigma_0 = m_prior1$Sigma,
                                          filename = "temp",
                                          warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                          file_refit = .file_refit))
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize",
                                                            Sigma_0 = m_prior2$Sigma,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})



context("fit_ideal_adaptor (known mu_0 and Sigma_0)")

test_that("test fitting NIX (one cue)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize",
                                                            mu_0 = m_prior1$mu,
                                                            Sigma_0 = m_prior1$Sigma,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting NIW (one cue)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            mu_0 = m_prior1$mu,
                                                            Sigma_0 = m_prior1$Sigma,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting NIW (two cues)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize",
                                                            mu_0 = m_prior2$mu,
                                                            Sigma_0 = m_prior2$Sigma,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})

test_that("test fitting MNIX (two cues)", {
  expect_warning(expect_no_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize",
                                                            mu_0 = m_prior2$mu,
                                                            Sigma_0 = m_prior2$Sigma,
                                                            filename = "temp",
                                                            warmup = .warmup, iter = .iter, chains = .chains, cores = .chains, silent = .silent, verbose = .verbose,
                                                            file_refit = .file_refit)))
})
