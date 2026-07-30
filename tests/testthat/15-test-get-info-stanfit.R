
source("../functions-to-make-or-load-models.R")

context("get information from stanfit NIW (one cue)")

fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", file_refit = "never")
test_that("Test is.ideal_adaptor_stanfit", {
  expect_false(is.ideal_adaptor_stanfit(NULL))
  expect_false(is.ideal_adaptor_stanfit(NA))
  expect_false(is.ideal_adaptor_stanfit(1))
  expect_false(is.ideal_adaptor_stanfit("1"))
  expect_false(is.ideal_adaptor_stanfit(TRUE))
  expect_false(is.ideal_adaptor_stanfit(list(1)))
  expect_false(is.ideal_adaptor_stanfit(example_exemplar_model(1)))
  expect_false(is.ideal_adaptor_stanfit(example_MVG_ideal_observer(1)))
  expect_false(is.ideal_adaptor_stanfit(example_NIW_ideal_adaptor(1)))
  expect_true(is.ideal_adaptor_stanfit(fit))
})

test_that("add draws - input check (1 cue)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "plus20")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "plus20"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (1 cue)", {
  expect_equal(nrow(get_draws(fit, groups = "prior", ndraws = 10, seed = 1, wide = F) %>% distinct(.draw)), 10)
  expect_equal(nrow(get_draws(fit, groups = "prior", wide = F, summarize = T)), 2)
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = F)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "cue", "cue2", "m", "S", "lapse_rate"))
})


context("get information from stanfit NIW (two cues)")

fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize", file_refit = "never")
test_that("add draws - input check (2 cues)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "plus20.20")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "plus20.20"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (2 cues)", {
  expect_equal(nrow(get_draws(fit, groups = "prior", ndraws = 10, seed = 1, wide = F) %>% distinct(.draw)), 10)
  expect_equal(nrow(get_draws(fit, groups = "prior", wide = F, summarize = T)), 2)
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = F)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "cue", "cue2", "m", "S", "lapse_rate"))
})

test_that("get exposure category statistic", {
  # Get error when *exposure* statistics is requested for *prior*
  expect_error(get_exposure_category_statistic(fit, groups = "prior"))
  expect_error(get_exposure_category_mean(fit, groups = "prior"))
  expect_error(get_exposure_category_cov(fit, groups = "prior"))
  # When prior is not requested
  # get_exposure_mean
  expect_true(is.vector(get_exposure_category_mean(fit, "/b/", "baseline")))
  expect_error(is.vector(get_exposure_category_mean(fit, "wrong", "baseline")))
  expect_error(is.vector(get_exposure_category_mean(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_mean(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_css
  expect_true(is.matrix(get_exposure_category_css(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_css(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_css(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_css(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_uss
  expect_true(is.matrix(get_exposure_category_uss(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_uss(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_uss(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_uss(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_cov
  expect_true(is.matrix(get_exposure_category_cov(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_cov(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_cov(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_cov(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get multiple exposure statistics
  expect_true(is_tibble(get_exposure_category_statistic(fit, c("/b/", "/p/"), c("baseline", "plus20.20"), c("n", "mean", "cov"))))
})

test_that("get expected category statistic", {
  expect_true(is.vector(get_expected_mu(fit, "/b/", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "wrong", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "/b/", "wrong")))
  expect_true(is.matrix(get_expected_sigma(fit, "/b/", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "wrong", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_expected_sigma(fit, c("/b/", "/p/"), c("prior", "plus20.20"))))
  expect_true(is_tibble(get_expected_category_statistic(fit, c("/b/", "/p/"), c("prior", "plus20.20"), c("mu", "Sigma"))))
})

test_that("get exposure category statistic", {
  # Get error when *exposure* statistics is requested for *prior*
  expect_error(get_exposure_category_statistic(fit, groups = "prior"))
  expect_error(get_exposure_category_mean(fit, groups = "prior"))
  expect_error(get_exposure_category_cov(fit, groups = "prior"))
  # When prior is not requested
  # get_exposure_mean
  expect_true(is.vector(get_exposure_category_mean(fit, "/b/", "baseline")))
  expect_error(is.vector(get_exposure_category_mean(fit, "wrong", "baseline")))
  expect_error(is.vector(get_exposure_category_mean(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_mean(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_css
  expect_true(is.matrix(get_exposure_category_css(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_css(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_css(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_css(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_uss
  expect_true(is.matrix(get_exposure_category_uss(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_uss(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_uss(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_uss(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_cov
  expect_true(is.matrix(get_exposure_category_cov(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_cov(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_cov(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_cov(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get multiple exposure statistics
  expect_true(is_tibble(get_exposure_category_statistic(fit, c("/b/", "/p/"), c("baseline", "plus20.20"), c("n", "mean", "cov"))))
})

test_that("get expected category statistic", {
  expect_true(is.vector(get_expected_mu(fit, "/b/", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "wrong", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "/b/", "wrong")))
  expect_true(is.matrix(get_expected_sigma(fit, "/b/", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "wrong", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_expected_sigma(fit, c("/b/", "/p/"), c("prior", "plus20.20"))))
  expect_true(is_tibble(get_expected_category_statistic(fit, c("/b/", "/p/"), c("prior", "plus20.20"), c("mu", "Sigma"))))
})


context("get information from stanfit NIW (three cues)")
fit <- get_example_stanfit(3, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize", file_refit = "never")
test_that("add ibbu draws - input check (3 cues)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "plus20.20.20")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "plus20.20.20"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (3 cues)", {
  expect_equal(nrow(get_draws(fit, groups = "prior", ndraws = 10, seed = 1, wide = F) %>% distinct(.draw)), 10)
  expect_equal(nrow(get_draws(fit, groups = "prior", wide = F, summarize = T)), 2)
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = F)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "cue", "cue2", "m", "S", "lapse_rate"))
})

# test_that("Add draws - check wide = T", {
#   expect_equal(nrow(get_draws(fit, groups = "prior", wide = T, summarize = T)), 1)
#   expect_equal(nrow(get_draws(fit, groups = "prior", wide = T, summarize = F)), 10)
# })




