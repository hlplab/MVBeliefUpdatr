context("Class information")

fit = readRDS("../test models/IBBU_mv_fit_2 exposure groups_2 categories_2 cues_Drouin et al 2016.rds")
g = add_ibbu_stanfit_draws(fit, which = "prior", summarize = T, nest = T)

test_that("Add ibbu draws - output check", {
  expect_true(is.NIW_ideal_adaptor_stanfit(fit))
  expect_false(is.NIW_ideal_adaptor_stanfit(g))
  expect_true(is.NIW_ideal_adaptor_MCMC(g))
  expect_false(is.NIW_ideal_adaptor_MCMC(fit))
})
