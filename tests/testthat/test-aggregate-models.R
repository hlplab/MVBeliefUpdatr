test_that("aggregate_models works for MVG_IdealObserver models", {
  m1 <- new_mvg_ideal_observer(
    category_template = new_category_representation_template(
      representations = list(
        c1 = new_mvg_category_representation("c1", c("cue1", "cue2"), mu = c(0, 0), Sigma = diag(2)),
        c2 = new_mvg_category_representation("c2", c("cue1", "cue2"), mu = c(2, 2), Sigma = diag(2))
      )
    ),
    category_prior = c(c1 = 0.6, c2 = 0.4),
    lapse_rate = 0.05
  )

  m2 <- new_mvg_ideal_observer(
    category_template = new_category_representation_template(
      representations = list(
        c1 = new_mvg_category_representation("c1", c("cue1", "cue2"), mu = c(2, 2), Sigma = 2 * diag(2)),
        c2 = new_mvg_category_representation("c2", c("cue1", "cue2"), mu = c(4, 4), Sigma = 2 * diag(2))
      )
    ),
    category_prior = c(c1 = 0.4, c2 = 0.6),
    lapse_rate = 0.15
  )

  agg <- aggregate_models(list(m1, m2))
  expect_true(S7::S7_inherits(agg, MVG_IdealObserver))
  expect_equal(get_lapse_rate(agg), 0.10)
  expect_equal(get_category_prior(agg), c(c1 = 0.5, c2 = 0.5))

  c1_rep <- agg@category_template@representations$c1
  expect_equal(c1_rep@mu, c(1, 1))
  expect_equal(c1_rep@Sigma, 1.5 * diag(2))

  # Test custom weights
  agg_w <- aggregate_models(list(m1, m2), weights = c(3, 1))
  expect_equal(get_lapse_rate(agg_w), 0.075)
  expect_equal(get_category_prior(agg_w), c(c1 = 0.55, c2 = 0.45))
  expect_equal(agg_w@category_template@representations$c1@mu, c(0.5, 0.5))
})

test_that("aggregate_models works for NIW_IdealAdaptor models", {
  m1 <- new_niw_ideal_adaptor(
    category_template = new_category_representation_template(
      representations = list(
        c1 = new_niw_category_representation("c1", c("cue1", "cue2"), m = c(0, 0), S = diag(2), kappa = 10, nu = 12),
        c2 = new_niw_category_representation("c2", c("cue1", "cue2"), m = c(2, 2), S = diag(2), kappa = 10, nu = 12)
      )
    ),
    category_prior = c(c1 = 0.5, c2 = 0.5),
    lapse_rate = 0.02
  )

  m2 <- new_niw_ideal_adaptor(
    category_template = new_category_representation_template(
      representations = list(
        c1 = new_niw_category_representation("c1", c("cue1", "cue2"), m = c(2, 4), S = 3 * diag(2), kappa = 20, nu = 16),
        c2 = new_niw_category_representation("c2", c("cue1", "cue2"), m = c(4, 6), S = 3 * diag(2), kappa = 20, nu = 16)
      )
    ),
    category_prior = c(c1 = 0.5, c2 = 0.5),
    lapse_rate = 0.04
  )

  agg <- aggregate_models(list(m1, m2))
  expect_true(S7::S7_inherits(agg, NIW_IdealAdaptor))
  expect_equal(get_lapse_rate(agg), 0.03)

  c1_rep <- agg@category_template@representations$c1

  expect_equal(c1_rep@m, c(1, 2))
  expect_equal(c1_rep@S, 2 * diag(2))
  expect_equal(c1_rep@kappa, 15)
  expect_equal(c1_rep@nu, 14)
})

test_that("aggregate_models validates inputs", {
  expect_error(aggregate_models(list()), "non-empty list")
  expect_error(aggregate_models("not_a_list"), "non-empty list")
  expect_error(aggregate_models(list(1, 2)), "CognitiveModel")

  m1 <- new_mvg_ideal_observer(
    category_template = new_category_representation_template(
      representations = list(
        c1 = new_mvg_category_representation("c1", c("cue1", "cue2"), mu = c(0, 0), Sigma = diag(2))
      )
    )
  )
  m_niw <- new_niw_ideal_adaptor(
    category_template = new_category_representation_template(
      representations = list(
        c1 = new_niw_category_representation("c1", c("cue1", "cue2"), m = c(0, 0), S = diag(2), kappa = 10, nu = 12)
      )
    )
  )

  expect_error(aggregate_models(list(m1, m_niw)), "same S7 class")
  expect_error(aggregate_models(list(m1, m1), weights = c(1, -1)), "positive weights")
  expect_error(aggregate_models(list(m1, m1), weights = c(1, 2, 3)), "same length")
})
