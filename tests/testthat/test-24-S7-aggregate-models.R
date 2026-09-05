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

test_that("aggregate S7 generic works on category representations", {
  r1 <- new_mvg_category_representation("c1", c("cue1", "cue2"), mu = c(0, 0), Sigma = diag(2))
  r2 <- new_mvg_category_representation("c1", c("cue1", "cue2"), mu = c(4, 8), Sigma = 3 * diag(2))

  # Variadic call
  agg_var <- aggregate(r1, r2, weights = c(1, 3))
  expect_true(S7::S7_inherits(agg_var, MVG_CategoryRepresentation))
  expect_equal(agg_var@mu, c(3, 6))
  expect_equal(agg_var@Sigma, 2.5 * diag(2))

  # List call
  agg_list <- aggregate(list(r1, r2), weights = c(1, 1))
  expect_equal(agg_list@mu, c(2, 4))
  expect_equal(agg_list@Sigma, 2 * diag(2))

  # UVG
  uvg1 <- new_uvg_category_representation("c1", "cue1", mu = 0, sigma2 = 1)
  uvg2 <- new_uvg_category_representation("c1", "cue1", mu = 10, sigma2 = 5)
  agg_uvg <- aggregate(uvg1, uvg2)
  expect_equal(agg_uvg@mu, 5)
  expect_equal(agg_uvg@sigma2, 3)

  # NIX
  nix1 <- new_nix_category_representation("c1", "cue1", m = 0, kappa = 4, nu = 8, sigma2 = 2)
  nix2 <- new_nix_category_representation("c1", "cue1", m = 10, kappa = 8, nu = 12, sigma2 = 6)
  agg_nix <- aggregate(nix1, nix2)
  expect_equal(agg_nix@m, 5)
  expect_equal(agg_nix@sigma2, 4)
  expect_equal(agg_nix@kappa, 6)
  expect_equal(agg_nix@nu, 10)
})

test_that("aggregate S7 generic works on category representation templates", {
  t1 <- new_category_representation_template(
    representations = list(
      c1 = new_mvg_category_representation("c1", c("cue1", "cue2"), mu = c(0, 0), Sigma = diag(2)),
      c2 = new_mvg_category_representation("c2", c("cue1", "cue2"), mu = c(2, 2), Sigma = diag(2))
    )
  )
  t2 <- new_category_representation_template(
    representations = list(
      c1 = new_mvg_category_representation("c1", c("cue1", "cue2"), mu = c(4, 4), Sigma = 3 * diag(2)),
      c2 = new_mvg_category_representation("c2", c("cue1", "cue2"), mu = c(6, 6), Sigma = 3 * diag(2))
    )
  )

  agg_t <- aggregate(t1, t2)
  expect_true(S7::S7_inherits(agg_t, MVBU_CategoryRepresentationTemplate))
  expect_equal(agg_t@representations$c1@mu, c(2, 2))
  expect_equal(agg_t@representations$c2@mu, c(4, 4))
  expect_equal(agg_t@representations$c1@Sigma, 2 * diag(2))
})

test_that("aggregate works for Exemplar models and representations", {
  ex1 <- matrix(c(1, 2, 3, 4), ncol = 2)
  ex2 <- matrix(c(5, 6, 7, 8), ncol = 2)

  r1 <- new_exemplar_category_representation("c1", c("cue1", "cue2"), exemplars = ex1, exemplar_weights = c(0.5, 0.5), c = 1)
  r2 <- new_exemplar_category_representation("c1", c("cue1", "cue2"), exemplars = ex2, exemplar_weights = c(0.4, 0.6), c = 3)

  agg_r <- aggregate(r1, r2, weights = c(1, 1))
  expect_true(S7::S7_inherits(agg_r, Exemplar_CategoryRepresentation))
  expect_equal(nrow(agg_r@exemplars), 4)
  expect_equal(agg_r@exemplars, rbind(ex1, ex2))
  expect_equal(sum(agg_r@exemplar_weights), 1)
  expect_equal(agg_r@c, 2)

  # Full exemplar model aggregation
  m1 <- new_exemplar_model(
    category_template = new_category_representation_template(representations = list(c1 = r1)),
    lapse_rate = 0.04
  )
  m2 <- new_exemplar_model(
    category_template = new_category_representation_template(representations = list(c1 = r2)),
    lapse_rate = 0.08
  )

  agg_m <- aggregate(m1, m2)
  expect_true(S7::S7_inherits(agg_m, Exemplar_Model))
  expect_equal(get_lapse_rate(agg_m), 0.06)
  expect_equal(nrow(agg_m@category_template@representations$c1@exemplars), 4)
})

test_that("aggregate falls back to stats::aggregate for non-S7 objects", {
  df <- data.frame(group = c("a", "a", "b", "b"), val = c(1, 3, 5, 7))
  res <- aggregate(val ~ group, data = df, FUN = mean)
  expect_equal(res$val, c(2, 6))
})
