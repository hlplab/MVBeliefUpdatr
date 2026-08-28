# Extracted from test-01-s7-core-classes.R:412

# test -------------------------------------------------------------------------
mvg_rep <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    Sigma = matrix(c(1, 0, 0, 2), nrow = 2)
  )
expect_true(S7::S7_inherits(mvg_rep, MVG_CategoryRepresentation))
mvg_params <- get_parameters(mvg_rep)
expect_equal(mvg_params$mu, c(0, 1))
niw_rep <- new_niw_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    m = c(0, 1),
    kappa = 1,
    nu = 3,
    S = matrix(c(1, 0, 0, 1), nrow = 2)
  )
expect_true(S7::S7_inherits(niw_rep, NIW_CategoryRepresentation))
niw_params <- get_parameters(niw_rep)
expect_equal(niw_params$kappa, 1)
uvg_rep <- new_uvg_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    mu = 0,
    sigma2 = 1
  )
expect_true(S7::S7_inherits(uvg_rep, UVG_CategoryRepresentation))
nix_rep <- new_nix_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    m = 0,
    kappa = 1,
    nu = 2,
    sigma2 = 1
  )
expect_true(S7::S7_inherits(nix_rep, NIX_CategoryRepresentation))
mnix_rep_lik <- new_mnix_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    component_m = c(0, 1),
    component_kappa = c(1, 2),
    component_nu = c(3, 4),
    component_sigma2 = c(1, 2)
  )
