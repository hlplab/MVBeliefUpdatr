representation <- example_mvg_category_representation(n_cues = 2)
uvg <- example_uvg_category_representation(n_cues = 1)
nix <- example_nix_category_representation(n_cues = 1)
muvg <- example_muvg_category_representation(n_cues = 2)
mnix <- example_mnix_category_representation(n_cues = 2, nu = 30)
mvg <- example_mvg_category_representation(n_cues = 2)
niw <- example_niw_category_representation(n_cues = 2, nu = 30)
exemplar <- example_exemplar_category_representation(n_cues = 2)

test_that("representation coercions cover supported family pairs", {
  expect_true(S7::S7_inherits(as_nix_category_representation(uvg, kappa = 10, nu = 30), NIX_CategoryRepresentation))
  expect_true(S7::S7_inherits(as_uvg_category_representation(nix), UVG_CategoryRepresentation))
  expect_true(S7::S7_inherits(as_mnix_category_representation(muvg, kappa = 10, nu = 30), MNIX_CategoryRepresentation))
  expect_true(S7::S7_inherits(as_muvg_category_representation(mnix), MUVG_CategoryRepresentation))
  expect_true(S7::S7_inherits(as_niw_category_representation(mvg, kappa = 10, nu = 30), NIW_CategoryRepresentation))
  expect_true(S7::S7_inherits(as_mvg_category_representation(niw), MVG_CategoryRepresentation))
  expect_true(S7::S7_inherits(as_exemplar_category_representation(representation, n = 20), Exemplar_CategoryRepresentation))
  expect_true(S7::S7_inherits(as_mvg_category_representation(exemplar), MVG_CategoryRepresentation))
})

test_that("template coercions preserve template structure", {
  template <- example_category_representation_template("MVG", n_cues = 2)
  expect_true(S7::S7_inherits(as_niw_category_representation_template(template, kappa = 10, nu = 30), MVBU_CategoryRepresentationTemplate))
  expect_true(S7::S7_inherits(as_exemplar_category_representation_template(template, n = 10), MVBU_CategoryRepresentationTemplate))
})

test_that("model coercions preserve model structure", {
  model <- example_mvg_ideal_observer(n_cues = 2)
  expect_true(S7::S7_inherits(as_niw_ideal_adaptor(model, kappa = 10, nu = 30), MVBU_CognitiveModel))
  expect_true(S7::S7_inherits(as_exemplar_model(model, n = 10), MVBU_CognitiveModel))
})

test_that("same-family coercions return the original object", {
  expect_identical(as_mvg_category_representation(mvg), mvg)
  expect_identical(as_exemplar_category_representation(exemplar), exemplar)
})

test_that("unsupported coercions fail clearly", {
  expect_error(as_nix_category_representation(mvg, kappa = 10, nu = 30), "not supported")
  expect_error(as_mvg_category_representation(uvg), "not supported")
})
