# Extracted from test-01-s7-core-classes.R:619

# test -------------------------------------------------------------------------
niw_dist <- new_niw_model_distribution(group_label = "fit-niw")
mvg_dist <- new_mvg_model_distribution(group_label = "fit-mvg")
muvg_dist <- new_muvg_model_distribution(group_label = "fit-muvg")
mnix_dist <- new_mnix_model_distribution(group_label = "fit-mnix")
uvg_dist <- new_uvg_model_distribution(group_label = "fit-uvg")
nix_dist <- new_nix_model_distribution(group_label = "fit-nix")
ex_dist <- new_exemplar_model_distribution(group_label = "fit-ex")
expect_true(S7::S7_inherits(niw_dist, NIW_IdealAdaptorDistribution))
expect_true(S7::S7_inherits(mvg_dist, MVG_IdealObserverDistribution))
expect_true(S7::S7_inherits(muvg_dist, MUVG_IdealObserverDistribution))
expect_true(S7::S7_inherits(mnix_dist, MNIX_IdealAdaptorDistribution))
expect_true(S7::S7_inherits(uvg_dist, UVG_IdealObserverDistribution))
expect_true(S7::S7_inherits(nix_dist, NIX_IdealAdaptorDistribution))
expect_true(S7::S7_inherits(ex_dist, Exemplar_ModelDistribution))
expect_equal(get_model_family(niw_dist), "NIW")
expect_equal(get_model_family(mvg_dist), "MVG")
expect_equal(get_model_family(muvg_dist), "MUVG")
expect_equal(get_model_family(mnix_dist), "MNIX")
expect_equal(get_model_family(uvg_dist), "UVG")
expect_equal(get_model_family(nix_dist), "NIX")
expect_equal(get_model_family(ex_dist), "EXEMPLAR")
baseline_families <- get_registered_model_families()
expect_true(all(c("EXEMPLAR", "MNIX", "MUVG", "MVG", "NIW", "NIX", "UVG") %in% baseline_families))
register_model_family(
    family = "custommix",
    category_representation_class = "CustomMix_CategoryRepresentation",
    cognitive_model_class = "CustomMix_Model",
    model_distribution_class = "CustomMix_ModelDistribution"
  )
