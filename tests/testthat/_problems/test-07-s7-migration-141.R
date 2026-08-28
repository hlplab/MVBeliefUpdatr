# Extracted from test-07-s7-migration.R:141

# test -------------------------------------------------------------------------
legacy_muvg <- data.frame(category = factor(c("A", "B")))
legacy_muvg$component_mu <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
legacy_muvg$component_sigma2 <- list(c(F1 = 1, F2 = 2), c(F1 = 1, F2 = 1))
muvg_template <- as_s7_category_representation_template(legacy_muvg, family = "MUVG")
expect_true(S7::S7_inherits(muvg_template, MVBU_CategoryRepresentationTemplate))
expect_true(all(vapply(muvg_template@representations, function(r) S7::S7_inherits(r, MUVG_CategoryRepresentation), logical(1))))
muvg_model <- as_s7_muvg_ideal_observer(legacy_muvg)
expect_true(S7::S7_inherits(muvg_model, MUVG_IdealObserver))
legacy_mnix <- data.frame(category = factor(c("A", "B")))
legacy_mnix$component_m <- list(c(0, 1), c(2, 3))
legacy_mnix$component_kappa <- list(c(1, 2), c(1, 2))
legacy_mnix$component_nu <- list(c(3, 4), c(3, 4))
legacy_mnix$component_sigma2 <- list(c(1, 2), c(1, 1))
mnix_template <- as_s7_category_representation_template(legacy_mnix, family = "MNIX")
