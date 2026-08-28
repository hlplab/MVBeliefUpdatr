# Extracted from test-05-fit-ideal-adaptor-minimal-examples.R:108

# test -------------------------------------------------------------------------
skip_if_not_installed("rstan")
skip_if_not_installed("tidybayes")
model_name <- "minimal-nix_ideal_adaptor-nix"
fit <- load_ideal_adaptor_fit_model(model_name)
recovered_stanfit <- tidybayes::recover_types(get_stanfit(fit))
recovered_fit <- set_stanfit(fit, recovered_stanfit)
expect_true(S7::S7_inherits(recovered_fit, IdealAdaptorStanfit))
expect_s4_class(get_stanfit(recovered_fit), "stanfit")
expect_true(!is.null(attr(get_stanfit(recovered_fit), "tidybayes_constructors")))
expect_true(is.function(get_constructor(recovered_fit, "group")))
