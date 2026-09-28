# Cross-family method parity test matrix
# Verifies that all concrete cognitive model families (UVG, NIX, MUVG, MNIX,
# MVG, NIW, Exemplar) adhere to the unified S7 method interface contract.

test_that("all cognitive model families satisfy the core generic method parity matrix", {
  families <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")
  
  for (fam in families) {
    nc <- if (fam %in% c("UVG", "NIX")) 1L else 2L
    m <- example_model(fam, n_cues = nc)
    
    # 1. Type and structure introspection
    expect_true(S7::S7_inherits(m, MVBU_CognitiveModel), info = fam)
    expect_true(is.character(get_model_type(m)), info = fam)
    expect_true(is.character(get_representation_type(m)), info = fam)
    
    cues <- get_cue_labels(m)
    expect_equal(length(cues), nc, info = fam)
    cats <- get_category_labels(m)
    expect_true(length(cats) >= 2L, info = fam)
    
    # 2. Parameters extraction
    pars <- get_parameters(m)
    expect_true(is.list(pars) || is.data.frame(pars), info = fam)
    
    # Generate test observation points
    set.seed(42)
    dat <- as.data.frame(matrix(rnorm(nc * 6L), nrow = 6L))
    colnames(dat) <- cues
    
    # 3. Likelihood computation
    lik <- likelihood(m, dat)
    expect_true(is.data.frame(lik) || is.matrix(lik), info = fam)
    expect_equal(nrow(lik), 6L, info = fam)
    expect_true(all(lik >= 0, na.rm = TRUE), info = fam)
    
    # 4. Posterior computation
    post <- posterior(m, dat)
    expect_true(is.data.frame(post) || is.matrix(post), info = fam)
    expect_equal(nrow(post), 6L, info = fam)
    # Probabilities should sum to 1 across categories for each observation
    row_sums <- rowSums(post[, cats, drop = FALSE])
    expect_equal(row_sums, rep(1, 6L), tolerance = 1e-6, info = fam)
    
    # 5. Categorization under decision rules
    # 5a. Deterministic MAP criterion
    cat_crit <- categorize(m, dat, decision_rule = "criterion")
    expect_true(is.character(cat_crit), info = fam)
    expect_equal(length(cat_crit), 6L, info = fam)
    expect_true(all(cat_crit %in% cats), info = fam)
    
    # 5b. Probability matching (proportional) with detailed return format
    cat_prop <- categorize(m, dat, decision_rule = "proportional", simplify = FALSE)
    expect_s3_class(cat_prop, "data.frame")
    expect_equal(nrow(cat_prop), 6L, info = fam)
    expect_true(all(c("response_category", "response_probability") %in% names(cat_prop)), info = fam)
    expect_true(all(cat_prop$response_probability >= 0 & cat_prop$response_probability <= 1), info = fam)
    
    # 5c. Posterior sampling
    set.seed(123)
    cat_samp <- categorize(m, dat, decision_rule = "sampling")
    expect_true(is.character(cat_samp), info = fam)
    expect_equal(length(cat_samp), 6L, info = fam)
    expect_true(all(cat_samp %in% cats), info = fam)
    
    # 6. Sample observations (generative synthesis)
    samp <- sample_observations(m, n = 8L)
    expect_s3_class(samp, "data.frame")
    expect_equal(nrow(samp), 8L, info = fam)
    expect_true("category" %in% names(samp), info = fam)
    expect_true(all(cues %in% names(samp)), info = fam)
    
    # 7. Model evaluation
    resp_vec <- sample(cats, size = 6L, replace = TRUE)
    ev_loglik <- evaluate_model(m, dat, response_category = resp_vec, method = "log_lik")
    expect_true(is.numeric(ev_loglik), info = fam)
    expect_true(is.finite(ev_loglik), info = fam)
    
    ev_acc <- evaluate_model(m, dat, response_category = resp_vec, method = "accuracy")
    expect_true(is.numeric(ev_acc), info = fam)
    expect_true(ev_acc >= 0 && ev_acc <= 1, info = fam)
    
    # 8. Printing and summary
    out_print <- capture.output(print(m))
    expect_true(length(out_print) > 0L, info = fam)
    out_sum <- capture.output(summary(m))
    expect_true(length(out_sum) > 0L, info = fam)
    
    # 9. Plotting methods
    p_cat <- plot_categories(m, interactive = FALSE)
    expect_s3_class(p_cat, "ggplot")
    
    p_dec <- plot_categorization_functions(m, interactive = FALSE)
    expect_s3_class(p_dec, "ggplot")
  }
})

test_that("cross-family template parity holds across representation families", {
  families <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")
  
  for (fam in families) {
    nc <- if (fam %in% c("UVG", "NIX")) 1L else 2L
    tpl <- example_category_representation_template(fam, n_cues = nc)
    
    expect_true(S7::S7_inherits(tpl, MVBU_CategoryRepresentationTemplate), info = fam)
    expect_true(is.character(get_representation_type(tpl)), info = fam)
    expect_equal(length(get_cue_labels(tpl)), nc, info = fam)
    expect_true(length(get_category_labels(tpl)) >= 2L, info = fam)
    
    # Check that individual representations inherit from MVBU_CategoryRepresentation
    for (cat_name in get_category_labels(tpl)) {
      rep_obj <- tpl@representations[[cat_name]]
      expect_true(S7::S7_inherits(rep_obj, MVBU_CategoryRepresentation), info = paste(fam, cat_name))
      expect_equal(get_cue_labels(rep_obj), get_cue_labels(tpl), info = paste(fam, cat_name))
    }
  }
})
