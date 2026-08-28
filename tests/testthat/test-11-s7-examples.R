test_that("representation examples cover all families and cue counts", {
  families <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")
  for (family in families) {
    cue_counts <- if (family %in% c("UVG", "NIX")) 1 else 1:3
    for (n_cues in cue_counts) {
      representation <- if (family == "UVG") {
        example_uvg_category_representation(n_cues = n_cues)
      } else if (family == "NIX") {
        example_nix_category_representation(n_cues = n_cues)
      } else {
        example_category_representation(family, n_cues = n_cues)
      }
      expect_true(S7::S7_inherits(representation, MVBU_CategoryRepresentation))
    }
  }
  expect_error(example_uvg_category_representation(n_cues = 2))
  expect_error(example_nix_category_representation(n_cues = 3))
})

test_that("representation example wrapper dispatches by type", {
  for (family in c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")) {
    representation <- example_category_representation(family, n_cues = 1)
    expect_true(S7::S7_inherits(representation, MVBU_CategoryRepresentation))
  }
  expect_error(example_category_representation("unknown"))
})

test_that("template examples cover all families and cue counts", {
  for (family in c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")) {
    cue_counts <- if (family %in% c("UVG", "NIX")) 1 else 1:3
    for (n_cues in cue_counts) {
      template <- example_category_representation_template(family, n_cues = n_cues)
      expect_true(S7::S7_inherits(template, MVBU_CategoryRepresentationTemplate))
    }
  }
  expect_error(example_category_representation_template("UVG", n_cues = 2))
  expect_error(example_category_representation_template("NIX", n_cues = 3))
})

test_that("model examples cover all families and cue counts", {
  for (family in c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")) {
    cue_counts <- if (family %in% c("UVG", "NIX")) 1 else 1:3
    for (n_cues in cue_counts) {
      model <- example_model(family, n_cues = n_cues)
      expect_true(S7::S7_inherits(model, MVBU_CognitiveModel))
    }
  }
  expect_error(example_model("UVG", n_cues = 2))
  expect_error(example_model("NIX", n_cues = 3))
  expect_error(example_model("unknown"))
})
