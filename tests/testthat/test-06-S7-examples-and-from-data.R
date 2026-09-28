test_that("representation examples cover all families and cue counts", {
  families <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")
  for (family in families) {
    cue_counts <- if (family %in% c("UVG", "NIX")) {
      1
    } else if (family %in% c("MUVG", "MNIX")) {
      2:3
    } else {
      1:3
    }
    for (n_cues in cue_counts) {
      representation <- if (family == "UVG") {
        example_uvg_category_representation(n_cues = n_cues)
      } else if (family == "NIX") {
        example_nix_category_representation(n_cues = n_cues)
      } else if (family == "MUVG") {
        example_muvg_category_representation(n_cues = n_cues)
      } else if (family == "MNIX") {
        example_mnix_category_representation(n_cues = n_cues)
      } else {
        example_category_representation(family, n_cues = n_cues)
      }
      expect_true(S7::S7_inherits(representation, MVBU_CategoryRepresentation))
    }
  }
  expect_error(example_uvg_category_representation(n_cues = 2))
  expect_error(example_nix_category_representation(n_cues = 3))
  expect_error(example_muvg_category_representation(n_cues = 1))
  expect_error(example_mnix_category_representation(n_cues = 1))
})

test_that("representation example wrapper dispatches by type", {
  for (family in c("UVG", "NIX", "MVG", "NIW", "EXEMPLAR")) {
    representation <- example_category_representation(family, n_cues = 1)
    expect_true(S7::S7_inherits(representation, MVBU_CategoryRepresentation))
  }
  for (family in c("MUVG", "MNIX")) {
    representation <- example_category_representation(family, n_cues = 2)
    expect_true(S7::S7_inherits(representation, MVBU_CategoryRepresentation))
  }
  expect_error(example_category_representation("unknown"))
})

test_that("template examples cover all families and cue counts", {
  for (family in c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")) {
    cue_counts <- if (family %in% c("UVG", "NIX")) {
      1
    } else if (family %in% c("MUVG", "MNIX")) {
      2:3
    } else {
      1:3
    }
    for (n_cues in cue_counts) {
      template <- example_category_representation_template(
        family,
        n_cues = n_cues
      )
      expect_true(S7::S7_inherits(template, MVBU_CategoryRepresentationTemplate))
    }
  }
  expect_error(example_category_representation_template("UVG", n_cues = 2))
  expect_error(example_category_representation_template("NIX", n_cues = 3))
  expect_error(example_category_representation_template("MUVG", n_cues = 1))
  expect_error(example_category_representation_template("MNIX", n_cues = 1))
})

test_that("model examples cover all families and cue counts", {
  for (family in c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")) {
    cue_counts <- if (family %in% c("UVG", "NIX")) {
      1
    } else if (family %in% c("MUVG", "MNIX")) {
      2:3
    } else {
      1:3
    }
    for (n_cues in cue_counts) {
      model <- example_model(family, n_cues = n_cues)
      expect_true(S7::S7_inherits(model, MVBU_CognitiveModel))
    }
  }
  expect_error(example_model("UVG", n_cues = 2))
  expect_error(example_model("NIX", n_cues = 3))
  expect_error(example_model("MUVG", n_cues = 1))
  expect_error(example_model("MNIX", n_cues = 1))
  expect_error(example_model("unknown"))
})

data("mixer6", package = "MVBeliefUpdatr")
stops_data <- mixer6[
  mixer6$speaker == "111138" &
    mixer6$stop %in% c("/b/", "/p/"),
  c("stop", "vot", "f0")
]
stops_data$stop <- factor(as.character(stops_data$stop))
stops_data <- stops_data[!is.na(stops_data$f0), ]

representation_types <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")
template_types <- representation_types
model_types <- representation_types

representation_cues <- function(type) {
  if (type %in% c("UVG", "NIX")) "vot" else c("vot", "f0")
}

representation_args <- function(type) {
  args <- list(
    data = stops_data,
    category = "stop",
    cues = representation_cues(type)
  )
  if (type %in% c("NIX", "MNIX", "NIW")) {
    args[c("kappa", "nu")] <- list(kappa = 10, nu = 30)
  }
  args
}

test_that("category representation from-data constructors cover all families", {
  for (type in representation_types) {
    args <- representation_args(type)
    args$data <- args$data[args$data$stop == "/b/", ]
    args$category <- "stop"
    result <- do.call(
      new_category_representation_from_data,
      c(args, type = type)
    )
    expect_true(S7::S7_inherits(result, MVBU_CategoryRepresentation))
  }
  expect_error(
    new_category_representation_from_data(stops_data, type = "unknown", cues = "vot")
  )
})

test_that("category representation from-data constructors require one category", {
  expect_true(S7::S7_inherits(
    new_mvg_category_representation_from_data(
      stops_data[stops_data$stop == "/b/", ],
      category = "stop", cues = "vot"
    ),
    MVG_CategoryRepresentation
  ))
  expect_error(
    new_mvg_category_representation_from_data(
      stops_data, category = "stop", cues = "vot"
    )
  )
})

test_that("category representation template from-data constructors cover all families", {
  for (type in template_types) {
    result <- do.call(
      new_category_representation_template_from_data,
      c(
        list(
          data = stops_data, category = "stop",
          cues = representation_cues(type)
        ),
        type = type,
        if (type %in% c("NIX", "MNIX", "NIW")) list(kappa = 10, nu = 30) else list()
      )
    )
    expect_true(S7::S7_inherits(result, MVBU_CategoryRepresentationTemplate))
  }
  expect_error(
    new_category_representation_template_from_data(
      stops_data, type = "unknown", cues = "vot"
    )
  )
})

test_that("model from-data constructors cover all families", {
  for (type in model_types) {
    result <- do.call(
      new_model_from_data,
      c(
        list(
          data = stops_data, category = "stop",
          cues = representation_cues(type)
        ),
        type = type,
        if (type %in% c("NIX", "MNIX", "NIW")) list(kappa = 10, nu = 30) else list()
      )
    )
    expect_true(S7::S7_inherits(result, MVBU_CognitiveModel))
  }
  expect_error(
    new_model_from_data(stops_data, type = "unknown", cues = "vot")
  )
})

test_that("one-cue family constraints are enforced", {
  expect_error(
    new_category_representation_from_data(
      stops_data, type = "UVG", cues = c("vot", "f0")
    )
  )
  expect_error(
    new_category_representation_from_data(
      stops_data, type = "NIX", cues = c("vot", "f0")
    )
  )
  expect_true(S7::S7_inherits(
    new_category_representation_from_data(
      stops_data[stops_data$stop == "/b/", ],
      type = "MUVG", category = "stop", cues = c("vot", "f0")
    ),
    MUVG_CategoryRepresentation
  ))
  expect_true(S7::S7_inherits(
    new_category_representation_from_data(
      stops_data[stops_data$stop == "/b/", ],
      type = "MNIX", category = "stop", cues = c("vot", "f0"),
      kappa = 10, nu = 30
    ),
    MNIX_CategoryRepresentation
  ))
})
