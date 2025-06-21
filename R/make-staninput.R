check_exposure_test_data <- function(data, cues, category, response, group, which.data = "the", verbose = F) {
  assert_that(is_tibble(data) | is.data.frame(data))
  assert_cols_in_data(data, cues, which.data, scalar = F)
  assert_cols_in_data(data, group, which.data, scalar = T)

  data %<>%
    drop_na(.env$cues, .env$group) %>%
    mutate(across(.env$group, as.factor))

  if(!is.null(category)) {
    assert_cols_in_data(data, category, which.data, scalar = T)
    data %<>%
      mutate(across(.env$category, as.factor)) %>%
      drop_na(.env$category)
  }

  if (!is.null(response)) {
    assert_cols_in_data(data, response, which.data, scalar = T)
    data %<>%
      mutate(across(.env$response, as.factor)) %>%
      drop_na(.env$response)
  }
  assert_that(nrow(data) > 0,
              msg = paste("There must be at least one observation in", which.data, "data."))

  if (verbose){
    print("In check_exposure_test_data():")
    print (data)
  }

  return(data)
}


get_test_counts <- function(test, cues, response, group, verbose = F) {
  test_counts <-
    test %>%
    as_tibble(.name_repair = "minimal") %>%
    group_by(
      !!! syms(cues),
      !! sym(response),
      !! sym(group)) %>%
    tally() %>%
    pivot_wider(
      names_from = !! response,
      values_from = n,
      values_fill = 0) %>%
    ungroup()

  if (verbose) {
    print("In get_test_counts():")
    print(test_counts)
  }

  return(test_counts)
}



#' Get category statistics from data as list of arrays
#'
#' Convenience function that obtains category statistics (aggregate functions) over `cues` in a data frame, and prepares them to be
#' used as input to a MVBeliefUpdatr Stan program. Each function is calculated for each unique combination of values in the columns
#' `category` and (optionally) `group`. Missing values will be automatically filled since Stan can't handle `NA` inputs.
#'
#' @param data `tibble` or `data.frame` with the data. Each row should be an observation of a category,
#'   and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param cues Names of columns with cue values.
#' @param category Name of column that contains the category label for the exposure data. This column must be a factor.
#'   Can be `NULL` for unsupervised updating (not yet implemented). (default: "category")
#' @param group Name of column that contains information about which observations form a group. This column must be
#'   a factor.
#' @param fill A value to fill in for missing combinations of `group` and `category` for which there is no data. (default: `0`)
#' @param simplify A list of logicals of the same length as `...` indicating whether the array resulting from each function
#'   should be simplified? See \code{\link{to_array}} for more detail. (default: `TRUE` for all functions)
#' @param ... A named list of functions (statistics) of `cues` to calculate for each unique combination of `group` and
#'   `category`.
#'
#' @return A named list of length `...` (names are the names of the functions that have been computed). The list will be
#'   sorted first by `group` and then by `category`, in ascending order of the levels of those variables. The elements of
#'   the list will be arrays.
#'   Missing values---resulting from combinations of `group` and `category` for which there is no data---will be filled
#'   with `fill` values coerced into the same structure as all other outputs for that function (e.g., if the function f results
#'   in a 2-element vector for all combinations of `group` and `category` for which there is data, then unobserved data
#'   results in a 2-element vector of `fill` values).
#'
#' @keywords TBD
#' @rdname get_category_statistics_as_list_of_arrays
#' @export
get_category_statistics_as_list_of_arrays <- function(
    data,
    group = "group",
    category = "category",
    cues,
    fill = as.list(rep(0, length(list(...)))),
    simplify = as.list(rep(T, length(list(...)))),
    verbose = F,
    ...
) {
  # The order of category and group is important since the Stan programs all expect category to be the outer-most array dimension,
  # followed by group
  result <- get_aggregates_from_grouped_data_as_list_of_arrays(data, groups = c(group, category), cols = cues, fill = fill, simplify = simplify, verbose = verbose, ...)

  return(result)
}


#' Prepare long data from incremental exposure-test design for input to Stan
#'
#' Takes \code{data.frame} or \code{tibble} that contains the exposure and test data from an incremental
#' exposure-test design in long format, and prepares it for input to the \code{\link{ideal_adaptor_stanfit}}
#' Stan programs. This is done by pretending that each incremental test block (and its preceding exposure)
#' constitute a separate between-participant condition. Note that this does not capture the dependency
#' between test responses of participants in the same between-participant conditions, but such dependencies
#' are not modeled by current `MVBeliefUpdatr` Stan programs anyway (which do not include random effects by
#' participants).
#'
#' @param data Data frame or tibble to be sliced. Each row should be a single exposure or test observation.
#' @param group Character string indicating the name of the column that contains the information about
#'   the between-participant condition. (default: "Group")
#' @param phase Character string indicating the name of the column that contains the information about
#'   whether an observation is part of "exposure" or "test". This column must contain the values "exposure"
#'   and "test". Observation with other values will be ignored. (default: "Phase")
#' @param block Character string indicating the name of the column that contains the information about the
#'   incremental exposure and test blocks. Must be a factor with the levels indicating the order of the blocks.
#'   (default: "Block")
#' @param join_adjacent_test_blocks Logical indicating whether adjacent test blocks without intervening
#'   exposure blocks should be joined into a single test block. This will speed up \code{\link{fit_ideal_adapor}}
#'   since there will be fewer conditions to iterate over but also means that the default plotting functions
#'   won't be able to plot the results of the different test blocks separately. (default: `FALSE`)
#' @param verbose Should verbose output be provided? (default: `FALSE`)
#'
#' @return A data frame or tibble in long format with a new column "ExposureGroup" that contains a unique
#'   label for each unique combination of `group` and `block`.
#'
#' @export
reshape_incremental_design_into_unique_exposure_test_combinations <- function(
    data,
    group = "Group",
    phase = "Phase",
    block = "Block",
    join_adjacent_test_blocks = FALSE,
    verbose = FALSE
) {
  stopifnot(is_tibble(data) || is.data.frame(data))
  data %<>% ungroup()
  stopifnot(all(c(phase, group, block) %in% names(data)))
  stopifnot(all(c("exposure", "test") %in% unique(data[[phase]])))
  stopifnot(is.factor(data[[block]]))
  if (
    any(
      unique(data %>% filter(!! sym(phase) == "exposure") %>% pull(!! sym(block))) %in%
      unique(data %>% filter(!! sym(phase) == "test") %>% pull(!! sym(block)))))
    stop2("The levels of the block variable in the exposure phase must not overlap with those in the test phase. Please check your data.")


  if (verbose && length(setdiff(unique(data[[phase]]), c("exposure", "test"))) > 0) {
    message(
      paste("The following values in the", phase, "column are not recognized as exposure or test and thus removed:",
            paste(setdiff(unique(data[[phase]]), c("exposure", "test")), collapse = ", ")))
  }
  data %<>%
    filter(!! sym(phase) %in% c("exposure", "test")) %>%
    mutate(..block_order = as.numeric(!! sym(block)))

  block_levels <- levels(data[[block]])
  phase_levels <- data %>% distinct(!! sym(phase), !! sym(block), ..block_order) %>% arrange(..block_order) %>% pull(!! sym(phase))
  if (join_adjacent_test_blocks)
    for (b in 1:(length(block_levels) - 1)) {
      if (all(phase_levels[b:(b+1)] == "test")) {
        if (verbose) message("Joining adjacent test blocks ", block_levels[b], " and ", block_levels[b+1], " into a single test block.")

        block_levels[b] <- paste(block_levels[b], block_levels[b+1], sep = "_")
        block_levels <- block_levels[-(b+1)]
        data %<>%
          mutate(
            ..block_order = ifelse(..block_order == b + 1, b, ..block_order),
            !! sym(block) := factor(ifelse(..block_order %in% c(b, b + 1), block_levels[b], as.character(!! sym(block))), levels = block_levels))
      }
    }
  if (verbose) message("Inferred block order: ", paste(block_levels, collapse = ", "))

  testblock_order <-
    data %>%
    distinct(!! sym(phase), ..block_order) %>%
    filter(!! sym(phase) == "test") %>%
    .[["..block_order"]]
  if (verbose) message("Inferred test block order: ", paste(testblock_order, collapse = ", "))

  df.new <- tibble()
  for (g in unique(data[[group]]))
    for (b in testblock_order) {
      df.new %<>%
        rbind(
          data %>%
            # Include only exposure blocks from the current exposure condition and the current test block
            # from the current exposure condition (but not earlier test blocks)
            filter(
              !! sym(group) == g, ..block_order <= b,
              ..block_order == b | ! (!! sym(phase) == "test")) %>%
            mutate(ExposureGroup = if (b == 1) "no exposure" else paste0("Group ", g, "_up to block ", block_levels[b])))
    }

  df.new  %>%
    select(ExposureGroup, !! group, !!phase, !! block, everything())
}


#' Specify control parameters for make_staninput()
#'
#' This function is used to specify control parameters for the `make_staninput()` function, and to provide
#' reasonable defaults for any of the unspecified parameters.
#'
#' @param tau_scale A vector of scales for the Cauchy priors for each cue's standard deviations. Used in
#'   both the prior for m_0 and the prior for S_0. (default: vector of `5`s, assuming that the data are standardized).
#' @param L_omega_eta A vector of etas of the LKJ prior for the correlations of the covariance matrix of \code{mu_0}. Only used for
#'   models with multivariate categories (e.g., NIW_ideal_adaptor). (default: `1`,
#'   which corresponds to a uniform prior of correlation matrices)
#' @param split_loglik_per_observation Optionally, split the log likelihood per observation. This can be helpful of leave-one-out
#'   estimation in order to avoid high Pareto k, but it also makes the stored stanfit object much larger. (default: `0`)
#' @param transform_type An affine transformation that can be applied to the data. See `type` in \code{\link{get_affine_transform}}
#'    for details. (default: "standardize", which standardizes each cue separately)
#'
#' @return A list of control parameters that can be passed to \code{\link{make_staninput}}.
#'
#' @export
control_staninput <- function(
    tau_scale = 5,
    L_omega_eta = 1,
    split_loglik_per_observation = 0,
    # THIS IS KEPT HERE JUST FOR NOW UNTIL I HAVE DETERMINED WHICH TRANSFORM IS BEST SUITED FOR FITTING.
    # (also remove the documentation for this once it's no longer needed AND remove transform_information
    # from the returned information AND change the documentation for the returned object above.)
    transform_type = c("identity", "center", "standardize", "PCA whiten", "ZCA whiten")[3]
) {
  list(
    tau_scale = tau_scale,
    L_omega_eta = L_omega_eta,
    split_loglik_per_observation = split_loglik_per_observation,
    transform_type = transform_type)
}

#' Prepare data to fit ideal_adaptor_stanfit via rstan
#'
#' Take exposure and test data as input, and prepare the data for input into an MVBeliefUpdatr Stan program.
#'
#' It is important to use \code{group} to identify individuals that had a specific exposure (or no exposure at all)
#' and specific test trials. You should \emph{not} use \code{group} to identify exposure conditions. Setting \code{group} to an exposure condition
#' results in an exposure that concatenates the exposure observations from all subjects in that condition. Typically, this
#' is not what users intend, as it models exposure to the combination of exposure tokens across all subjects, rather than
#' exposure to one set of those exposure tokens. To achieve this intended outcome, use
#' \code{group.unique} to identify groups with identical exposure. This will correctly use only one unique instance of the
#' observations that any level of \code{group} receives during exposure.
#'
#' @param exposure `tibble` or `data.frame` with the exposure data. Each row should be an observation of a category,
#'   and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param test `tibble` or `data.frame` with the test data. Each row should be an observation, and contain information
#'   about the cue values of the test stimulus and the participant's response.
#' @param cues Names of columns with cue values. Must exist in both exposure and test data.
#' @param category Name of column in exposure data that contains the category label. Can be \code{NULL} for unsupervised updating
#'   (not yet implemented). (default: "category")
#' @param response Name of column in test data that contains participants' responses. (default: "response")
#' @param group Name of column that contains information about which observations form a group. Typically, this is
#'   a variable identifying subjects/participants. Must exist in both exposure and test data. (default: "group")
#' @param group.unique Name of column that uniquely identifies each group with identical exposure. This could be a
#'   variable indicating the different conditions in an experiment. Using group.unique is optional, but can be
#'   substantially more efficient if many groups share the same exposure. To ignore, set to \code{NULL}. (default: \code{NULL})
#' @param lapse_rate,mu_0,Sigma_0 Optionally, lapse rate, prior expected category means (\code{mu_0}) and/or prior expected
#'   category covariance matrices (\code{Sigma_0}) for all categories. Lapse rate should be a number between 0 and 1. For \code{mu_0}
#'   and \code{Sigma_0}, each should be a list, with each element being the expected mean/covariance matrix for a specific
#'   category prior to updating. Elements of \code{mu_0} and \code{Sigma_0} should be ordered in the same order as the levels of the
#'   category variable in \code{exposure} and \code{test}. These prior expected means and covariance matrices could be
#'   estimated, for example, from phonetically annotated speech recordings (see \code{\link{make_MVG_from_data}}
#'   for a convenient way to do so). Internally, m_0 is then set to \code{mu_0} (so that the expected value of the prior
#'   distribution of means is mu_0) and S_0 is set so that the expected value of the inverse-Wishart is \code{Sigma_0} given nu_0.
#'   Importantly, \strong{Sigma_0 should be convolved with perceptual noise (i.e., add perceptual noise covariance matrix to
#'   the category variability covariance matrices when you specify \code{Sigma_0})} since the stancode for the inference of the
#'   NIW ideal adaptor does \emph{not} infer category and noise variability separately.
#' @param control A list of control parameters that only experienced users should change since it can change the fitting and
#'   interpretation of the model. See \code{\link{control_staninput}} for details.
#'
#' @return A list consisting of:
#' \itemize{
#'   \item{`data`: }{A data.frame with the exposure and test data after exclusion of NAs and other checks.}
#'   \item{`staninput`: }{A named list of variables and values to be handed to Stan.}
#'   \item{`transform_information`: }{A list with information about the transformation that was applied to the data.}
#'}
#'
#' @seealso \code{\link{is.ideal_adaptor_staninput}}
#' @keywords TBD
#'
#' @importFrom purrr map_lgl map_int
#' @importFrom rstan nlist
#' @rdname make_staninput
#' @export
make_staninput <- function(
    exposure, test,
    cues, category = "category", response = "response",
    group = "group", group.unique = NULL,
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    control = control_staninput(),
    stanmodel = "NIW_ideal_adaptor",
    verbose = F
) {
  stopifnot(is.list(control))
  stopifnot(all(c("tau_scale", "L_omega_eta", "split_loglik_per_observation", "transform_type") %in% names(control)))
  tau_scale <- control$tau_scale
  if (length(tau_scale) == 1) tau_scale <- rep(tau_scale, length(cues))
  L_omega_eta <- control$L_omega_eta
  split_loglik_per_observation <- control$split_loglik_per_observation
  transform_type <- control$transform_type

  if (!is.null(lapse_rate)) {
    assert_that(is.number(lapse_rate), msg = "If not NULL, lapse_rate must be a number.")
    assert_that(between(lapse_rate, 0, 1), msg = "If not NULL, lapse rate must be a number between 0 and 1.")
  }

  assert_that(
    transform_type %in% c("identity", "center", "standardize", "PCA whiten", "ZCA whiten"),
    msg = paste("transform_type must be one of the following: identity, center, standardize, PCA whiten, ZCA whiten."))

  cues <- unique(cues)
  n.cues <- length(cues)
  if (stanmodel == "NIX_ideal_adaptor") {
    assert_that(n.cues == 1, msg = paste0("NIX_ideal_adaptor requires univariate data (!= ", n.cues, " cues found)."))
  } else if (stanmodel == "MNIX_ideal_adaptor") {
    assert_that(n.cues >= 2, msg = paste0("MNIX_ideal_adaptor requires multivariate data (!= ", n.cues, " cue found)."))
  }

  assert_that(
    length(tau_scale) == n.cues,
    msg = paste0("tau_scale must be a vector with the same number of elements as cues (", n.cues, ")."))

  exposure <-
    check_exposure_test_data(
      data = exposure,
      cues = cues,
      category = category,
      response = NULL,
      group = group,
      which.data = "exposure",
      verbose = verbose)

  if (!is.null(group.unique)) {
    assert_that(group.unique %in% names(exposure),
                msg = paste("Column for group.unique ", group.unique, "not found in exposure data."))
    if (verbose)
      message(paste0("Collapsing *exposure* observations to unique values of group.unique (", group.unique, ") by
      discarding the data from all but the first group member. This means that each unique exposure condition will
      only be counted once. All test observations are still counted, but aggregated for each unique value of group.unique."))

    exposure %<>%
      mutate(across(all_of(group.unique), as.factor)) %>%
      group_by(!! sym(group.unique), !! sym(category), !!! syms(cues)) %>%
      filter(!! sym(group) == unique(!! sym(group))[1])

    # Keep information about original group around so that it can be stored in the data object below
    group.original <- group
    group <- group.unique
  }

  # Keep only necessary columns but store original exposure data to attach it below
  exposure %<>% ungroup()
  exposure_original <- exposure %>% select(all_of(c(group.unique, group.original, category, cues)))
  exposure %<>% select(all_of(c(group, category, cues)))

  test <-
    check_exposure_test_data(
      data = test,
      cues = cues,
      category = NULL,
      response = response,
      group = group,
      which.data = "test",
      verbose = verbose)

  # Store original test data
  test_original <- test %>% select(all_of(c(group.unique, group.original, cues, response)))
  test %<>% select(all_of(c(group, cues, response)))
  exposure_original %<>% mutate(across(all_of(group.unique), ~ factor(.x, levels = levels(test_original[[group.unique]]))))
  exposure_original %<>% mutate(across(all_of(group), ~ factor(.x, levels = levels(test_original[[group]]))))

  # -----------------------------------------------------------------
  # Check whether exposure and test data are aligned in terms of factor levels
  # -----------------------------------------------------------------
  assert_that(all(levels(exposure[[category]]) == levels(test[[response]])),
              msg = paste("category variable", category, "in exposure data and response variable", response, "in test data must be factors with the same levels in the same order. Either the levels do not match, or they are not in the same order."))
  assert_that(all(levels(exposure[[group]]) %in% levels(test[[group]])),
              msg = paste("All levels of the grouping variable", group, "found in exposure must also be present in test."))
  if (verbose && !all(levels(test[[group]]) %in% levels(exposure[[group]])))
    message(paste("Not all levels of the grouping variable", group, "that are present in test were found in exposure.
    This is expected if and only if the data contained a test prior to (or without any) exposure.
    Creating 0 exposure data for these groups and aligning factor levels for group across exposure and test data."))
  exposure %<>% mutate(across(all_of(group), ~ factor(.x, levels = levels(test[[group]]))))

  if (!is.null(mu_0)) {
    if (nlevels(exposure[[category]]) == 1) {
      assert_that(is.vector(mu_0),
                  msg = "If mu_0 is not NULL and there is only one category, mu_0 must be a vector.")
    } else {
      assert_that(is.list(mu_0) & length(mu_0) == nlevels(exposure[[category]]),
                  msg = "If mu_0 is not NULL, mu_0 must be a list of vectors with as many elements as there are categories.")
    }
    assert_that(all(map_lgl(mu_0, is.numeric), map_lgl(mu_0, ~ is.null(dim(.x)) | length(dim(.x)) == 1)),
                msg = "If mu_0 is a list, each element must be a vector.")
    assert_that(all((map_int(mu_0, length)) == n.cues),
                msg = paste0(
                  "At least one element of mu_0 does not have the correct dimensionality. Observations have ",
                  n.cues,
                  " dimensions. Dimensionality of mu_0 ranges from ",
                  paste(map_int(mu_0, length) %>% range(), collapse = " to "),
                  "."))
  }
  if (!is.null(Sigma_0)) {
    if (nlevels(exposure[[category]]) == 1) {
      assert_that(is.array(Sigma_0),
                  msg = "If Sigma_0 is not NULL and there is only one category, Sigma_0 must be a positive-definite  matrix.")
    } else {
      assert_that(is.list(Sigma_0) & length(Sigma_0) == nlevels(exposure[[category]]),
                  msg = "If Sigma_0 not NULL, Sigma_0 must be a list of positive-definite matrices with as many elements as there are categories.")
    }
    assert_that(all(map_lgl(Sigma_0, is.numeric), map_lgl(Sigma_0, ~ length(dim(.x)) == 2)),
                msg = "If Sigma_0 is a list, each element must be a k x k matrix.")
    assert_that(all(map_lgl(Sigma_0, ~ all(dim(.x) == n.cues))),
                msg = paste0(
                  "At least one element of Sigma_0 does not have the correct dimensionality. Observations have ",
                  n.cues,
                  " dimensions. Sigma_0 includes matrices of dimension ",
                  paste(paste(map(Sigma_0, ~ dim(.x) %>% paste(collapse = " x "))) %>% unique(), collapse = ", "),
                  "."))
  }

  # -----------------------------------------------------------------
  # Transform data and category representations (if provided)
  # (for now, store original untransformed data and parameters so that they can also be attached to the stanfit object)
  # -----------------------------------------------------------------
  exposure_untransformed <- exposure
  test_untransformed <- test
  mu_0_untransformed <- mu_0
  Sigma_0_untransformed <- Sigma_0
  test_counts_untransformed <-
    get_test_counts(
      test = test,
      cues = cues,
      response = response,
      group = group,
      verbose = verbose)

  transform <- get_affine_transform(exposure, cues, transform_type)
  exposure <- transform[["transform.function"]](exposure)
  test <- transform[["transform.function"]](test)
  if (!is.null(mu_0)) mu_0 %<>% map(~ transform_category_mean(m = .x, transform))
  if (!is.null(Sigma_0)) Sigma_0 %<>% map(~ transform_category_cov(S = .x, transform))
  test_counts <-
    get_test_counts(
      test = test,
      cues = cues,
      response = response,
      group = group,
      verbose = verbose)

  # -----------------------------------------------------------------
  # Prepare nlist for Stan input
  # -----------------------------------------------------------------
  # Only for the NIX model: simplify some parameters to scalars
  K <- length(cues)
  M <- nlevels(exposure[[category]])
  L <- nlevels(exposure[[group]])

  staninput <-
    nlist(
      K, M, L,

      tau_scale = tau_scale %>% to_array(inner_dims = K, simplify = if (stanmodel == "NIX_ideal_adaptor") T else F),
      L_omega_eta,
      split_loglik_per_observation,

      lapse_rate_known = if (is.null(lapse_rate)) 0 else 1,
      lapse_rate_data = lapse_rate %>% to_array(outer_dims = if (is.null(lapse_rate)) NULL else 1),

      p_cat = rep(1/M, M) %>% to_array(), # Currently only used by MNIX

      mu_0_known = if (is.null(mu_0)) 0 else 1,
      Sigma_0_known = if (is.null(Sigma_0)) 0 else 1,
      mu_0_data = mu_0 %>% to_array(inner_dims = K, outer_dims = M, simplify = if (stanmodel == "NIX_ideal_adaptor") T else F),
      Sigma_0_data = Sigma_0 %>% to_array(inner_dims = c(K, K), outer_dims = M, simplify = if (stanmodel == "NIX_ideal_adaptor") T else F),

      shift = transform$transform.parameters[["shift"]] %>% to_array(inner_dims = K, simplify = if (stanmodel == "NIX_ideal_adaptor") T else F),
      INV_SCALE = transform$transform.parameters[["INV_SCALE"]] %>% to_array(inner_dims = c(K, K), simplify = if (stanmodel == "NIX_ideal_adaptor") T else F)
    )

  staninput_untransformed <- staninput
  staninput_untransformed[["mu_0_data"]] <- mu_0_untransformed %>% to_array(inner_dims = K, outer_dims = M, simplify = if (stanmodel == "NIX_ideal_adaptor") T else F)
  staninput_untransformed[["Sigma_0_data"]] <- Sigma_0_untransformed %>% to_array(inner_dims = c(K, K), outer_dims = M, simplify = if (stanmodel == "NIX_ideal_adaptor") T else F)

  # -----------------------------------------------------------------
  # Add the more model-specific parts of the staninput
  # (could be further consolidated later)
  # -----------------------------------------------------------------
  if (stanmodel == "NIX_ideal_adaptor") {
    staninput %<>%
      append(
        make_staninput_for_NIX_ideal_adaptor(
          exposure = exposure, test_counts = test_counts,
          cues = cues, category = category, group = group,
          verbose = verbose))

    staninput_untransformed %<>%
      append(
        make_staninput_for_NIX_ideal_adaptor(
          exposure = exposure_untransformed, test_counts = test_counts_untransformed,
          cues = cues, category = category, group = group,
          verbose = verbose))
  } else if (stanmodel == "NIW_ideal_adaptor") {
    staninput %<>%
      append(
        make_staninput_for_NIW_ideal_adaptor(
          exposure = exposure, test_counts = test_counts,
          cues = cues, category = category, group = group,
          verbose = verbose))

    staninput_untransformed %<>%
      append(
        make_staninput_for_NIW_ideal_adaptor(
          exposure = exposure_untransformed, test_counts = test_counts_untransformed,
          cues = cues, category = category, group = group,
          verbose = verbose))
  } else if (stanmodel == "MNIX_ideal_adaptor") {
    staninput %<>%
      append(
        make_staninput_for_MNIX_ideal_adaptor(
          exposure = exposure, test_counts = test_counts,
          cues = cues, category = category, group = group,
          verbose = verbose))

    staninput_untransformed %<>%
      append(
        make_staninput_for_MNIX_ideal_adaptor(
          exposure = exposure_untransformed, test_counts = test_counts_untransformed,
          cues = cues, category = category, group = group,
          verbose = verbose))
  } else {
    stop2("Model type ", stanmodel, " not recognized.")
  }
  # (For now): combine transformed and untransformed versions of staninput
  staninput <-
    list(
      transformed = staninput,
      untransformed = staninput_untransformed)

  # Bind exposure and test data as processed, as long with factor level information used above.
  # Also attach information about variable name mapping as attribute to data.frame
  data <-
    bind_rows(
      exposure_original[, c(group.unique, group.original, category, cues)] %>%
        drop_na() %>%
        mutate(Phase = "exposure"),
      test_original[, c(group.unique, group.original, response, cues)] %>%
        drop_na() %>%
        mutate(Phase = "test")) %>%
    mutate(
      Phase = factor(Phase, levels = c("exposure", "test")),
      !! sym(group.unique) := factor(!! sym(group.unique), levels = levels(exposure_original[[.env$group.unique]])),
      !! sym(group.original) := factor(!! sym(group.original), levels = levels(exposure_original[[.env$group.original]])),
      !! sym(category) := factor(!! sym(category), levels = levels(exposure_original[[.env$category]])),
      !! sym(response) := factor(!! sym(response), levels = levels(exposure_original[[.env$response]]))) %>%
    relocate(Phase, all_of(c(group.unique, group, category, cues, response)))

  attr(data, "group.unique") <- group.unique
  attr(data, "group") <- group.original
  attr(data, "category") <- category
  attr(data, "cues") <- cues
  attr(data, "response") <- response

  return(
    list(
      staninput = staninput,
      data = data,
      transform_information = transform))
}

make_staninput_for_NIX_ideal_adaptor <- function(
    exposure, test_counts,
    cues, category, group,
    verbose = F
) {
  staninput <-
    exposure %>%
    get_category_statistics_as_list_of_arrays(
      cues = cues, category = category, group = group,
      simplify = list(T, T, T),
      N_exposure = length, x_mean_exposure = mean, x_sd_exposure = sd) %>%
    within({
      # x_test is different from definition in make_staninput_for_NIW_ideal_adaptor (since vector, rather than matrix, is expected)
      # could be changed if stan program instead is changed to accept matrix and then turn it into vector.

      # It might be possible to simplify these steps through use of get_category_statistics_as_list_of_arrays
      x_test <-
        test_counts[[cues]] %>%
        as.numeric() %T>%
        { attr(., which = "cues") <- cues }
      y_test <-
        test_counts[[group]] %>%
        as.numeric() %T>%
        # This should work since we check above that exposure and test have the same levels for group
        { attr(., which = "levels") <- levels(exposure[[group]]) }
      z_test_counts <-
        test_counts %>%
        mutate(rownames = paste0("group=", !! sym(group), "; ", paste(cues, collapse = "-"), "=", paste(!!! syms(cues), sep = ","))) %>%
        column_to_rownames("rownames") %>%
        # This should work since we check above that exposure categories and test responses have the same levels
        select(levels(exposure[[category]])) %>%
        as.matrix()
      N_test <- length(x_test)
    })

  return(staninput)
}

make_staninput_for_NIW_ideal_adaptor <- function(
    exposure, test_counts,
    cues, category, group,
    verbose = F
) {
  # It might be possible to collapse the different make_staninput functions even further since they seem to only differ in
  # a) the agregate functions, b) that the multivariate models need simplify = F in a few places where to_array() is used,
  # and c) some specifics of how the test data is entered (again the multivariate and univariate models differ in the
  # dimensionality of the input they expect)
  staninput <-
    exposure %>%
    get_category_statistics_as_list_of_arrays(
      cues = cues, category = category, group = group,
      # Different from make_staninput_for_NIX_ideal_adaptor (since mean and cov are vector and matrix):
      fill = list(0, rep(0, length(cues)), diag(length(cues))),
      simplify = list(T, F, F),
      verbose = verbose,
      N_exposure = nrow, x_mean_exposure = colMeans, x_ss_exposure = get_sum_of_uncentered_squares_from_df) %>%
    within({
      # It might be possible to simplify these steps through use of get_category_statistics_as_list_of_arrays
      x_test <-
        test_counts %>%
        select(all_of(cues)) %>%
        as.matrix() %T>%
        { attr(., which = "cues") <- cues }
      y_test <-
        test_counts[[group]] %>%
        as.numeric() %T>%
        { attr(., which = "levels") <- levels(exposure[[group]]) }
      z_test_counts <-
        test_counts %>%
        mutate(rownames = paste0("group=", !! sym(group), "; ", paste(cues, collapse = "-"), "=", paste(!!! syms(cues), sep = ","))) %>%
        column_to_rownames("rownames") %>%
        select(levels(exposure[[category]])) %>%
        as.matrix()
      N_test <- nrow(x_test)
    })

  return(staninput)
}

make_staninput_for_MNIX_ideal_adaptor <- function(
    exposure, test_counts,
    cues, category, group,
    verbose = F
) {
  staninput <-
    get_category_statistics_as_list_of_arrays(
      data = exposure,
      cues = cues, category = category, group = group,
      fill = list(0, rep(0, length(cues)), diag(length(cues))),
      # Different from make_staninput_for_NIX_ideal_adaptor (since mean and cov are vector and matrix):
      simplify = list(T, F, F),
      verbose = verbose,
      N_exposure = nrow, x_mean_exposure = colMeans, x_cov_exposure = cov) %>%
    within({
      # It might be possible to simplify these steps through use of get_category_statistics_as_list_of_arrays
      x_test <-
        test_counts %>%
        select(all_of(cues)) %>%
        as.matrix() %T>%
        { attr(., which = "cues") <- cues }
      y_test <-
        test_counts[[group]] %>%
        as.numeric() %T>%
        { attr(., which = "levels") <- levels(exposure[[group]]) }
      z_test_counts <-
        test_counts %>%
        mutate(rownames = paste0("group=", !! sym(group), "; ", paste(cues, collapse = "-"), "=", paste(!!! syms(cues), sep = ","))) %>%
        column_to_rownames("rownames") %>%
        select(levels(exposure[[category]])) %>%
        as.matrix()
      N_test <- nrow(x_test)
    })

  return(staninput)
}


#' Compose data to fit ideal_adaptor_stanfit via rstan
#'
#' DEPRECATED: Use \code{make_staninput(stanmodel = "NIW_ideal_adaptor")} instead.
#'
#' @export
compose_data_to_infer_NIW_ideal_adaptor <- function(stanmodel = "NIW_ideal_adaptor", ...)
  make_staninput(..., stanmodel = "NIW_ideal_adaptor")
