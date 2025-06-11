#' @export
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


#' @export
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



#' Get category statistics (functions) from data as list of lists or list of arrays
#'
#' Get category statistics (aggregate functions) over `cues` in a data frame. Each function is calculated for each unique combination
#' of values in the columns `category` and (optionally) `group`.
#'
#' @param data `tibble` or `data.frame` with the data. Each row should be an observation of a category,
#'   and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param cues Names of columns with cue values.
#' @param category Name of column that contains the category label for the exposure data. This column must be a factor.
#'   Can be `NULL` for unsupervised updating (not yet implemented). (default: "category")
#' @param group Name of column that contains information about which observations form a group. This column must be
#'   a factor.
#' @param simplify A list of logicals of the same length as `...` indicating whether the array resulting from each function
#'   should be simplified? See \code{\link{to_array}} for more detail. (default: `TRUE` for all functions)
#' @param ... A named list of functions (statistics) of `cues` to calculate for each unique combination of `group` and
#'   `category`. Functions must return `NA` for any combination `category` and `group` for which `data` does not contain
#'   any observations.
#'
#' @return A named list of length `...` (names are the names of the functions that have been computed). The list will be
#'   sorted first by `group` and then by `category`, in ascending order of the levels of those variables. The elements of
#'   the list will either be lists (for `get_category_statistics_as_list_of_lists`) or arrays (for `get_category_statistics_as_list_of_arrays`).
#'   Missing values---resulting from combinations of `group` and `category` for which there is no data---will be filled
#'   with `NA`s coerced into the same structure as all other outputs for that function (e.g., if the function f results
#'   in a 2-element vector for all combinations of `group` and `category` for which there is data, then unobserved data
#'   results in a 2-element vector of `NA`s).
#'
#' @keywords TBD
#' @rdname get_category_statistics_as_list
#' @export
get_category_statistics_as_list_of_lists <- function(
    data,
    group,
    category,
    cues,
    verbose = F,
    ...
) {
  # Capture the named functions
  fn_list <- list(...)

  if (!is.factor(data[[group]])) {
    stop2("Group variable must be a factor.")
  }

  if (!is.factor(data[[category]])) {
    stop2("Category variable must be a factor.")
  }

  data %<>%
    select(all_of(group), all_of(category), all_of(cues)) %>%
    # Remove "." from group or category names to prevent it from leading to issues below
    # (since split uses "." as a separator)
    mutate(
      across(
        all_of(c(group, category)),
        ~ factor(
          gsub("\\.", "_", .x),
          levels = gsub("\\.", "_", levels(.x)))))

  # Get sorted levels
  group_levels <- levels(data[[group]])
  category_levels <- levels(data[[category]])

  # Split the data by group and category
  split_data <- split(data, list(data[[group]], data[[category]]), drop = FALSE)

    # Sort the split list by factor levels of group (first) and category (second)
  split_keys <- names(split_data)
  split_order <- order(
    match(sapply(strsplit(split_keys, "\\."), `[`, 1), group_levels),
    match(sapply(strsplit(split_keys, "\\."), `[`, 2), category_levels))
  split_data <- split_data[split_order]

  # Initialize result list
  result <- lapply(fn_list, function(f) vector("list", length(split_data)))
  names(result) <- names(fn_list)

  # Apply functions
  for (i in seq_along(split_data)) {
    x <- do.call(rbind, lapply(seq_len(nrow(split_data[[i]])), function(j) {
      as.numeric(split_data[[i]][j, cues, drop = FALSE])
    }))
    for (fname in names(fn_list)) {
      # Call function and make sure that NA is returned instead if x is empty data frame
      result[[fname]][[i]] <- if (is.null(x)) NA else fn_list[[fname]](x)
    }
  }

  # For each function, get the highest dimensionality of any of its outputs, and then coerce
  # NA outputs into that same dimensionality.
  dim_target <- result %>% map(~ reduce((map(.x, dim2)), pmax))
  result %<>%
    map2(
      .y = dim_target,
      ~ map(
        .x,
        function(x) {
          if (any(is.na(x))) {
            x <-
              if (length(.y) == 1) {
                # Coerce NA into vector
                rep(NA, .y)
              } else if (length(.y) == 2) {
                # Coerce NA into matrix
                matrix(rep(NA, prod(.y)), nrow = .y[1], ncol = .y[2])
              } else if (length(.y) > 2) {
                # Coerce NA into array
                array(rep(NA, prod(.y)), dim = .y)
              } else NA

            return(x)
          } else return(x)
        }))

  return(result)
}

#' @rdname get_category_statistics_as_list
#' @export
get_category_statistics_as_list_of_arrays <- function(
    data,
    group,
    category,
    cues,
    simplify = as.list(rep(T, length(list(...)))),
    verbose = F,
    ...
) {
  if (!is.list(simplify) || !all(map_lgl(simplify, is.logical))) stop2("Argument simplify must be a list of logicals.")
  if (length(simplify) != length(list(...))) stop2("Simplify must be a list of equal length as the number of functions provided in `...`.")

  result <-
    get_category_statistics_as_list_of_lists(
      data = data,
      group = group,
      category = category,
      cues = cues,
      verbose = verbose,
      ...)

  # (for now, only) name outer dimensions of array
  dimnames <- list(levels(data[[category]]), levels(data[[group]]))
  result %<>%
    map2(
      .y = simplify,
      .f =
        function(.x, .y) {
          to_array(
          .x,
          outer_dims = c(nlevels(data[[category]]), nlevels(data[[group]])),
          dimnames = dimnames,
          simplify = .y)
      })

  return(result)
}


#' Turn (lists of) numerics into arrays (for Stan input)
#'
#' Takes `NULL`, numeric atomics (single scalars, vectors, matrices) or lists of numeric elements as inputs and turns them into numeric arrays.
#'
#' @param x The input to be turned into an array.
#' @param inner_dims Integer vector of intended dimensions of x.
#'   If `NULL`, `inner_dims` will be inferred from the input.
#'   If not `NULL`, this will be used for checks and to enforce the correct dimensions of `x`. (default: `NULL`)
#' @param outer_dims An integer vector of outer dimensions for the array. If `x` is not a list, an array of `x` will be repeated for each
#'   outer dimension. If `x` is list, the outer dimensions will each index one element of `x`.
#'   The total length of `x` must match the product of outer dimensions. The first dimension will be iterated
#'   over first (and thus alternating fastest), the second dimension will be iterated over second, etc. See details.
#' @param dimnames A list of character vectors with names for each dimension of the array. If not `NULL`, the length
#'   of this list must match the total number of dimensions, and each character vector must match the length of the
#'   corresponding dimension. Dimensions are named from the outside-in (i.e. first outer than inner dimensions).
#'   (default: `NULL`)
#' @param simplify Should the array be simplified by removing inner dimensions whenever the inner element(s) are just scalars?
#'   (e.g., if the inner elements are a one-element vector or a 1x1 matrix). Note that atomic values that can be simplified
#'   will be returned as scalars (and thus without dimnames) if there are no outer dimensions, since the returned object cannot
#'   be an array in that case. (default: `TRUE`)
#'
#' @return An array.
#'
#' @details
#' \itemize{
#'   \item{`NULL` will be turned in a `array(numeric(0), dim = c(.outer_dims, .inner_dims))`. This is sometimes required for
#'   optional Stan inputs.}
#'   \item{Atomics (scalars, vectors, and matrices) will be turned into arrays of appropriate dimensions.}
#'   \item{Lists will be turned into arrays of appropriate dimensions.}
#' }
#'
#' E.g., if the input is a list of six 2x2 covariance matrices, and `...` is `3, 2`, then output will
#' be an array of dimensionality `c(3, 2, 2, 2)`.
#'
#'
#' @keywords TBD
#' @export
to_array <- function(
    x,
    inner_dims = NULL,
    outer_dims = NULL,
    dimnames = NULL,
    simplify = T
) {
  stopifnot(is.null(inner_dims) || all(inner_dims == round(inner_dims)))
  stopifnot(is.null(outer_dims) || all(outer_dims == round(outer_dims)))

  if (is.null(x)) {
    # optionally simplify by removing inner_dims if all inner dimension have length 1 (i.e., if all elements are intended to be scalars)
    if (simplify) { inner_dims <- if (all(inner_dims == 1)) NULL else inner_dims }

    total_dim_length <- length(inner_dims) + length(outer_dims)
    arr <- array(numeric(), dim = if (total_dim_length == 0) 0 else rep(0, total_dim_length))
  } else
    # Handle atomic input
    if (is.atomic(x)) {
    # NOTE: no simplification is currently being applied for atomic elements.
    found_inner_dims <- dim2(x)

    if (is.null(inner_dims)) inner_dims <- found_inner_dims
    if (any(found_inner_dims != inner_dims))
      stop2(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ", "), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))

    # If desired inner_dims don't match found inner_dims, see whether the input can be coerced into the desired inner_dims format
    if (length(found_inner_dims) != length(inner_dims)) {
      if (prod(inner_dims) == prod(found_inner_dims)) {
        if (length(inner_dims) == 1) {
          # Coerce into scalar or vector format
          x <- as.vector(x)
        } else if (length(inner_dims) > 1) {
          # Coerce into matrix format
          x <- matrix(x, nrow = inner_dims[1], ncol = inner_dims[2])
        }
      } else {
        stop2(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ", "), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))
      }
    }

    # optionally simplify by removing inner_dims if all inner dimension have length 1 (i.e., if all elements are essentially scalars)
    # Since array with no dimensions are not allowed, special handling is required when there are no outer dimensions and we still
    # want to simplify. In that case, x as a scalar without dimnames.
    if (simplify) { inner_dims <- if (all(inner_dims == 1)) { if (!is.null(outer_dims)) NULL else return(as.numeric(x)) } else inner_dims }

    # First cast an array with outer_dims on the outside (necessary to get correct handling of matrices when outer_dims are requested)
    total_dim_length <- length(inner_dims) + length(outer_dims)
    arr <- rep(x, prod(outer_dims)) %>% array(dim = if (total_dim_length == 0) 0 else c(inner_dims, outer_dims))
    # Move outer dimensions to the front
    perm <- c(if (is.null(outer_dims)) NULL else seq(length(inner_dims) + 1, total_dim_length), if (is.null(inner_dims)) NULL else 1:length(inner_dims))
    arr <- aperm(a = arr, perm = perm)
  } else
    # Handle list input
    if (is.list(x)) {
    # Use first list element to infer detected dimensions
    found_inner_dims <- dim2(x[[1]])

    if (is.null(inner_dims)) inner_dims <- found_inner_dims
    if (any(found_inner_dims != inner_dims)) stop2(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ","), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))

    # Check that x has the right number of elements
    expected_len <- if (length(outer_dims) == 0) 1 else prod(outer_dims)
    if (length(x) != expected_len) {
      stop2(paste0("Length of list (", length(x), ") input does not match product of provided outer dimensions (", expected_len, ")."))
    }

    # Combine list elements into a matrix
    arr <- simplify2array(x, except = NULL)

    # Move last dimensions to the front
    total_dims <- dim(arr)
    perm <- c(length(total_dims), 1:(length(total_dims) - 1))
    arr <- aperm(a = arr, perm = perm)

    # optionally simplify by removing inner_dims if all inner dimension have length 1 (i.e., if all elements are essentially scalars)
    if (simplify) { inner_dims <- if (all(inner_dims == 1)) NULL else inner_dims }
    # split first dimension into all outer dimensions
    arr <- array(arr, dim = c(outer_dims, inner_dims))
  } else
    stop2("Unsupported input type.")

  # Name dimensions of array
  # (no special handling and checking necessary since dimnames already does that: it exhaustively names from the outer to
  # inner dimensions until it hits a mismatch)
  dimnames(arr) <- dimnames

  return(arr)
}


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
#'   interpretation of the model:
#' \itemize{
#' \item {`tau_scale`: A vector of scales for the Cauchy priors for each cue's standard deviations. Used in
#'   both the prior for m_0 and the prior for S_0. (default: vector of `5`s, assuming that the data are standardized).}
#' \item {`L_omega_eta`: A vector of etas of the LKJ prior for the correlations of the covariance matrix of \code{mu_0}. Only used for
#'   models with multivariate categories (e.g., NIW_ideal_adaptor). (default: `1`,
#'   which corresponds to a uniform prior of correlation matrices)}
#' \item{`split_loglik_per_observation`: Optionally, split the log likelihood per observation. This can be helpful of leave-one-out
#'   estimation in order to avoid high Pareto k, but it also makes the stored stanfit object much larger. (default: `0`)}
#' \item{`transform_type`: An affine transformation that can be applied to the data. See `type` in \code{\link{get_affine_transform}}
#'    for details. (default: "standardize", which standardizes each cue separately)}
#' }
#'
#' @return A list consisting of:
#' \itemize{
#'   \item{`data`: A data.frame with the exposure and test data after exclusion of NAs and other checks.}
#'   \item{`staninput`: A named list of variables and values to be handed to Stan.}
#'   \item{`transform_information`: A list with information about the transformation that was applied to the data.}
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

    group <- group.unique
  }

  # Make sure data is ungrouped so that transform_cues works correctly, and keep only the necessary columns
  exposure %<>%
    ungroup() %>%
    select(all_of(c(group, cues, category)))

  test <-
    check_exposure_test_data(
      data = test,
      cues = cues,
      category = NULL,
      response = response,
      group = group,
      which.data = "test",
      verbose = verbose) %>%
    select(c(!! group, !!! cues, !! response))

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
  exposure %<>%
    mutate(across(all_of(group), ~ factor(.x, levels = levels(test[[!! group]]))))

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
  data <-
    bind_rows(
      exposure[, c(group.unique, if (group == group.unique) NULL else group, category, cues)] %>%
        drop_na() %>%
        mutate(Phase = "exposure"),
      test[, c(group.unique, if (group == group.unique) NULL else group, response, cues)] %>%
        drop_na() %>%
        mutate(Phase = "test")) %>%
    mutate(
      !! sym(group.unique) := factor(!! sym(group.unique), levels = levels(exposure[[.env$group.unique]])),
      !! sym(group) := factor(!! sym(group), levels = levels(exposure[[.env$group]])),
      !! sym(category) := factor(!! sym(category), levels = levels(exposure[[.env$category]])),
      !! sym(response) := factor(!! sym(response), levels = levels(test[[.env$response]]))) %>%
    relocate(Phase, all_of(c(group.unique, group, category, cues, response)))

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

  # Deal with empty exposure condition (pre-exposure tests)
  staninput$N_exposure %<>% replace_na_in_array(fill = 0)
  staninput$x_mean_exposure %<>% replace_na_in_array(fill = 0)
  staninput$x_sd_exposure %<>% replace_na_in_array(fill = 0)

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

  # Deal with empty exposure condition (pre-exposure tests)
  staninput$N_exposure %<>% replace_na_in_array(fill = 0)
  staninput$x_mean_exposure %<>% replace_na_in_array(fill = 0)
  staninput$x_ss_exposure %<>% replace_na_in_array(fill = 0)

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

  # Deal with empty exposure condition (pre-exposure tests)
  staninput$N_exposure %<>% replace_na_in_array(fill = 0)
  staninput$x_mean_exposure %<>% replace_na_in_array(fill = 0)
  staninput$x_cov_exposure %<>% replace_na_in_array(fill = 0)

  return(staninput)
}


#' Compose data to fit ideal_adaptor_stanfit via rstan
#'
#' DEPRECATED: Use \code{make_staninput(model_type = "NIW_ideal_adaptor")} instead.
#'
#' @export
compose_data_to_infer_NIW_ideal_adaptor <- function(model_type = "NIW_ideal_adaptor", ...)
  make_staninput(..., model_type = "NIW_ideal_adaptor")
