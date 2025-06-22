#' @importFrom tidybayes spread_draws recover_types
#' @importFrom tidyr spread gather unite
#' @importFrom dplyr relocate
NULL

# Functions general to stanfit objects ----------------------------------------

#' Get parameter names of stanfit
#'
#' Get the names for all parameters in `fit`.
#'
#' @param fit A \code{\link{stanfit}} object.
#' @param original_pars Should the original parameter names used in the Stan code be returned?
#'   This includes parameters for which no samples were stored. Alternatively, only the
#'   (renamed) parameters for which samples were stored are returned. (default: `FALSE`)
#' @export
get_params <- function(fit, original_pars = F) {
  stanfit <- get_stanfit(fit)

  names <- if (original_pars) stanfit@model_pars else names(stanfit)
  return(names)
}

#' Get number of post-warmup MCMC samples from stanfit
#'
#' Get the total number of post-warmup MCMC samples in `fit`.
#'
#' @param fit A \code{\link{stanfit}} object.
#' @export
get_number_of_draws <- function(fit) {
  stanfit <- get_stanfit(fit)
  return(length(stanfit@sim$samples[[1]][[1]]))
}

#' Get indices for random MCMC draws from stanfit
#'
#' Returns `ndraws` indices for random post-warmup MCMC draws (without replacement) from
#' `fit`.
#'
#' @param fit A \code{\link{stanfit}} object.
#' @param ndraws Number of indices to be returned. Can't be larger than total number of
#' post-warmup samples across all MCMC chains in `fit`.
#'
#' @return A numeric vector.
get_random_draw_indices <- function(fit, ndraws)
{
  n.all.draws <- get_number_of_draws(fit)
  assert_that(ndraws <= n.all.draws,
              msg = paste0("Cannot return ", ndraws, " draws because there are only ", n.all.draws, " in the object."))

  draws <- sample(1:n.all.draws, size = ndraws)
  return(draws)
}


# Functions general to NIW ideal adaptor stanfit objects ----------------------------------------

# Method exported for compatibility with rstan
#' @export
get_stanmodel.ideal_adaptor_stanfit <- function(x) {
  stanfit <- get_stanfit(x)

  return(get_stanmodel(stanfit))
}

#' Get the name of the stanmodel from an ideal adaptor stanfit
#'
#' Returns the name of the stanfit model.
#'
#' @param x \code{\link{ideal_adaptor_stanfit}} object.
#'
#' @return A character string.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
get_stanmodelname <- function(x, ...) {
  UseMethod("get_stanmodelname")
}

#' @rdname get_stanmodelname
#' @export
get_stanmodelname.ideal_adaptor_stanfit <- function(x) {
  stanfit <- get_stanfit(x)

  return(stanfit@model_name)
}

#' Get the stanfit from an ideal adaptor stanfit
#'
#' Returns the stanfit created by \code{stan} during the creation of the \code{\link{ideal_adaptor_stanfit}}
#' object.
#'
#' @param x \code{\link{ideal_adaptor_stanfit}} object.
#'
#' @return An \code{\link{rstan::stanfit}} object.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
get_stanfit <- function(x, ...) {
  UseMethod("get_stanfit")
}

#' @rdname get_stanfit
#' @export
get_stanfit.ideal_adaptor_stanfit <- function(x) {
  # Check if the stanfit is older than version 0.0.1.0015 when we introduced the version slot
  if (hasSlot(x, "version")) {
    return(x$stanfit)
  } else {
    stop("It appears that the model was fit with an old version of MVBeliefUpdatr (< 0.0.1.0015). Please refit the model.")
  }
}

#' Set the stanfit of an ideal adaptor stanfit
#'
#' Sets the stanfit of an the \code{\link{ideal_adaptor_stanfit}} object.
#'
#' @param x \code{\link{ideal_adaptor_stanfit}} object.
#' @param stanfit An \code{\link[rstan]{stanfit}} object of adequate structure.
#'
#' @return An \code{\link{ideal_adaptor_stanfit}} object with the updated stanfit.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
set_stanfit <- function(x, ...) {
  UseMethod("set_stanfit")
}

#' @rdname set_stanfit
#' @export
set_stanfit.ideal_adaptor_stanfit <- function(x, stanfit = NULL) {
  if (!is.null(stanfit)) {
    assert_that(is.stanfit(stanfit))
    assert_that(
      stanfit@model_name %in% names(MVBeliefUpdatr:::stanmodels),
      msg = paste0("stanfit object was not created by one of the accepted stancodes:\n\t",
                   paste(names(MVBeliefUpdatr:::stanmodels), collapse = "\n\t"),
                   "\n(you can get the name of your model from your_stanfit@model_name)."))
  }

  # Check if the stanfit is older than version 0.0.1.0015 when we introduced the version slot
  if (hasSlot(x, "version")) {
    x$stanfit <- stanfit
    return(x)
  } else {
    stop2("It appears that the model was fit with an old version of MVBeliefUpdatr (< 0.0.1.0015). Please refit the model.")
  }
}

#' Get the transform/untransform information from an ideal adaptor stanfit
#'
#' Returns the transform/untransform information handed to \code{stan} or \code{sampling} during the creation of the \code{stanfit}
#' object.
#'
#' @param x An \code{\link{ideal_adaptor_stanfit}} object.
#'
#' @return A function.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
get_transform_information <- function(x, ...) {
  UseMethod("get_transform_information")
}

#' @export
get_transform_function <- function(x, ...) {
  UseMethod("get_transform_function")
}

#' @export
get_untransform_function <- function(x, ...) {
  UseMethod("get_untransform_function")
}

#' @rdname get_transform_information
#' @export
get_transform_information.ideal_adaptor_stanfit <- function(x) {
  return(x$transform_information)
}

#' @rdname get_transform_information
#' @export
get_transform_function.ideal_adaptor_stanfit <- function(x) {
  return(get_transform_information(x)$transform.function)
}

#' @rdname get_transform_information
#' @export
get_untransform_function.ideal_adaptor_stanfit <- function(x) {
  return(get_transform_information(x)$untransform.function)
}

#' Get the input data from an ideal adaptor stanfit
#'
#' Returns the inputs handed to \code{stan} or \code{sampling} during the creation of the \code{\link{ideal_adaptor_stanfit}}
#' object.
#'
#' @param x \code{\link{ideal_adaptor_stanfit}} object.
#' @param which Should "transformed" or "untransformed" staninput be returned or "both"? (default: `"untransformed"`)
#'
#' @return A list with element names and structures determined by the type of stanfit model.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
get_staninput <- function(x, ...) {
  UseMethod("get_staninput")
}

#' @rdname get_staninput
#' @export
get_staninput.ideal_adaptor_stanfit <- function(x, which = c("untransformed", "transformed", "both")[1]) {
  if (which == "untransformed") {
    return(x$staninput$untransformed)
  } else if (which == "transformed") {
    return(x$staninput$transformed)
  } else if (which == "both") {
    return(x$staninput)
  } else {
    stop2('which must be one of "untransformed", "transformed", or "both".')
  }
}

#' Set the input data from an ideal adaptor stanfit
#'
#' Sets the inputs handed to \code{stan} or \code{sampling} during the creation of the \code{\link{ideal_adaptor_stanfit}}
#' object.
#'
#' @param x \code{\link{ideal_adaptor_stanfit}} object.
#' @param staninput A list with element names and structures determined by the type of stanfit model.
#' @param which Should "transformed" or "untransformed" staninput be set or "both"? (default: `"both"`)
#'
#' @return A list with element names and structures determined by the type of stanfit model.
#'
#' @return An \code{\link{ideal_adaptor_stanfit}} object with the updated staninput.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
set_staninput <- function(x, ...) {
  UseMethod("set_staninput")
}

#' @rdname set_staninput
#' @export
set_staninput.ideal_adaptor_stanfit <- function(x, staninput, which = c("untransformed", "transformed", "both")[3]) {
  if (which == "both") {
    x$staninput <- staninput
  } else if (which == "transformed") {
    x$staninput$transformed <- staninput
  } else if (which == "untransformed") {
    x$staninput$untransformed <- staninput
  } else {
    stop2('which must be one of "untransformed", "transformed", or "both".')
  }

  return(x)
}

#' Get exposure category statistic from ideal adaptor stanfit
#'
#' Returns the category means mu and/or category covariance matrix Sigma for the exposure data for an
#' \code{\link{ideal_adaptor_stanfit}}.
#'
#' @param x An \code{\link{ideal_adaptor_stanfit}}.
#' @param categories Character vector with categories for which category statistics are to be
#' returned. (default: all categories)
#' @param groups Character vector with groups for which category statistics are to be
#' returned. (default: all groups except `"prior"`)
#' @param statistic Which exposure statistic should be returned? `n` for number of observations, `mean` for
#' category mean, `css` or `uss` for centered or uncentered category sum-of-square matrix, `cov` for the
#' category covariance matrix, or a character vector with any combination thereof. (default: all)
#' @param untransform_cues DEPRECATED. Should statistics be transformed back into the original cue space? (default: `FALSE`)
#' @param ... additional arguments to \code{\link{get_staninput}}.
#'
#' @return If just one group and category was requested, a vector (for the mean) or matrix (for the covariance
#' matrix). If more than one group or category was requested, a tibble with one row for each unique combination
#' of group and category.
#'
#' @seealso TBD
#' @keywords TBD
#' @importFrom magrittr %<>%
#' @export
get_exposure_category_statistic <- function(x, ...) {
  UseMethod("get_exposure_category_statistic")
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_mean <- function(x, ...) {
  UseMethod("get_exposure_category_mean")
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_css <- function(x, ...) {
  UseMethod("get_exposure_category_css")
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_uss <- function(x, ...) {
  UseMethod("get_exposure_category_uss")
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_cov <- function(x, ...) {
  UseMethod("get_exposure_category_cov")
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_statistic.ideal_adaptor_stanfit <- function(
  x,
  categories = get_category_levels(x),
  groups = get_group_levels(x, include_prior = FALSE),
  statistic = c("n", "mean", "css", "uss", "cov"),
  untransform_cues = FALSE,
  ...
) {
  assert_that(all(statistic %in% c("n", "mean", "css", "uss", "cov")),
              msg = "statistic must be one of 'n', mean', 'css', 'uss', or 'cov'.")
  assert_that(any(is.factor(categories), is.character(categories), is.numeric(categories)))
  assert_that(any(is.factor(groups), is.character(groups), is.numeric(groups)))
  assert_that(all(categories %in% get_category_levels(x)),
              msg = paste("Some categories not found in the exposure data:",
                          paste(setdiff(categories, get_category_levels(x)), collapse = ", ")))
  assert_that(all(groups %in% get_group_levels(x)),
              msg = paste("Some groups not found in the exposure data:",
                          paste(setdiff(groups, get_group_levels(x, include_prior = FALSE)), collapse = ", ")))
  staninput <- get_staninput(x, ...)

  # Get names for requested dimensions from model object
  category_names <- get_category_levels(x)
  group_names <- get_group_levels(x)
  cue_names<- get_cue_levels(x)

  stanmodelname <- get_stanmodelname(x)
  df <- NULL
  # (If untransform_cues is TRUE, all statistics first need to be calculated because of the approach
  # taken below to untransform the statistics.)

  # Get counts n
  if (any(untransform_cues, c("n", "css", "cov") %in% statistic)) {
    n <- staninput$N_exposure
    d <- dim(n)
    if (!length(category_names)) category_names <- paste0("category_", seq_len(d[1]))
    if (!length(group_names)) group_names <- paste0("group_", seq_len(d[2]))
    dn <- list(category = category_names, group = group_names)

    df.n <- tibble()
    for (c in 1:d[1]) { # category
      for (g in 1:d[2]) { # group/condition
        df.n %<>%
          bind_rows(
            tibble(
              group = dn[[2]][g],
              category = dn[[1]][c],
              n = n[c, g]))
      }
    }

    df <- if (!is.null(df)) df %<>% left_join(df.n, by = c("group", "category")) else df.n
  }

  # Get central tendencies m
  if (any(untransform_cues, c("mean", "css", "cov") %in% statistic)) {
    m <- staninput$x_mean_exposure
    d <- dim(m)
    if (!length(category_names)) category_names <- paste0("category_", seq_len(d[1]))
    if (!length(group_names)) group_names <- paste0("group_", seq_len(d[2]))
    if (stanmodelname == "NIX_ideal_adaptor") {
      dn <- list(category = category_names, group = group_names)
    } else {
      if (!length(cue_names)) cue_names <- paste0("cues", seq_len(d[3]))
      dn <- list(category = category_names, group = group_names, cue = cue_names)
    }

    df.m <- tibble()
    # Get dimnames (skips if they are all null)
    for (c in 1:d[1]) { # category
      for (g in 1:d[2]) { # group/condition
        if (stanmodelname == "NIX_ideal_adaptor") {
          df.m %<>%
            bind_rows(
              tibble(
                group = dn[[2]][g],
                category = dn[[1]][c],
                value = m[c, g]))
        } else {
          for (f in 1:d[3]) { # cue
            df.m %<>%
              bind_rows(
                tibble(
                  group = dn[[2]][g],
                  category = dn[[1]][c],
                  cue = dn[[3]][f],
                  value = m[c, g, f]))
          }
        }
      }
    }

    # Pivot wide to have mean for all cues in the same row for any unique combination of group and category
    # (no additional grouping necessary)
    df.m %<>%
      pivot_wider(names_from = "cue", values_from = "value") %>%
      make_vector_column(cols = dn[[3]], vector_col = "mean", .keep = "unused")

    df <- if (!is.null(df)) df %<>% left_join(df.m, by = c("group", "category")) else df.m
  }

  # Get scatter or covariance matrices s
  if (any(untransform_cues, c("uss", "css", "cov") %in% statistic)) {
    if (stanmodelname == "NIX_ideal_adaptor") {
      s <- staninput$x_sd_exposure
      # For now: stopping here, but note that some conditionals have already been added below
      # (within this block for uss, css, and cov) to extract the currently adequate quantity
      # form the model. But it's probably easier to standardize the models and store e.g., the cov
      # for all types of models (or all stats), rather than storing different stats for each model.
      # The stancode could then transform the input data to the correct quantities. This would make
      # the handling here a lot easier.
      stop2("Extraction of uss, css, or cov not yet implemented for NIX ideal adaptor stanfit.")
    } else if (stanmodelname == "NIW_ideal_adaptor") {
      s <- staninput$x_ss_exposure
    } else if (stanmodelname == "MNIX_ideal_adaptor") {
      s <- staninput$x_cov_exposure
      stop2("Extraction of uss, css, or cov not yet implemented for MNIX ideal adaptor stanfit.")
    } else {
      stop2("Unrecognized stanmodel. No method available to extract category variance.")
    }

    d <- dim(s)
    if (!length(category_names)) category_names <- paste0("category_", seq_len(d[1]))
    if (!length(group_names)) group_names <- paste0("group_", seq_len(d[2]))
    if (stanmodelname == "NIX_ideal_adaptor") {
      dn <- list(category = category_names, group = group_names)
    } else {
      if (!length(cue_names)) cue_names <- paste0("cues", seq_len(d[3]))
      dn <- list(category = category_names, group = group_names, cue = cue_names, cue2 = cue_names)
    }

    df.s <- tibble()
    for (c in 1:d[1]) { # category
      for (g in 1:d[2]) { # group/condition
        if (stanmodelname == "NIX_ideal_adaptor") {
          df.s %<>%
            bind_rows(
              tibble(
                group = dn[[2]][g],
                category = dn[[1]][c],
                value = s[c, g, f1, f2]))
        } else {
          for (f1 in 1:d[3]) { # cue1
            for (f2 in 1:d[4]) { # cue2
              df.s %<>%
                bind_rows(
                  tibble(
                    group = dn[[2]][g],
                    category = dn[[1]][c],
                    cue = dn[[3]][f1],
                    cue2 = dn[[4]][f2],
                    value = s[c, g, f1, f2]))
            }
          }
        }
      }
    }

    df.s %<>%
      group_by(category, group) %>%
      summarise(uss = list(matrix(value, nrow = sqrt(length(value)))))

    df <- if (!is.null(df)) df %<>% left_join(df.s, by = c("group", "category")) else df.s
  }

  if (any(untransform_cues, c("css", "cov") %in% statistic)) {
    df %<>% mutate(css = pmap(.l = list(.data$uss, .data$n, .data$mean), uss2css))
  }

  if (any(untransform_cues, c("cov") %in% statistic)) {
    df %<>%
      mutate(cov = map2(.data$css, .data$n, css2cov))
  }

  df %<>%
    # Obtain untransformed uss and/or css from cov by walking back/undoing above steps
    # (direct transform of uss)
    # BE VERY CAREFUL IN CHANGING THIS ORDER OR MOVING PARTS OF THIS UP, BECAUSE OF THE DEPENDENCIES BETWEEN THE OPERATIONS
    { if (untransform_cues & "cov" %in% statistic) mutate(., cov = map(.data$cov, ~ untransform_category_cov(.x, get_transform_information(x)))) else . } %>%
    { if (untransform_cues & any(c("css", "uss") %in% statistic)) mutate(., css = map2(.data$cov, n, cov2css)) else . } %>%
    { if (untransform_cues & any(c("uss", "mean") %in% statistic)) mutate(., mean = map(.data$mean, ~ untransform_category_mean(.x, get_transform_information(x)))) else . } %>%
    { if (untransform_cues & "uss" %in% statistic) mutate(., uss = pmap(.l = list(.data$cov, .data$n, .data$mean), css2uss)) else . } %>%
    select(group, category, !!! syms(statistic)) %>%
    filter(., .data[["group"]] %in% .env[["groups"]]) %>%
    filter(., .data[["category"]] %in% .env[["categories"]]) %>%
    mutate(
      category = factor(category, levels = .env[["categories"]]),
      group = factor(group, levels = .env[["groups"]]))

  # If just one category, group, and statistic was requested, just return that object
  # (rather than the tibble)
  if (nrow(df) == 1 & length(statistic) == 1) df <- df[, statistic][[1]][[1]]
  return(df)
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_mean.ideal_adaptor_stanfit <- function(...) {
  return(get_exposure_category_statistic.ideal_adaptor_stanfit(..., statistic = "mean"))
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_css.ideal_adaptor_stanfit <- function(...) {
  return(get_exposure_category_statistic.ideal_adaptor_stanfit(..., statistic = "css"))
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_uss.ideal_adaptor_stanfit <- function(...) {
  return(get_exposure_category_statistic.ideal_adaptor_stanfit(..., statistic = "uss"))
}

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_cov.ideal_adaptor_stanfit <- function(...) {
  return(get_exposure_category_statistic.ideal_adaptor_stanfit(..., statistic = "cov"))
}


#' Get the exposure and/or test data from an ideal adaptor stanfit.
#'
#' Returns the exposure and/or test data used during the creation of the \code{\link[rstan]{stanfit}}.
#' object.
#'
#' @param x \code{\link{ideal_adaptor_stanfit}} object.
#' @param groups Character vector of groups for which test data is requested. Typically, the levels of these factors
#'   are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#'   \code{\link[tidybayes]{recover_types}} on the stanfit object to add or change these levels later.
#'   (default: all categories/groups will be selected)
#' @param .rename_to_MVB_default Should the data columns be renamed to the default names used internally in
#'   `MVBeliefUpdatr`? his option is included primarily for internal use. (default: `FALSE`)
#' @param .from_staninput Should the data be extracted from the staninput (rather than data) object stored in the
#'   stanfit object? This returns the data in a somewhat different format. This option is included primarily for
#'   internal use. (default: `FALSE`)
#' @param ... additional arguments to \code{\link{get_staninput}}.
#'
#' @return A \code{tibble} in which each row is a test token. Columns include the cues
#'   and the response counts (one column per category) for all test tokens and all groups.
#'
#' @seealso TBD
#' @keywords TBD
#' @rdname get_exposure_test_data
#' @export
get_data <- function(x, ...) {
  UseMethod("get_data")
}

#' @rdname get_exposure_test_data
#' @export
get_data.ideal_adaptor_stanfit <- function(
    x,
    groups = get_group_levels(x, include_prior = FALSE),
    .rename_to_MVB_default = FALSE,
    ...
) {
  data <- x$data

  group.unique <- attr(data, "group.unique")
  if (.rename_to_MVB_default) {
    group <- attr(data, "group")
    category <- attr(data, "category")
    cues <- attr(data, "cues")
    response <- attr(data, "response")

    data <-
      data %>%
      rename(
        group.unique = !! sym(group.unique),
        group = !! sym(group),
        category = !! sym(category),
        response = !! sym(response))

    for (c in 1:length(cues)) {
      data <- data %>% rename(!! sym(paste0("cue", c)) := !! sym(cues[c]))
    }

    group.unique <- "group.unique"
  }

  data %>% filter(!! sym(group.unique) %in% groups)
}

#' @rdname get_exposure_test_data
#' @export
get_exposure_data <- function(x, ...) {
  UseMethod("get_exposure_data")
}

#' @rdname get_exposure_test_data
#' @export
get_exposure_data.ideal_adaptor_stanfit <- function(
    x,
    groups = get_group_levels(x, include_prior = FALSE),
    ...
) {
  get_data(x, ...) %>%
    filter(Phase == "exposure")
}

#' @rdname get_exposure_test_data
#' @export
get_test_data <- function(x, ...) {
  UseMethod("get_test_data")
}

#' @rdname get_exposure_test_data
#' @export
get_test_data.ideal_adaptor_stanfit <- function(
  x,
  groups = get_group_levels(x, include_prior = FALSE),
  .from_staninput = FALSE, # for now included so as to allow the previous way of extracting test data from the staninput
  ...
) {
  if (.from_staninput) {
    data <- get_staninput(x, ...)
    data <-
      data[["x_test"]] %>%
      cbind(data[["z_test_counts"]]) %>%
      as_tibble(.name_repair = "minimal") %>%
      mutate(
        group.id = data[["y_test"]],
        group = factor(attr(data[["y_test"]], "levels")[group.id],
                       levels = attr(data[["y_test"]], "levels"))) %>%
      filter(group %in% groups)
  } else {
    data <-
      get_data(x, ...) %>%
      filter(Phase == "test")
  }

  return(data)
}


#' Get or restore the original group or category levels from an ideal adaptor stanfit.
#'
#' Checks if information is available about the original values and order of the factor levels
#' for the category variable (for which beliefs about means and covariances are inferred) or
#' group variable (e.g., subject or exposure group), respectively. If available,
#' that information is returned. `get_category_levels()` and `get_group_levels()` are
#' convenience functions, calling `get_staninput_variable_levels()`.
#'
#' @param x \code{\link{ideal_adaptor_stanfit}} object.
#' @param variable Either "category" or "group".
#' @param indeces A vector of category or group indices that should be turned into the original
#' category levels, or `NULL` if only the unique levels in their original order (as vector of characters)
#' should be returned. (default: `NULL`)
#'
#' @return If no category or group indices are provided, the levels of the category/group are returned (in the
#' original order). Otherwise a vector of the same length as \code{indices}
#'
#' @seealso \code{\link[tidybayes]{recover_types}} from tidybayes, \code{\link{get_constructor}}
#' @keywords TBD
#' @export
get_staninput_variable_levels <- function(x, ...) {
  UseMethod("get_staninput_variable_levels")
}

#' @rdname get_staninput_variable_levels
#' @export
get_category_levels <- function(x, ...) {
  UseMethod("get_category_levels")
}

#' @rdname get_staninput_variable_levels
#' @export
get_group_levels <- function(x, ...) {
  UseMethod("get_group_levels")
}

#' @rdname get_staninput_variable_levels
#' @export
get_cue_levels <- function(x, ...) {
  UseMethod("get_cue_levels")
}

#' @rdname get_staninput_variable_levels
#' @export
get_staninput_variable_levels.ideal_adaptor_stanfit <- function(x, variable = c("category", "group", "cue"), indices = NULL) {
  assert_that(is.null(indices) | all(indices > 0))
  f <- get_constructor(x, variable)

  if (is.null(indices)) return(levels(f(c()))) else return(f(indices))
}

#' @rdname get_staninput_variable_levels
#' @export
get_category_levels.ideal_adaptor_stanfit <- function(x, indices = NULL) {
  return(get_staninput_variable_levels(x, "category", indices))
}

#' @rdname get_staninput_variable_levels
#' @export
get_group_levels.ideal_adaptor_stanfit <- function(x, indices = NULL, include_prior = F) {
  groups <- get_staninput_variable_levels(x, "group", indices)
  if (include_prior) groups <- append("prior", groups)

  return(groups)
}

#' @rdname get_staninput_variable_levels
#' @export
get_cue_levels.ideal_adaptor_stanfit <- function(x, indices = NULL) {
  cues <- get_staninput_variable_levels(x, "cue", indices)

  return(cues)
}


#' Get tidybayes constructor from an ideal adaptor stanfit.
#'
#' Gets the tidybayes constructor function from the stanfit object. `get_category_constructor()` and
#' `get_group_constructor()` are convenience functions, calling `get_constructor()`. See \code{
#' \link[tidybayes]{recover_types}}. If variable is
#'
#' @param x \code{\link{ideal_adaptor_stanfit}} object.
#' @param variable Either "category" or "group". If set to `NULL` then a list of all constructors is
#' returned. That list is `NULL` if not tidybayes constructors are found in fit. (default: c("category", "group"))
#'
#' @return A constructor function, a list of constructor functions, or `NULL`. If a specific constructor
#' function is requested but not found, a warning is shown.
#'
#' @seealso \code{\link[tidybayes]{recover_types}} from tidybayes, \code{\link{get_staninput_variable_levels}}
#' @keywords TBD
#' @rdname get_constructor
#' @export
#' @importFrom rlang sym
get_constructor <- function(x, variable = NULL) {
  assert_that(class(x) %in% c("stanfit", "ideal_adaptor_stanfit"))
  if (class(x) == "ideal_adaptor_stanfit") stanfit <- get_stanfit(x) else stanfit <- x

  available_constructors <- c("category", "group", "cue", "cue2")
  if (is.null(variable)) return(attr(stanfit, "tidybayes_constructors"))

  assert_that(variable %in% available_constructors,
              msg = paste0("Variable name must be one of ", paste(available_constructors, collapse = "or"), "."))

  if (is.null(attr(stanfit, "tidybayes_constructors")[[rlang::sym(variable)]])) {
    warning(paste0(class(fit), " object ", deparse(substitute(fit)), " does not contain type information about the variable ", variable,
                   ". Applying recover_types() to the object might fix this."))
    return(NULL)
  }

  f <- attr(stanfit, "tidybayes_constructors")[[rlang::sym(variable)]]

  return(f)
}

#' @rdname get_constructor
#' @export
get_category_constructor <- function(x) {
  return(get_constructor(x, "category"))
}

#' @rdname get_constructor
#' @export
get_group_constructor <- function(x) {
  return(get_constructor(x, "group"))
}

#' @rdname get_constructor
#' @export
get_cue_constructor <- function(x) {
  return(get_constructor(x, "cue"))
}

#' @rdname get_constructor
#' @export
get_cue2_constructor <- function(x) {
  return(get_constructor(x, "cue2"))
}


#' Get expected category mean mu or covariance matrix sigma
#'
#' Returns the expected value of posterior marginal distribution over category means mu and/or
#' category covariance matrix Sigma, marginalized over all MCMC samples.
#'
#' Each MCMC samples' expected value for the category mean \code{E[mu] = m_n}
#' (i.e, the posterior/updated mean of the multivariate Normal over category means \code{mu}).
#' Marginalizing across all MCMC samples (representing uncertainty in the true value of
#' \code{m_n}), we get \code{E[E[mu]] = mean(m_n)}.
#'
#' Each MCMC samples' expected value for the category covariance matrix
#' \code{E[Sigma] = S_n / (nu_n - D - 1)}, where \code{S_n} is the posterior/updated scatter matrix,
#' \code{nu_n} is the posterior/updated pseudocount representing the strength of the posterior/updated
#' beliefs over category covariance matrices sigma (i.e., the inverse-Wishart), and \code{D} is
#' the dimension of the multivariate Normal. Marginalizing across all MCMC samples
#' (representing uncertainty in the true value of \code{S_n}), we get
#' \code{E[E[Sigma]] = mean(S_n / (nu_n - D - 1))}.
#'
#' @param x An \code{\link[=is.ideal_adaptor_stanfit]{mv_ibbu_stanfit}} object.
#' @param categories Character vector with categories for which category statistics are to be
#' returned. (default: all categories)
#' @param groups Character vector with groups for which category statistics are to be returned.
#' (default: all groups)
#' @param statistic Which category statistic should be returned? `mu` for category mean or `Sigma` for category
#' covariance matrix, or `c("mu", "Sigma")` for both. (default: both)
#' @param ... additional arguments to \code{\link{get_draws}}.
#'
#' @return If just one group and category was requested, a vector (for the mean) or matrix (for the covariance
#' matrix). If more than one group or category was requested, a tibble with one row for each unique combination
#' of group and category.
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
get_expected_category_statistic <- function(x, ...) {
  UseMethod("get_expected_category_statistic")
}

get_expected_mu <- function(x, ...) {
  UseMethod("get_expected_mu")
}

get_expected_sigma <- function(x, ...) {
  UseMethod("get_expected_sigma")
}

#' @rdname get_expected_category_statistic
#' @export
get_expected_category_statistic.ideal_adaptor_stanfit <- function(
  x,
  categories = get_category_levels(x),
  groups = get_group_levels(x, include_prior = TRUE),
  statistic = c("mu", "Sigma"),
  ...
) {
  assert_that(all(statistic %in% c("mu", "Sigma")))
  assert_that(any(is.factor(categories), is.character(categories), is.numeric(categories)))
  assert_that(any(is.factor(groups), is.character(groups), is.numeric(groups)))
  assert_that(all(categories %in% get_category_levels(x)),
              msg = paste("Some categories not found in model:",
                          paste(setdiff(categories, get_category_levels(x)), collapse = ", ")))
  assert_that(all(groups %in% get_group_levels(x, include_prior = T)),
              msg = paste("Some groups not found in model:",
                          paste(setdiff(groups, get_group_levels(x, include_prior = T)), collapse = ", ")))

  x <-
    get_draws(x, categories = categories, groups = groups, wide = F, nest = T, summarize = FALSE, ...) %>%
    mutate(Sigma = get_expected_Sigma_from_S(S, nu)) %>%
    group_by(group, category) %>%
    summarise(
      mu.mean = list(m %>% reduce(`+`) / length(m)),
      Sigma.mean = list(Sigma %>% reduce(`+`) / length(Sigma))) %>%
    ungroup() %>%
    select(group, category, !!! rlang::syms(paste0(statistic, ".mean"))) %>%
    mutate(
      category = factor(category, levels = .env[["categories"]]),
      group = factor(group, levels = .env[["groups"]]))

  # If just one category, group, and statistic was requested, just return that object
  # (rather than the tibble)
  if (nrow(x) == 1 & length(statistic) == 1) x <- x[, paste0(statistic, ".mean")][[1]][[1]]
  return(x)
}

#' @rdname get_expected_category_statistic
#' @export
get_expected_mu.ideal_adaptor_stanfit <- function(x, ...) {
  return(get_expected_category_statistic(x, statistic = "mu", ...))
}

#' @rdname get_expected_category_statistic
#' @export
get_expected_sigma.ideal_adaptor_stanfit <- function(x, ...) {
  return(get_expected_category_statistic(x, statistic = "Sigma", ...))
}


#' Get categorization function from a model object
#'
#' Returns a categorization function that can be used to categorize a set of cues into categories.
#'
#' @param x A model object. Currently, only \code{\link{ideal_adaptor_stanfit}} objects are supported.
#' @param lapse_treatment Should the consequences of attentional lapses be included in the categorization function
#'   ("marginalize") or not ("no_lapses")? (default: "marginalize")
#' @param ... Optionally, additional arguments handed to \code{\link{get_draws}}.
#'
#' @return A function that takes as input cue values and returns posterior probabilities of the first category,
#'   based on the posterior predictive of the cues given the ideal adaptor categories' parameters. The function
#'   will accept the following arguments:
#'   \itemize{
#'     \item{`target_category`:}{The index of the category for which categorization should be shown. (default: `1`)}
#'     \item{`logit`:}{Should the function return log-odds (TRUE) or probabilities (FALSE)? (default: FALSE)}
#'   }
#'
#' @rdname get_categorization_function
#' @export
get_categorization_function <- function(x, ...) {
  UseMethod("get_categorization_function")
}

#' @rdname get_categorization_function
#' @export
get_categorization_function.ideal_adaptor_stanfit <- function(
    x,
    lapse_treatment = c("no_lapses", "sample", "marginalize")[3],
    groups = get_group_levels(x, include_prior = TRUE),
    ...
) {
  d.pars <-
    get_draws(
      x,
      groups = groups,
      summarize = F,
      wide = F,
      ...)

  d.pars %<>%
    group_by(group, .draw) %>%
    do(f =
         get_categorization_function_from_stanfit_draws(
           .,
           # Set noise treatment to no noise for now since ideal adaptor stanfits
           # include the noise in the inferred category variability
           noise_treatment = "no_noise",
           lapse_treatment = lapse_treatment,
         ))

  return(d.pars)
}

get_categorization_function_from_stanfit_draws <- function(x, ...) {
  get_NIW_categorization_function(
    ms = x$m,
    Ss = x$S,
    kappas = x$kappa,
    nus = x$nu,
    lapse_rate = unlist(x$lapse_rate)[1],
    ...
  )
}


#' Get MCMC draws from an ideal adaptor stanfit
#'
#' Get MCMC draws of all parameters from incremental Bayesian belief-updating (IBBU) as a tibble. Both wide
#' (`wide=TRUE`) or long format (`wide=FALSE`) can be chosen as output. By default all post-warmup draws are
#' returned, but if `summarize=TRUE` then just the mean of each parameter is returned instead.
#'
#' By default, the category means and scatter matrices are nested, rather than each of their elements being
#' stored separately (`nest=TRUE`).
#'
#' @aliases add_draws
#'
#' @param fit \code{\link{ideal_adaptor_stanfit}} object.
#' @param which DEPRECATED. Use `groups` instead. Should parameters for the prior, posterior, or both be added? (default: `"posterior"`)
#' @param ndraws Number of random draws or `NULL` if all draws are to be returned. Only `draws` or `ndraws` should be non-zero. (default: `NULL`)
#' @param untransform_cues DEPRECATED Should m_0 and S_0 be transformed back into the original cue space? (default: `FALSE`)
#' @param summarize Should the mean of the draws be returned instead of all of the draws? (default: `FALSE`)
#' @param wide Should all parameters be returned in one row? (default: `FALSE`)
#' @param nest Should the category mean vectors and scatter matrices be nested into one cell each, or should each element
#' be stored in a separate cell? (default: `TRUE`)
#' @param categories Character vector of categories for which draws are to be returned. (default: all categories)
#' @param groups Character vector of groups for which draws are to be returned. (default: all groups)
#' @param seed A seed to use when subsampling draws (i.e. when ndraws is not NULL).
#'
#' @return tibble with post-warmup (posterior) MCMC draws of the prior/posterior parameters of the IBBU model
#' (\code{kappa, nu, m, S, lapse_rate}). \code{kappa} and \code{nu} are the pseudocounts that determine the strength of the beliefs
#' into the mean and covariance matrix, respectively. \code{m} is the mean of the multivariate normal distribution over category
#' means mu. \code{S} is the scatter matrix that determines both the covariance of the category means mu, and the
#' Inverse Wishart distribution over category covariance matrices Sigma.
#'
#' The expected value of the category mean mu is \code{m}. The expected value of the category covariance matrix Sigma
#' is \code{S / (nu - D - 1)}, where \code{D} is the dimension of the multivariate Gaussian category. For details,
#' \insertCite{@see @murphy2012 p. 134;textual}{MVBeliefUpdatr}.
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#'
#'
#' @importFrom assertthat assert_that
#' @importFrom tidyselect ends_with
#' @export
get_draws <- function(fit, ...) {
  UseMethod("get_draws")
}

#' @rdname get_draws
#' @export
get_draws.ideal_adaptor_stanfit <- function(
  fit,
  categories = get_category_levels(fit),
  groups = get_group_levels(fit, include_prior = TRUE),
  ##### SPECIAL HANDLING OF WHICH, WHICH IS NOW DEPRECATED.
  which = if ("prior" %in% groups) { if (length(groups) > 1) "both" else "prior" } else "posterior",
  ##### END OF SPECIAL HANDLING
  ndraws = NULL,
  ##### UNTRANSFORM_CUES (HERE AND ELSEWHERE) IS NO LONGER NEEDED, NOW THAT CUES ARE UNTRANSFORMED
  # WITHIN STAN. BUT WE CAN CHANGE THIS ARGUMENT TO TRANSFORM, AND ALLOW USERS TO HAND AN AFFINE
  # TRANSFORM THAT IS THEN APPLIED TO THE CUES?
  untransform_cues = FALSE,
  #####
  summarize = FALSE,
  wide = FALSE,
  nest = TRUE,
  seed = if (!is.null(ndraws)) runif(1, -1e6, 1e6) else NULL
) {
  # Binding variables that RMD Check gets confused about otherwise
  # (since they are in non-standard evaluations)
  .chain <- .iteration <- .draw <- group <- category <- kappa <- nu <- m <- S <- lapse_rate <- NULL

  assert_contains_draws(fit)
  assert_that(any(is.factor(categories), is.character(categories), is.numeric(categories)))
  assert_that(any(is.factor(groups), is.character(groups), is.numeric(groups)))
  assert_that(all(categories %in% get_category_levels(fit)),
              msg = paste("Some categories not found in model:",
                          paste(setdiff(categories, get_category_levels(fit)), collapse = ", ")))
  assert_that(all(groups %in% get_group_levels(fit, include_prior = T)),
              msg = paste("Some groups not found in model:",
                          paste(setdiff(groups, get_group_levels(fit, include_prior = T)), collapse = ", ")))

  ##### SPECIAL HANDLING OF WHICH, WHICH IS NOW DEPRECATED.
  assert_that(which %in% c("prior", "posterior", "both"),
              msg = "which must be one of 'prior', 'posterior', or 'both'.")
  ##### END OF SPECIAL HANDLING

  assert_that(any(is.null(ndraws), is.count(ndraws)),
              msg = "If not NULL, ndraw must be a count.")
  assert_that(any(is.null(ndraws), !is.null(seed)),
              msg = "If not ndraws is not NULL, seed must be specified.")
  assert_that(is.flag(summarize))
  assert_that(is.flag(wide))
  assert_that(!all(wide, !nest),
              msg = "Wide format is currently not implemented without nesting.")

  if ("prior" %in% groups & length(groups) > 1) {
    d.prior <-
      get_draws(
        fit = fit,
        categories = categories, groups = "prior",
        ndraws = ndraws,
        untransform_cues = untransform_cues,
        summarize = summarize, wide = wide, nest = nest, seed = seed)
    d.posterior <-
      get_draws(
        fit = fit, categories = categories, groups = setdiff(groups, "prior"),
        ndraws = ndraws,
        untransform_cues = untransform_cues,
        summarize = summarize, wide = wide, nest = nest, seed = seed)
    d.pars <-
      rbind(d.prior, d.posterior) %>%
      mutate(group = factor(.data$group, levels = c(levels(d.prior$group), levels(d.posterior$group))))

    return(d.pars)
  } else {
    # Parameters' names depend on whether prior or posterior is to be extracted.
    postfix <- if ("prior" %in% groups) "_0" else "_n"
    kappa <- paste0("kappa", postfix)
    nu <- paste0("nu", postfix)
    m <- paste0("m", postfix)
    S <- paste0("S", postfix)

    # Variables by which parameters are indexed
    pars.index <- if ("prior" %in% groups) "category" else c("category", "group")

    stanfit <- get_stanfit(fit)
    # Get non-nested draws
    if ("prior" %in% groups) {
      d.pars <-
        stanfit %>%
        spread_draws(
          !! rlang::sym(kappa),
          !! rlang::sym(nu),
          (!! rlang::sym(m))[!!! rlang::syms(pars.index), cue],
          (!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2],
          lapse_rate,
          ndraws = ndraws,
          seed = seed)
    } else {
      d.pars <-
        stanfit %>%
        spread_draws(
          (!! rlang::sym(kappa))[!!! rlang::syms(pars.index)],
          (!! rlang::sym(nu))[!!! rlang::syms(pars.index)],
          (!! rlang::sym(m))[!!! rlang::syms(pars.index), cue],
          (!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2],
          lapse_rate,
          ndraws = ndraws,
          seed = seed)
    }

    d.pars %<>%
      # If group is prior, then add the group variable with value "prior" to d.pars first.
      rename_at(vars(ends_with(postfix)), ~ sub(postfix, "", .)) %>%
      { if (summarize) {
        group_by(., !!! syms(pars.index), cue, cue2) %>%
          summarise(
            .,
            across(
              c(kappa, nu, m, S, lapse_rate),
              mean)) %>%
          mutate(.chain = "all", .iteration = "all", .draw = "all")
      } else . } %>%
      { if ("prior" %in% groups) mutate(., group = "prior") else . } %>%
      filter(.data$group %in% .env$groups, .data$category %in% .env$categories) %>%
      relocate(.chain, .iteration, .draw, group, category, kappa, nu, cue, cue2, m, S, lapse_rate)

    # Check that all samples are well-formed
    d.pars %>%
      ungroup() %>%
      summarise(
        across(
          c(kappa, nu, m, S, lapse_rate),
          ~ all(!is.na(.x) & !is.nan(.x) & !is.infinite(.x)))) %>%
      { if (!all(.)) stop2("Some draws are not well-formed (NA, NaN, infinite values). This likely means that there was an issue during the fitting of the stanfit object") }

    if (untransform_cues) {
      d.pars %<>%
        nest_cue_information_in_model() %>%
        untransform_model(transform = fit$transform_information)

      if (!nest) d.pars %<>% unnest_cue_information_in_model()
    } else if (nest) {
      d.pars %<>%
        nest_cue_information_in_model()
    }

    # Clean-up
    d.pars %<>%
      ungroup() %>%
      # Make sure that group and category are factors (even if group or category ids are just numbers)
      mutate(
        category = factor(category, levels = .env$categories),
        group = factor(group, levels = .env$groups))

    if (wide) {
      if ("prior" %in% groups)
        d.pars %<>%
        pivot_longer(
          cols = c(m, S),
          names_to = "variable",
          values_to = "value")
      else
        d.pars %<>%
        pivot_longer(
          cols = c(kappa, nu, m, S),
          names_to = "variable",
          values_to = "value")

      d.pars %<>%
        unite(temp, !!! rlang::syms(pars.index), variable) %>%
        pivot_wider(names_from = temp, values_from = value)
    }

    return(d.pars)
  }
}
