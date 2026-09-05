#' @include S7-core-classes.R
#' @include S7-generics.R
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr mutate filter select transmute arrange slice rename pull relocate
#' @importFrom lifecycle deprecate_warn
NULL

# Deprecation warning helper
.warn_legacy_tibble_conversion <- function(what = "as_tibble()") {
  lifecycle::deprecate_warn(
    when = "0.0.9",
    what = what,
    details = paste0(
      "Converting an S7 model/template to a legacy tibble. Note that legacy tibbles ",
      "do not capture all information from the new S7 objects (such as decision_rule, ",
      "noise_treatment, lapse_treatment, and structured validation). ",
      "Please see vignette(\"s7-class-structure-and-workflows\") for the modern S7 workflow."
    )
  )
}

# -------------------------------------------------------------------------
# as_tibble methods for S7 objects
# -------------------------------------------------------------------------

#' Convert S7 objects to legacy tibble format
#'
#' Converts modern S7 cognitive models, category representation templates, or
#' individual category representations into legacy tibbles matching the schemas
#' from pre-S7 versions of \pkg{MVBeliefUpdatr}.
#'
#' @param x An S7 object inheriting from \code{\link{MVBU_Object}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\link[tibble]{tibble}} with columns and classes adhering to the legacy format.
#'
#' @name as_tibble.MVBU_Object
#' @rdname as_tibble.MVBU_Object
#' @exportS3Method tibble::as_tibble
as_tibble.MVBU_Object <- function(x, ...) {
  .warn_legacy_tibble_conversion("as_tibble()")
  .stop("as_tibble() is not implemented for class: ", paste(class(x), collapse = ", "))
}

#' @rdname as_tibble.MVBU_Object
#' @exportS3Method tibble::as_tibble
as_tibble.MVBU_CategoryRepresentation <- function(x, ...) {
  .warn_legacy_tibble_conversion("as_tibble()")
  # Wrap single representation into a 1-category template and convert
  cat_name <- get_category_labels(x)
  if (length(cat_name) == 0L || is.na(cat_name[1])) cat_name <- "category1"
  rep_list <- list(x)
  names(rep_list) <- cat_name[1]
  tpl <- new_category_representation_template(rep_list)
  as_tibble_from_template(tpl)
}

#' @rdname as_tibble.MVBU_Object
#' @exportS3Method tibble::as_tibble
as_tibble.MVBU_CategoryRepresentationTemplate <- function(x, ...) {
  .warn_legacy_tibble_conversion("as_tibble()")
  as_tibble_from_template(x)
}

#' @rdname as_tibble.MVBU_Object
#' @exportS3Method tibble::as_tibble
as_tibble.MVBU_CognitiveModel <- function(x, ...) {
  .warn_legacy_tibble_conversion("as_tibble()")
  as_tibble_from_model(x)
}

#' @exportS3Method base::as.data.frame
as.data.frame.MVBU_Object <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(tibble::as_tibble(x, ...), row.names = row.names, optional = optional, ...)
}

# -------------------------------------------------------------------------
# Internal conversion engines
# -------------------------------------------------------------------------

as_tibble_from_template <- function(tpl) {
  reps <- tpl@representations
  if (length(reps) == 0L) {
    df <- tibble::tibble(category = factor())
    class(df) <- c("MVBU_representation", class(df))
    return(df)
  }

  first_rep <- reps[[1]]
  cat_names <- names(reps)
  if (is.null(cat_names) || any(cat_names == "")) {
    cat_names <- vapply(reps, function(r) {
      lbls <- get_category_labels(r)
      if (length(lbls) > 0L) lbls[1] else "category"
    }, character(1))
  }
  cat_factor <- factor(cat_names, levels = cat_names)

  if (S7::S7_inherits(first_rep, MVG_CategoryRepresentation)) {
    mu_list <- lapply(reps, function(r) {
      v <- as.numeric(r@mu)
      names(v) <- get_cue_labels(r)
      v
    })
    Sigma_list <- lapply(reps, function(r) {
      m <- as.matrix(r@Sigma)
      cues <- get_cue_labels(r)
      dimnames(m) <- list(cues, cues)
      m
    })
    df <- tibble::tibble(
      category = cat_factor,
      mu = mu_list,
      Sigma = Sigma_list
    )
    class(df) <- c("MVG", "MVBU_representation", class(df))
    return(df)
  }

  if (S7::S7_inherits(first_rep, NIW_CategoryRepresentation)) {
    m_list <- lapply(reps, function(r) {
      v <- as.numeric(r@m)
      names(v) <- get_cue_labels(r)
      v
    })
    S_list <- lapply(reps, function(r) {
      m <- as.matrix(r@S)
      cues <- get_cue_labels(r)
      dimnames(m) <- list(cues, cues)
      m
    })
    kappa_vec <- vapply(reps, function(r) as.numeric(r@kappa), numeric(1))
    nu_vec <- vapply(reps, function(r) as.numeric(r@nu), numeric(1))
    df <- tibble::tibble(
      category = cat_factor,
      m = m_list,
      S = S_list,
      kappa = kappa_vec,
      nu = nu_vec
    )
    class(df) <- c("NIW_belief", "MVBU_representation", class(df))
    return(df)
  }

  if (S7::S7_inherits(first_rep, Exemplar_CategoryRepresentation)) {
    ex_list <- lapply(reps, function(r) {
      m <- as.matrix(r@exemplars)
      colnames(m) <- get_cue_labels(r)
      m
    })
    sim_fn <- function(x, y, weights = rep(1, length(x)), distance_metric = 2, distance_decay_factor = 1) {
      distance <- sum(weights * abs(x - y)^distance_metric)^(1 / distance_metric)
      exp(-distance^distance_decay_factor)
    }
    sim_list <- replicate(length(cat_names), sim_fn, simplify = FALSE)
    df <- tibble::tibble(
      category = cat_factor,
      exemplars = ex_list,
      sim_function = sim_list
    )
    class(df) <- c("exemplars", "MVBU_representation", class(df))
    return(df)
  }

  if (S7::S7_inherits(first_rep, UVG_CategoryRepresentation)) {
    mu_vec <- vapply(reps, function(r) as.numeric(r@mu), numeric(1))
    sigma_vec <- vapply(reps, function(r) sqrt(as.numeric(r@sigma2)), numeric(1))
    df <- tibble::tibble(
      category = cat_factor,
      mu = mu_vec,
      sigma = sigma_vec
    )
    class(df) <- c("UVG", "MVBU_representation", class(df))
    return(df)
  }

  if (S7::S7_inherits(first_rep, NIX_CategoryRepresentation)) {
    m_vec <- vapply(reps, function(r) as.numeric(r@m), numeric(1))
    s2_vec <- vapply(reps, function(r) as.numeric(r@sigma2), numeric(1))
    kappa_vec <- vapply(reps, function(r) as.numeric(r@kappa), numeric(1))
    nu_vec <- vapply(reps, function(r) as.numeric(r@nu), numeric(1))
    df <- tibble::tibble(
      category = cat_factor,
      m = m_vec,
      S = s2_vec,
      kappa = kappa_vec,
      nu = nu_vec
    )
    class(df) <- c("NIX", "MVBU_representation", class(df))
    return(df)
  }

  if (S7::S7_inherits(first_rep, MUVG_CategoryRepresentation)) {
    mu_list <- lapply(reps, function(r) as.numeric(r@mu))
    sigma_list <- lapply(reps, function(r) as.numeric(r@sigma))
    df <- tibble::tibble(
      category = cat_factor,
      mu = mu_list,
      sigma = sigma_list
    )
    class(df) <- c("MUVG", "MVBU_representation", class(df))
    return(df)
  }

  if (S7::S7_inherits(first_rep, MNIX_CategoryRepresentation)) {
    m_list <- lapply(reps, function(r) as.numeric(r@m))
    S_list <- lapply(reps, function(r) as.numeric(r@S))
    kappa_list <- lapply(reps, function(r) as.numeric(r@kappa))
    nu_list <- lapply(reps, function(r) as.numeric(r@nu))
    df <- tibble::tibble(
      category = cat_factor,
      m = m_list,
      S = S_list,
      kappa = kappa_list,
      nu = nu_list
    )
    class(df) <- c("MNIX", "MVBU_representation", class(df))
    return(df)
  }

  # Fallback generic representation
  df <- tibble::tibble(category = cat_factor)
  class(df) <- c("MVBU_representation", class(df))
  df
}

as_tibble_from_model <- function(model) {
  df <- as_tibble_from_template(model@category_template)
  cat_names <- levels(df$category)
  n_cat <- length(cat_names)

  # Prior
  prior <- get_category_prior(model)
  if (is.null(prior) || length(prior) == 0L) {
    prior_vals <- rep(1 / n_cat, n_cat)
  } else if (!is.null(names(prior))) {
    prior_vals <- as.numeric(prior[cat_names])
    if (any(is.na(prior_vals))) {
      prior_vals <- as.numeric(prior)
    }
  } else {
    prior_vals <- as.numeric(prior)
  }
  if (length(prior_vals) != n_cat) {
    prior_vals <- rep(prior_vals, length.out = n_cat)
  }
  names(prior_vals) <- cat_names

  # Lapse rate
  lr <- get_lapse_rate(model)
  if (is.null(lr) || length(lr) == 0L || is.na(lr)) {
    lr <- 0
  }
  lapse_rate_vals <- rep(as.numeric(lr)[1], n_cat)

  # Lapse bias
  lb <- get_lapse_bias(model)
  if (is.null(lb) || length(lb) == 0L) {
    lb_vals <- prior_vals
  } else if (!is.null(names(lb))) {
    lb_vals <- as.numeric(lb[cat_names])
    if (any(is.na(lb_vals))) {
      lb_vals <- as.numeric(lb)
    }
  } else {
    lb_vals <- as.numeric(lb)
  }
  if (length(lb_vals) != n_cat) {
    lb_vals <- rep(lb_vals, length.out = n_cat)
  }
  names(lb_vals) <- cat_names

  # Sigma_noise
  sn <- get_noise(model)
  cue_names <- get_cue_labels(model)
  if (!is.null(sn) && is.matrix(sn)) {
    dimnames(sn) <- list(cue_names, cue_names)
    noise_list <- replicate(n_cat, sn, simplify = FALSE)
  } else {
    noise_list <- vector("list", n_cat)
  }

  df$prior <- prior_vals
  df$lapse_rate <- lapse_rate_vals
  df$lapse_bias <- lb_vals
  df$Sigma_noise <- noise_list

  # Model class tagging
  if (S7::S7_inherits(model, MVG_IdealObserver)) {
    class(df) <- c("MVG_ideal_observer", "MVBU_model", class(df))
  } else if (S7::S7_inherits(model, NIW_IdealAdaptor)) {
    class(df) <- c("NIW_ideal_adaptor", "MVBU_model", class(df))
  } else if (S7::S7_inherits(model, Exemplar_Model)) {
    class(df) <- c("exemplar_model", "MVBU_model", class(df))
  } else if (S7::S7_inherits(model, UVG_IdealObserver)) {
    class(df) <- c("UVG_ideal_observer", "MVBU_model", class(df))
  } else if (S7::S7_inherits(model, NIX_IdealAdaptor)) {
    class(df) <- c("NIX_ideal_adaptor", "MVBU_model", class(df))
  } else {
    class(df) <- c("MVBU_model", class(df))
  }

  df
}

# -------------------------------------------------------------------------
# Deprecated dplyr methods forwarding to as_tibble()
# -------------------------------------------------------------------------

#' @exportS3Method dplyr::mutate
mutate.MVBU_Object <- function(.data, ...) {
  dplyr::mutate(tibble::as_tibble(.data), ...)
}

#' @exportS3Method dplyr::transmute
transmute.MVBU_Object <- function(.data, ...) {
  dplyr::transmute(tibble::as_tibble(.data), ...)
}

#' @exportS3Method dplyr::filter
filter.MVBU_Object <- function(.data, ...) {
  dplyr::filter(tibble::as_tibble(.data), ...)
}

#' @exportS3Method dplyr::select
select.MVBU_Object <- function(.data, ...) {
  dplyr::select(tibble::as_tibble(.data), ...)
}

#' @exportS3Method dplyr::arrange
arrange.MVBU_Object <- function(.data, ...) {
  dplyr::arrange(tibble::as_tibble(.data), ...)
}

#' @exportS3Method dplyr::slice
slice.MVBU_Object <- function(.data, ...) {
  dplyr::slice(tibble::as_tibble(.data), ...)
}

#' @exportS3Method dplyr::rename
rename.MVBU_Object <- function(.data, ...) {
  dplyr::rename(tibble::as_tibble(.data), ...)
}

#' @exportS3Method dplyr::pull
pull.MVBU_Object <- function(.data, ...) {
  dplyr::pull(tibble::as_tibble(.data), ...)
}

#' @exportS3Method dplyr::relocate
relocate.MVBU_Object <- function(.data, ...) {
  dplyr::relocate(tibble::as_tibble(.data), ...)
}
