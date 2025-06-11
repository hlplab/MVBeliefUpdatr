get_current_versions <- function()
  list(
    MVBeliefUpdatr = utils::packageVersion("MVBeliefUpdatr"),
    rstan = utils::packageVersion("rstan"),
    stanHeaders = utils::packageVersion("StanHeaders"))

# Placeholder functions since we're not using S4 classes
hasSlot <- function(object, slot) {
  slot_names <- slotNames(class(object))

  if (length(slot_names) == 0) {
    slot_names <- names(object)
  }

  slot %in% slot_names
}


# from brms
# get pattern matches in text as vector
# @param simplify return an atomic vector of matches?
# @param first only return the first match in each string?
# @return character vector containing matches
get_matches <- function(pattern, text, simplify = TRUE,
                        first = FALSE, ...) {
  x <- regmatches(text, gregexpr(pattern, text, ...))
  if (first) {
    x <- lapply(x, function(t) if (length(t)) t[1] else t)
  }
  if (simplify) {
    if (first) {
      x <- lapply(x, function(t) if (length(t)) t else "")
    }
    x <- unlist(x)
  }
  x
}

# from brms
# check if x is a try-error resulting from try()
is_try_error <- function(x) {
  inherits(x, "try-error")
}


is_atomic_or_null <- function(x) {
  is.atomic(x) || is.null(x)
}

isNA <- function(x) {
  length(x) == 1L && is.na(x)
}

replace_na_in_array <- function(x, fill = 0) {
  stopifnot(is.array(x))
  stopifnot(is.scalar(fill))

  x[is.na(x)] <- fill
  return(x)
}

is_equal <- function(x, y, check.attributes = FALSE, ...) {
  isTRUE(all.equal(x, y, check.attributes = check.attributes, ...))
}

# extract factor levels from an arbitrary variable
extract_levels <- function(x) {
  # do not check for NAs according to #1355
  if (!is.factor(x)) {
    x <- factor(x)
  }
  levels(x)
}

# check if 'x' will behave like a factor in design matrices
is_like_factor <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

# from brms
# as.factor but allows to pass levels
as_factor <- function(x, levels = NULL) {
  if (is.null(levels)) {
    out <- as.factor(x)
  } else {
    out <- factor(x, levels = levels)
  }
  out
}

# from brms
# coerce 'x' to a single logical value
as_one_logical <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.logical(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse1(s, width.cutoff = 100L)
    stop2("Cannot coerce '", s, "' to a single logical value.")
  }
  x
}

# from brms
# coerce 'x' to a single integer value
as_one_integer <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- SW(as.integer(x))
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse1(s, width.cutoff = 100L)
    stop2("Cannot coerce '", s, "' to a single integer value.")
  }
  x
}

# from brms
# coerce 'x' to a single numeric value
as_one_numeric <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- SW(as.numeric(x))
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse1(s, width.cutoff = 100L)
    stop2("Cannot coerce '", s, "' to a single numeric value.")
  }
  x
}

# from brms
# coerce 'x' to a single character string
as_one_character <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.character(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse1(s, width.cutoff = 100L)
    stop2("Cannot coerce '", s, "' to a single character value.")
  }
  x
}

# from brms
# coerce 'x' to a single character variable name
as_one_variable <- function(x, allow_na = TRUE) {
  x <- as_one_character(x)
  if (x == "NA" && allow_na) {
    return(x)
  }
  if (!nzchar(x) || !is_equal(x, all_vars(x))) {
    stop2("Cannot coerce '", x, "' to a single variable name.")
  }
  x
}

# dim that returns length of vector for vector (not from brms)
dim2 <- function(x) {
  if (is.null(dim(x))) return(length(x))
  return(dim(x))
}


warning2 <- function(...) {
  warning(..., call. = FALSE)
}

stop2 <- function(...) {
  stop(..., call. = FALSE)
}

# From brms
#' Execute a Function Call
#'
#' Execute a function call similar to \code{\link{do.call}}, but without
#' deparsing function arguments. For large number of arguments (i.e., more
#' than a few thousand) this function currently is somewhat inefficient
#' and should be used with care in this case.
#'
#' @param what Either a function or a non-empty character string naming the
#'   function to be called.
#' @param args A list of arguments to the function call. The names attribute of
#'   \code{args} gives the argument names.
#' @param pkg Optional name of the package in which to search for the
#'   function if \code{what} is a character string.
#' @param envir An environment within which to evaluate the call.
#'
#' @return The result of the (evaluated) function call.
#'
#' @keywords internal
#' @export
do_call <- function(what, args, pkg = NULL, envir = parent.frame()) {
  call <- ""
  if (length(args)) {
    if (!is.list(args)) {
      stop2("'args' must be a list.")
    }
    fun_args <- names(args)
    if (is.null(fun_args)) {
      fun_args <- rep("", length(args))
    } else {
      nzc <- nzchar(fun_args)
      fun_args[nzc] <- paste0("`", fun_args[nzc], "` = ")
    }
    names(args) <- paste0(".x", seq_along(args))
    call <- paste0(fun_args, names(args), collapse = ",")
  } else {
    args <- list()
  }
  if (is.function(what)) {
    args$.fun <- what
    what <- ".fun"
  } else {
    what <- paste0("`", as_one_character(what), "`")
    if (!is.null(pkg)) {
      what <- paste0(as_one_character(pkg), "::", what)
    }
  }
  call <- paste0(what, "(", call, ")")
  eval2(call, envir = args, enclos = envir)
}


# like 'eval' but parses characters before evaluation
eval2 <- function(expr, envir = parent.frame(), ...) {
  if (is.character(expr)) {
    expr <- str2expression(expr)
  }
  eval(expr, envir, ...)
}

# evaluate an expression without printing output or messages
# @param expr expression to be evaluated
# @param type type of output to be suppressed (see ?sink)
# @param try wrap evaluation of expr in 'try' and
#   not suppress outputs if evaluation fails?
# @param silent actually evaluate silently?
eval_silent <- function(expr, type = "output", try = FALSE,
                        silent = TRUE, ...) {
  try <- as_one_logical(try)
  silent <- as_one_logical(silent)
  type <- match.arg(type, c("output", "message"))
  expr <- substitute(expr)
  envir <- parent.frame()
  if (silent) {
    if (try && type == "message") {
      try_out <- try(utils::capture.output(
        out <- eval(expr, envir), type = type, ...
      ))
      if (is_try_error(try_out)) {
        # try again without suppressing error messages
        out <- eval(expr, envir)
      }
    } else {
      utils::capture.output(out <- eval(expr, envir), type = type, ...)
    }
  } else {
    out <- eval(expr, envir)
  }
  out
}

# find the name that 'x' had in a specific environment
substitute_name <- function(x, envir = parent.frame(), nchar = 50) {
  out <- substitute(x)
  out <- eval2(paste0("substitute(", out, ")"), envir = envir)
  if (missing(out)) {
    return(NULL)
  }
  substr(collapse(deparse(out)), 1, nchar)
}
