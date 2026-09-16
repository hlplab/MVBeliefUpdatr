# Copied and modified from brms

#' Compile a Stan model with the requested backend.
#'
#' @param model Stan model code.
#' @param backend Backend name to use for compilation.
#' @param ... Additional arguments passed to the backend-specific compiler.
#' @return A compiled Stan model object.
#' @keywords internal
.compile_model <- function(model, backend, ...) {
  backend <- as_one_character(backend)
  .compile_model <- get(paste0(".compile_model_", backend), mode = "function")
  .compile_model(model, ...)
}

#' Compile a Stan model with rstan.
#'
#' @param model Stan model code.
#' @param threads Number of threads to use for compilation.
#' @param silent Verbosity level for compilation output.
#' @param ... Additional arguments passed to rstan::stan_model().
#' @return A compiled Stan model object from rstan.
#' @keywords internal
.compile_model_rstan <- function(
  model,
  threads,
  # opencl,
  silent = 1,
  ...
) {
  args <- list(...)
  args$model_code <- model
  if (silent < 2) {
    message("Compiling Stan program...")
  }
  # if (use_threading(threads, force = TRUE)) {
  #   if (utils::packageVersion("rstan") >= "2.26") {
  #     threads_per_chain_def <- rstan::rstan_options("threads_per_chain")
  #     on.exit(rstan::rstan_options(threads_per_chain = threads_per_chain_def))
  #     rstan::rstan_options(threads_per_chain = threads$threads)
  #   } else {
  #     .stop("Threading is not supported by backend 'rstan' version ",
  #           utils::packageVersion("rstan"), ".")
  #   }
  # }
  # if (use_opencl(opencl)) {
  #   .stop("OpenCL is not supported by backend 'rstan' version ",
  #         utils::packageVersion("rstan"), ".")
  # }
  eval_silent(
    do_call(rstan::stan_model, args),
    type = "message", try = TRUE, silent = silent >= 2
  )
}

#' Compile a Stan model with cmdstanr.
#'
#' @param model Stan model code.
#' @param threads Number of threads to use for compilation.
#' @param silent Verbosity level for compilation output.
#' @param ... Additional arguments passed to cmdstanr::cmdstan_model().
#' @return A compiled Stan model object from cmdstanr.
#' @keywords internal
.compile_model_cmdstanr <- function(
  model,
  threads,
  # opencl,
  silent = 1,
  ...
) {
  require_package("cmdstanr")
  args <- list(...)
  args$stan_file <- cmdstanr::write_stan_file(model)
  # if (cmdstanr::cmdstan_version() >= "2.29.0") {
  #   .canonicalize_stan_model(args$stan_file, overwrite_file = TRUE)
  # }
  # if (use_threading(threads, force = TRUE)) {
  #   args$cpp_options$stan_threads <- TRUE
  # }
  # if (use_opencl(opencl)) {
  #   args$cpp_options$stan_opencl <- TRUE
  # }
  eval_silent(
    do_call(cmdstanr::cmdstan_model, args),
    type = "message", try = TRUE, silent = silent >= 2
  )
}

#' Normalize Stan code to avoid unnecessary recompilation after whitespace
#' or comment changes.
#'
#' @param x A string containing the Stan code.
#' @return A normalized character string.
#' @keywords internal
.normalize_stancode <- function(x) {
  x <- as_one_character(x)
  # Remove single-line comments
  x <- gsub("//[^\n\r]*[\n\r]", " ", x)
  x <- gsub("//[^\n\r]*$", " ", x)
  # Remove multi-line comments
  x <- gsub("/\\*([^*]*(\\*[^/])?)*\\*/", " ", x)
  # Standardize whitespace (including newlines)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}
