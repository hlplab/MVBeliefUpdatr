# Copied and modified from brms

# compile Stan model
# @param model Stan model code
# @return validated Stan model code
compile_model <- function(model, backend, ...) {
  backend <- as_one_character(backend)
  .compile_model <- get(paste0(".compile_model_", backend), mode = "function")
  .compile_model(model, ...)
}

# compile Stan model with rstan
# @param model Stan model code
# @return model compiled with rstan
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
  #     stop2("Threading is not supported by backend 'rstan' version ",
  #           utils::packageVersion("rstan"), ".")
  #   }
  # }
  # if (use_opencl(opencl)) {
  #   stop2("OpenCL is not supported by backend 'rstan' version ",
  #         utils::packageVersion("rstan"), ".")
  # }
  eval_silent(
    do_call(rstan::stan_model, args),
    type = "message", try = TRUE, silent = silent >= 2
  )
}

# compile Stan model with cmdstanr
# @param model Stan model code
# @return model compiled with cmdstanr
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

# Normalizes Stan code to avoid triggering refit after whitespace and
# comment changes in the generated code.
# In some distant future, StanC3 may provide its own normalizing functions,
# until then this is a set of regex hacks.
# @param x a string containing the Stan code
normalize_stancode <- function(x) {
  x <- as_one_character(x)
  # Remove single-line comments
  x <- gsub("//[^\n\r]*[\n\r]", " ", x)
  x <- gsub("//[^\n\r]*$", " ", x)
  # Remove multi-line comments
  x <- gsub("/\\*([^*]*(\\*[^/])?)*\\*/", " ", x)
  # Standardize whitespace (including newlines)
  x <- gsub("[[:space:]]+"," ", x)
  trimws(x)
}
