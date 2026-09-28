if (!exists("is_tibble", inherits = TRUE)) {
  suppressPackageStartupMessages(library("tibble", character.only = TRUE, quietly = TRUE))
}
if (!exists(".nlist", inherits = TRUE)) {
  suppressPackageStartupMessages(library("rlang", character.only = TRUE, quietly = TRUE))
}
if (!exists("sampling", inherits = TRUE)) {
  suppressPackageStartupMessages(library("rstan", character.only = TRUE, quietly = TRUE))
}
if (!exists("recover_types", inherits = TRUE)) {
  suppressPackageStartupMessages(library("tidybayes", character.only = TRUE, quietly = TRUE))
}

# Ensure lifecycle deprecation warnings are always emitted during testing
# (rather than throttled to once per session)
options(lifecycle_verbosity = "warning")

r_dir_candidates <- c(
  "R",
  file.path("..", "..", "R"),
  tryCatch(
    testthat::test_path("..", "..", "R"),
    error = function(e) character(0)
  )
)
r_dir_hits <- vapply(
  r_dir_candidates,
  function(path) {
    length(path) == 1L && !is.na(path) &&
      file.exists(file.path(path, "internal-globals.R"))
  },
  logical(1)
)
if (any(r_dir_hits)) {
  r_dir <- r_dir_candidates[r_dir_hits][1]
  pkg_root <- normalizePath(
    file.path(r_dir, ".."),
    winslash = "/",
    mustWork = TRUE
  )
  old_wd <- getwd()
  # Package code comes from the loaded namespace. Re-sourcing R/ into the global
  # environment would shadow the S7 generics with method-less copies.
  # rstan resolves the model files in inst/stan relative to the package root.
  setwd(pkg_root)
  reg.finalizer(environment(), function(e) {
    if (is.character(old_wd) && length(old_wd) == 1L && nzchar(old_wd)) {
      setwd(old_wd)
    }
  }, onexit = TRUE)
}

