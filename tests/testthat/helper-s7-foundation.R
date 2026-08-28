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

r_dir_candidates <- c(
  "R",
  file.path("..", "..", "R"),
  testthat::test_path("..", "..", "R")
)
r_dir_hits <- vapply(r_dir_candidates, function(path) file.exists(file.path(path, "internal-globals.R")), logical(1))
r_dir <- r_dir_candidates[r_dir_hits]
if (length(r_dir) == 0L || all(is.na(r_dir))) {
  stop("Could not locate the package R directory for test helper sourcing.")
}
r_dir <- r_dir[1]

pkg_root <- normalizePath(file.path(r_dir, ".."), winslash = "/", mustWork = TRUE)
pkg_root <- pkg_root[1]
r_dir_abs <- normalizePath(file.path(pkg_root, "R"), winslash = "/", mustWork = TRUE)
old_wd <- getwd()

# Package code comes from the loaded namespace. Re-sourcing R/ into the global
# environment would shadow the S7 generics with method-less copies, breaking dispatch
# for any code whose enclosing environment is the global environment.

# rstan resolves the model files in inst/stan relative to the package root.
setwd(pkg_root)
on.exit <- function() {
  if (is.character(old_wd) && length(old_wd) == 1L && nzchar(old_wd)) {
    setwd(old_wd)
  }
}
reg.finalizer(environment(), function(e) on.exit(), onexit = TRUE)

