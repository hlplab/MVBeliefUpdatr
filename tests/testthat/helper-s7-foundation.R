if (!requireNamespace("assertthat", quietly = TRUE)) {
  stop("The assertthat package is required for these S7 tests.")
}
if (!exists(".assert_that", inherits = TRUE)) {
  suppressPackageStartupMessages(library("assertthat", character.only = TRUE, quietly = TRUE))
}
if (!exists("is_tibble", inherits = TRUE)) {
  suppressPackageStartupMessages(library("tibble", character.only = TRUE, quietly = TRUE))
}
if (!exists("nlist", inherits = TRUE)) {
  suppressPackageStartupMessages(library("rlang", character.only = TRUE, quietly = TRUE))
}
if (!exists("sampling", inherits = TRUE)) {
  suppressPackageStartupMessages(library("rstan", character.only = TRUE, quietly = TRUE))
}

source_with_env <- function(path) {
  sys.source(path, envir = globalenv())
}

r_dir_candidates <- c("R", file.path("..", "..", "R"))
r_dir <- r_dir_candidates[vapply(r_dir_candidates, function(path) file.exists(file.path(path, "globals.R")), logical(1))][1]
if (is.na(r_dir)) {
  stop("Could not locate the package R directory for test helper sourcing.")
}

pkg_root <- normalizePath(file.path(r_dir, ".."), winslash = "/", mustWork = TRUE)
r_dir_abs <- normalizePath(file.path(pkg_root, "R"), winslash = "/", mustWork = TRUE)
old_wd <- getwd()

if (!exists("MVBU_PROB_TOL", inherits = TRUE)) {
  source_with_env(file.path(r_dir_abs, "globals.R"))
}
source_with_env(file.path(r_dir_abs, "asserts.R"))
source_with_env(file.path(r_dir_abs, "utils.R"))
source_with_env(file.path(r_dir_abs, "misc_imported.R"))
source_with_env(file.path(r_dir_abs, "to-array.R"))
source_with_env(file.path(r_dir_abs, "basics.R"))
source_with_env(file.path(r_dir_abs, "MVBeliefUpdatr-package.R"))
source_with_env(file.path(r_dir_abs, "get-info-from-NIW-IA-stanfit.R"))
source_with_env(file.path(r_dir_abs, "S7-core-classes.R"))
source_with_env(file.path(r_dir_abs, "S7-generics.R"))
source_with_env(file.path(r_dir_abs, "S7-transform-information.R"))
source_with_env(file.path(r_dir_abs, "S7-staninput.R"))
source_with_env(file.path(r_dir_abs, "S7-stanfit-input.R"))
source_with_env(file.path(r_dir_abs, "deprecated-make-staninput.R"))
source_with_env(file.path(r_dir_abs, "S7-stanfit.R"))
source_with_env(file.path(r_dir_abs, "S7-core-methods.R"))
source_with_env(file.path(r_dir_abs, "S7-stanfit-methods.R"))
source_with_env(file.path(r_dir_abs, "S7-stanfit-input-methods.R"))

# Source Stan model definitions from the package root so rstan can locate the
# model files in inst/stan when this helper is loaded from the test directory.
setwd(pkg_root)
on.exit <- function() {
  if (is.character(old_wd) && length(old_wd) == 1L && nzchar(old_wd)) {
    setwd(old_wd)
  }
}
reg.finalizer(environment(), function(e) on.exit(), onexit = TRUE)

source_with_env(file.path(r_dir_abs, "stanmodels.R"))
source_with_env(file.path(r_dir_abs, "fit-IA-stanfit.R"))
source_with_env(file.path(r_dir_abs, "S7-phase2-migration.R"))
