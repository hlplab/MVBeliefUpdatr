# S7 registration hook for MVBeliefUpdatr.

# This file should only register methods; the implementation of the methods
# lives in R/S7-core-methods.R.
.onLoad <- function(libname, pkgname) {
  # 1️⃣ Set package-wide defaults (options, env vars)
  #    e.g. options(myPkg.verbose = FALSE)
  options(
    lifecycle_verbosity = getOption("lifecycle_verbosity", "warning"),
    MVBU_deprecation_behavior = getOption("MVBU_deprecation_behavior", "warn")
  )

  # 2️⃣ Register S3 methods that need the namespace early
  #    e.g. S3method("print", "myClass")
  S7::methods_register()

  # 3️⃣ Load compiled code (if any)
  #    e.g. library.dynam("myPkg", pkgname, libname, .registration = TRUE)

  # 4️⃣ Initialize any package-level caches or lazy-load objects
  #    e.g. assign("cache", new.env(parent = emptyenv()), envir = parent.env(environment()))
}


