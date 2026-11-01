# S7 registration hook for MVBeliefUpdatr.

# This file should only register methods; the implementation of the methods
# lives in R/S7-methods.R.
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
  .init_family_registry()

  # 5️⃣ Deprecated dataset active binding
  ns <- asNamespace(pkgname)
  if (exists("ChodroffWilson2018", envir = ns, inherits = FALSE)) {
    rm("ChodroffWilson2018", envir = ns)
  }
  makeActiveBinding(
    "ChodroffWilson2018",
    function() {
      if (!nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
        lifecycle::deprecate_warn(
          "0.1.0",
          I("The 'ChodroffWilson2018' dataset"),
          with = I("the 'mixer6' dataset"),
          details = paste(
            "The legacy ChodroffWilson2018.RData was replaced by",
            "the 'mixer6' dataset."
          ),
          env = ns
        )
      }
      MVBeliefUpdatr::mixer6
    },
    env = ns
  )
}
