# Agents MUST ignore everything below this line

+ many of the internal functions use tryCatch calls with function that should be sufficient to handle any errors themselves. Isn't that unnecessary? e.g.

format_fixed_parameters_table <- function(fixed_params, object = NULL) {
  if (is.null(fixed_params) || length(fixed_params) == 0L) {
    return(tibble::tibble())
  }
  if (is.data.frame(fixed_params)) {
    return(tibble::as_tibble(fixed_params))
  }

  cue_labels <- tryCatch(get_cue_labels(object), error = function(e) NULL)
  category_labels <- tryCatch(
    get_category_labels(object),
    error = function(e) NULL
  )
  fam <- tryCatch(
    sub("_.*$", "", get_model_family(object)),
    error = function(e) ""
  )


# To do after Phase 9

+ check docs/to-do/possible-removal-of-typed-category-representations.md for further simplifications/standardization of class structure.

+ think about whether evaluate_model should be renamed. it's job is similar to the add information criteria (IC) functions in brms, just that it's frequentist log-likelihood or accuracy vs. Bayesian ICs. perhaps this could be unified by changing the evaluate_models function into several separate add_* functions that add the IC to the model (similar to brms), allowing also a common compare models function that would extend that compare_models function from non-stanfit to stanfit model (allowing comparison on any IC)?

+ also check whether evaluate_model applied to stanfit or stanfit posterior objects , which currently uses method = "loglik", makes sense 
