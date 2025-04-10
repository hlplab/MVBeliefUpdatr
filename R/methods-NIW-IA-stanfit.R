#' loo for NIW ideal adaptor stanfit
#'
#' \code{loo} method for \code{\link{NIW_ideal_adaptor_stanfit}} objects.
#'
#' @param x An \code{\link{NIW_ideal_adaptor_stanfit}} object.
#'
#' @return A \code{loo} object.
#'
#' @method loo NIW_ideal_adaptor_stanfit
#'
#' @importFrom loo loo extract_log_lik relative_eff loo.array
#' @export
loo.NIW_ideal_adaptor_stanfit <- function(
    x,
    pars = "log_lik",
    ...,
    save_psis = FALSE,
    cores = getOption("mc.cores", 1)
) {
  stopifnot(length(pars) == 1L)
  stanfit <- get_stanfit(x)
  LLarray <- loo::extract_log_lik(stanfit = stanfit,
                                  parameter_name = pars,
                                  merge_chains = FALSE)
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
  loo::loo.array(LLarray,
                 r_eff = r_eff,
                 cores = cores,
                 save_psis = save_psis)
}


#' Summarize NIW ideal adaptor stanfit
#'
#' \code{summary} method for \code{\link{NIW_ideal_adaptor_stanfit}} objects. Specifies reasonable defaults
#' for the parameters to be summarized for the stanfit object.
#'
#' @param x An \code{\link{NIW_ideal_adaptor_stanfit}} object.
#' @param pars A character vector of parameter names to be summarized. If `NULL`, all parameters of the
#' NIW_ideal_adaptor are summarized. (default: `NULL`)
#' @param prior_only Should only the priors and other sufficient parameters be summarized? (default: `FALSE`)
#' @param include_transformed_pars Should transformed parameters be included in the summary? (default: `FALSE`)
#' @param ... Additional arguments passed to \code{\link[rstan]{summary}}.
#'
#' @method summary NIW_ideal_adaptor_stanfit
#'
#' @importFrom rstan summary
#' @importFrom tibble rownames_to_column
#' @export
summary.NIW_ideal_adaptor_stanfit <- function(x, pars = NULL, prior_only = FALSE, include_transformed_pars = F, ...) {
  stanfit <- get_stanfit(x)
  if (is.null(pars)) {
    pars <- names(stanfit)
    pars <- grep("^((kappa|nu|m|S)_|lapse_rate)", pars, value = T)
    pars <- grep("^m_0_(tau|L_omega|cov)", pars, value = T, invert = T)
    pars <- grep("^((m|S)_0|lapse_rate)_param", pars, value = T, invert = T)
    if (!include_transformed_pars) pars <- grep("^_transformed", pars, value = T, invert = T)
  }

  # Sort and filter output
  full_summary <-
    rstan::summary(stanfit, pars = pars, ...)$summary %>%
    as.data.frame() %>%
    rownames_to_column("Parameter") %>%
    mutate(
      name = factor(gsub("^(.*)_(0|n).*$", "\\1", Parameter), levels = c("kappa", "nu", "m", "S", "lapse_rate")),
      distribution = gsub("^(.*)_(0|n).*$", "\\2", Parameter),
      distribution = factor(ifelse(distribution == Parameter, "0", distribution), levels = c("0", "n")),
      index = gsub("^.*_(0|n)\\[(.*)\\]$", "\\2", Parameter),
      index = ifelse(index == Parameter, 1, index)) %>%
    { if (prior_only) filter(., distribution == "0") else . } %>%
    separate(index, into = c("i1", "i2", "i3", "i4"), sep = ",", fill = "right") %>%
    arrange(distribution, name, i1, i2, i3, i4) %>%
    select(-c(name, distribution, i1, i2, i3, i4))


  Rhats <- full_summary[, "Rhat"]
  if (any(Rhats > 1.05, na.rm = TRUE)) {
    warning2(
      "Parts of the model have not converged (some Rhats are > 1.05). ",
      "Be careful when analysing the results! We recommend running ",
      "more iterations and/or setting stronger priors."
    )
  }
  div_trans <- sum(nuts_params(x, pars = "divergent__")$Value)
  adapt_delta <- control_params(x)$adapt_delta
  if (div_trans > 0) {
    warning2(
      "There were ", div_trans, " divergent transitions after warmup. ",
      "Increasing adapt_delta above ", adapt_delta, " may help. See ",
      "http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"
    )
  }

  return(full_summary)
}

#' #' loo moment matching for NIW ideal adaptor stanfit
#' #'
#' #' \code{loo_moment_match} method for \code{\link{NIW_ideal_adaptor_stanfit}} objects. For detaisl,
#' #' see \code{\link{loo::loo_moment_match}}.
#' #'
#' #' @param x An \code{\link{NIW_ideal_adaptor_stanfit}} object.
#' #' @param loo A `loo` object.
#' #' @param post_draws,log_lik_i,unconstrain_pars,log_prob_upars,log_lik_i_upars User-defined functions.
#' #' @method loo NIW_ideal_adaptor_stanfit
#' #'
#' #' @importFrom loo loo_moment_match.default
#' #' @export
#' # wrapper around loo_moment_match.default from loo package
#' loo_moment_match.stanfit <- function(x, loo = loo, ...) {
#'   loo::loo_moment_match.default(
#'     x = x, loo = loo,
#'     post_draws = rstan:::post_draws_stanfit,
#'     log_lik_i = rstan:::log_lik_i_stanfit,
#'     unconstrain_pars = unconstrain_pars_stanfit,
#'     log_prob_upars = log_prob_upars_stanfit,
#'     log_lik_i_upars = log_lik_i_upars_stanfit,
#'     ...)
#' }
#'
#' # extract original posterior draws
#' post_draws_stanfit <- function(x, ...) {
#'   as.matrix(x)
#' }
#'
#' # compute a matrix of log-likelihood values for the ith observation
#' # matrix contains information about the number of MCMC chains
#' #' @importFrom loo extract_log_lik
#' log_lik_i_stanfit <- function(x, i, parameter_name = "log_lik", ...) {
#'   loo::extract_log_lik(x, parameter_name, merge_chains = FALSE)[, , i]
#' }
#'
#' # transform parameters to the unconstraint space
#' unconstrain_pars_stanfit <- function(x, pars, ...) {
#'   skeleton <- rstan:::.create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
#'   upars <- apply(pars, 1, FUN = function(theta) {
#'     rstan::unconstrain_pars(x, .rstan_relist(theta, skeleton))
#'   })
#'   # for one parameter models
#'   if (is.null(dim(upars))) {
#'     dim(upars) <- c(1, length(upars))
#'   }
#'   t(upars)
#' }
#'
#'
#'
#' # Copied and modified from https://github.com/paul-buerkner/brms/blob/master/R/loo_moment_match.R
#' loo_moment_match.NIW_ideal_adaptor_stanfit <- function(
#'     model,
#'     loo = NULL,
#'     k_threshold = 0.7,
#'     newdata = NULL,
#'     recompile = FALSE,
#'     ...
#' ) {
#'   stopifnot(is.NIW_ideal_adaptor_stanfit(model))
#'   loo <- loo %||% x$criteria[["loo"]]
#'
#'   if (is.null(loo)) {
#'     stop2("No 'loo' object was provided and none is stored within the model.")
#'   } else if (!is.loo(loo)) {
#'     stop2("Inputs to the 'loo' argument must be of class 'loo'.")
#'   }
#'
#'   if (is.null(newdata)) {
#'     newdata <- get_staninput(model)
#'   } else {
#'     newdata <- as.data.frame(newdata)
#'   }
#'
#'   # check <- as_one_logical(check)
#'   # if (check) {
#'   #   yhash_loo <- attr(loo, "yhash")
#'   #   yhash_fit <- hash_response(x, newdata = newdata)
#'   #   if (!is_equal(yhash_loo, yhash_fit)) {
#'   #     stop2(
#'   #       "Response values used in 'loo' and 'x' do not match. ",
#'   #       "If this is a false positive, please set 'check' to FALSE."
#'   #     )
#'   #   }
#'   # }
#'
#'   # otherwise loo_moment_match may fail in a new R session or on another machine
#'   # x <- update_misc_env(x, recompile = recompile)
#'   out <- try(
#'     loo::loo_moment_match.default(
#'       x,
#'       loo = loo,
#'       post_draws = as.matrix,
#'       log_lik_i = .log_lik_i,
#'       unconstrain_pars = .unconstrain_pars,
#'       log_prob_upars = .log_prob_upars,
#'       log_lik_i_upars = .log_lik_i_upars,
#'       k_threshold = k_threshold,
#'       # newdata = newdata,
#'       # resp = resp,
#'       ...))
#'   if (is_try_error(out)) {
#'     stop2(
#'       "Moment matching failed. Perhaps you did not set ",
#'       "'save_pars = save_pars(all = TRUE)' when fitting your model? "
#'       # "If you are running moment matching on another machine than the one ",
#'       # "used to fit the model, you may need to set recompile = TRUE."
#'     )
#'   }
#'   out
#' }
#'
#' # compute a vector of log-likelihood values for the ith observation
#' .log_lik_i <- function(x, i, ...) {
#'   as.vector(log_lik(x, newdata = newdata[i, , drop = FALSE], ...))
#' }
#'
#' # transform parameters to the unconstrained space
#' .unconstrain_pars <- function(x, pars, ...) {
#'   unconstrain_pars_stanfit(x$fit, pars = pars, ...)
#' }
#'
#' # compute log_prob for each posterior draws on the unconstrained space
#' .log_prob_upars <- function(x, upars, ...) {
#'   x <- update_misc_env(x, only_windows = TRUE)
#'   log_prob_upars_stanfit(x$fit, upars = upars, ...)
#' }
#'
#'
#' loo_moment_match.default <- function(x, loo, post_draws, log_lik_i,
#'                                      unconstrain_pars, log_prob_upars,
#'                                      log_lik_i_upars, max_iters = 30L,
#'                                      k_threshold = NULL, split = TRUE,
#'                                      cov = TRUE, cores = getOption("mc.cores", 1),
#'                                      ...) {
#'   assert_that(is.NIW_ideal_adaptor_stanfit(x))
#'   assert_that(class(loo) == "loo")
#'
#'   # input checks
#'   checkmate::assertClass(loo, classes = "loo")
#'   checkmate::assertFunction(post_draws)
#'   checkmate::assertFunction(log_lik_i)
#'   checkmate::assertFunction(unconstrain_pars)
#'   checkmate::assertFunction(log_prob_upars)
#'   checkmate::assertFunction(log_lik_i_upars)
#'   checkmate::assertNumber(max_iters)
#'   checkmate::assertNumber(k_threshold, null.ok=TRUE)
#'   checkmate::assertLogical(split)
#'   checkmate::assertLogical(cov)
#'   checkmate::assertNumber(cores)
#'
#'
#'   if ("psis_loo" %in% class(loo)) {
#'     is_method <- "psis"
#'   } else {
#'     stop("loo_moment_match currently supports only the \"psis\" importance sampling class.")
#'   }
#'
#'
#'   S <- dim(loo)[1]
#'   N <- dim(loo)[2]
#'   if (is.null(k_threshold)) {
#'     k_threshold <- ps_khat_threshold(S)
#'   }
#'   pars <- post_draws(x, ...)
#'   # transform the model parameters to unconstrained space
#'   upars <- unconstrain_pars(x, pars = pars, ...)
#'   # number of parameters in the **parameters** block only
#'   npars <- dim(upars)[2]
#'   # if more parameters than samples, do not do Cholesky transformation
#'   cov <- cov && S >= 10 * npars
#'   # compute log-probabilities of the original parameter values
#'   orig_log_prob <- log_prob_upars(x, upars = upars, ...)
#'
#'   # loop over all observations whose Pareto k is high
#'   ks <- loo$diagnostics$pareto_k
#'   kfs <- rep(0,N)
#'   I <- which(ks > k_threshold)
#'
#'   loo_moment_match_i_fun <- function(i) {
#'     loo_moment_match_i(i = i, x = x, log_lik_i = log_lik_i,
#'                        unconstrain_pars = unconstrain_pars,
#'                        log_prob_upars = log_prob_upars,
#'                        log_lik_i_upars = log_lik_i_upars,
#'                        max_iters = max_iters, k_threshold = k_threshold,
#'                        split = split, cov = cov, N = N, S = S, upars = upars,
#'                        orig_log_prob = orig_log_prob, k = ks[i],
#'                        is_method = is_method, npars = npars, ...)
#'   }
#'
#'   if (cores == 1) {
#'     mm_list <- lapply(X = I, FUN = function(i) loo_moment_match_i_fun(i))
#'   }
#'   else {
#'     if (!os_is_windows()) {
#'       mm_list <- parallel::mclapply(X = I, mc.cores = cores,
#'                                     FUN = function(i) loo_moment_match_i_fun(i))
#'     }
#'     else {
#'       cl <- parallel::makePSOCKcluster(cores)
#'       on.exit(parallel::stopCluster(cl))
#'       mm_list <- parallel::parLapply(cl = cl, X = I,
#'                                      fun = function(i) loo_moment_match_i_fun(i))
#'     }
#'   }
#'
#'   # update results
#'   for (ii in seq_along(I)) {
#'     i <- mm_list[[ii]]$i
#'     loo$pointwise[i, "elpd_loo"] <- mm_list[[ii]]$elpd_loo_i
#'     loo$pointwise[i, "p_loo"] <- mm_list[[ii]]$p_loo
#'     loo$pointwise[i, "mcse_elpd_loo"] <- mm_list[[ii]]$mcse_elpd_loo
#'     loo$pointwise[i, "looic"] <- mm_list[[ii]]$looic
#'
#'     loo$diagnostics$pareto_k[i] <- mm_list[[ii]]$k
#'     loo$diagnostics$n_eff[i] <- mm_list[[ii]]$n_eff
#'     kfs[i] <- mm_list[[ii]]$kf
#'
#'     if (!is.null(loo$psis_object)) {
#'       loo$psis_object$log_weights[, i] <- mm_list[[ii]]$lwi
#'     }
#'   }
#'   if (!is.null(loo$psis_object)) {
#'     attr(loo$psis_object, "norm_const_log") <- matrixStats::colLogSumExps(loo$psis_object$log_weights)
#'     loo$psis_object$diagnostics <- loo$diagnostics
#'   }
#'
#'   # combined estimates
#'   cols_to_summarize <- !(colnames(loo$pointwise) %in% c("mcse_elpd_loo",
#'                                                         "influence_pareto_k"))
#'   loo$estimates <- table_of_estimates(loo$pointwise[, cols_to_summarize,
#'                                                     drop = FALSE])
#'
#'   # these will be deprecated at some point
#'   loo$elpd_loo <- loo$estimates["elpd_loo","Estimate"]
#'   loo$p_loo <- loo$estimates["p_loo","Estimate"]
#'   loo$looic <- loo$estimates["looic","Estimate"]
#'   loo$se_elpd_loo <- loo$estimates["elpd_loo","SE"]
#'   loo$se_p_loo <- loo$estimates["p_loo","SE"]
#'   loo$se_looic <- loo$estimates["looic","SE"]
#'
#'   # Warn if some Pareto ks are still high
#'   throw_pareto_warnings(loo$diagnostics$pareto_k, k_threshold)
#'   # if we don't split, accuracy may be compromised
#'   if (!split) {
#'     throw_large_kf_warning(kfs, k_threshold)
#'   }
#'
#'   loo
#' }
