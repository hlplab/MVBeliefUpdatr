# Extract diagnostic information from MVBU_Stanfit object

Extract quantities that can be used to diagnose sampling behavior of the
algorithms applied by Stan at the back-end of MVBeliefUpdatr. These
diagnostic functions are all copied and modified from brms.

## Usage

``` r
log_posterior(x, ...)

nuts_params(x, pars = NULL, ...)

rhat(x, pars = NULL, ...)

neff_ratio(x, pars = NULL, ...)

## S7 method for class <MVBU_Stanfit>
log_posterior(x, ...)

## S7 method for class <MVBU_Stanfit>
nuts_params(x, pars = NULL, ...)

## S7 method for class <MVBU_Stanfit>
rhat(x, pars = NULL, ...)

## S7 method for class <MVBU_Stanfit>
neff_ratio(x, pars = NULL, ...)

# S3 method for class 'MVBU_Stanfit'
as_draws(x, ...)

# S3 method for class 'MVBU_Stanfit'
as_draws_df(x, ...)

# S3 method for class 'MVBU_Stanfit'
as_draws_array(x, ...)

# S3 method for class 'MVBU_Stanfit'
as_draws_matrix(x, ...)

# S3 method for class 'MVBU_Stanfit'
as_draws_list(x, ...)

# S3 method for class 'MVBU_Stanfit'
as_draws_rvars(x, ...)
```

## Arguments

- x:

  An `MVBeliefUpdatr` object.

- ...:

  Arguments passed to individual methods.

- pars:

  An optional character vector of parameter names. For `nuts_params`
  these will be NUTS sampler parameter names rather than model
  parameters. If pars is omitted all parameters are included.

## Value

Log posterior draws, NUTS sampler parameters, Rhat, or Neff ratio.

The exact form of the output depends on the method.

## Details

For more details see
[bayesplot::bayesplot-extractors](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html).

## Examples

``` r
# \donttest{
fit <- example_ideal_adaptor_stanfit()
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
rhat(fit)
#>            kappa_0_shared_param[1]               nu_0_shared_param[1] 
#>                          0.9903867                          0.9896015 
#>          category_prior_param[1,1]          category_prior_param[1,2] 
#>                          0.9610481                          0.9610481 
#>              lapse_bias_param[1,1]              lapse_bias_param[1,2] 
#>                          0.9882127                          0.9882127 
#>             m_0_unknown_param[1,1]             m_0_unknown_param[2,1] 
#>                          1.1495576                          1.0317271 
#>             m_0_tau_transformed[1]       m_0_L_omega_transformed[1,1] 
#>                          1.0920547                                 NA 
#>       tau_0_param_transformed[1,1]       tau_0_param_transformed[2,1] 
#>                          0.9993326                          0.9606640 
#> L_omega_0_param_transformed[1,1,1] L_omega_0_param_transformed[2,1,1] 
#>                                 NA                                 NA 
#>                  category_prior[1]                  category_prior[2] 
#>                          0.9610481                          0.9610481 
#>                         lapse_rate                      lapse_bias[1] 
#>                          1.0898848                          0.9882127 
#>                      lapse_bias[2]                         kappa_0[1] 
#>                          0.9882127                          0.9903867 
#>                         kappa_0[2]                            nu_0[1] 
#>                          0.9903867                          0.9896015 
#>                            nu_0[2]               m_0_transformed[1,1] 
#>                          0.9896015                          1.1495576 
#>               m_0_transformed[2,1]             S_0_transformed[1,1,1] 
#>                          1.0317271                          1.0145663 
#>             S_0_transformed[2,1,1]                       kappa_n[1,1] 
#>                          0.9606640                          0.9903867 
#>                       kappa_n[2,1]                       kappa_n[1,2] 
#>                          0.9903867                          0.9903867 
#>                       kappa_n[2,2]                       kappa_n[1,3] 
#>                          0.9903867                          0.9903867 
#>                       kappa_n[2,3]                          nu_n[1,1] 
#>                          0.9903867                          0.9896015 
#>                          nu_n[2,1]                          nu_n[1,2] 
#>                          0.9896015                          0.9896015 
#>                          nu_n[2,2]                          nu_n[1,3] 
#>                          0.9896015                          0.9896015 
#>                          nu_n[2,3]             m_n_transformed[1,1,1] 
#>                          0.9896015                          1.1495576 
#>             m_n_transformed[2,1,1]             m_n_transformed[1,2,1] 
#>                          1.0317271                          1.0204904 
#>             m_n_transformed[2,2,1]             m_n_transformed[1,3,1] 
#>                          1.0103325                          0.9974461 
#>             m_n_transformed[2,3,1]           S_n_transformed[1,1,1,1] 
#>                          1.0295546                          1.0145663 
#>           S_n_transformed[2,1,1,1]           S_n_transformed[1,2,1,1] 
#>                          0.9606640                          1.0086900 
#>           S_n_transformed[2,2,1,1]           S_n_transformed[1,3,1,1] 
#>                          0.9609168                          1.0194280 
#>           S_n_transformed[2,3,1,1]       t_scale_transformed[1,1,1,1] 
#>                          0.9597241                          0.9633621 
#>       t_scale_transformed[2,1,1,1]       t_scale_transformed[1,2,1,1] 
#>                          0.9625473                          0.9754434 
#>       t_scale_transformed[2,2,1,1]       t_scale_transformed[1,3,1,1] 
#>                          0.9633092                          0.9863899 
#>       t_scale_transformed[2,3,1,1]                        log_lik_sum 
#>                          0.9756763                          1.0169554 
#>                         log_lik[1]                         log_lik[2] 
#>                          0.9627076                          0.9719042 
#>                         log_lik[3]                         log_lik[4] 
#>                          0.9902845                          0.9640097 
#>                         log_lik[5]                         log_lik[6] 
#>                          1.0332323                          0.9817348 
#>                         log_lik[7]                         log_lik[8] 
#>                          0.9845119                          0.9687474 
#>                         log_lik[9]                        log_lik[10] 
#>                          0.9634242                          1.0108312 
#>                        log_lik[11]                        log_lik[12] 
#>                          0.9735873                          1.0513459 
#>                        log_lik[13]                        log_lik[14] 
#>                          0.9968547                          1.1269671 
#>                        log_lik[15]                           m_0[1,1] 
#>                          1.0111466                          1.1495576 
#>                           m_0[2,1]                         S_0[1,1,1] 
#>                          1.0317271                          1.0145663 
#>                         S_0[2,1,1]                         m_n[1,1,1] 
#>                          0.9606640                          1.1495576 
#>                         m_n[2,1,1]                         m_n[1,2,1] 
#>                          1.0317271                          1.0204904 
#>                         m_n[2,2,1]                         m_n[1,3,1] 
#>                          1.0103325                          0.9974461 
#>                         m_n[2,3,1]                       S_n[1,1,1,1] 
#>                          1.0295546                          1.0145663 
#>                       S_n[2,1,1,1]                       S_n[1,2,1,1] 
#>                          0.9606640                          1.0086900 
#>                       S_n[2,2,1,1]                       S_n[1,3,1,1] 
#>                          0.9609168                          1.0194280 
#>                       S_n[2,3,1,1]                               lp__ 
#>                          0.9597241                          1.2986253 
neff_ratio(fit)
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#>            kappa_0_shared_param[1]               nu_0_shared_param[1] 
#>                          0.5272677                          0.7065037 
#>          category_prior_param[1,1]          category_prior_param[1,2] 
#>                          0.4123210                          0.4123210 
#>              lapse_bias_param[1,1]              lapse_bias_param[1,2] 
#>                          0.5272677                          0.5272677 
#>             m_0_unknown_param[1,1]             m_0_unknown_param[2,1] 
#>                          0.4434381                          1.3250028 
#>             m_0_tau_transformed[1]       m_0_L_omega_transformed[1,1] 
#>                          0.6211765                                 NA 
#>       tau_0_param_transformed[1,1]       tau_0_param_transformed[2,1] 
#>                          0.4626955                          0.5750953 
#> L_omega_0_param_transformed[1,1,1] L_omega_0_param_transformed[2,1,1] 
#>                                 NA                                 NA 
#>                  category_prior[1]                  category_prior[2] 
#>                          0.4123210                          0.4123210 
#>                         lapse_rate                      lapse_bias[1] 
#>                          0.7227376                          0.5272677 
#>                      lapse_bias[2]                         kappa_0[1] 
#>                          0.5272677                          0.5272677 
#>                         kappa_0[2]                            nu_0[1] 
#>                          0.5272677                          0.7065037 
#>                            nu_0[2]               m_0_transformed[1,1] 
#>                          0.7065037                          0.4434381 
#>               m_0_transformed[2,1]             S_0_transformed[1,1,1] 
#>                          1.3250028                          0.4626955 
#>             S_0_transformed[2,1,1]                       kappa_n[1,1] 
#>                          0.5750953                          0.5272677 
#>                       kappa_n[2,1]                       kappa_n[1,2] 
#>                          0.5272677                          0.5272677 
#>                       kappa_n[2,2]                       kappa_n[1,3] 
#>                          0.5272677                          0.5272677 
#>                       kappa_n[2,3]                          nu_n[1,1] 
#>                          0.5272677                          0.7065037 
#>                          nu_n[2,1]                          nu_n[1,2] 
#>                          0.7065037                          0.7065037 
#>                          nu_n[2,2]                          nu_n[1,3] 
#>                          0.7065037                          0.7065037 
#>                          nu_n[2,3]             m_n_transformed[1,1,1] 
#>                          0.7065037                          0.4434381 
#>             m_n_transformed[2,1,1]             m_n_transformed[1,2,1] 
#>                          1.3250028                          0.4619834 
#>             m_n_transformed[2,2,1]             m_n_transformed[1,3,1] 
#>                          1.3250028                          0.7696966 
#>             m_n_transformed[2,3,1]           S_n_transformed[1,1,1,1] 
#>                          1.3250028                          0.4626955 
#>           S_n_transformed[2,1,1,1]           S_n_transformed[1,2,1,1] 
#>                          0.5750953                          0.4987560 
#>           S_n_transformed[2,2,1,1]           S_n_transformed[1,3,1,1] 
#>                          0.5273649                          0.4253552 
#>           S_n_transformed[2,3,1,1]       t_scale_transformed[1,1,1,1] 
#>                          0.5869040                          0.5270021 
#>       t_scale_transformed[2,1,1,1]       t_scale_transformed[1,2,1,1] 
#>                          0.5272677                          0.6971519 
#>       t_scale_transformed[2,2,1,1]       t_scale_transformed[1,3,1,1] 
#>                          0.5272677                          1.0264243 
#>       t_scale_transformed[2,3,1,1]                        log_lik_sum 
#>                          0.4441121                          0.9735211 
#>                         log_lik[1]                         log_lik[2] 
#>                          0.9735211                          0.4123210 
#>                         log_lik[3]                         log_lik[4] 
#>                          0.6211765                          0.9735211 
#>                         log_lik[5]                         log_lik[6] 
#>                          0.6421622                          0.4950000 
#>                         log_lik[7]                         log_lik[8] 
#>                          0.9735211                          1.3250028 
#>                         log_lik[9]                        log_lik[10] 
#>                          0.8815803                          0.9722762 
#>                        log_lik[11]                        log_lik[12] 
#>                          1.3250028                          0.5272677 
#>                        log_lik[13]                        log_lik[14] 
#>                          0.5272677                          0.4100459 
#>                        log_lik[15]                           m_0[1,1] 
#>                          0.8840930                          0.4434381 
#>                           m_0[2,1]                         S_0[1,1,1] 
#>                          1.3250028                          0.4626955 
#>                         S_0[2,1,1]                         m_n[1,1,1] 
#>                          0.5750953                          0.4434381 
#>                         m_n[2,1,1]                         m_n[1,2,1] 
#>                          1.3250028                          0.4619834 
#>                         m_n[2,2,1]                         m_n[1,3,1] 
#>                          1.3250028                          0.7696966 
#>                         m_n[2,3,1]                       S_n[1,1,1,1] 
#>                          1.3250028                          0.4626955 
#>                       S_n[2,1,1,1]                       S_n[1,2,1,1] 
#>                          0.5750953                          0.4987560 
#>                       S_n[2,2,1,1]                       S_n[1,3,1,1] 
#>                          0.5273649                          0.4253552 
#>                       S_n[2,3,1,1]                               lp__ 
#>                          0.5869040                          0.1493471 
# }
# \donttest{
fit <- example_ideal_adaptor_stanfit()
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
lp <- log_posterior(fit)
head(lp)
#>   Chain Iteration     Value
#> 1     1         1 -21.48362
#> 2     1         2 -24.47470
#> 3     1         3 -27.48986
#> 4     1         4 -25.44730
#> 5     1         5 -28.66702
#> 6     1         6 -27.05637
r <- rhat(fit)
head(r)
#>   kappa_0_shared_param[1]      nu_0_shared_param[1] category_prior_param[1,1] 
#>                 0.9903867                 0.9896015                 0.9610481 
#> category_prior_param[1,2]     lapse_bias_param[1,1]     lapse_bias_param[1,2] 
#>                 0.9610481                 0.9882127                 0.9882127 
# }
```
