# Evaluate a cognitive model or fitted model

Evaluates the predictions of a cognitive model or fitted Stan model
against observed categorization responses. When applied to fitted model
objects with attached test data (such as
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
or
[`MVBU_StanfitPosterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitPosterior.md)),
leaving `x = NULL` and `response_category = NULL` automatically
evaluates the model against the model's attached test data.

## Usage

``` r
evaluate_model(
  model,
  x = NULL,
  response_category = NULL,
  method = "log_lik",
  decision_rule = if (identical(method, "accuracy")) "criterion" else "proportional",
  return_by_x = FALSE,
  ...
)

## S7 method for class <MVBU_CognitiveModel>
evaluate_model(
  model,
  x = NULL,
  response_category = NULL,
  method = "log_lik",
  decision_rule = if (identical(method, "accuracy")) "criterion" else "proportional",
  return_by_x = FALSE,
  ...
)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
evaluate_model(
  model,
  x = NULL,
  response_category = NULL,
  method = "log_lik",
  decision_rule = if (identical(method, "accuracy")) "criterion" else "proportional",
  return_by_x = FALSE,
  ...
)
```

## Arguments

- model:

  A cognitive model object (e.g.,
  [`MVG_IdealObserver`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
  [`NIW_IdealAdaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md))
  or
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md).

- x:

  Observations/stimuli (matrix, data frame, list of numeric vectors, or
  numeric vector). If `NULL` and `model` has attached test data (e.g.,
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  or
  [`MVBU_StanfitPosterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitPosterior.md)),
  stimuli `x` are automatically extracted from the model's attached test
  data.

- response_category:

  Vector of observed category responses for each observation. If `NULL`
  and `model` has attached test data (e.g.,
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  or
  [`MVBU_StanfitPosterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitPosterior.md)),
  responses are automatically extracted from the model's attached test
  data column `"response"`.

- method:

  Evaluation method(s): `"log_lik"` (default),
  `"log_lik_permutation_constant"`, `"accuracy"`, or
  `"likelihood-up-to-constant"` (deprecated alias). Can be a character
  vector to return multiple metrics.

- decision_rule:

  Decision rule to apply: `"proportional"`, `"criterion"`, or
  `"sampling"`. Defaults to `"criterion"` when `method == "accuracy"`,
  and `"proportional"` otherwise.

- return_by_x:

  Logical; if `TRUE`, returns evaluation metrics grouped by unique
  stimulus values \\x\\. Defaults to `FALSE`.

- ...:

  Additional arguments passed to methods.

## Value

A numeric scalar, a tibble (if `return_by_x = TRUE`), or a named list of
these if multiple methods were requested.

## Details

Supported evaluation methods include:

- `"log_lik"`: Trial-by-trial categorical log-likelihood (default):
  \$\$\log L\_{\text{cat}} = \sum\_{i=1}^N \log P(Y = y_i \mid X = x_i)
  = \sum\_{u=1}^U \sum\_{j=1}^K n\_{u,j} \log p\_{u,j}\$\$ where \\u\\
  indexes unique stimulus locations, \\j\\ indexes categories,
  \\n\_{u,j}\\ is the count of responses for category \\j\\ at stimulus
  \\u\\, and \\p\_{u,j}\\ is the model's predicted probability of
  category \\j\\ at \\u\\. This represents the log probability of
  observing the exact sequence of individual responses \\(y_1, \dots,
  y_N)\\. This is the standard likelihood metric in GLMs, logistic
  regression, and Bayesian Stan models.

- `"accuracy"`: Proportion of correct / matching responses.

- `"likelihood-up-to-constant"`: Deprecated alias for `"log_lik"`
  (supported until version 0.1).

Additionally, the following method is available to support calculation
of order invariant log-likelihoods:

- `"log_lik_permutation_constant"`: The log multinomial combinatorial
  permutation factor: \$\$C = \sum\_{u=1}^U \left( \log(N_u!) -
  \sum\_{j=1}^K \log(n\_{u,j}!) \right) = \sum\_{u=1}^U \log
  \binom{N_u}{n\_{u,1}, \dots, n\_{u,K}}\$\$ where \\N_u = \sum_j
  n\_{u,j}\\ is the total number of presentations at stimulus location
  \\u\\. This factor accounts for all permutations of response order at
  each stimulus location. Crucially, \\C\\ depends purely on the
  observed data and is completely independent of the cognitive model.

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
x_stim <- matrix(c(10, -1, 55, 0.5), nrow = 2, byrow = TRUE)
evaluate_model(
  mod,
  x = x_stim,
  response_category = c("/b/", "/k/"),
  method = "log_lik"
)
#> [1] -1.624713
```
