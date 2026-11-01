# S7 class for fixed cognitive parameter specifications

Represents fixed parameter specifications in a dual-space (original vs.
transformed cue space) framework, supporting declarative parameter
fixing across global, category-specific, and cue-specific cognitive
parameters.

Constructs a fixed parameter specification for
[`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md).
Supports fixing parameters globally, per-category, or per-cue, in either
original or transformed cue space.

## Usage

``` r
MVBU_FixedParameter(
  par_class = character(),
  value = list(),
  category = character(),
  cue = character(),
  cue_space = character()
)

fix_parameter(
  value,
  class = c("m_0", "S_0", "kappa_0", "nu_0", "expected_mu_0", "expected_Sigma_0", "mu_0",
    "Sigma_0", "lapse_rate", "lapse_bias", "category_prior", "weights"),
  category = NULL,
  cue = NULL,
  cue_space = c("original", "transformed")
)

# S3 method for class 'MVBU_FixedParameter'
x[i, ...]
```

## Arguments

- par_class:

  Character vector of parameter classes (e.g., `"m_0"`, `"S_0"`,
  `"kappa_0"`, `"nu_0"`, `"expected_mu_0"`, `"expected_Sigma_0"`,
  `"mu_0"`, `"Sigma_0"`, `"lapse_rate"`, `"lapse_bias"`,
  `"category_prior"`, `"weights"`).

- value:

  Numeric value, vector, matrix, or character assignment (e.g.,
  `"lapse_rate = 0.05"`).

- category:

  Character string specifying category name (default: `""` or `NULL` for
  all categories).

- cue:

  Character string specifying cue name (default: `""` or `NULL` for all
  cues).

- cue_space:

  Character string specifying cue space: `"original"` (default) or
  `"transformed"`.

- class:

  Character string specifying parameter class: `"m_0"`, `"S_0"`,
  `"kappa_0"`, `"nu_0"`, `"expected_mu_0"`, `"expected_Sigma_0"`,
  `"mu_0"`, `"Sigma_0"`, `"lapse_rate"`, `"lapse_bias"`,
  `"category_prior"`, or `"weights"`.

- x:

  An `MVBU_FixedParameter` object.

- i:

  Indices specifying elements to extract.

- ...:

  Additional arguments (ignored).

## Value

An object of class `MVBU_FixedParameter`.

A subsetted `MVBU_FixedParameter` object.

## Category Moment Parameter Interpretations

When fixing category representation parameters, note the distinction
between prior hyperparameters, expected prior moments, and known
category parameters:

- `"m_0"` and `"S_0"`:

  Prior hyperparameters of the conjugate prior distribution
  (Normal-Inverse-Wishart or Normal-Inverse-Chi-squared). `"m_0"` is the
  prior category mean vector (center of the normal prior on \\\mu_0\\),
  and `"S_0"` is the prior scatter matrix.

- `"expected_mu_0"` and `"expected_Sigma_0"`:

  Expected moments of category representations under the prior. For
  category means, \\E\[\mu_0\] = m_0\\, so `"expected_mu_0"` is
  canonicalized directly to `"m_0"`. For covariance matrices,
  \\E\[\Sigma_0\] = S_0 / (\nu_0 - K - 1)\\ for NIW (or \\E\[\Sigma_0\]
  = S_0 / (\nu_0 - 2)\\ for NIX/MNIX). Specifying `"expected_Sigma_0"`
  allows users to fix the intuitive expected covariance. Stan input
  preparation automatically rescales this value to the conjugate scatter
  matrix `"S_0"` given \\\nu_0\\.

- `"mu_0"` and `"Sigma_0"`:

  Truly known category representation parameters entering the model as
  fixed, non-random constants without latent prior parameters. Fixing
  either or both moments determines the conjugate belief updating and
  the predictive likelihood behavior:

  Case 1 (both unknown)

  :   Standard NIW or NIX conjugate updating for both mean and
      covariance; the predictive likelihood is Student-\\t\\.

  Case 2 (`"Sigma_0"` known, `"mu_0"` unknown)

  :   Category covariance is fixed to \\\Sigma_0\\. Category mean is
      updated via Normal-mean updating (\\\kappa_n = \kappa_0 + n\\,
      \\\mu_n = (\kappa_0 m_0 + n \bar{x}) / \kappa_n\\). Because the
      covariance is known without uncertainty, the predictive likelihood
      is Normal with covariance \\(1 + 1/\kappa_n) \Sigma_0\\.

  Case 3 (`"mu_0"` known, `"Sigma_0"` unknown)

  :   Category mean is fixed to \\\mu_0\\. Category covariance is
      updated via Inverse-Wishart (NIW) or Inverse-Chi-squared (NIX)
      around the known mean (\\\nu_n = \nu_0 + n\\, \\S_n = S_0 + \sum
      (x_i - \mu_0)(x_i - \mu_0)^T\\). Because the mean is known, the
      predictive likelihood is Student-\\t\\ with \\\nu_n\\ degrees of
      freedom and scale \\S_n / (\nu_n - K + 1)\\ (or \\S_n / \nu_n\\
      for NIX).

  Case 4 (both `"mu_0"` and `"Sigma_0"` known)

  :   Both category moments are known and fixed. No belief updating
      occurs with exposure, and the predictive likelihood is Normal with
      mean \\\mu_0\\ and covariance \\\Sigma_0\\.

## See also

[`parameter_prior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Prior.md),
[`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md)

## Examples

``` r
# Fix cognitive parameters
fix_kappa <- fix_parameter(value = 10, class = "kappa_0")
fix_nu <- fix_parameter(value = 5, class = "nu_0")
fix_lapse <- fix_parameter(value = 0.02, class = "lapse_rate")

# Combine fixed parameters
c(fix_kappa, fix_nu, fix_lapse)
#> MVBU_FixedParameter with 3 specification(s):
#> # A tibble: 3 × 5
#>   class      category cue   cue_space value
#>   <chr>      <chr>    <chr> <chr>     <chr>
#> 1 kappa_0    (all)    (all) original  10   
#> 2 nu_0       (all)    (all) original  5    
#> 3 lapse_rate (all)    (all) original  0.02 
```
