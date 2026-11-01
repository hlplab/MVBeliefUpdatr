# S7 class for prior distribution specifications

Represents prior distribution specifications in a dual-space (original
vs. transformed cue space) framework, supporting flexible brms-style
declarations across parameter classes.

Constructs a prior specification in the style of brms.

## Usage

``` r
MVBU_Prior(
  declaration = character(),
  par_class = character(),
  category = character(),
  cue = character(),
  par = character(),
  cue_space = character()
)

parameter_prior(
  declaration,
  class = c("tau_0", "cor_0", "m_0", "weight", "weights", "lapse_rate", "lapse_bias",
    "category_prior", "kappa_0", "nu_0"),
  category = NULL,
  cue = NULL,
  par = NULL,
  cue_space = c("transformed", "original")
)

# S3 method for class 'MVBU_Prior'
c(..., recursive = FALSE)

# S3 method for class 'MVBU_Prior'
x[i, ...]

# S3 method for class 'MVBU_Prior'
length(x)

# S3 method for class 'MVBU_Prior'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'MVBU_Prior'
as_tibble(x, ...)

# S3 method for class 'MVBU_Prior'
print(x, ...)
```

## Arguments

- declaration:

  Character string or numeric scalar specifying the prior distribution
  (e.g., `"cauchy(0, 5)"`, `"lkj(1)"`, `"normal(0, 1)"`, or numeric
  `5`).

- par_class:

  Character vector of parameter classes (e.g., `"tau_0"`, `"cor_0"`,
  `"m_0"`, `"weight"`, `"lapse_rate"`, `"lapse_bias"`,
  `"category_prior"`, `"kappa_0"`, `"nu_0"`).

- category:

  Character string specifying category name (default: `""` or `NULL` for
  class-wide priors).

- cue:

  Character string specifying cue name (default: `""` or `NULL` for all
  cues).

- par:

  Character string specifying an optional coefficient, cue, or category
  name (default: `NULL` for class-wide priors).

- cue_space:

  Character string specifying whether the prior is formulated in
  transformed cue space (`"transformed"`, the default) or original
  physical cue space (`"original"`).

- class:

  Character string specifying the parameter class. Supported classes
  include:

  - `"tau_0"`: Scale / standard deviation for category covariance.

  - `"cor_0"`: Correlation matrix / Cholesky factor `L_omega`.

  - `"m_0"`: Category prior mean `m_0`.

  - `"weight"`: Cue integration weights.

  - `"lapse_rate"`: Lapse probability.

  - `"lapse_bias"`: Lapse response distribution simplex.

  - `"category_prior"`: Category prior probability simplex.

  - `"kappa_0"`: Prior pseudocount for mean.

  - `"nu_0"`: Prior pseudocount for covariance.

- ...:

  Additional arguments passed to print.

- recursive:

  Logical (default: FALSE).

- x:

  An `MVBU_Prior` object.

- i:

  Indices specifying elements to extract.

- row.names:

  Optional row names (default: NULL).

- optional:

  Logical (default: FALSE).

## Value

An object of class `MVBU_Prior`.

A single combined `MVBU_Prior` object.

A subsetted `MVBU_Prior` object.

An integer length.

A `data.frame` of prior specifications.

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) of
prior specifications.

The object `x`, invisibly.

## See also

[`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md),
[`fix_parameter`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_FixedParameter.md)

## Examples

``` r
# Specify class-wide priors
p1 <- parameter_prior("cauchy(0, 5)", class = "tau_0")
p2 <- parameter_prior("lkj(1)", class = "cor_0")
p3 <- parameter_prior("normal(0, 1)", class = "m_0")

# Combine priors into a list
c(p1, p2, p3)
#> MVBeliefUpdatr Prior Specifications:
#>   declaration class category   cue par   cue_space
#>  cauchy(0, 5) tau_0    (all) (all)     transformed
#>        lkj(1) cor_0    (all) (all)     transformed
#>  normal(0, 1)   m_0    (all) (all)     transformed
```
