# loo for MVBU stanfit objects

`loo` method for
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
objects. For details, see
[loo::loo](https://mc-stan.org/loo/reference/loo.html).

## Usage

``` r
# S3 method for class 'MVBU_Stanfit'
loo(
  x,
  pars = "log_lik",
  ...,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1)
)
```

## Arguments

- x:

  An
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object.

- pars:

  Parameter name containing the pointwise log likelihood. (default:
  `"log_lik"`)

- ...:

  Additional arguments passed to
  [`loo::loo.array`](https://mc-stan.org/loo/reference/loo.html).

- save_psis:

  Logical scalar passed to
  [`loo::loo.array`](https://mc-stan.org/loo/reference/loo.html).
  (default: `FALSE`)

- cores:

  Number of cores for parallel execution. (default:
  `getOption("mc.cores", 1)`)

## Value

A `loo` object.
