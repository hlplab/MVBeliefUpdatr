# Extract a category-posterior function from a cognitive model

Extract a category-posterior function from a cognitive model

## Usage

``` r
get_category_posterior_function(x, noise_treatment, lapse_treatment, ...)
```

## Arguments

- x:

  A cognitive model object.

- noise_treatment:

  Optional character string specifying treatment of perceptual noise:
  `"no_noise"`, `"sample"`, or `"marginalize"`.

  - `"no_noise"`: No perceptual noise is applied to input cues, and no
    noise is assumed during categorization.

  - `"sample"`: Simulates sensory noise on a specific observation. Noise
    drawn from \\\mathcal{N}(0, \Sigma\_{\text{noise}})\\ is added to
    input cues (and projected into composite cue space in early
    integration), and the observer accounts for their own noise
    uncertainty during categorization by adding
    \\\Sigma\_{\text{noise}}\\ to their likelihood scale.

  - `"marginalize"`: Evaluates categorization in the limit across
    infinitely many noisy encounters (expected consequences). No noise
    is sampled onto the input stimulus, but \\\Sigma\_{\text{noise}}\\
    is added to the category covariance scale.

- lapse_treatment:

  Optional character string specifying treatment of attentional lapses:
  `"no_lapses"`, `"sample"`, or `"marginalize"`.

  - `"no_lapses"`: Lapses are ignored; pure normative posterior
    probabilities are returned (even if `lapse_rate > 0`).

  - `"sample"`: Simulates stochastic attentional lapses on individual
    observations. With probability `lapse_rate`, a lapse occurs and the
    posterior is replaced by `lapse_bias`.

  - `"marginalize"`: Evaluates categorization in the limit across
    infinitely many encounters (expected response probabilities), mixing
    normative posteriors with `lapse_bias`: \\(1 - \lambda) p(c
    \mid x) + \lambda b_c\\.

  Defaults to `NULL` (uses the model's configured lapse treatment).

- ...:

  Additional arguments passed to methods.

## Value

A function that computes category posterior probabilities for one or
more observations.

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
post_fn <- get_category_posterior_function(mod)
post_fn(matrix(c(0, 0), nrow = 1))
#>            /b/       /d/       /g/        /k/        /p/         /t/
#> [1,] 0.4315091 0.3092369 0.2205673 0.01204114 0.02026484 0.006380768
```
