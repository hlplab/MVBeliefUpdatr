# Transform and untransform means, covariances, and cognitive models

Applies an affine (linear) transformation specified in an
[`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
object (from
[`get_affine_transform`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_affine_transform.md))
to a mean vector, covariance matrix, or within an S7 cognitive model.

## Usage

``` r
transform_model(model, transform)

untransform_model(model, transform)

transform_mean(m, transform)

untransform_mean(m, transform)

transform_cov(S, transform)

untransform_cov(S, transform)
```

## Arguments

- model:

  An
  [`MVBU_CognitiveModel`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
  object.

- transform:

  An
  [`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
  object.

- m:

  A mean vector (or scalar for univariate representations).

- S:

  A covariance matrix (or scalar variance for univariate
  representations).

## Value

A model, mean vector, or covariance matrix of the same class as the
input.

## Details

Gaussian and Normal-Inverse-Wishart distributions transform analytically
under affine maps \\y = A(x + b)\\: \$\$\mu' = A(\mu + b), \quad \Sigma'
= A \Sigma A^T\$\$

## Examples

``` r
df <- data.frame(vot = rnorm(30, 20, 10), f0 = rnorm(30, 220, 25))
trans <- get_affine_transform(df, cues = c("vot", "f0"),
                              type = "standardize")

# Transform mean vector and covariance matrix
m <- c("vot" = 20, "f0" = 220)
S <- matrix(c(100, 20, 20, 625), nrow = 2,
            dimnames = list(c("vot", "f0"), c("vot", "f0")))
m_trans <- transform_mean(m, trans)
S_trans <- transform_cov(S, trans)

# Untransform back to original space
untransform_mean(m_trans, trans)
#> vot  f0 
#>  20 220 
untransform_cov(S_trans, trans)
#>     vot  f0
#> vot 100  20
#> f0   20 625

# Transform an entire cognitive model
mod <- example_mvg_ideal_observer(n_cues = 2)
mod_trans <- transform_model(mod, trans)
untransform_model(mod_trans, trans)
#> <MVBeliefUpdatr::MVG_IdealObserver>
#>   Categories (6): /b/, /d/, /g/, /k/, /p/, /t/
#>   Cues (2): VOT, f0_semitones
#>   Decision rule: sampling
#>   Category prior: /b/: 0.167, /d/: 0.167, /g/: 0.167, /k/: 0.167, /p/: 0.167, /t/: 0.167
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (6):
#>     /b/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /d/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /g/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /k/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /p/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /t/: MVG(mu = vector(2), Sigma = matrix(2, 2))
```
