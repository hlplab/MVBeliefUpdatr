# Get perceptual noise covariance matrix from a cognitive model

Extract perceptual noise covariance from a cognitive model.

## Usage

``` r
get_noise(x, ...)

## S7 method for class <MVBU_CognitiveModel>
get_noise(x, ...)
```

## Arguments

- x:

  A cognitive model.

## Value

A matrix or numeric value representing perceptual noise covariance
(\\\Sigma\_{\text{noise}}\\), or `NULL` if no noise is specified.
