# Get perceptual noise treatment from a cognitive model

Extract perceptual noise treatment from a cognitive model.

## Usage

``` r
get_noise_treatment(x, ...)

## S7 method for class <MVBU_CognitiveModel>
get_noise_treatment(x, ...)
```

## Arguments

- x:

  A cognitive model.

## Value

A character string representing the noise treatment (`"no_noise"`,
`"sample"`, or `"marginalize"`).
