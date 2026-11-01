# Add Gaussian noise to cues in a data set

Adds Gaussian noise to cues in a data set.

## Usage

``` r
add_noise(data, cues, method = c("sample", "marginalize"), Sigma)
```

## Arguments

- data:

  `data.frame` or `tibble`.

- cues:

  Column names of cues in `data`.

- method:

  If `"sample"` (default) then a sample of noise is added to the input.
  If `"marginalize"` then each observation is transformed into the
  marginal distribution that results from convolving the input with
  noise. Currently, only `"sample"` is implemented.

- Sigma:

  Covariance matrix of noise.

## Value

Same as `data` with noise added to cues.

## Examples

``` r
df <- data.frame(vot = c(10, 20, 30), f0 = c(150, 180, 210))
noise_cov <- matrix(c(4, 0, 0, 9), nrow = 2)
add_noise(df, cues = c("vot", "f0"), Sigma = noise_cov)
#>   vot  f0                 cues
#> 1  10 150 8.914342, 151.742989
#> 2  20 180  21.53636, 181.39130
#> 3  30 210  28.22845, 206.70066
```
