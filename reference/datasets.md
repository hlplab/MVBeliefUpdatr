# Phonetic Datasets Bundled with MVBeliefUpdatr

Bundled phonetic corpora and experimental datasets for modeling speech
perception, category representation, and incremental belief updating.

## See also

[`h95`](https://hlplab.github.io/MVBeliefUpdatr/reference/h95.md),
[`pb52`](https://hlplab.github.io/MVBeliefUpdatr/reference/pb52.md),
[`swehvd`](https://hlplab.github.io/MVBeliefUpdatr/reference/swehvd.md),
[`mixer6`](https://hlplab.github.io/MVBeliefUpdatr/reference/mixer6.md)

## Examples

``` r
# Inspect bundled datasets
data(h95)
data(pb52)
head(h95)
#> # A tibble: 6 × 9
#>   speaker sex   type  vowel    f0    f1    f2    f3 duration
#>   <fct>   <fct> <fct> <fct> <dbl> <dbl> <dbl> <dbl>    <dbl>
#> 1 1       male  boys  /æ/     238   630  2423  3166      257
#> 2 1       male  boys  /ɑ/     241   831  1676  2602      212
#> 3 1       male  boys  /ɔ/     247   725  1384  2642      242
#> 4 1       male  boys  /ɛ/     214   713  2095  3129      184
#> 5 1       male  boys  /eɪ/    230   534  2690  3335      222
#> 6 1       male  boys  /ɝ/     240   608  1733  2159      227
head(pb52)
#> # A tibble: 6 × 9
#>   speaker sex   type  vowel repetition    f0    f1    f2    f3
#>   <fct>   <fct> <fct> <fct>      <int> <dbl> <dbl> <dbl> <dbl>
#> 1 1       male  men   /i/            1   160   240  2280  2850
#> 2 1       male  men   /i/            2   186   280  2400  2790
#> 3 1       male  men   /ɪ/            1   203   390  2030  2640
#> 4 1       male  men   /ɪ/            2   192   310  1980  2550
#> 5 1       male  men   /ɛ/            1   161   490  1870  2420
#> 6 1       male  men   /ɛ/            2   155   570  1700  2600
```
