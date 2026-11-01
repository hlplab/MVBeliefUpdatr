# SwehVd: Swedish Vowels in /h-V-d/ Context

Phonetic measurements from 24 native Central Swedish female speakers
producing 21 vowel categories in carrier sentences in /h-V-d/ context
(Persson et al. 2021) .

## Usage

``` r
data(swehvd)
```

## Format

A tibble with 23,685 rows and 17 variables:

- speaker:

  Talker identifier as a factor (e.g., `"SW_001"`).

- sex:

  Speaker sex as a factor with levels `"female"` and `"male"` (all 24
  talkers in SwehVd are female).

- age:

  Speaker age in years as an integer.

- vowel:

  Standardized Unicode IPA vowel label in phonemic slashes as a factor
  with 21 levels (see `levels(swehvd$vowel)`).

- transcribed_vowel:

  Transcription description character string.

- quantity:

  Phonological vowel quantity as a factor with levels `"long"` and
  `"short"`.

- word:

  Produced stimulus word.

- token:

  Token index within word and speaker.

- trial:

  Experimental trial index.

- location:

  Relative measurement timepoint along vowel duration (20, 35, 50, 65,
  80 percent).

- f0:

  Fundamental frequency in Hertz (Hz).

- f1:

  First formant frequency in Hertz (Hz).

- f2:

  Second formant frequency in Hertz (Hz).

- f3:

  Third formant frequency in Hertz (Hz).

- duration:

  Total vowel duration in seconds.

- speech_rate:

  Local speech rate in syllables per minute.

- unreliable_measurement:

  Logical flag indicating whether measurement was flagged as potentially
  unreliable.

## Source

Swedish Vowel Database (SwehVd) project (Persson et al. 2021) .

## References

Persson C, Kurumada C, Jaeger TF (2021). “SwehVd: A Swedish database of
vowels in /hVd/ context.” *OSF*.
<https://github.com/hlplab/SwehVd-article>.

## See also

[`h95`](https://hlplab.github.io/MVBeliefUpdatr/reference/h95.md),
[`pb52`](https://hlplab.github.io/MVBeliefUpdatr/reference/pb52.md),
[`mixer6`](https://hlplab.github.io/MVBeliefUpdatr/reference/mixer6.md)

## Examples

``` r
data(swehvd)
head(swehvd)
#> # A tibble: 6 × 17
#>   speaker sex      age vowel transcribed_vowel quantity word  token trial
#>   <fct>   <fct>  <int> <fct> <chr>             <fct>    <chr> <int> <int>
#> 1 SW_001  female    30 /ɑː/  targeted          long     had       1     1
#> 2 SW_001  female    30 /ɑː/  targeted          long     had       1     1
#> 3 SW_001  female    30 /ɑː/  targeted          long     had       1     1
#> 4 SW_001  female    30 /ɑː/  targeted          long     had       1     1
#> 5 SW_001  female    30 /ɑː/  targeted          long     had       1     1
#> 6 SW_001  female    30 /ɑː/  targeted          long     had       2     2
#> # ℹ 8 more variables: location <dbl>, f0 <dbl>, f1 <dbl>, f2 <dbl>, f3 <dbl>,
#> #   duration <dbl>, speech_rate <dbl>, unreliable_measurement <lgl>
table(swehvd$quantity, swehvd$vowel)
#>        
#>         /ɑː/  /a/ /eː/  /ɛ/ /iː/  /ɪ/ /uː/  /ʊ/ /ʉː/  /ɵ/ /yː/  /ʏ/ /ɛː/ /æː/
#>   long  1090    0 1075    0 1055    0 1110    0 1070    0 1075    0 1075 1065
#>   short    0 1070    0 2145    0 1095    0 1105    0 1075    0 1075    0    0
#>        
#>          /æ/ /oː/  /ɔ/ /øː/  /ø/ /œː/  /œ/
#>   long     0 1065    0 1090    0 1055    0
#>   short 1075    0 1070    0 1075    0 1075
```
