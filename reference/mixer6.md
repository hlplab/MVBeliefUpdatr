# Mixer 6 Connected Speech Word-Initial Stop Consonant Data

Word-initial stop consonant acoustic measurements from connected speech
productions across 180 native US English speakers in the Mixer 6 corpus
(Chodroff and Wilson 2018) . Contains annotations for voice onset time
(VOT), fundamental frequency (f0 in Hz and talker-mean-centered
semitones), vowel duration, word duration, and speaking rate.

## Usage

``` r
data(mixer6)
```

## Format

A tibble with 96,357 rows and 21 variables:

- speaker:

  Speaker identifier as a factor.

- sex:

  Speaker sex as a factor with levels `"female"` and `"male"`.

- stop:

  Standardized Unicode IPA stop consonant category in phonemic slashes
  as a factor with 6 levels: `/b/`, `/p/`, `/d/`, `/t/`, `/g/`, `/k/`.

- stop_poa:

  Place of articulation as a factor with levels `"labial"`, `"coronal"`,
  and `"dorsal"`.

- voicing:

  Phonological voicing as a factor with levels `"voiced"` and
  `"voiceless"`.

- word:

  Produced word string.

- vowel:

  Following vowel context in ARPABET notation (e.g., `"EH1"`, `"IY1"`).

- trial:

  Trial number within recording session.

- filename:

  Audio recording filename.

- start:

  Start timestamp in seconds.

- end:

  End timestamp in seconds.

- f0:

  Fundamental frequency at vowel onset in Hertz (Hz).

- f0_semitones:

  Talker-mean-centered fundamental frequency in semitones: \\12 \times
  \log_2(f_0 / \bar{f}\_{0,\text{speaker}})\\.

- vot:

  Voice onset time in milliseconds (ms).

- vowel_duration:

  Following vowel duration in milliseconds (ms).

- word_duration:

  Word duration in milliseconds (ms).

- speech_rate:

  Speaking rate in syllables per second.

- pos:

  Utterance position (`"utt_init"`, `"utt_mid"`).

- syll:

  Syllable count category (`"one"`, `"two"`, `"more"`).

- type:

  Word type category (`"lex"` for lexical/content, `"func"` for
  function).

- session:

  Session number.

## Source

Mixer 6 corpus via Eleanor Chodroff's OSF repository
(<https://osf.io/jt5mc/>).

## References

Chodroff E, Wilson C (2018). “Predictability of stop consonant phonetics
across talkers: Between-category and within-category dependencies among
cues for place and voice.” *Linguistics Vanguard*, **4**(s2), 20170047.

## See also

[`h95`](https://hlplab.github.io/MVBeliefUpdatr/reference/h95.md),
[`pb52`](https://hlplab.github.io/MVBeliefUpdatr/reference/pb52.md),
[`swehvd`](https://hlplab.github.io/MVBeliefUpdatr/reference/swehvd.md)

## Examples

``` r
data(mixer6)
head(mixer6)
#> # A tibble: 6 × 21
#>   speaker sex    stop  stop_poa voicing word  vowel trial filename  start    end
#>   <fct>   <fct>  <fct> <fct>    <fct>   <chr> <chr> <int> <chr>     <dbl>  <dbl>
#> 1 111138  female /b/   labial   voiced  BENE… EH1       1 111138_…   4.43   4.44
#> 2 111138  female /b/   labial   voiced  BEING IY1       2 111138_…   4.95   4.96
#> 3 111138  female /b/   labial   voiced  BY    AY1      37 111138_… 154.   154.  
#> 4 111138  female /b/   labial   voiced  BUSH  UH1     129 111138_… 563.   563.  
#> 5 111138  female /b/   labial   voiced  BACK  AE1     130 111138_… 575.   575.  
#> 6 111138  female /b/   labial   voiced  BE    IY0     218 111138_… 871.   871.  
#> # ℹ 10 more variables: f0 <dbl>, f0_semitones <dbl>, vot <dbl>,
#> #   vowel_duration <dbl>, word_duration <dbl>, speech_rate <dbl>, pos <chr>,
#> #   syll <chr>, type <chr>, session <int>
table(mixer6$stop, mixer6$voicing)
#>      
#>       voiced voiceless
#>   /b/  17659         0
#>   /p/      0     13522
#>   /d/  24369         0
#>   /t/      0      8219
#>   /g/  16179         0
#>   /k/      0     16409
```
