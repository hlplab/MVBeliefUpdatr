#' Phonetic Datasets Bundled with MVBeliefUpdatr
#'
#' Bundled phonetic corpora and experimental datasets for modeling speech perception,
#' category representation, and incremental belief updating.
#'
#' @docType data
#' @name datasets
NULL

#' Hillenbrand et al. (1995) Vowel Acoustic Data
#'
#' Acoustic measurements of 12 American English vowels produced by men, women,
#' boys, and girls in /h-V-d/ syllables. This dataset is a standardized copy of
#' the dataset bundled in the \pkg{phonTools} package \insertCite{barreda2015}{MVBeliefUpdatr},
#' originally published by \insertCite{hillenbrand1995;textual}{MVBeliefUpdatr}.
#'
#' @format A tibble with 1,668 rows and 9 variables:
#' \describe{
#'   \item{speaker}{Speaker identifier as a factor.}
#'   \item{sex}{Speaker sex as a factor with levels \code{"female"} and \code{"male"}.}
#'   \item{type}{Speaker demographic group as a factor with levels \code{"men"},
#'     \code{"women"}, \code{"boys"}, and \code{"girls"}.}
#'   \item{vowel}{Standardized Unicode IPA vowel label in phonemic slashes as a factor
#'     with 12 levels: \code{/æ/}, \code{/ɝ/}, \code{/ɑ/}, \code{/eɪ/}, \code{/ɛ/},
#'     \code{/i/}, \code{/ɪ/}, \code{/oʊ/}, \code{/ɔ/}, \code{/u/}, \code{/ʊ/}, \code{/ʌ/}.}
#'   \item{f0}{Fundamental frequency in Hertz (Hz).}
#'   \item{f1}{First formant frequency in Hertz (Hz).}
#'   \item{f2}{Second formant frequency in Hertz (Hz).}
#'   \item{f3}{Third formant frequency in Hertz (Hz).}
#'   \item{duration}{Vowel duration in milliseconds (ms).}
#' }
#'
#' @source Standardized from \code{phonTools::h95} \insertCite{barreda2015}{MVBeliefUpdatr}.
#' @references
#' \insertRef{hillenbrand1995}{MVBeliefUpdatr}
#'
#' \insertRef{barreda2015}{MVBeliefUpdatr}
#'
#' @seealso \code{\link{pb52}}, \code{\link{swehvd}}, \code{\link{mixer6}}
#' @docType data
#' @keywords data
#' @usage data(h95)
"h95"

#' Peterson & Barney (1952) Vowel Acoustic Data
#'
#' Acoustic measurements of 10 American English vowels produced by 76 speakers
#' (men, women, and children) in /h-V-d/ words. This dataset is a standardized copy
#' of the dataset bundled in the \pkg{phonTools} package \insertCite{barreda2015}{MVBeliefUpdatr},
#' originally published by \insertCite{peterson-barney1952;textual}{MVBeliefUpdatr}.
#'
#' @format A tibble with 1,520 rows and 9 variables:
#' \describe{
#'   \item{speaker}{Speaker identifier as a factor.}
#'   \item{sex}{Speaker sex as a factor with levels \code{"female"} and \code{"male"}.}
#'   \item{type}{Speaker demographic group as a factor with levels \code{"men"},
#'     \code{"women"}, and \code{"children"}.}
#'   \item{vowel}{Standardized Unicode IPA vowel label in phonemic slashes as a factor
#'     with 10 levels: \code{/i/}, \code{/ɪ/}, \code{/ɛ/}, \code{/æ/}, \code{/ʌ/},
#'     \code{/ɑ/}, \code{/ɔ/}, \code{/ʊ/}, \code{/u/}, \code{/ɝ/}.}
#'   \item{repetition}{Repetition index (1 or 2).}
#'   \item{f0}{Fundamental frequency in Hertz (Hz).}
#'   \item{f1}{First formant frequency in Hertz (Hz).}
#'   \item{f2}{Second formant frequency in Hertz (Hz).}
#'   \item{f3}{Third formant frequency in Hertz (Hz).}
#' }
#'
#' @source Standardized from \code{phonTools::pb52} \insertCite{barreda2015}{MVBeliefUpdatr}.
#' @references
#' \insertRef{peterson-barney1952}{MVBeliefUpdatr}
#'
#' \insertRef{barreda2015}{MVBeliefUpdatr}
#'
#' @seealso \code{\link{h95}}, \code{\link{swehvd}}, \code{\link{mixer6}}
#' @docType data
#' @keywords data
#' @usage data(pb52)
"pb52"

#' SwehVd: Swedish Vowels in /h-V-d/ Context
#'
#' Phonetic measurements from 24 native Central Swedish female speakers producing
#' 21 vowel categories in carrier sentences in /h-V-d/ context \insertCite{persson2021}{MVBeliefUpdatr}.
#'
#' @format A tibble with 23,685 rows and 17 variables:
#' \describe{
#'   \item{speaker}{Talker identifier as a factor (e.g., \code{"SW_001"}).}
#'   \item{sex}{Speaker sex as a factor with levels \code{"female"} and \code{"male"}
#'     (all 24 talkers in SwehVd are female).}
#'   \item{age}{Speaker age in years as an integer.}
#'   \item{vowel}{Standardized Unicode IPA vowel label in phonemic slashes as a factor
#'     with 21 levels: \code{/ɑː/}, \code{/a/}, \code{/eː/}, \code{/ɛ/}, \code{/iː/},
#'     \code{/ɪ/}, \code{/uː/}, \code{/ʊ/}, \code{/ʉː/}, \code{/ɵ/}, \code{/yː/},
#'     \code{/ʏ/}, \code{/ɛː/}, \code{/æː/}, \code{/æ/}, \code{/oː/}, \code{/ɔ/},
#'     \code{/øː/}, \code{/ø/}, \code{/œː/}, \code{/œ/}.}
#'   \item{transcribed_vowel}{Transcription description character string.}
#'   \item{quantity}{Phonological vowel quantity as a factor with levels \code{"long"}
#'     and \code{"short"}.}
#'   \item{word}{Produced stimulus word.}
#'   \item{token}{Token index within word and speaker.}
#'   \item{trial}{Experimental trial index.}
#'   \item{location}{Relative measurement timepoint along vowel duration (20, 35, 50, 65, 80 percent).}
#'   \item{f0}{Fundamental frequency in Hertz (Hz).}
#'   \item{f1}{First formant frequency in Hertz (Hz).}
#'   \item{f2}{Second formant frequency in Hertz (Hz).}
#'   \item{f3}{Third formant frequency in Hertz (Hz).}
#'   \item{duration}{Total vowel duration in seconds.}
#'   \item{speech_rate}{Local speech rate in syllables per minute.}
#'   \item{unreliable_measurement}{Logical flag indicating whether measurement was flagged
#'     as potentially unreliable.}
#' }
#'
#' @source Swedish Vowel Database (SwehVd) project \insertCite{persson2021}{MVBeliefUpdatr}.
#' @references
#' \insertRef{persson2021}{MVBeliefUpdatr}
#'
#' @seealso \code{\link{h95}}, \code{\link{pb52}}, \code{\link{mixer6}}
#' @docType data
#' @keywords data
#' @usage data(swehvd)
"swehvd"

#' Mixer 6 Connected Speech Word-Initial Stop Consonant Data
#'
#' Word-initial stop consonant acoustic measurements from connected speech productions
#' across 180 native US English speakers in the Mixer 6 corpus \insertCite{chodroff-wilson2018}{MVBeliefUpdatr}.
#' Contains annotations for voice onset time (VOT), fundamental frequency (f0 in Hz
#' and talker-mean-centered semitones), vowel duration, word duration, and speaking rate.
#'
#' @format A tibble with 96,357 rows and 21 variables:
#' \describe{
#'   \item{speaker}{Speaker identifier as a factor.}
#'   \item{sex}{Speaker sex as a factor with levels \code{"female"} and \code{"male"}.}
#'   \item{stop}{Standardized Unicode IPA stop consonant category in phonemic slashes
#'     as a factor with 6 levels: \code{/b/}, \code{/p/}, \code{/d/}, \code{/t/},
#'     \code{/g/}, \code{/k/}.}
#'   \item{stop_poa}{Place of articulation as a factor with levels \code{"labial"},
#'     \code{"coronal"}, and \code{"dorsal"}.}
#'   \item{voicing}{Phonological voicing as a factor with levels \code{"voiced"}
#'     and \code{"voiceless"}.}
#'   \item{word}{Produced word string.}
#'   \item{vowel}{Following vowel context in ARPABET notation (e.g., \code{"EH1"}, \code{"IY1"}).}
#'   \item{trial}{Trial number within recording session.}
#'   \item{filename}{Audio recording filename.}
#'   \item{start}{Start timestamp in seconds.}
#'   \item{end}{End timestamp in seconds.}
#'   \item{f0}{Fundamental frequency at vowel onset in Hertz (Hz).}
#'   \item{f0_semitones}{Talker-mean-centered fundamental frequency in semitones:
#'     \eqn{12 \times \log_2(f_0 / \bar{f}_{0,\text{speaker}})}.}
#'   \item{vot}{Voice onset time in milliseconds (ms).}
#'   \item{vowel_duration}{Following vowel duration in milliseconds (ms).}
#'   \item{word_duration}{Word duration in milliseconds (ms).}
#'   \item{speech_rate}{Speaking rate in syllables per second.}
#'   \item{pos}{Utterance position (\code{"utt_init"}, \code{"utt_mid"}).}
#'   \item{syll}{Syllable count category (\code{"one"}, \code{"two"}, \code{"more"}).}
#'   \item{type}{Word type category (\code{"lex"} for lexical/content, \code{"func"} for function).}
#'   \item{session}{Session number.}
#' }
#'
#' @source Mixer 6 corpus via Eleanor Chodroff's OSF repository (\url{https://osf.io/jt5mc/}).
#' @references
#' \insertRef{chodroff-wilson2018}{MVBeliefUpdatr}
#'
#' @seealso \code{\link{h95}}, \code{\link{pb52}}, \code{\link{swehvd}}
#' @docType data
#' @keywords data
#' @usage data(mixer6)
"mixer6"


# deprecated --------------------------------------------------------------

#' Chodroff & Wilson (2018) data on word-initial stop-voicing in native US English
#'
#' `r lifecycle::badge("deprecated")`
#'
#' `ChodroffWilson2018` is deprecated. Please use [mixer6] instead, which provides
#' the full Mixer 6 connected speech dataset with standardized column names and
#' Unicode IPA notation.
#'
#' @docType data
#' @name ChodroffWilson2018
#' @usage data(ChodroffWilson2018)
#' @keywords data
#' @references
#' \insertRef{chodroff-wilson2018}{MVBeliefUpdatr}
#' @seealso \code{\link{mixer6}}
"ChodroffWilson2018"
