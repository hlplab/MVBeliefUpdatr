# data-raw/import_phonetic_datasets.R
# Prepares bundled phonetic datasets for MVBeliefUpdatr:
# 1. h95 (Hillenbrand et al., 1995 via phonTools)
# 2. pb52 (Peterson & Barney, 1952 via phonTools)
# 3. swehvd (Persson et al., 2021 / SwehVd project)
# 4. mixer6 (Chodroff & Wilson, 2018 / Mixer 6 corpus via OSF jt5mc)
#
# Consistent column ordering:
#   speaker -> demographics (sex, type, age) -> category (vowel or stop) ->
#   sub-features/context (stop_poa, voicing, word, etc.) -> task/trial metadata ->
#   acoustic cues (f0, formants, vot) -> duration/rate -> quality flags.

set.seed(42)
options(timeout = 600)

cat("=== 1. Processing h95 ===\n")
data(h95, package = "phonTools", envir = environment())

# Type: m = men, w = women, b = boys, g = girls
# Sex: m/b = male, w/g = female
h95_type_map <- c(m = "men", w = "women", b = "boys", g = "girls")
h95_sex_map <- c(m = "male", w = "female", b = "male", g = "female")

# Vowel mapping to standardized Unicode IPA notation:
h95_vowel_map <- c(
  "{"  = "/æ/",
  "3'" = "/ɝ/",
  "A"  = "/ɑ/",
  "e"  = "/eɪ/",
  "E"  = "/ɛ/",
  "i"  = "/i/",
  "I"  = "/ɪ/",
  "o"  = "/oʊ/",
  "O"  = "/ɔ/",
  "u"  = "/u/",
  "U"  = "/ʊ/",
  "V"  = "/ʌ/"
)

h95_clean <- tibble::tibble(
  speaker = factor(as.character(h95$speaker)),
  sex = factor(h95_sex_map[as.character(h95$type)], levels = c("female", "male")),
  type = factor(h95_type_map[as.character(h95$type)], levels = c("men", "women", "boys", "girls")),
  vowel = factor(h95_vowel_map[as.character(h95$vowel)], levels = unname(h95_vowel_map)),
  f0 = as.numeric(h95$f0),
  f1 = as.numeric(h95$f1),
  f2 = as.numeric(h95$f2),
  f3 = as.numeric(h95$f3),
  duration = as.numeric(h95$dur)
)

h95 <- h95_clean
save(h95, file = "data/h95.rda", compress = "xz")
cat("Saved data/h95.rda:", nrow(h95), "rows\n")


cat("=== 2. Processing pb52 ===\n")
data(pb52, package = "phonTools", envir = environment())

pb52_type_map <- c(m = "men", w = "women", c = "children")
pb52_sex_map <- c(f = "female", m = "male")

pb52_vowel_map <- c(
  "{"  = "/æ/",
  "3'" = "/ɝ/",
  "A"  = "/ɑ/",
  "E"  = "/ɛ/",
  "i"  = "/i/",
  "I"  = "/ɪ/",
  "O"  = "/ɔ/",
  "u"  = "/u/",
  "U"  = "/ʊ/",
  "V"  = "/ʌ/"
)

pb52_clean <- tibble::tibble(
  speaker = factor(as.character(pb52$speaker)),
  sex = factor(pb52_sex_map[as.character(pb52$sex)], levels = c("female", "male")),
  type = factor(pb52_type_map[as.character(pb52$type)], levels = c("men", "women", "children")),
  vowel = factor(pb52_vowel_map[as.character(pb52$vowel)], levels = unname(pb52_vowel_map)),
  repetition = as.integer(pb52$repetition),
  f0 = as.numeric(pb52$f0),
  f1 = as.numeric(pb52$f1),
  f2 = as.numeric(pb52$f2),
  f3 = as.numeric(pb52$f3)
)

pb52 <- pb52_clean
save(pb52, file = "data/pb52.rda", compress = "xz")
cat("Saved data/pb52.rda:", nrow(pb52), "rows\n")


cat("=== 3. Processing swehvd ===\n")
swehvd_url <- "https://raw.githubusercontent.com/hlplab/SwehVd-article/master/data/phonetic%20vowel%20statistics/Swedish/Persson_2021_L1_vowels_wDistrOutliers.csv"
swehvd_raw <- readr::read_csv(swehvd_url, show_col_types = FALSE)

# Swedish vowels: replace square brackets with standard phonemic slashes
swehvd_vowels_raw <- as.character(swehvd_raw$category)
swehvd_vowels_ipa <- gsub("^\\[(.*)\\]$", "/\\1/", swehvd_vowels_raw)
swehvd_unique_vowels <- unique(swehvd_vowels_ipa)

swehvd_clean <- tibble::tibble(
  speaker = factor(as.character(swehvd_raw$Talker)),
  sex = factor(rep("female", nrow(swehvd_raw)), levels = c("female", "male")),
  age = as.integer(swehvd_raw$Age),
  vowel = factor(swehvd_vowels_ipa, levels = swehvd_unique_vowels),
  transcribed_vowel = as.character(swehvd_raw$Transcribed_vowel),
  quantity = factor(as.character(swehvd_raw$Quantity), levels = c("long", "short")),
  word = as.character(swehvd_raw$Word),
  token = as.integer(swehvd_raw$Token),
  trial = as.integer(swehvd_raw$Trial),
  location = as.numeric(swehvd_raw$Location),
  f0 = as.numeric(swehvd_raw$F0),
  f1 = as.numeric(swehvd_raw$F1),
  f2 = as.numeric(swehvd_raw$F2),
  f3 = as.numeric(swehvd_raw$F3),
  duration = as.numeric(swehvd_raw$Duration),
  speech_rate = as.numeric(swehvd_raw$SR),
  unreliable_measurement = as.logical(swehvd_raw$Unreliable_measurement)
)

swehvd <- swehvd_clean
save(swehvd, file = "data/swehvd.rda", compress = "xz")
cat("Saved data/swehvd.rda:", nrow(swehvd), "rows\n")


cat("=== 4. Processing mixer6 ===\n")
mixer6_url <- "https://osf.io/download/7xc8z/"
# Read using readr to avoid base R interpreting 'F' as FALSE:
mixer6_raw <- readr::read_csv(
  mixer6_url,
  col_types = readr::cols(
    subj = readr::col_character(),
    gender = readr::col_character(),
    stop = readr::col_character(),
    poa = readr::col_character(),
    word = readr::col_character(),
    vowel = readr::col_character(),
    filename = readr::col_character(),
    pos = readr::col_character(),
    syll = readr::col_character(),
    type = readr::col_character()
  )
)

# Standardize stop consonants to IPA slashes:
stop_map <- c(
  B = "/b/",
  P = "/p/",
  D = "/d/",
  T = "/t/",
  G = "/g/",
  K = "/k/"
)

poa_map <- c(
  lab = "labial",
  cor = "coronal",
  dor = "dorsal"
)

sex_map <- c(
  F = "female",
  M = "male"
)

voicing_map <- c(
  "/b/" = "voiced",
  "/d/" = "voiced",
  "/g/" = "voiced",
  "/p/" = "voiceless",
  "/t/" = "voiceless",
  "/k/" = "voiceless"
)

stop_ipa <- stop_map[toupper(mixer6_raw$stop)]
voicing <- factor(voicing_map[stop_ipa], levels = c("voiced", "voiceless"))
stop_poa <- factor(poa_map[mixer6_raw$poa], levels = c("labial", "coronal", "dorsal"))
sex <- factor(sex_map[toupper(mixer6_raw$gender)], levels = c("female", "male"))

# Compute talker-mean-centered f0 in semitones:
f0_hz <- as.numeric(mixer6_raw$usef0)
spk_chr <- as.character(mixer6_raw$subj)
mean_f0_by_spk <- ave(f0_hz, spk_chr, FUN = function(x) mean(x, na.rm = TRUE))
f0_semitones <- 12 * log2(f0_hz / mean_f0_by_spk)

mixer6_clean <- tibble::tibble(
  speaker = factor(spk_chr),
  sex = sex,
  stop = factor(stop_ipa, levels = c("/b/", "/p/", "/d/", "/t/", "/g/", "/k/")),
  stop_poa = stop_poa,
  voicing = voicing,
  word = as.character(mixer6_raw$word),
  vowel = as.character(mixer6_raw$vowel),
  trial = as.integer(mixer6_raw$trial),
  filename = as.character(mixer6_raw$filename),
  start = as.numeric(mixer6_raw$start),
  end = as.numeric(mixer6_raw$end),
  f0 = f0_hz,
  f0_semitones = f0_semitones,
  vot = as.numeric(mixer6_raw$vot),
  vowel_duration = as.numeric(mixer6_raw$vdur),
  word_duration = as.numeric(mixer6_raw$wdur),
  speech_rate = as.numeric(mixer6_raw$spk_rate),
  pos = as.character(mixer6_raw$pos),
  syll = as.character(mixer6_raw$syll),
  type = as.character(mixer6_raw$type),
  session = as.integer(mixer6_raw$session)
)

mixer6 <- mixer6_clean
save(mixer6, file = "data/mixer6.rda", compress = "xz")
cat("Saved data/mixer6.rda:", nrow(mixer6), "rows\n")

cat("=== Phonetic datasets successfully processed and saved! ===\n")
