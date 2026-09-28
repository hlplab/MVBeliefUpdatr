# test-00-datasets.R
# Tests for bundled phonetic datasets: h95, pb52, swehvd, mixer6, and legacy ChodroffWilson2018.

test_that("h95 loads with standardized schema, Unicode IPA vowels, and valid values", {
  data("h95", package = "MVBeliefUpdatr")
  expect_s3_class(h95, "tbl_df")
  expect_equal(nrow(h95), 1668L)
  expect_equal(
    colnames(h95),
    c("speaker", "sex", "type", "vowel", "f0", "f1", "f2", "f3", "duration")
  )
  expect_equal(levels(h95$sex), c("female", "male"))
  expect_equal(levels(h95$type), c("men", "women", "boys", "girls"))
  expect_equal(nlevels(h95$vowel), 12L)
  expect_true(all(grepl("^/.+/$", levels(h95$vowel))))
  expect_true(all(c("/i/", "/ɪ/", "/ɛ/", "/æ/", "/ɑ/", "/ɔ/", "/ʊ/", "/u/", "/ʌ/", "/ɝ/", "/eɪ/", "/oʊ/") %in% levels(h95$vowel)))
  expect_true(all(h95$f0 > 0))
  expect_true(all(h95$f1 > 0))
  expect_true(all(h95$duration > 0))
})

test_that("pb52 loads with standardized schema, Unicode IPA vowels, and valid values", {
  data("pb52", package = "MVBeliefUpdatr")
  expect_s3_class(pb52, "tbl_df")
  expect_equal(nrow(pb52), 1520L)
  expect_equal(
    colnames(pb52),
    c("speaker", "sex", "type", "vowel", "repetition", "f0", "f1", "f2", "f3")
  )
  expect_equal(levels(pb52$sex), c("female", "male"))
  expect_equal(levels(pb52$type), c("men", "women", "children"))
  expect_equal(nlevels(pb52$vowel), 10L)
  expect_true(all(grepl("^/.+/$", levels(pb52$vowel))))
  expect_true(all(c("/i/", "/ɪ/", "/ɛ/", "/æ/", "/ʌ/", "/ɑ/", "/ɔ/", "/ʊ/", "/u/", "/ɝ/") %in% levels(pb52$vowel)))
  expect_true(all(pb52$f0 > 0))
  expect_true(all(pb52$f1 > 0))
  expect_true(all(pb52$repetition %in% c(1L, 2L)))
})

test_that("swehvd loads with standardized schema, Swedish IPA vowels, and valid values", {
  data("swehvd", package = "MVBeliefUpdatr")
  expect_s3_class(swehvd, "tbl_df")
  expect_equal(nrow(swehvd), 23685L)
  expect_equal(
    colnames(swehvd),
    c("speaker", "sex", "age", "vowel", "transcribed_vowel", "quantity",
      "word", "token", "trial", "location", "f0", "f1", "f2", "f3",
      "duration", "speech_rate", "unreliable_measurement")
  )
  expect_equal(levels(swehvd$sex), c("female", "male"))
  expect_true(all(swehvd$sex == "female"))
  expect_equal(nlevels(swehvd$speaker), 24L)
  expect_equal(nlevels(swehvd$vowel), 21L)
  expect_true(all(grepl("^/.+/$", levels(swehvd$vowel))))
  expect_equal(levels(swehvd$quantity), c("long", "short"))
  expect_true(all(swehvd$f0[!is.na(swehvd$f0)] > 0))
  expect_true(all(swehvd$duration > 0))
})

test_that("mixer6 loads with standardized schema, stop categories, and valid values", {
  data("mixer6", package = "MVBeliefUpdatr")
  expect_s3_class(mixer6, "tbl_df")
  expect_equal(nrow(mixer6), 96357L)
  expect_equal(
    colnames(mixer6),
    c("speaker", "sex", "stop", "stop_poa", "voicing", "word", "vowel",
      "trial", "filename", "start", "end", "f0", "f0_semitones", "vot",
      "vowel_duration", "word_duration", "speech_rate", "pos", "syll",
      "type", "session")
  )
  expect_equal(levels(mixer6$sex), c("female", "male"))
  expect_equal(levels(mixer6$stop), c("/b/", "/p/", "/d/", "/t/", "/g/", "/k/"))
  expect_equal(levels(mixer6$stop_poa), c("labial", "coronal", "dorsal"))
  expect_equal(levels(mixer6$voicing), c("voiced", "voiceless"))
  expect_equal(nlevels(mixer6$speaker), 180L)
  expect_true(all(!is.na(mixer6$vot)))
})

test_that("deprecated ChodroffWilson2018 dataset loads for backward compatibility", {
  data("ChodroffWilson2018", package = "MVBeliefUpdatr")
  expect_s3_class(ChodroffWilson2018, "tbl_df")
  expect_equal(nrow(ChodroffWilson2018), 65800L)
  expect_true(all(c("category", "VOT", "f0", "vowel_duration") %in% colnames(ChodroffWilson2018)))
})
