make_vowel_test_data <- function() {
  category_labels <- c("AA", "AE", "AH", "AO", "EH", "ER", "IH", "IY", "OW", "UH", "UW")
  data <- expand.grid(
    vowel = category_labels,
    token = seq_len(8L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  category_index <- match(data$vowel, category_labels)
  data$vowel <- factor(data$vowel, levels = category_labels)
  data$F1 <- 400 + 35 * category_index + 3 * data$token
  data$F2 <- 1800 - 25 * category_index + data$token^2
  data$F3 <- 2500 + 15 * category_index + data$token^3 / 10
  tibble::as_tibble(data)
}