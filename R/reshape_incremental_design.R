#' Prepare long data from incremental exposure-test design for input to Stan
#'
#' Takes \code{data.frame} or \code{tibble} that contains the exposure and test data from an incremental
#' exposure-test design in long format, and prepares it for input to the \code{\link{ideal_adaptor_stanfit}}
#' Stan programs. This is done by pretending that each incremental test block (and its preceding exposure)
#' constitute a separate between-participant condition. Note that this does not capture the dependency
#' between test responses of participants in the same between-participant conditions, but such dependencies
#' are not modeled by current `MVBeliefUpdatr` Stan programs anyway (which do not include random effects by
#' participants).
#'
#' @param data Data frame or tibble to be sliced. Each row should be a single exposure or test observation.
#' @param group Character string indicating the name of the column that contains the information about
#'   the between-participant condition. (default: "Group")
#' @param phase Character string indicating the name of the column that contains the information about
#'   whether an observation is part of "exposure" or "test". This column must contain the values "exposure"
#'   and "test". Observation with other values will be ignored. (default: "Phase")
#' @param block Character string indicating the name of the column that contains the information about the
#'   incremental exposure and test blocks. Must be a factor with the levels indicating the order of the blocks.
#'   (default: "Block")
#' @param join_adjacent_test_blocks Logical indicating whether adjacent test blocks without intervening
#'   exposure blocks should be joined into a single test block. This will speed up \code{\link{fit_ideal_adaptor}}
#'   since there will be fewer conditions to iterate over but also means that the default plotting functions
#'   won't be able to plot the results of the different test blocks separately. (default: `FALSE`)
#' @param verbose Should verbose output be provided? (default: `FALSE`)
#'
#' @return A data frame or tibble in long format with a new column "ExposureGroup" that contains a unique
#'   label for each unique combination of `group` and `block`.
#'
#' @export
reshape_incremental_design_into_unique_exposure_test_combinations <- function(
    data,
    group = "Group",
    phase = "Phase",
    block = "Block",
    join_adjacent_test_blocks = FALSE,
    verbose = FALSE
) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  .assert_true(
    all(c(phase, group, block) %in% names(data)),
    msg = "The data must contain phase, group, and block columns."
  )
  .assert_true(
    all(c("exposure", "test") %in% unique(as.character(data[[phase]]))),
    msg = "The phase column must contain both 'exposure' and 'test' values."
  )
  if (!is.factor(data[[block]])) {
    data[[block]] <- factor(data[[block]])
  }

  exposure_blocks <- unique(as.character(data[[block]][data[[phase]] == "exposure"]))
  test_blocks <- unique(as.character(data[[block]][data[[phase]] == "test"]))
  .assert_true(
    !any(exposure_blocks %in% test_blocks),
    msg = "The levels of the block variable in the exposure phase must not overlap with those in the test phase. Please check your data."
  )

  if (verbose && length(setdiff(unique(as.character(data[[phase]])), c("exposure", "test"))) > 0) {
    message(
      paste("The following values in the", phase, "column are not recognized as exposure or test and thus removed:",
            paste(setdiff(unique(as.character(data[[phase]])), c("exposure", "test")), collapse = ", ")))
  }

  keep_rows <- data[[phase]] %in% c("exposure", "test")
  data <- data[keep_rows, , drop = FALSE]
  data[["..block_order"]] <- as.numeric(data[[block]])

  block_levels <- levels(data[[block]])
  phase_table <- unique(data[, c(phase, block, "..block_order"), drop = FALSE])
  phase_table <- phase_table[order(phase_table[["..block_order"]]), , drop = FALSE]
  phase_levels <- as.character(phase_table[[phase]])

  if (join_adjacent_test_blocks) {
    for (b in seq_len(length(block_levels) - 1)) {
      if (all(phase_levels[b:(b + 1)] == "test")) {
        if (verbose) {
          message("Joining adjacent test blocks ", block_levels[b], " and ", block_levels[b + 1], " into a single test block.")
        }

        block_levels[b] <- paste(block_levels[b], block_levels[b + 1], sep = "_")
        block_levels <- block_levels[-(b + 1)]

        data[["..block_order"]] <- ifelse(data[["..block_order"]] == b + 1, b, data[["..block_order"]])
        data[[block]] <- factor(
          ifelse(data[["..block_order"]] %in% c(b, b + 1), block_levels[b], as.character(data[[block]])),
          levels = block_levels
        )
      }
    }
  }

  if (verbose) {
    message("Inferred block order: ", paste(block_levels, collapse = ", "))
  }

  testblock_order <- sort(unique(as.numeric(data[["..block_order"]][data[[phase]] == "test"])))
  if (verbose) {
    message("Inferred test block order: ", paste(testblock_order, collapse = ", "))
  }

  rows <- list()
  group_levels <- unique(as.character(data[[group]]))
  idx <- 1L
  for (g in group_levels) {
    for (b in testblock_order) {
      subset <- data[
        data[[group]] == g & data[["..block_order"]] <= b & (data[["..block_order"]] == b | data[[phase]] != "test"),
        , drop = FALSE
      ]
      subset[["ExposureGroup"]] <- if (b == 1) "no exposure" else paste0("Group ", g, "_up to block ", block_levels[b])
      rows[[idx]] <- subset
      idx <- idx + 1L
    }
  }

  if (length(rows) == 0L) {
    df.new <- data.frame(ExposureGroup = character(0), stringsAsFactors = FALSE)
  } else {
    df.new <- do.call(rbind, rows)
  }

  keep_cols <- c("ExposureGroup", group, phase, block, setdiff(names(df.new), c("ExposureGroup", group, phase, block)))
  df.new[, keep_cols, drop = FALSE]
}
