# Prepare long data from incremental exposure-test design for input to Stan

Takes a `data.frame` or `tibble` that contains exposure and test data
from an incremental exposure-test design in long format, and prepares it
for input to the
[`ideal_adaptor_stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/deprecated-functions.md)
Stan programs. This is done by treating each incremental test block (and
its preceding exposure) as a separate between-participant condition.

## Usage

``` r
reshape_incremental_design_into_unique_exposure_test_combinations(
  data,
  group = "Group",
  phase = "Phase",
  block = "Block",
  join_adjacent_test_blocks = FALSE,
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame or tibble to be sliced. Each row should be a single
  exposure or test observation.

- group:

  Character string indicating the name of the column that contains the
  information about the between-participant condition. (default:
  `"Group"`)

- phase:

  Character string indicating the name of the column that contains the
  information about whether an observation is part of `"exposure"` or
  `"test"`. This column must contain the values `"exposure"` and
  `"test"`. Observations with other values will be ignored. (default:
  `"Phase"`)

- block:

  Character string indicating the name of the column that contains the
  information about the incremental exposure and test blocks. Must be a
  factor with levels indicating the order of the blocks. (default:
  `"Block"`)

- join_adjacent_test_blocks:

  Logical indicating whether adjacent test blocks without intervening
  exposure blocks should be joined into a single test block. (default:
  `FALSE`)

- verbose:

  Should verbose output be provided? (default: `FALSE`)

## Value

A data frame or tibble in long format with a new column
`"ExposureGroup"` that contains a unique label for each unique
combination of `group` and `block`.

## Examples

``` r
df <- data.frame(
  Group = "G1",
  Phase = c("exposure", "test", "exposure", "test"),
  Block = factor(c("E1", "T1", "E2", "T2"),
                 levels = c("E1", "T1", "E2", "T2")),
  vot = c(10, 20, 15, 25),
  response = c(NA, "b", NA, "p")
)
reshape_incremental_design_into_unique_exposure_test_combinations(df)
#>              ExposureGroup Group    Phase Block vot response ..block_order
#> 1  Group G1_up to block T1    G1 exposure    E1  10     <NA>             1
#> 2  Group G1_up to block T1    G1     test    T1  20        b             2
#> 11 Group G1_up to block T2    G1 exposure    E1  10     <NA>             1
#> 3  Group G1_up to block T2    G1 exposure    E2  15     <NA>             3
#> 4  Group G1_up to block T2    G1     test    T2  25        p             4
```
