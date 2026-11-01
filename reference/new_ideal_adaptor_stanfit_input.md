# Construct a streamlined ideal-adaptor Stanfit input object

This function composes the exposure and test data in a compact form that
can be used as input to the ideal-adaptor Stan models. It accepts empty
exposure data as long as test data are present, and returns an S7
fit-input object with the prepared data (in both original and
transformed spaces) and transformation metadata.

## Usage

``` r
new_ideal_adaptor_stanfit_input(
  exposure,
  test,
  cues,
  category = "category",
  response = "response_category",
  group = "group",
  group_unique = NULL,
  control = list(transform_type = "standardize"),
  ...
)
```

## Arguments

- exposure:

  A data frame or tibble with exposure data. Each row is assumed to
  contain one observation.

- test:

  A data frame or tibble with test data. Each row is assumed to contain
  one observation.

- cues:

  Names of columns with cue values. These columns must be present in
  both `exposure` and `test`.

- category:

  Name of the column in `exposure` that stores the category information.
  (default: "category")

- response:

  Name of the column in `test` that stores the response information.
  (default: "response_category")

- group:

  Name of the column in `exposure` and `test` that stores the grouping
  information. (default: "group")

- group_unique:

  Optional grouping column used to collapse repeated exposure
  conditions.

- control:

  A list of control parameters for the constructor, such as
  `list(transform_type = "standardize")` or produced by
  [`control_staninput()`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitInput.md).

- ...:

  Additional arguments passed to methods or deprecated arguments.

## Value

An object of class `IdealAdaptorStanfitInput` containing the prepared
data (original and transformed) and transform metadata.

## Details

Exposure and test data are checked for the required cue, group,
category, and response columns. The exposure and test data are
transformed to facilitate model-fitting. This transformation is undone
in the stanfit object returned by the Stan model, so that the model
parameters are returned in the original cue space.

It is important to use `group` to identify individuals that had a
specific exposure (or no exposure at all) and specific test trials. You
should not use `group` to identify exposure conditions. Setting `group`
to an exposure condition results in exposure that concatenates the
exposure observations from all subjects in that condition. Typically,
this is not what users intend, as it models exposure to the combination
of exposure tokens across all subjects, rather than exposure to one set
of those exposure tokens. To achieve this intended outcome, use
`group_unique` to identify groups with identical exposure. This will
correctly use only one unique instance of the observations that any
level of `group` receives during exposure.

## Examples

``` r
exp_df <- data.frame(
  group = "g1",
  category = c("b", "p"),
  vot = c(0, 50)
)
test_df <- data.frame(
  group = "g1",
  response_category = c("b", "p"),
  vot = c(10, 40)
)
inp <- new_ideal_adaptor_stanfit_input(
  exposure = exp_df,
  test = test_df,
  cues = "vot"
)
inp
#> <IdealAdaptorStanfitInput>
#>   Categories (2): b, p
#>   Cues (1): vot
#>   Groups (1): g1
#>   Data: 4 rows
```
