# Model Family Class

Encapsulates the mathematical primitives, moment extraction, and naming
conventions for a distributional model family in `MVBeliefUpdatr`.

Registers a model family in the package registry. Accepts either an
MVBU_Family object directly (recommended) or a family name with
associated class names and Stan model configuration.

## Usage

``` r
MVBU_Family(
  family_names = list(),
  density_function = function() NULL,
  sampling_function = function() NULL,
  expected_moments = function() NULL,
  marginal_moments = NULL,
  sufficient_statistics_function = NULL,
  stanmodel = NULL
)

new_mvbu_family(
  name_short,
  name_long = NULL,
  name_details = NULL,
  density_function,
  sampling_function,
  expected_moments,
  marginal_moments = NULL,
  sufficient_statistics_function = NULL,
  stanmodel = NULL
)

register_model_family(
  family,
  category_representation_class = NULL,
  cognitive_model_class = NULL,
  stanmodel = NULL,
  overwrite = TRUE
)
```

## Arguments

- family_names:

  A list containing `name_short`, `name_long`, and `name_details`.
  Default: [`list()`](https://rdrr.io/r/base/list.html).

- density_function:

  Density function `function(x, rep, log = FALSE)`.

- sampling_function:

  Sampling function
  `function(rep, n = 1L, with_replacement = TRUE, ...)`.

- expected_moments:

  Function `function(rep)` returning `list(mu, Sigma, weights)`.

- marginal_moments:

  Optional function returning marginal predictive moments. Default:
  `NULL`.

- sufficient_statistics_function:

  Optional function returning sufficient statistics. Default: `NULL`.

- stanmodel:

  Optional list bundling `stanmodel_name` and `stanfit_class` for Stan
  model families, or `NULL`. Used when `family` is a character string.
  Default: `NULL`.

- name_short:

  Short uppercase identifier for the family (e.g., "UVG").

- name_long:

  Long descriptive name. Default: `name_short`.

- name_details:

  Optional list specifying tier class/type names:
  `CategoryRepresentation_name`, `CategoryTemplate_name`,
  `CognitiveModel_name`, `CognitiveModelList_name`. Default: `NULL`.

- family:

  An MVBU_Family object, or character string specifying the model family
  short name.

- category_representation_class:

  Name of category-representation class. Used when `family` is a
  character string. Default: `NULL`.

- cognitive_model_class:

  Name of cognitive-model class. Used when `family` is a character
  string. Default: `NULL`.

- overwrite:

  Logical indicating whether to overwrite an existing registration for
  this family. Default: `TRUE`.

## See also

`new_mvbu_family()`, `register_model_family()`,
[`get_model_family()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_model_family.md)

## Examples

``` r
# Inspect a registered model family
reg <- get_model_family_registration("UVG")
reg$names$name_long
#> [1] "Univariate Gaussian Ideal Observer"

# Create a custom minimal family object
dummy_fam <- new_mvbu_family(
  name_short = "CUSTOM",
  density_function = function(x, rep, log = FALSE) rep(0, nrow(x)),
  sampling_function = function(rep, n = 1L, ...) matrix(0, n, 1),
  expected_moments = function(rep) list(mu = 0, Sigma = matrix(1, 1, 1))
)
dummy_fam
#> <MVBeliefUpdatr::MVBU_Family>
#>  @ family_names                  :List of 3
#>  .. $ name_short  : chr "CUSTOM"
#>  .. $ name_long   : chr "CUSTOM"
#>  .. $ name_details:List of 4
#>  ..  ..$ CategoryRepresentation_name: chr "CUSTOM_CategoryRepresentation"
#>  ..  ..$ CategoryTemplate_name      : chr "CUSTOM_CategoryRepresentationTemplate"
#>  ..  ..$ CognitiveModel_name        : chr "CUSTOM_CognitiveModel"
#>  ..  ..$ CognitiveModelList_name    : chr "CUSTOM_CognitiveModelList"
#>  @ density_function              : function (x, rep, log = FALSE)  
#>  @ sampling_function             : function (rep, n = 1L, ...)  
#>  @ expected_moments              : function (rep)  
#>  @ marginal_moments              : NULL
#>  @ sufficient_statistics_function: NULL
#>  @ stanmodel                     : NULL
```
