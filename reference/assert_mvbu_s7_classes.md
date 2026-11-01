# Assert That an Object Inherits from an MVBeliefUpdatr S7 Class

These helpers validate whether an object inherits from the S7 classes
introduced for category representations, templates, cognitive models,
Stan input objects, and Stan fit objects.

## Usage

``` r
assert_MVBU_CategoryRepresentation(x, msg = NULL)

assert_MVBU_CategoryRepresentationTemplate(x, msg = NULL)

assert_MVBU_CognitiveModel(x, msg = NULL)

assert_MVBU_Stanfit(x, msg = NULL)

assert_IdealAdaptorStanfit(x, msg = NULL)

assert_IdealAdaptorStanfitInput(x, msg = NULL)

assert_MVBU_Staninput(x, msg = NULL)

assert_IdealAdaptorStaninput(x, msg = NULL)
```

## Arguments

- x:

  Object to test.

- msg:

  Optional error message.

## Value

Invisibly `TRUE` if the assertion passes.

## Examples

``` r
# Assert that an object inherits from MVBU_CategoryRepresentation
rep <- example_niw_category_representation(n_cues = 2)
assert_MVBU_CategoryRepresentation(rep)

# Assert that an object inherits from MVBU_CognitiveModel
mod <- example_mvg_ideal_observer(n_cues = 2)
assert_MVBU_CognitiveModel(mod)

# Assert that an object inherits from IdealAdaptorStanfitInput
sfi <- example_ideal_adaptor_stanfit_input()
assert_IdealAdaptorStanfitInput(sfi)
```
