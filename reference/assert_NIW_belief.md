# Deprecated: assert_NIW_belief

**\[deprecated\]** `assert_NIW_belief()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[NIW_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
or S7 assertions instead.

## Usage

``` r
assert_NIW_belief(x, category = "category", verbose = F, strict = F)
```

## Arguments

- x:

  Object to test.

- category:

  Name of the category field to inspect.

- verbose:

  Logical. Whether to emit additional diagnostics while testing the
  object.

- strict:

  Logical. Whether to require the object to be an NIW belief strictly,
  rather than allowing an ideal adaptor.

## Value

Invisibly TRUE if the assertion passes.

## See also

[NIW_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
