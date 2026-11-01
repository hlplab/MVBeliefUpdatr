# Singular and Plural Naming Conventions in MVBeliefUpdatr

## 1. Principles and Core Conventions

To ensure transparent, predictable, and clean interfaces across `MVBeliefUpdatr`, function and argument names adhere to systematic singular vs. plural semantic rules:

### A. Arguments: Single-Entity Specifiers (Singular)
Use **singular** nouns when an argument identifies a single column name, a single scalar setting, or an individual attribute:
- **Dataframe Column Identifiers**:
  - `category`: Single column name storing category identity (e.g., `category = "category"`).
  - `group`: Single column name storing grouping identity (e.g., `group = "group"`).
  - `group_unique`: Single column name storing unique group identity (e.g., `group_unique = "group_unique"`).
  - `response`: Single column name storing test response category (e.g., `response = "response_category"`).
  - `slice_cue`: Single cue name along which to slice higher-dimensional spaces (e.g., `slice_cue = "F3"`).
- **Configuration & Behavioral Options**:
  - `decision_rule`: Single decision rule algorithm (`"criterion"`, `"proportional"`, `"sampling"`).
  - `noise_treatment`: Single treatment of perceptual noise (`"no_noise"`, `"sample"`, `"marginalize"`).
  - `lapse_treatment`: Single treatment of lapses (`"no_lapses"`, `"lapses"`).
  - `statistic`: Single or multiple statistic selector (`statistic = c("n", "mean", "css", "uss", "cov")`).

### B. Arguments: Sets, Vectors & Dimension Lists (Plural)
Use **plural** nouns when an argument accepts multiple values, dimensions, or filters:
- **Dimension Names & Vectors**:
  - `cues`: Character vector naming continuous cue dimensions (e.g., `cues = c("F1", "F2")`). Even if only 1 cue is provided, the argument remains `cues`.
  - `categories`: Character vector of category labels to filter or evaluate (e.g., `categories = c("b", "d", "g")`).
  - `groups`: Character vector of group labels to filter or evaluate (e.g., `groups = c("native", "nonnative")`).
  - `response_categories`: Character vector of response category labels to filter on test data.
  - `pars`: Character vector of parameter names to filter (e.g., `pars = c("mu", "Sigma")`).
  - `slices`: Numeric vector of values along the third cue to slice at.
  - `levels`: Numeric vector of probability / density contour levels.
  - `limits`: List or numeric vector of boundary limits across dimensions.
  - `step_labels`: Character vector of custom labels for sequential update steps.
  - `observations`: Matrix or data frame of multiple empirical cue observations.

### C. Quantity & Count Parameters
- **Observation Counts**: `n_samples` (integer count of data rows/points).
- **Exemplar Counts**: `n_exemplars` (integer count of exemplars to sample for plotting).
- **MCMC Draw Counts**: `ndraws` (established Stan/posterior/brms standard naming).
- **Sample Size / Multiplicity**: `n` (standard R scalar count in generative sampling, e.g., `sample_observations(x, n = 100)`).

### D. Function & Generic Naming
- **Extractor Functions Returning Collections**:
  - Plural when returning vectors/lists of attributes:
    - `get_cue_labels()`, `get_category_labels()`, `get_group_labels()`, `get_labels()`
    - `get_parameter_names()`, `get_parameters()`
    - `get_category_representations()`
    - `get_draws()`
- **Extractor Functions Returning Single Structures or Scalars**:
  - Singular when returning an individual object, scalar, or single closure:
    - `get_category_template()`
    - `get_category_prior()`
    - `get_lapse_rate()`
    - `get_lapse_bias()`
    - `get_noise()`
    - `get_noise_treatment()`, `get_lapse_treatment()`
    - `get_category_likelihood_function()`, `get_category_posterior_function()`
    - `get_transform_function()`, `get_untransform_function()`
    - `get_stanfit()`, `get_staninput()`
    - `get_metadata()`, `get_model_family()`, `get_model_type()`, `get_representation_type()`
- **Generative Sampling**:
  - `sample_observations()`: Plural (samples one or more observations). Deprecates legacy singular `sample_observation()`.
- **Visualization Generics**:
  - `plot_categories()`: Plural (plots category distributions).
  - `plot_categorization_function()`: Singular function (plots the categorization decision surface).
  - `plot_parameters()`, `plot_parameter_correlations()`, `plot_parameters_pairwise()`: Plural (plots parameter estimates and their joint relations).
  - `plot_sample()`, `plot_exposure_sample()`, `plot_test_sample()`: Singular mass noun (`sample` referring to empirical sample data).
  - `plot_model_updates()`: Plural (plots sequential updates).
  - `plot_diagnostics()`: Plural (plots Stan MCMC diagnostics).

---

## 2. Audit Table of Exported Functions & Arguments

| Function / Generic | Argument Name | Form | Rationale & Consistency Notes |
| :--- | :--- | :--- | :--- |
| `get_data` / `get_exposure_data` / `get_test_data` | `groups` | Plural | Filters subset of groups |
| `get_data` / `get_exposure_data` | `categories` | Plural | Filters subset of categories |
| `get_data` / `get_test_data` | `response_categories` | Plural | Filters subset of response categories |
| `get_data` / `get_exposure_data` / `get_test_data` | `n_samples` | Plural | Subsample row count |
| `get_data` / `get_exposure_data` / `get_test_data` | `original_names` | Plural | Logical flag whether to restore original variable names |
| `new_ideal_adaptor_stanfit_input` | `cues` | Plural | Vector of column names for cue dimensions |
| `new_ideal_adaptor_stanfit_input` | `category` | Singular | Column name for category in exposure data |
| `new_ideal_adaptor_stanfit_input` | `response` | Singular | Column name for response in test data |
| `new_ideal_adaptor_stanfit_input` | `group` | Singular | Column name for group in data |
| `new_ideal_adaptor_stanfit_input` | `group_unique` | Singular | Column name for unique group |
| `plot_categories` | `cues` | Plural | Vector of 1-3 cue dimensions to plot |
| `plot_categories` | `categories` | Plural | Vector of categories to plot |
| `plot_categories` | `groups` | Plural | Vector of groups to plot |
| `plot_categories` | `levels` | Plural | Vector of probability/density levels |
| `plot_categories` | `limits` | Plural | Axis limits |
| `plot_categories` | `slices` | Plural | Coordinate values along 3rd cue |
| `plot_categories` | `slice_cue` | Singular | Single cue name along which to slice |
| `plot_categories` | `n_exemplars` | Plural | Max number of exemplars to sample |
| `plot_categories` | `ndraws` | Plural | Number of MCMC posterior draws |
| `plot_sample` | `sample` | Singular | Subsets to `"exposure"`, `"test"`, or `c("exposure", "test")` |
| `plot_sample` | `densities` | Plural | *Note*: Specifies density estimation method (`"gaussian"` or `"kernel"`). Named plural for historical parity with density layers. |
| `plot_categorization_function` | `categories` | Plural | Categories for which response surface is evaluated |
| `plot_categorization_function` | `decision_rule` | Singular | Decision rule algorithm |
| `plot_parameters` | `pars` | Plural | Vector of parameter names |
| `plot_parameters` | `groups` | Plural | Vector of groups |
| `evaluate_model` | `response_category` | Singular | Vector of trial-by-trial response categories (each observation has 1 category) |
| `evaluate_model` | `method` | Singular | Evaluation metric (or vector of metrics) |
| `evaluate_model` | `decision_rule` | Singular | Decision rule |
| `sample_observations` | `n` | Singular | Count of rows to sample |
| `sample_observations` | `with_replacement` | Singular | Logical flag |
| `sample_observations` | `randomize_order` | Singular | Logical flag |
| `get_original_variable_names` | `variable` | Singular | Selected variable name to query |
| `get_expected_category_statistic` | `statistic` | Singular | Statistic name(s) to compute (`"n"`, `"mean"`, `"cov"`, etc.) |
| `get_exposure_category_statistic` | `statistic` | Singular | Statistic name(s) to compute |

---

## 3. Exceptions and Intentional Distinctions

1. **`densities` in `plot_sample`**:
   - `densities` is plural (`densities = c("gaussian", "kernel")`), whereas `decision_rule` and `method` are singular. This is intentional: it specifies the family of empirical density estimators to fit across categories.
2. **`response` vs `response_category` vs `response_categories`**:
   - In `new_ideal_adaptor_stanfit_input()`, `response` (singular) names the raw column in the user's test dataframe.
   - In `evaluate_model()`, `response_category` (singular) refers to the vector of actual observed trial responses (one response per trial).
   - In `get_data()` and `get_test_data()`, `response_categories` (plural) refers to the filter vector of category levels to keep.
3. **`ndraws` vs `n_samples`**:
   - `ndraws` is preserved without an underscore to adhere to the `posterior`, `brms`, and `bayesplot` ecosystem standard.
   - `n_samples` and `n_exemplars` use standard snake_case for observations and exemplars within `MVBeliefUpdatr`.
