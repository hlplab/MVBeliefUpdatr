---
trigger: glob
globs: *.R, *.Rmd, *.qmd, *.py, *.ipynb, *.stan
---

## Coding guidelines

Accessibility & Verifiability: Produce clearly structured and well-documented code that is easy to read and understand. Use consistent naming conventions (e.g., snake_case for functions and variables), modular design, and clear function interfaces. Include comments and docstrings to explain purpose, logic, inputs, outputs, and any assumptions made in the code. Ensure that all code is verifiable by providing test cases, example inputs/outputs, and clear instructions for running the code. 

Often some light (or more developed) class structure (e.g., S3/S4 in R, or Python classes) is appropriate for clarity and accessibility, but avoid over-engineering. Use classes when they provide clear benefits in terms of organization, reusability, and maintainability of the code.

Elegance & clarity: Collect repeated code into functions or methods to improve readability and maintainability. Declare global variables and constants in a single location to avoid scattering them throughout the codebase. This includes global settings for visualizations (e.g., set_theme() in ggplot2) and other global environment options. Avoid hardcoding values in functions or elsewhere deeply hidden in the code. Avoid using global state or hidden dependencies that can lead to unpredictable behavior.

Correctness: It is better to fail with informative messages / warnings / errors. Avoid code that fails quietly or produces incorrect results without notifying the user. Validate inputs and outputs, and include unit tests to ensure code correctness. Use assertions to check for expected conditions and raise informative errors when assumptions are violated.

No 'hacks': Do not quietly 'fix' failures; instead, raise informative exceptions that clearly indicate the nature of the problem and how to resolve it (e.g., when encountering overflow or other numerical issues, do not quietly enforce non-zero values; instead, allow users to specify them and document that behavior). Never make figures, tables or other data summaries work by silently dropping data points or changing the underlying data (or its display) without explicit user consent. Hardcoded 'fallback' options similarly should be avoided, and only ever be used after first explaining the rationale to the user and getting their consent.

When the user states that something "doesn't make sense"  (or similar statements), this should prompt careful investigation and debugging, not a quick fix that hides the underlying issue. 

Reproducibility: Every script must be deterministic. Always set seeds for random number generators (e.g., set.seed(42)), but do so transparently (e.g., not hard-coded within functions but by handing functions seed parameters). Document all dependencies and package versions in a `requirements.txt` file, and include a README section on how to set up the development environment.

Efficiency: Optimize code for performance unless this conflicts with correctness, accessibility or reproductibility, especially for large datasets or computationally intensive tasks. Use vectorized operations and efficient data structures where possible. Consider parallelization or (local) distributed computing for tasks that can be parallelized.

Coding format: Follow established coding conventions for the chosen programming language (e.g., the Tidyverse style guide for R, PEP 8 for Python).

Clean up before commit: Before committing staged changes, review and refactor code to remove unused functions, variables, and imports. Ensure that the codebase that is about to be committed remains clean, organized, and free of clutter (no need to check unstages changes). Avoid leaving commented-out code or stale temporary scripts in the codebase; instead, use version control to track changes and revert if necessary. Explicitly propose changes and get the user's consent before removing or changing any code that is in use.

Analysis & Modeling: Use appropriate statistical and computational modeling techniques for the research questions at hand. Ensure that models are well-specified, validated, and interpreted correctly. Provide clear explanations of model assumptions, limitations, and implications. Use simulation studies to validate model performance and robustness. Use Stan to write novel generative models, and use `brms` or `rstanarm` for Bayesian modeling when appropriate. Use `lme4` for frequentist mixed-effects models. Avoid using black-box machine learning models without clear interpretability and justification.

Visualization: Use clear and informative visualizations to communicate results. Ensure that plots are properly labeled, with appropriate titles, axis labels, legends, and units. Avoid cluttered or misleading visualizations. Use `ggplot2` (preferred) or `plotly` for all plots. 

Equations: Write out mathematical equations in LaTeX format.

Accessible output format: Output results in standard, portable formats like .csv or .tab. The only exception are outputs that can be reproduced from the code in the repo (e.g., for complex computations). Those can be stored in binary formats (e.g., .rds, .pkl) with clear documentation on how to regenerate them."

Output locations: Figures should be stored in a `figures/` directory, and tables in a `tables/` directory. Use descriptive names for files and directories to facilitate navigation and understanding. All raw and augmented data should be stored in a `data/` directory. Models should be stored in a `models/` directory, and scripts in a `code/` or `scripts/` directory (use only one). Use subdirectories to organize files by analysis or experiment when appropriate.

All temporary or intermediate files should be stored in a `tmp/` directory, with subfolders mimicking the organization of the repo, so that temporary scripts pertaining to code in e.g., subfolder code/a/b/c will be found in tmp/a/b/c/. 
