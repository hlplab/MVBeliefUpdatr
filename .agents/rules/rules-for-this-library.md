---
trigger: always_on
---

1. Document R code using `roxygen`. This is how package documentation is generated, and NAMESPACE and other derivative information is updated. A few conventions:

1a. Think about what help pages can be grouped together on a single help page (i.e., which functions share rdnames).

1b. Link relevant related functions using the seealso tag.

1c. Dependencies between different code files should be indicated using the @include tag at the top of the R file.

1d. Dependencies on external packages should be clearly marked for each function using @import tags.


2. Aim for transparency. Code should be human readable to make it easier to maintain. Avoid unnecessary complexity or convoluted code.


3. Aim for consistency in your code across functions (e.g., in function and argument naming). You can propose revisions to previously proposed code, if that improves consistency. A few conventions:

3a. Names for internal functions should start with ".". Their roxygen documentation should mark them as internal, and not list them in help overview pages.

3b. Deprecated functions should be marked as deprecated in a consistent way in the roxygen documentation. They should also evoke call the `lifecycle` package's deprecation warning with informative messages.

3c. Within an R file, functions should be grouped together based on the semantics (based on what they do and which object types they relate to), and ordered in the same way across files. 


4. Avoid redundancy. In particular, check what functions are already available and think about whether they could be reused---if necessary, with small changes or extensions---to achieve your goal before you create additional functions. 


5. Always develop tests for new code using the `testthat` package. 


