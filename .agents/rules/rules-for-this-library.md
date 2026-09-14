---
trigger: always_on
---

1. Document R code using `roxygen`. This is how package documentation is generated, and NAMESPACE and other derivative information is updated. A few conventions:

1a. Think about what help pages can be grouped together on a single help page (i.e., which functions share rdnames).

1b. Link relevant related functions using the seealso tag.

1c. Dependencies between different code files should be indicated using the @include tag at the top of the R file.

1d. Dependencies on external packages should be clearly marked for each function using @import tags.

1e. Documentation of class objects should be marked following roxygen recommendations (which might differ for S7 and S3 classes).

1f. All functions/methods/classes should have roxygen documentation for all of their params/slots/etc. Defaults (explicit in the argument statements of the function or implicit in the code) should always be stated as part of the param description.


2. Aim for transparency. Code should be human readable to make it easier to maintain. Avoid unnecessary complexity or convoluted code.

2a. Keep non-roxygen comments contained in the original code. Do not remove them unless explicitly instructed to do so.



3. Aim for consistency in your code across functions (e.g., in function and argument naming). You can propose revisions to previously proposed code, if that improves consistency. A few conventions:

3a. Names for internal functions should start with ".". Their roxygen documentation should mark them as internal, and not list them in help overview pages.

3b. Within an R file, functions should be grouped together based on the semantics (based on what they do and which object types they relate to), and ordered in the same way across files. 

3c. Deprecated functions should be marked as deprecated in a consistent way in the roxygen documentation. They should also evoke call the `lifecycle` package's deprecation warning with informative messages.

3d. Always collect all deprecated functions in a file at the bottom of the file below a commented line "deprecated". Once all  functions in an R file are deprecated, rename that file to "deprecated-{original name}, following examples that already exist in the R folder.

3e. Deprecated S3 methods that have name-identical S7 methods should be completely removed, rather than to keep them as wrappers. If that means that an R file has not meaningful content anymore, delete that R file. Also make sure to remove any tests of those removed S3 methods.

3f. S7 functions should be defined in existing R files for S7-related methods / objects, or---when such files don't exist---in new files starting with "S7-...". Each R file should only contain semantically related code that is related to the name of the file (making future maintenance of the code easier).

3g. S7 functions should aim to reduce dependency on external packages. Some dependencies are ok, if it makes the code substantially more transparent without *relevant* loss of computational efficiency. In particular, the following are ok: dependency on rstan, loo, posterior.



4. Avoid redundancy. In particular, check what functions are already available and think about whether they could be reused---if necessary, with small changes or extensions---to achieve your goal before you create additional functions. 

4a. Whenever starting a new topic in a conversation, make sure that you read in existing functionality, including internal .is_, .assert_* and util functions. Use these functions if it keeps code lean. If changes to those functions would allow more elegant code, alert me to it.

4b. Do not silently introduce aliases without my explicit instruction. You can propose aliases to be reviewed by me. 

4c. When you note that highly similar code has been created in multiple places in the library, alert me to it and propose a unification of that functionality, as long as it does not conflict with rule 2 (transparency).



5. Aim for *relevant* computational efficiency. It is not important that functions that are typically applied once (e.g., print, summary, aping and transforming draws from a stanfit object, or updating models, should be kept efficient. 

5a. When you make choices for efficiency motivate them by comments in the code (non-roxygen comments).



6. Always develop tests for new code using the `testthat` package.

6a. Order tests by naming them test-XX-{test name}.R, where XX is a number. 

6b. For test names follow the naming scheme used for the main R files. E.g., if a test assess S7 functionality, {test name} should start with "S7-"; if a test assesses deprecated functions, {test name} should start with "deprecated-", followed by the old test name.

6c. Remove tests that refer to functions/methods/objects that do no longer exist.

6d. Unless there is a good reason, try to avoid duplicating existing functionality from the package in helper functions for the tests (see rule 5). E.g., don't write new read/write or fit functions for stanfit objects, when those already exist in the code---unless there is a good reason for it. In that case, include that reason on the comments in the test code, and make me aware of those reasons.