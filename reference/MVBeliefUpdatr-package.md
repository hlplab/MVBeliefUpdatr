# MVBeliefUpdatr: Fitting, Summarizing, and Visualizing of Gaussian Ideal Observers and Adaptors

Model, fit, and visualize incremental Bayesian categorization and belief
updating for uni- and multivariate Gaussian, cue integration, and
exemplar-based cognitive models using 'S7' class architecture. Provides
tools for simulating normative ideal observers and ideal adaptors under
conjugate priors including Normal-Inverse-Chi-squared (NIX),
Multivariate Normal-Inverse-Chi-squared (MNIX), and
Normal-Inverse-Wishart (NIW) distributions. Supports Bayesian parameter
inference via Markov Chain Monte Carlo (MCMC) sampling using 'Stan',
model comparison via 'loo', and rich visualizations with 'ggplot2' and
'patchwork'.

## Overview

`MVBeliefUpdatr` provides a unified, object-oriented framework for
Bayesian ideal observers, ideal adaptors, incremental belief updating,
and hierarchical Bayesian cognitive modeling.

The package supports continuous distributional representations
(univariate Gaussian, multivariate Gaussian, multi-univariate Gaussian
cue integration), conjugate prior distributions over category parameters
(Normal-Inverse- Chisquare, Normal-Inverse-Wishart), and exemplar-based
representations.

Using analytical conjugate updating or hierarchical Bayesian inference
in Stan via rstan, `MVBeliefUpdatr` facilitates simulation, parameter
estimation, categorization prediction, and visualization of speech
perception and perceptual adaptation experiments.

## Package Vignettes

For detailed guides and worked examples, see the package vignettes:

- [`vignette("s7-class-structure-and-workflows")`](https://hlplab.github.io/MVBeliefUpdatr/articles/s7-class-structure-and-workflows.md):

  S7 Class Architecture, Constructors, and Workflows

- [`vignette("visualizing-models-and-categories")`](https://hlplab.github.io/MVBeliefUpdatr/articles/visualizing-models-and-categories.md):

  Visualizing Models and Categories

- [`vignette("fitting-and-working-with-stanfit-models")`](https://hlplab.github.io/MVBeliefUpdatr/articles/fitting-and-working-with-stanfit-models.md):

  Fitting and Working with MVBeliefUpdatr Stanfit Models

## S7 Core Class Architecture

All model objects are organized under an explicit, compositional S7
class hierarchy:

- `MVBU_CategoryRepresentation`:

  Abstract base class for individual category representations,
  encapsulating category parameters, cue labels, and likelihood
  functions.

- `MVBU_CategoryRepresentationTemplate`:

  Validated collection of category representations forming a perceptual
  inventory or template.

- `MVBU_CognitiveModel`:

  Decision model pairing a category template with decision behavior
  (decision rule, category priors, lapse rate, lapse bias, and
  perceptual noise).

- `MVBU_IdealObserver`:

  Subclass of `MVBU_CognitiveModel` for observers with fixed generative
  category distributions.

- `MVBU_IdealAdaptor`:

  Subclass of `MVBU_CognitiveModel` for adaptors with conjugate prior
  beliefs that update incrementally from exposure.

- `MVBU_Stanfit`:

  Wrapper for `rstan` `stanfit` objects inferring prior beliefs, lapse
  parameters, and noise from behavioral test data.

- `MVBU_StanfitInput`:

  Container for exposure data, test data, priors, and transformation
  functions used to prepare Stan inputs.

- `MVBU_Staninput`:

  Validated input data list passed to Stan models.

## Supported Model Families

`MVBeliefUpdatr` provides concrete S7 classes and constructors across
seven families:

- `UVG` / `NIX`:

  Univariate Gaussian ideal observers
  ([UVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-uvg.md))
  and Normal-Inverse-Chisquare ideal adaptors
  ([NIX_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md)).

- `MUVG` / `MNIX`:

  Multi-univariate Gaussian ideal observers
  ([MUVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-muvg.md))
  and independent Normal-Inverse-Chisquare cue integration adaptors
  ([MNIX_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md)).

- `MVG` / `NIW`:

  Multivariate Gaussian ideal observers
  ([MVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md))
  and Normal-Inverse-Wishart ideal adaptors
  ([NIW_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)).

- `EXEMPLAR`:

  Exemplar-based categorization models
  ([Exemplar_Model](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md))
  with kernel density estimation over stored exemplars.

## Acknowledgments

Belief-updating formulations build upon conjugate Bayesian theory
(Murphy 2012) , phonetic sliding template models (Nearey and Assmann
2007) , Dave Kleinschmidt's `BeliefUpdatr` (Kleinschmidt and Jaeger
2011) , (Kleinschmidt and Jaeger 2012) , (Kleinschmidt and Jaeger 2015)
, (Kleinschmidt and Jaeger 2016) , and unsupervised adaptation modeling
(Yan and Jaeger 2018) .

Contributions and feedback from Zach Burchill, Anna Persson, and Xin Xie
are gratefully acknowledged.

## References

Kleinschmidt DF, Jaeger TF (2011). “A Bayesian belief updating model of
phonetic recalibration and selective adaptation.” In *ACL Workshop on
Cognitive Modeling and Computational Linguistics*.  
  
Kleinschmidt DF, Jaeger TF (2012). “A continuum of phonetic adaptation :
Evaluating an incremental belief-updating model of recalibration and
selective adaptation.” In *Proceedings of the 34th Annual Meeting of the
Cognitive Science Society (CogSci12)*, 605-610.  
  
Kleinschmidt DF, Jaeger TF (2015). “Robust speech perception: Recognize
the familiar, generalize to the similar, and adapt to the novel.”
*Psychological Review*, **122**(2), 148-203. ISSN 0033-295X.
[doi:10.1037/a0038695](https://doi.org/10.1037/a0038695) .
<http://www.ncbi.nlm.nih.gov/pubmed/25844873>.  
  
Kleinschmidt DF, Jaeger TF (2016). “What do you expect from an
unfamiliar talker?” In Trueswell JC, Papafragou A, Grodner D, Mirman D
(eds.), *Proceedings of the 38th Annual Meeting of the Cognitive Science
Society*.  
  
Murphy KP (2012). *Machine learning: a probabilistic perspective*. MIT
Press.  
  
Nearey TM, Assmann PF (2007). “Probabilistic 'sliding template' models
for indirect vowel normalization.” In Solé M, Beddor PS, Ohala M (eds.),
*Experimental Approaches to Phonology*, 246–269. Oxford University
Press, Oxford.  
  
Yan S, Jaeger TF (2018). “Comparing models of unsupervised adaptation in
speech perception.” In *Proceedings of the AMLaP conference*.

## See also

Useful links:

- <https://github.com/hlplab/MVBeliefUpdatr>

- Report bugs at <https://github.com/hlplab/MVBeliefUpdatr/issues>

## Author

**Maintainer**: T. Florian Jaeger <fjaeger@ur.rochester.edu>
([ORCID](https://orcid.org/0000-0002-1158-7308))

Authors:

- T. Florian Jaeger <fjaeger@ur.rochester.edu>
  ([ORCID](https://orcid.org/0000-0002-1158-7308))
