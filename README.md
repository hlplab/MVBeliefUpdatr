# MVBeliefUpdatr
`MVBeliefUpdatr` is a package that aims to facilitate the fitting, summarizing, and visualizing of univariate, multivariate, and cue-integration ideal observers and adaptors (as well as basic support for exemplar models). These models can be used, for example, to model speech perception and changes in speech perception as a response to recently experienced speech input. Stan programs are available that are designed to infer the categorization parameters (e.g., category priors, lapse rates) and the adaptation parameters (e.g., prior beliefs about category representations, incl. confidence parameters) based on behavioral data from participants in exposure-test experiments (e.g., visually- or lexically-guided perceptual learning experiments, distributional learning experiments, accent adaptation experiments, etc.). 


# Installation
You can install `MVBeliefUpdatr` from R:

```
devtools::install_github("hlplab/MVBeliefUpdatr")
```

That should install `MVBeliefUpdatr` and all of its dependencies. 


# How to cite this package
If you find this package useful, please cite the works in:

```
citation(package = "MVBeliefUpdatr")
```

Thank you!


# Examples
`MVBeliefUpdatr` code has now been used in a number of published studies. You can find some examples in: 

  * [Tan, Xie, & Jaeger (2021) *Frontiers in Psychology*](https://www.frontiersin.org/articles/10.3389/fpsyg.2021.676271/full)
  * [Xie, Buxo-Lugo, & Kurumada (2021) *Cognition*](https://www.sciencedirect.com/science/article/pii/S001002772100038X)
  * [Cummings & Theodore (2023), *Cognition*](https://www.sciencedirect.com/science/article/pii/S0010027723000380)
  * [Xie, Jaeger, & Kurumada (2023) *Cortex*](https://www.sciencedirect.com/science/article/pii/S0010945223001107) [[OSF](https://osf.io/q7gjp/)]
  * [Persson, Barreda, & Jaeger (2025) *JASA*](https://pubs.aip.org/asa/jasa/article/157/2/1458/3337250)
  * [Tan and Jaeger (2025) *Cognition*](https://www.sciencedirect.com/science/article/pii/S0010027725001350) [[OSF](https://osf.io/hxcy4)]

Several of these articles come with public OSF repositories with an R Markdown document that illustrates some of the ways in which `MVBeliefUpdatr` can be used in research on speech perception. 


# Acknowledgments
The stancode and some of the R code to this library were originally based on Dave Kleinschmidt's `beliefupdatr`, which implements Bayesian belief-updating for univariate Gaussian categories. Those models were used in Kleinschmidt & Jaeger (2011, 2012, 2015, 2016) and Kleinschmidt (2020). The stancode has further benefited tremendously from help by Paul Buerkner, Niko Huurre, and others amazing contributors on the [stan discourse forum](https://discourse.mc-stan.org/). Thanks for helpful pull requests go to Zach Burchill, Shawn Cummings, and Xin Xie.


# Reporting bugs
The package is still very much under development, and would benefit from your feedback. 

You can help improve this library by [Reporting a bug](https://github.com/hlplab/MVBeliefUpdatr/issues). Or even better, you can initiate a pull request, e.g., through GitHub's pull request feature](https://github.com/hlplab/MVBeliefUpdatr/pulls) or by using RStudio. To use RStudio, simply open the project with the cloned repository on your computer, go to the Git tab (usually one of the panels in the top right corner). Click "New Branch". 

When submitting pull requests, please choose a name for the new branch (e.g., `YOUR-NAME-plotting-patch`). Switch to that branch if it didn't happen automatically (by clicking on the pull down menu right next to the "New Branch" button; it usually says "master", but not should say the name of your new branch). Then edit the code the way you think is best (don't worry, you're not destroying anything). Then press push.

Your help is much appreciated.



