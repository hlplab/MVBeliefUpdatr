# MVBeliefUpdatr
`MVBeliefUpdatr` is a package that aims to facilitate the fitting, summarizing, and visualizing of multivariate Gaussian ideal observers and adaptors. This models can be used, for example, to model speech perception and changes in speech perception as a response to recently experienced speech input. The package is still very much under development, and would benefit from your feedback. 


# Installation
You can install `MVBeliefUpdatr` from R:

```
devtools::install_github("hlplab/MVBeliefUpdatr")
```

That should install `MVBeliefUpdatr` and all of its dependencies. If some dependencies are *not* automatically installed, please let email
Florian Jaeger (<fjaeger@ur.rochester.edu>), so that I can add them to the dependencies. 


# Examples
`MVBeliefUpdatr` code was used in [Tan, Xie, & Jaeger (2021)](https://www.frontiersin.org/articles/10.3389/fpsyg.2021.676271/full), [Xie, Buxo-Lugo, & Kurumada (2021)](https://www.sciencedirect.com/science/article/pii/S001002772100038X), [Xie, Jaeger, & Kurumada (2023)](https://osf.io/q7gjp/), and [Tan and Jaeger (2025)](https://osf.io/hxcy4). The latter two of these articles come with public OSF repositories with an R Markdown document that illustrates some of the ways in which `MVBeliefUpdatr` can be used in research on speech perception.


# Acknowledgments
The stancode and some of the R code to this library were originally based on Dave Kleinschmidt's `beliefupdatr`, which implements Bayesian belief-updating for univariate Gaussian categories. Those models were used in Kleinschmidt & Jaeger (2011, 2012, 2015, 2016) and Kleinschmidt (2020). The stancode further benefitted tremendously from help by Paul Buerkner, Niko Huurre, and others amazing contributors on the [stan discourse forum](https://discourse.mc-stan.org/). Thanks for helpful pull requests go to Zach Burchill and Xin Xie.


# Reporting bugs
You can help improve this library by opening issues (click Issues tab above). Or even better, you can initiate a pull request either through the web interface or by using RStudio. To use RStudio, simply open the project with the cloned repository on your computer, go to the Git tab (usually one of the panels in the top right corner). Click "New Branch", choose a name for the new branch (e.g., "YOUR-NAME-plotting-patch"; you don't have to select the remote or any other option--- just stick with the defaults). Switch to that branch if it didn't happen automatically (by clicking on the pull down menu right next to the "New Branch" button; it usually says "master", but not should say the name of your
new branch). Then edit the code the way you think is best (don't worry, you're not destroying anything). Then press push.

Your help is much appreciated.



