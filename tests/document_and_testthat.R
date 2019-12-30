# use_test("compose-input-for-stan.R")
# use_test("get-information-from-stanfit.R")
document()
devtools::test()


source("./R/class.R")
source("./R/NIW.R")
source("./R/get-information-from-stanfit.R")
source("./R/visualize-fit.R")

# Load data
fit = readRDS(file = "./tests/test models/IBBU_mv_fit_3 exposure groups_2 categories_2 cues_Drouin et al 2016.rds")
fit.input = readRDS(file = "./tests/test models/DATA_mv_fit_3 exposure groups_2 categories_2 cues_Drouin et al 2016.rds")
  
# class(fit) = "mv_ibbu_stanfit"
categories = c("s", "sh")
# groups = c("SH-BIAS", "SS-BIAS")
groups = c("SH-BIAS", "SS-BIAS", "Control")
cues = c("PC1", "PC2")
d = tidyr::crossing(
  category = factor(categories, levels = categories),
  group = factor(groups, levels = groups), 
  cue = cues,
  cue2 = cues
)

fit %<>%
  recover_types(d) 
fit %>%
  add_ibbu_draws(wide = F, which = "both", draws = c(10, 5), nest = F) -> f
fit %>%
  add_ibbu_draws(wide = F, which = "both", draws = c(10, 5), nest = T) -> g



source("./R/NIW.R")
get_expected_mu(g, "sh", "prior")
get_expected_sigma(g, "sh", "prior")
get_expected_sigma(g, c("s","sh"), c("prior", "control"))
get_expected_sigma(g, c("s","sh"), c("prior", "Control"))
get_expected_category_statistic(g, c("s","sh"), c("prior", "Control"), c("mu", "Sigma"))

source("./R/visualize-fit.R")
group.colors = c("darkgray", "blue", "red", "black")
plot_ibbu_parameters(fit, which = "both", n.draws = 2)
plot_ibbu_parameters(fit, which = "both", n.draws = 5,
                     group.colors = group.colors)
plot_ibbu_test_categorization(fit, fit.input,
                              which = "both", summarize = T, n.draws = 10,
                              group.colors = group.colors)
# should throw error: 
# fit %>%
#   add_ibbu_draws(wide = F, which = "both", draws = c(10, 5), nest = F) -> f
# plot_expected_ibbu_categories_2D(f) 
fit %>%
  add_ibbu_draws(wide = F, which = "both", draws = c(10, 5), nest = T) -> g
get_expected_category_statistic(g)
plot_expected_ibbu_categories_2D(g, type = "contour")
plot_expected_ibbu_categories_contour2D(g, fit.input)
plot_expected_ibbu_categories_2D(g, fit.input, type = "density", resolution = 5)









# Not for test (slow! (at least the 2D density one))
fit %>%
  add_ibbu_draws(wide = F, which = "both", nest = T) -> g
plot_expected_ibbu_categories_2D(g, fit.input)
plot_expected_ibbu_categories_density2D(g, fit.input, resolution = 100)
