plot_ibbu_parameters(fit, which = "prior", n.draws = 5)
# use_test("compose-input-for-stan.R")
# use_test("get-information-from-stanfit.R")
document()
devtools::test()


source("./R/class.R")
source("./R/get-information-from-stanfit.R")
source("./R/visualize-fit.R")

# Load data
fit = readRDS(file = "./tests/test models/IBBU_mv_fit_2 exposure groups_2 categories_2 cues_Drouin et al 2016.rds")

# class(fit) = "mv_ibbu_stanfit"
categories = c("s", "sh")
groups = c("SH-BIAS", "SS-BIAS")
# groups = c("SH-BIAS", "SS-BIAS", "Control")
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

plot_ibbu_parameters(fit, which = "prior", n.draws = 5)
plot_ibbu_parameters(fit, which = "both", n.draws = 5,
                     group.colors = c("darkgray", "red", "blue"))
