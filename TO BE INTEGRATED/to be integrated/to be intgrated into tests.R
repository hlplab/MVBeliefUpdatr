library(dplyr)
library(tidybayes)
install_github("git@github.com:tfjaeger/MVbeliefupdatr.git")
library(MVBeliefUpdatr)
# devtools::reload(MVBeliefUpdatr)
library(magrittr)

# Load data
fit = readRDS(file = "../models/TEMP-ibbu_mv_pca_by bias-all 3 exposure conditions.rds")
fit.input = readRDS(file = "../models/TEMP-data_mv_pca_by bias-all 3 exposure conditions.rds")

# class(fit) = "mv_ibbu_stanfit"
categories = c("s", "sh")
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
  add_ibbu_draws(wide = F, which = "prior", draws = 10) -> f
unique(f$group)

plot_ibbu_test_categorization(fit, fit.input, n.draws = 500, sort.by = "prior", group.colors = c("darkgray", "blue", "red", "black")) +
  theme_bw()
last_plot() + facet_wrap(~group)



fit = readRDS(file = "../models/TEMP-ibbu_mv_pca_by bias-SS & SH exposure conditions.rds")
fit.input = readRDS(file = "../models/TEMP-data_mv_pca_by bias-SS & SH exposure conditions.rds")
fit %<>%
  recover_types(d) 
plot_ibbu_test_categorization(fit, fit.input, n.draws = 500, sort.by = "prior", group.colors = c("darkgray", "blue", "red"),
                              group.ids = c("prior", "SS-BIAS", "SH-BIAS")) +
  theme_bw()
last_plot() + facet_wrap(~group)



fit %>%
  spread_draws(mu_0)