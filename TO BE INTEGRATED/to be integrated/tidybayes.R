#' @import assertthat
NULL

compose_data = function()
  # in composing the data and fitting the model make sure that the model inherits
  # variable names and values for e.g., the categories and cues, so that they can
  # can be used in spread_draws and alike.

recover_types = function()
  # add some function that adds in the variable names to the model
