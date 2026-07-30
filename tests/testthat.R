library(tidyverse)
library(testthat)
library(MVBeliefUpdatr)

files <- list.files("tests/testthat", pattern = "^\\d{2}-test.*\\.R$", full.names = TRUE)
for (file in files) {
  testthat::test_file(file)
}
