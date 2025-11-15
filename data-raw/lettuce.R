## code to prepare `lettuce` dataset goes here
library(tidyverse)
lettuce_phenotypes <- read_delim("data-raw/lettuce-data/phenotypes_lettuce.txt",
                                 col_types = cols(
                                   loc = col_factor(),
                                   gen = col_factor(),
                                   rep = col_factor(),
                                   y = col_double()
                                 )) |>
  # get rid of attributes
  as.data.frame() |>
  as_tibble()
str(lettuce_phenotypes)

lettuce_markers <- read_delim("data-raw/lettuce-data/markers_lettuce.txt") |>
  # get rid of attributes
  as.data.frame() |>
  as_tibble()

str(lettuce_markers)

lettuce_GRM <- tcrossprod(as.matrix(lettuce_markers[, -1]))
dimnames(lettuce_GRM) <- list(lettuce_markers$gen, lettuce_markers$gen)

usethis::use_data(lettuce_phenotypes, overwrite = TRUE)
usethis::use_data(lettuce_markers, overwrite = TRUE)
usethis::use_data(lettuce_GRM, overwrite = TRUE)
