## code to prepare `lettuce` dataset goes here
library(tidyverse)
lettuce_phenotypes <- read_delim("data-raw/lettuce-data/phenotypes_lettuce.txt")
lettuce_markers <- read_delim("data-raw/lettuce-data/markers_lettuce.txt")

usethis::use_data(lettuce_phenotypes, overwrite = TRUE)
usethis::use_data(lettuce_markers, overwrite = TRUE)
