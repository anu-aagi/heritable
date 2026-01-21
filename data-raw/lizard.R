library(tidyverse)
library(janitor)

# install.packages("BiocManager")
# BiocManager::install("impute")
# remotes::install_github("italo-granato/snpReady")
library(snpReady)

growth_data <- read.csv("data-raw/lizard-data/Ldeli_quangen_growth.csv", stringsAsFactors = F)

lizard_phenotypes <- read_csv("data-raw/lizard-data/Ldeli_quangen_growth_DA.csv")

# Identifying key to match sample name and lizard/genotype
key <- lizard_phenotypes |>
  select(liz_id, F1_Genotype) |>
  distinct() |>
  rename(genotype = F1_Genotype)

lizard_phenotypes <- lizard_phenotypes |>
  select(liz_id, treatment, dam_id, sire_id, mass, lnMass, days_since_hatch, z_days_since_hatch, z_days_since_hatch_I2) |>
  clean_names()

lizard_markers <- read_csv("data-raw/lizard-data/LD_SNP_high.csv") |>
    rename(genotype = ...1)

growth_f1_ids <- unique(key$genotype)

# Subset to genotypes in study
markers_f1 <- lizard_markers |> semi_join(key)

# Left join liz_id
markers_f1 <- left_join(markers_f1, key) |>
  relocate(liz_id, .before = genotype) |>
  select(-genotype) |>
  rename(genotype = liz_id)

#Setting row names and checking colnames
colnames(markers_f1)

#Removing genotype column and just have the SNP data left
hq_snps <- markers_f1[,-1] |> as.matrix()
rownames(hq_snps) <- markers_f1$genotype

# Clean and impute
lizard_snpready <- snpReady::raw.data(data = hq_snps, call.rate = 0.60, maf = 0.0001,
                                frame="wide", base = FALSE, plot = TRUE, imput.type="wright")

dim(lizard_snpready$M.clean)
rownames(lizard_snpready$M.clean)
lizard_markers <- as_tibble(lizard_snpready$M.clean)
lizard_markers$gen <- rownames(lizard_snpready$M.clean)
lizard_markers <- lizard_markers |>
  relocate(gen, 1)

# G matrix
lizard_GRM <- snpReady::G.matrix(M = lizard_markers[-1], format = "wide",
                        method = "VanRaden",  plot = FALSE)$Ga

dimnames(lizard_GRM) <- list(rownames(lizard_snpready$M.clean), rownames(lizard_snpready$M.clean))

usethis::use_data(lizard_phenotypes, overwrite = TRUE)
usethis::use_data(lizard_markers, overwrite = TRUE)
usethis::use_data(lizard_GRM, overwrite = TRUE)


