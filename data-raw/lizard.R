library(tidyverse)

# install.packages("BiocManager")
# BiocManager::install("impute")
# remotes::install_github("italo-granato/snpReady")

library(snpReady)

growth_data <- read.csv("data-raw/lizard-data/Ldeli_quangen_growth.csv", stringsAsFactors = F)
str(growth_data)

lizard_phenotypes <- read_csv("data-raw/lizard-data/Ldeli_quangen_growth_DA.csv")
dim(lizard_phenotypes)
names(lizard_phenotypes)

lizard_phenotypes <- lizard_phenotypes |>
  select(liz_id, treatment, dam_id, sire_id, mass, lnMass, days_since_hatch, z_days_since_hatch, z_days_since_hatch_I2)

lizard_markers <- read_csv("data-raw/lizard-data/LD_SNP_high.csv") |>
    rename(genotype = ...1)
dim(lizard_markers)

growth_f1_ids <- unique(growth_data$F1_Genotype)

key <- growth_data |>
  select(liz_id, F1_Genotype) |>
  distinct() |>
  rename(genotype = F1_Genotype)

markers_f1 <- lizard_markers[lizard_markers$genotype %in% growth_f1_ids ,]
dim(markers_f1)

# Left join liz_id
markers_f1 <- left_join(markers_f1, key) |>
  relocate(liz_id, .before = genotype) |>
  select(-genotype) |>
  rename(genotype = liz_id)

#Setting row names and checking colnames
colnames(markers_f1)
rownames(markers_f1) <- markers_f1$genotype

#Removing genotype column and just have the SNP data left
hq_snps <- markers_f1[,-1] |> as.matrix()
rownames(hq_snps) <- markers_f1$genotype

# Clean and impute
lizard_snpready <- snpReady::raw.data(data = hq_snps, call.rate = 0.60, maf = 0.0001,
                                frame="wide", base = FALSE, plot = TRUE, imput.type="wright")

dim(lizard_snpready$M.clean)
rownames(lizard_snpready$M.clean)
lizard_markers <- lizard_snpready$M.clean

# G matrix
lizard_GRM <- snpReady::G.matrix(M=lizard_markers, format = "wide",
                        method = "VanRaden",  plot = TRUE)

usethis::use_data(lizard_phenotypes, overwrite = TRUE)
usethis::use_data(lizard_markers, overwrite = TRUE)
usethis::use_data(lizard_GRM, overwrite = TRUE)


