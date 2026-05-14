# # This will be removed in the future
# library(dplyr)
# library(ggplot2)
# library(reshape2)
# library(pheatmap)
# library(ggrepel)
# library(lme4)
# library(asreml)
# library(heritable)
#
# data <- read.csv("../support-2026-04-usyd-stomata/data/data-raw/2025-GRDC-Data-20251027-subset.csv", check.names = FALSE)
#
# # Get variables
# meta_variable <- c("Genotype", "Row", "Range", "Row_Range", "Rep", "Day", "TOS",
#                    "Leaf Number", "Development", "Date", "Zadok Score")
# auxiliary_pheno <- colnames(data)[35:61]
# main_pheno <- colnames(data)[62:73]
#
# # Get hybrid lines
# hybrid <- readLines("../support-2026-04-usyd-stomata/data/data-raw/hybrid.txt")
#
# # Subset data
# data_subset <- dplyr::select(data, all_of(c(meta_variable, auxiliary_pheno, main_pheno)))
# data_subset$Date <- factor(data_subset$Date)
# data_subset$Development <- factor(data_subset$Development, levels = c("Booting", "Anthesis", "Fill"))
# data_subset$Genotype <- factor(data_subset$Genotype)
# data_subset$TOS <- factor(data_subset$TOS)
# data_subset$Row_Range <- factor(data_subset$Row_Range)
# data_subset$`Leaf Number` <- factor(data_subset$`Leaf Number`)
# data_subset$Leaf_Number <- data_subset$`Leaf Number`
# data_subset$Hybrid <- ifelse(data_subset$Genotype %in% hybrid, "Hybrid", "Non-hybrid") %>% factor
#
#
# model <- lmer(
#   `Gs` ~ TOS * Development +
#     (Development || Genotype / TOS) +
#     (1 | Date) +
#     (1 | Row_Range / `Leaf Number`),
#   data = data_subset
# )
# H2(model, "Genotype", marginal = TRUE)
# H2(model, "Genotype", stratification = data.frame(TOS = "TOS2", Development = "Fill"))
#
# model <- asreml(fixed = Gs ~ TOS * Development,
#                 random = ~  Genotype/TOS +
#                   Genotype:ar1(Development) +
#                   Date +
#                   Row_Range / Leaf_Number,
#                 data = data_subset)
# summary(model)
# H2(model, "Genotype", marginal = TRUE)
# H2(model, "Genotype", stratification = data.frame(TOS = "TOS2", Development = "Fill"))
#
#
# model <- asreml(fixed = y ~ rep,
#                 random = ~  loc*gen,
#                 data = lettuce_phenotypes)
# H2(model, "gen", solver = "LMM")
#
# table(data_subset$Genotype, data_subset$Development, data_subset$TOS)
#
# model <- lme4::lmer(
#   `Gs`^(1/3) ~
#     Hybrid*TOS * Development +
#     (Development | Genotype) +
#     (Development | Genotype : TOS) +
#     (1 | Date) +
#     (1 | Row_Range / `Leaf Number`),
#   data = data_subset
# )
# H2(model, "Genotype")
# H2(model, "Genotype", stratification = data.frame(TOS = "TOS2", Development = "Fill"))
