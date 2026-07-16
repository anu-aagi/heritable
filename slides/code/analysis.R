library(dplyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ggrepel)
library(lme4)
library(asreml)
library(heritable)
library(stringr)
library(ggsci)
library(patchwork)

data <- read.csv("data/2025-GRDC-Data-20251027-subset.csv",
                 check.names = FALSE, stringsAsFactors = TRUE
)
colnames(data) <- colnames(data) %>%
  str_remove_all("\\s*\\([^\\)]*\\)") %>% # remove (...) including content
  str_replace_all(" ", "_")

# Get variables
meta_variable <- c(
  "Genotype", "Row", "Range", "Row_Range", "Rep", "Day", "TOS",
  "Leaf_Number", "Development", "Date", "Zadok_Score"
)
auxiliary_pheno <- colnames(data)[35:61]
main_pheno <- colnames(data)[62:73]

# Get hybrid lines
hybrid <- readLines("data/hybrid.txt")

# Subset data
data_subset <- dplyr::select(data, all_of(c(meta_variable, auxiliary_pheno, main_pheno))) %>%
  mutate(
    Leaf_Number = factor(Leaf_Number),
    Date = factor(Date),
    Development = factor(Development, levels = c("Booting", "Anthesis", "Fill")),
    Hybrid = factor(ifelse(Genotype %in% hybrid, "Hybrid", "Non-hybrid")),
    across(main_pheno, log),
    Gs = log(Gs)
  )

# Strata
strata <- with(
  data,
  expand.grid(levels(TOS), levels(Development))
) %>%
  `colnames<-`(c("TOS", "Development")) %>%
  data.frame()

# Response
responses <- c("Gs", main_pheno)
responses <- responses[-grep("Stdev", responses)]

data_subset$Treatment <- data_subset$TOS
data_subset$Stage <- data_subset$Development
data_subset$Plot <- data_subset$Row_Range
data_subset$Replicate <- data_subset$Leaf_Number
data_subset$Type <- data_subset$Hybrid
data_subset$y <- data_subset$Gs

frm <- formula(
  paste0(r, "~ Hybrid * TOS * Development + (1 | Genotype)",
         " + (Development | Genotype)",
         " + (Development | Genotype:TOS)",
         " + (1 | Date) ",
         " + (1 | Row_Range / Leaf_Number)")
)

model <- lmer(
  y ~ Type * Treatment * Stage +
    (1 | Genotype) +
    (Stage | Genotype) +
    (Stage | Genotype:Treatment) +
    (1 | Date) +
    (1 | Plot / Replicate),
  data = data,
  REML = FALSE
)

my_H2 <- H2(model, "Genotype")

saveRDS(data_subset, "code/data1.rds")
saveRDS(model, "code/model1.rds")
saveRDS(my_H2, "code/H21.rds")

df <- readRDS("data/merredin-2025-traits_valid.rds")
df$y <- df$`Growth Score (Maturity)`
model <- asreml::asreml(
  fixed = y ~ 1,
  random = ~ Line + Range + Row + RowRep:ColRep,
  residual = ~ ar1(Range):ar1(Row),
  data = df
)
my_H2 <- H2(model, "Line")

saveRDS(my_H2, "code/H22.rds")
