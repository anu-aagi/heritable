# library(here)
# library(tidyverse)
# library(asreml)
#
# # ------- Description of data --------
#
# # Field trial phenotype data columns:
#
# # dataset: trial location + year combination
# # column, row: spatial location of plot in each trial
# # genotype: 584 breeding lines & check varieties
# # blcklegsur: percentage survival against blackleg
# # blcklgsur1: plant count at emergence
# # blcklgsur2: plant count at maturity
# # Genomic marker data consists of:
#
# # 584 rows named by the genotypes
# # 8,179 columns named by the post-filtered markers
# # Allele encodings are (AA, AB, BB) = (1, 0, -1)
# # Missing marker data was imputed using marker means (Montesinos-López et al. 2016; Rutkoski et al. 2016) of lines from within each breeding cohort
#
# # df is a bivariate single-site data-frame with the following columns:
#
# # pheno - phenotypic values of each trait (2n x 1 vector of observations, where n = c x r)
# # Trait - factor of 2 levels (BLeme, BLmat)
# # column - factor of c levels, where c = number of columns
# # row - factor of r levels, where r = number of rows
# # gkeep - factor of g levels, where g represents number of genotypes with genomic marker data
# # Kmat - g x g genomic relationship matrix ordered according to levels in gkeep
# # gdrop - factor of d levels, where d represents the number of genotypes without genomic marker data
#
# # df is ordered as rows within cols within traits
#
# # For genotype and phenotype data, please refer to Figshare link as stated in paper
#
# ## Load data
# bl_df <- read_csv(here("ignore/29095535/blackleg_survival_data.csv"))
# markers <- read_csv(here("ignore/29095535/genomic_marker_matrix_M.csv")) |> as.data.frame() |> as_tibble()
#
# ## Wrangle data
# ### Rename according to paper first
# # bl_df <- bl_df |>
# #  rename(
# #     blsur = blcklegsur,
# #     bleme = blcklgsur1,
# #     blmat = blcklgsur2
# #  )
#
# ### Calculate n
#  unique(bl_df$column)
#  unique(bl_df$row)
#
# # bl_df |>
# #     group_by(dataset) |>
# #     summarise(
# #         n_unique_c = n_distinct(column),
# #         n_unique_r = n_distinct(row),
# #         max_c = max(column),
# #         max_r = max(row),
# #         n = max(column) * max(row)
# #         )
#
# bl_df_n <- bl_df |>
#     group_by(dataset) |>
#     mutate(n = n_distinct(column) * n_distinct(row)) |>
#     ungroup()
#
# ### Long format
# bl_df_long <- bl_df_n |>
#     pivot_longer(
#         cols = starts_with("blcklgsur"),
#         names_to = "Trait",
#         values_to = "pheno"
#     )
#
#
# # bl_df_long |>
# #   filter(dataset == "RS23") |>
# #   arrange(Trait)
#
# # bl_df_long |>
# #     group_by(dataset, Trait) |>
# #     count()
#
# # G matrix and scaled G matrix (From paper Gs = MM'/r, r = trace G/m where m = number of genotypes))
# # M: individuals x markers
# # Tim uses pedicure package here
# M <- markers[,-1] |> as.matrix()
#
# marker_var <- apply(M, 2, var, na.rm = TRUE)
# scaling <- mean(marker_var)
#
# G <- tcrossprod(M) / scaling
#
# Gscaled <- tcrossprod(M) / (diag(G) / nrow(G)) # Is this the same as pedicure??
#
# # gkeep and gdrop
# # Are gkeep the genotypes where their marker data was NOT imputed? (-1,0,1)
# # and are gdrop the genotypes where their marker data was imputed? (e.g. fractional values?)
# # Is this across all markers? If one SNP was imputed then you count it for gdrop?
# gdrop_ids <- setdiff(unique(bl_df$genotype), markers$genotype) # These are gdrop
# gkeep_ids <- intersect(unique(bl_df$genotype), markers$genotype)  # These are gkeeps
#
# bl_df_valid <- bl_df_long |>
#   mutate(
#     gkeep = case_when(genotype %in% gkeep_ids ~ genotype,
#                       TRUE ~ NA),
#     gdrop = case_when(genotype %in% gdrop_ids ~ genotype,
#                       TRUE ~ NA)
#   ) |>
#   mutate(
#     Trait = fct(Trait),
#     genotype = fct(genotype),
#     gdrop = fct(gdrop),
#     gkeep = fct(gkeep),
#     column = as_factor(column),
#     row = as_factor(row)
#   )
#
# levels(bl_df_valid$genotype) |>  length()
# levels(bl_df_valid$gkeep) |>  length()
# levels(bl_df_valid$gdrop) |>  length()
#
# bl_df
# str(bl_df_valid)
#
# # ------- Univariate baseline --------
#
# # Genetic model : diag
# # Residual model : a separate spatial model for each Trait
# #########################################################
# # gkeep/gdrop a trick to use to maximise data usage, rather than throwing away missing data
# #TODO: There is a paper, animals
# bl_ss.df <- bl_df_valid |> filter(dataset == "RS23") |> arrange(Trait)
#
# bl_uni.b.asr <- asreml(pheno ~ Trait + Trait:gdrop,
#                     random = ~ diag(Trait):gkeep +
#                                at(Trait):column +
#                                at(Trait):row,
#                     residual = ~ dsum( ~ ar1(column):ar1(row) | Trait),
#                     na.action=na.method(x='include', y='include'),
#                     data = bl_ss.df,
#                     maxit = 20,
#                     workspace = "1000mb")
#
#
# # H2_Oakey(bl_uni.b.asr, "gkeep")
# n_g <- bl_uni.b.asr$noeff[["Trait:gkeep"]]
# vc_g <- bl_uni.b.asr$vparameters[["Trait:gkeep!Trait_blcklgsur1"]] * bl_uni.b.asr$sigma2
# vcov_g <- predict(bl_uni.b.asr, # I think I'm getting it wrong here. What am I predicting??
#                   classify = "gkeep",
#                   only = "gkeep",
#                   vcov = TRUE,
#                   trace = FALSE
# )$vcov
#
# Gg_inv <- diag(1 / vc_g, nrow = n_g, ncol = n_g)
#
# M <- diag(n_g) - (Gg_inv %*% vcov_g)
# eM <- eigen(M)
#
# ##########################################################
#
# # ------- Bivariate baseline --------
#
# # Genetic model : corgh
# # Residual model : corgh three-way separable
# #########################################################
#
# bl_biv.b.asr <- asreml(pheno ~ Trait + Trait:gdrop,
#                     random = ~ corgh(Trait):gkeep +
#                       at(Trait):column +
#                       at(Trait):row,
#                     residual =  ~ corgh(Trait):ar1(column):ar1(row),
#                     na.action=na.method(x='include', y='include'),
#                     data = bl_ss.df,
#                     maxit = 20,
#                     workspace = "1400mb")
#
# ##### Looking at Tim's Rdata
#
# # load(here("ignore/FONTI/timtam.RData"))
# #
# #
# # glimpse(df)
# # df |> count(Trait)
# # df |> count(genotype)
# # df |> count(gkeep)
# # df |> count(gdrop)
# #
# # M_df <- M |> as_tibble()
# # M_df$genotype <- rownames(M)
# #
# # M_df |> filter(str_detect(genotype, c("G585|G586"))) # They don't exist in M
# #
# # df |> filter(str_detect(genotype, c("G585|G586")))
# #
# # setdiff(levels(df$genotype), M_df$genotype) # These are gdrop
# # intersect(levels(df$genotype), M_df$genotype) |> length() # These are gkeeps
# # levels(df$gkeep)
# # levels(df$gdrop)
#
#
