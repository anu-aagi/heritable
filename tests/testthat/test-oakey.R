library(here)
library(tidyverse)
library(asreml)

# ------- Description of data --------

# Field trial phenotype data columns:

# dataset: trial location + year combination
# column, row: spatial location of plot in each trial
# genotype: 584 breeding lines & check varieties
# blcklegsur: percentage survival against blackleg
# blcklgsur1: plant count at emergence
# blcklgsur2: plant count at maturity
# Genomic marker data consists of:

# 584 rows named by the genotypes
# 8,179 columns named by the post-filtered markers
# Allele encodings are (AA, AB, BB) = (1, 0, -1)
# Missing marker data was imputed using marker means (Montesinos-López et al. 2016; Rutkoski et al. 2016) of lines from within each breeding cohort

# df is a bivariate single-site data-frame with the following columns:

# pheno - phenotypic values of each trait (2n x 1 vector of observations, where n = c x r)
# Trait - factor of 2 levels (BLeme, BLmat)
# column - factor of c levels, where c = number of columns
# row - factor of r levels, where r = number of rows
# gkeep - factor of g levels, where g represents number of genotypes with genomic marker data
# Kmat - g x g genomic relationship matrix ordered according to levels in gkeep
# gdrop - factor of d levels, where d represents the number of genotypes without genomic marker data

# df is ordered as rows within cols within traits

# For genotype and phenotype data, please refer to Figshare link as stated in paper

## Load data
df <- read_csv(here("ignore/29095535/blackleg_survival_data.csv")) 
markers <- read_csv(here("ignore/29095535/genomic_marker_matrix_M.csv")) |> as.data.frame() |> as_tibble()

## Wrangle data
### Rename according to paper first
df <- df |> 
 rename(
    blsur = blcklegsur,
    bleme = blcklgsur1,
    blmat = blcklgsur2
 )

### Calculate n 
 unique(df$column)
 unique(df$row)

# df |> 
#     group_by(dataset) |> 
#     summarise(
#         n_unique_c = n_distinct(column),
#         n_unique_r = n_distinct(row),
#         max_c = max(column),
#         max_r = max(row),
#         n = max(column) * max(row)
#         )

df_n <- df |> 
    group_by(dataset) |> 
    mutate(n = n_distinct(column) * n_distinct(row)) |> 
    ungroup() 

### Long format
df_long <- df_n |> 
    pivot_longer(
        cols = c(blsur, bleme, blmat),
        names_to = "Trait",
        values_to = "pheno"
    )

# df_long |> 
#     group_by(dataset, Trait) |> 
#     count()

# G matrix and scaled G matrix (From paper Gs = MM'/r, r = trace G/m where m = number of genotypes))
# M: individuals x markers
M <- markers[,-1] |> as.matrix()

marker_var <- apply(M, 2, var, na.rm = TRUE)
scaling <- mean(marker_var)

G <- tcrossprod(M) / scaling

Gscaled <- tcrossprod(M) / (diag(G) / nrow(G)) 

# gkeep and gdrop
# Are gkeep the genotypes where their marker data was NOT imputed? (-1,0,1)
# and are gdrop the genotypes where their marker data was imputed? (e.g. fractional values?)
# Is this across all markers? If one SNP was imputed then you count it for gdrop? 


# ------- Univariate baseline --------

# Genetic model : diag
# Residual model : a separate spatial model for each Trait 
#########################################################

uni.b.asr <- asreml(pheno ~ Trait + Trait:gdrop,
                    random = ~ diag(Trait):gkeep + 
                               at(Trait):column +
                               at(Trait):row,
                    residual = ~ dsum( ~ ar1(prange):ar1(prow) | Trait),
                    na.action=na.method(x='include', y='include'), 
                    data = df, 
                    maxit = 20,
                    workspace = "1000mb")
