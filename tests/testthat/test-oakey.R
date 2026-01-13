test_that("OAKEY heritability estimation works", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  skip_on_ci()

  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))
  lettuce_asreml <- readRDS(test_path("fixtures/lettuce_asreml.rds"))

  # Method implemented by ET as per Oakey et al. 2006
  H2_Oakey(asreml_model_random, target = "gen") 

  # Method implemented by YD as per Schmid et al. 2009
  H2_Oakey_YD <- function(model, target) {
  n_g <- model$noeff[[target]]
  vc_g <- heritable:::get_vc_g_asreml(model, target)
  vcov_g <- predict(model,
    classify = target,
    only = target,
    vcov = TRUE,
    trace = FALSE
  )$vcov

  Gg_inv <- diag(1 / vc_g, nrow = n_g, ncol = n_g)
  svds <- svd(Gg_inv)
  Gg_inv_sqrt <- sweep(svds$u, 2, sqrt(svds$d), "*") %*% t(svds$v)
  M <- diag(n_g) - (Gg_inv_sqrt %*% vcov_g %*% Gg_inv_sqrt)
  eM <- eigen(M)
  thres <- 1e-5

  mean(eM$values[eM$values > thres])
  }

  expect_equal(
    H2_Oakey(asreml_model_random, target = "gen"),
    H2_Oakey_YD(asreml_model_random, target = "gen"),
    tolerance = 1e-5
  )

  # lettuce model
  expect_equal(
    H2_Oakey(lettuce_asreml, target = "gen"), 
    H2_Oakey_YD(lettuce_asreml, target = "gen"),
    tolerance = 1e-5
  )
})
library(here)
library(tidyverse)
library(asreml)

# -------- Lourenco et al 2020 --------
 
## The R code below allows the user to generate phenotypic data from an alpha-design
## as well as generate SNP marker data and kinship matrix from it;
## This is because we cannot share the true SNP data and kinship matrix
set.seed(9999999)
library(dplyr) 
# install.packages("agricolae")
library(agricolae) # for using the design.alpha() function; needs the right R version installed
# enables the use of operator %>%
## generating an alpha-design with 702 genotypes and 39 blocks each with size 18
## NOTE that in our case, the simulated datasets had 698 genotypes
## in which case positions 699--702 in the design were set to missing values
# defining the parameters for the design
trt <- 1:702 # number of treatments, i.e., genotypes
k <- 18 # block size
r <- 2 # nubmer of replicates
# s <- trt[702]/k # number of blocks; don’t need to specify this
# creating the design, which is only ~ 90% efficient
# the user can create more efficient alpha-desgins by using the program CycDesignN
outdesign <- design.alpha(trt, k, r, serie = 2, seed = 0, kinds = "Super-Duper",
randomization=TRUE)
# getting the dataset obtained from the design
ourdata<-outdesign$book[,c(2,5,3,4)]
# changing colnames
colnames(ourdata)<-c("plot","rep","block","geno")
# changing blocks 40:78 to 1:39 (this is actually not needed)
ourdata[703:1404,3]<-as.integer(ourdata[703:1404,3])-39
## simulating maize yield
# fixing the variances for the synthetic data
var_s <- 0.005892 # marker effect variance
var_b <- 6.3148 # block effect variance
var_e <- 53.8715 # residual variance
# simulationg the SNP effects from N(0,var_s)
set.seed(79753)
u<-rnorm(11646,0,sqrt(var_s))
# generating SNP marker 702x111646 matrix (in this case a matrix of 1s and -1s)
# as well as the kisnhip matrix
# generating the SNP markers as below is just an had-doc way of doing it;
# it assumes the SNPs are independent and does not account for population structure
Z<-matrix(rbinom(n=702*11646,size=1,prob=0.5),702,11646)
for(i in 1:702){for(j in 1:11646){Z[i,j]<-2*Z[i,j]-1}}
K<-Z%*%t(Z) # or K<-Z%*%t(Z)/11646
# computing the true 702 breeding values
g<-Z%*%u
1
# computing plot and block effects for the 2 replicates
# block effect is the same within each of the 39*2 blocks
# but different across all the 39*2 blocks
block_eff <- vector()
for(i in 1:78){block_eff<-c(block_eff,rep(rnorm(1,0,sqrt(var_b)),18))}
res_eff <- rnorm(702*2,0,sqrt(var_e))
plot_eff <- block_eff+res_eff
ourdata <- cbind(ourdata,block_eff,res_eff,plot_eff)
# computing yield values for the 2 replicates; first 702 pos are rep1; last 702 are rep 2
# the TBVs are the same for the 702 genotypes on both replicates
dataset <- (ourdata[order(ourdata$rep,ourdata$geno),])
dataset <- cbind(dataset,rep(g,2))
colnames(dataset)[8]<-c("TBVs")
yield <- dataset[,7]+dataset[,8]
dataset <- cbind(dataset,yield)
# # UNCOMMENT THE CODE BELOW IF YOU WANT TO CONTAMINATE THE DATA
# # RANDOM contamination of a % of observations on replicate 1
# # define the % of contamination, e.g.,
# perc.cont<-0.05
# # compute the number of outliers referring to that %
# n.outliers<-round(702*perc.cont)
# # randomly chose the positions/genotypes to contaminate
# # in this case, because we ordered the dataset, the positions to contaminate
# # match the genotypes of the first replicate with the same index
# outlier.positions<- sample(1:702,n.outliers)
# # replace the observations in the assigned outlier.positions
# # by their value plus, e.g., 5 times the standard deviation
# dataset[outlier.positions,9]<-dataset[outlier.positions,9]+5*sqrt(var_e)
# # whole BLOCK contamination on replicate 1
# # define the number of blocks to contaminate, e.g.,
# n.blocks<-3
# # randomly chose the blocks on the first replicate that should be contaminated
# # the number of contaminated observations = n.blocks*18
# outlier.blocks<- sort(sample(1:39,n.blocks))
# # replace the observations in the assigned outlier.blocks
# # by their value plus, e.g., 5 times the standard deviation
# for(i in 1:n.blocks){
# dataset[dataset$block==outlier.blocks[i]&dataset$rep==1,9]<-
# dataset[dataset$block==outlier.blocks[i]&dataset$rep==1,9]+5*sqrt(var_e)
# }
# sorting the dataset by geno so that the robust 2-stage R code can be applied directly
toydata<-dataset[order(dataset$geno),]
toyG<-K
# saving the dataset
save(toydata,toyG,file=test_path("fixtures/toydata.Rdata"))

## loading the data
load(test_path("fixtures/toydata.Rdata"))
## loading the libraries
library(lme4) # for classical LMM
# install.packages("robustlmm")
library(robustlmm) # for robust LMM
library(psych) # to being able of using the tr() function
library(asreml) # for the classical LMM with kinship matrix
# if tr() does not work load library(matrixcalc) and use matrix.trace() function
## CLASSICAL first-stage fit below is commented;
## uncomment when you wish to run it, in which case you should comment the robust fit
# ## fitting the first-stage classical model
# fit<- lmer(yield ~ -1 + geno + (1|rep)+(1|rep:block), toydata)
# # getting the lsmeans and var-covar structure
# R <- summary(fit)$vcov
# mu <- summary(fit)$coefficients[,1]
# # computing the Smith’s and Standard weights
# w <-(1/diag(R)) #Standard weights
# wsmith<-diag(solve(R)) #Smith’s weights
# rm(R)
# # keeping also the estimated random effects variances
# STDs<-matrix(0,3,1)
# STDs[1,1]<- attr(VarCorr(fit)$’rep:block’, "stddev")
# STDs[2,1]<- attr(VarCorr(fit)$’rep’, "stddev")
# STDs[3,1]<- attr(VarCorr(fit),"sc")
# colnames(STDs)<-"std"
# rownames(STDs)<-c("REP:BLOCK","REP","Residual")
# stage1.vars<-STDs^2
# rm(STDs)
## fitting the first-stage robust model
## Takes a while to run
fit<- robustlmm::rlmer(yield ~ -1 + geno + (1|rep)+(1|rep:block), toydata,
rho.sigma.e = psi2propII(smoothPsi, k = 2.28)) 

save(fit,file=test_path("fixtures/robust_fit_oakey.Rdata"))

# getting the lsmeans and var-covar structure
R <- summary(fit)$vcov
mu <- summary(fit)$coefficients[,1]
# getting the robust weights
# do not confuse these with the Smith’s and Standard weights
# Note that if your data has missing values of yield, no robust weights are estimated
# and therefore the process of getting the robust weights for the 2nd-stage
# will not be as straighforward as it is in this case
rob.weights<-getME(fit,name="w_e")
# the robust weights need not the same for genos in rep1 and rep2
# but we want only 1 robust weight per-genotype
1
# thus we will choose the min between the 2 robust weights from the 2 replicates
# the next computations need to be adpated for each dataset because the order of
# the weights matches the one of the dataset as also do the order of the residuals
n<-length(mu)
aux<-vector()
for(k in seq(1,(n*2-1),by=2)){aux<-c(aux,min(rob.weights[k],rob.weights[k+1]))}
rob.weights<-aux
rm(aux,k,n)
# computing the Smith’s and Standard weights, which incorporate the robust weights
w <-(1/diag(R))*rob.weights #Standard weights
wsmith <-diag(solve(R))*rob.weights #Smith’s weights
rm(rob.weights,R)
# keeping also the estimated random effects variances
STDs<-matrix(0,3,1)
STDs[1,1]<- attr(VarCorr(fit)$’rep:block’, "stddev")
STDs[2,1]<- attr(VarCorr(fit)$’rep’, "stddev")
STDs[3,1]<- attr(VarCorr(fit),"sc")
colnames(STDs)<-"std"
rownames(STDs)<-c("REP:BLOCK","REP","Residual")
stage1.vars<-STDs^2
rm(STDs)
## fitting the second-stage model -- classical approach used
# # try out G=I to see how H2.M5 and H2.Oakey match
# toyG<-diag(dim(toyG)[1])
# colnames(toyG)<-names(mu)
# rownames(toyG)<-names(mu)
# preparing the data
plantid<-names(mu)
colnames(toyG)<-names(mu)
rownames(toyG)<-names(mu)
inv.toyG<-solve(toyG)
ourdata<-data.frame(plantid = as.factor(plantid),
mu = mu,
wsmith = wsmith,
w = w)
# fitting the model
# one can change the Smith’s weights (wsmith) for the Standard weights (w) below
fit.cls <-asreml( data =ourdata,
fixed =mu ~ 1,
random =~ vm(plantid, inv.toyG),
rcov =~ units, 
na.action = na.method("include"),
weights = wsmith,
family = asr_gaussian(dispersion=1.0),
workspace=16e7,
maxiter=1000
)
# computing the eBLUPs, estimated genetic variance and C22 matrix
gBLUP <- fit.cls$coefficients$random
s.var <- summary(fit.cls)$varcomp[’giv(plantid).giv’,’component’]
C22 <- predict(fit.cls, classify="giv(plantid)", only="giv(plantid)", vcov=T)$pred$vcov
2
# removing stuff from memory
rm(fit,fit.cls,ourdata,plantid)
rm(inv.toyG)
## third-stage -- heritability and predictive accuracy estimation
# preparing the matrices and auxiliary variables as in the paper notation
G <-toyG
n <-dim(G)[1]
G.tilde<-G*s.var
R.tilde<-solve(diag(n)*(wsmith))
rm(toyG)
# computing heritability and predictive accuracy via METHOD 5
V <-G.tilde+R.tilde
P <-(1/(n-1))*(diag(n)-matrix(1,n,n)/n)
one <-as.matrix(rep(1,n))
Q <-diag(n)-one %*% solve(t(one)%*%solve(V)%*%one) %*% t(one) %*% solve(V)
C <-G.tilde%*%solve(V)%*%Q
PA.est.m5<-tr(P%*%C%*%G.tilde)/sqrt(tr(P%*%G.tilde)*tr(t(C)%*%P%*%C%*%V))
H2.est.m5<-PA.est.m5^2
rm(V,P,Q,C,one)
# computing reliability and predictive accuracy via METHOD 7
v1<-G.tilde
v2<-G.tilde-C22
rho2<-vector()
for(j in 1:n){rho2[j]<-(v2[j,j])^2/(v1[j,j]*v2[j,j])}
rm(j,v1,v2)
RL.est.m7 <-mean(rho2)
PA.est.m7 <-mean(sapply(rho2,sqrt))
rm(rho2)
# computing heritability via OAKEY’s METHOD
D <-diag(n)-solve(G.tilde)%*%C22
eival <-eigen(D)$values
s <-length(eival[eival<0.0001])
H2.OAKEY<-tr(D)/(n-s)
rm(D,eival,s,G.tilde,R.tilde)
rm(n,w,wsmith)
rm(G,C22)
# printing out the results
cbind(t(stage1.vars), s.var)
cbind(H2.est.m5,H2.OAKEY,PA.est.m5,PA.est.m7)

# ------- Tim's data --------

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
bl_df <- read_csv(here("ignore/29095535/blackleg_survival_data.csv"))
markers <- read_csv(here("ignore/29095535/genomic_marker_matrix_M.csv")) |> as.data.frame() |> as_tibble()

## Wrangle data
### Rename according to paper first
# bl_df <- bl_df |>
#  rename(
#     blsur = blcklegsur,
#     bleme = blcklgsur1,
#     blmat = blcklgsur2
#  )

### Calculate n
 unique(bl_df$column)
 unique(bl_df$row)

# bl_df |>
#     group_by(dataset) |>
#     summarise(
#         n_unique_c = n_distinct(column),
#         n_unique_r = n_distinct(row),
#         max_c = max(column),
#         max_r = max(row),
#         n = max(column) * max(row)
#         )

bl_df_n <- bl_df |>
    group_by(dataset) |>
    mutate(n = n_distinct(column) * n_distinct(row)) |>
    ungroup()

### Long format
bl_df_long <- bl_df_n |>
    pivot_longer(
        cols = starts_with("blcklgsur"),
        names_to = "Trait",
        values_to = "pheno"
    )


# bl_df_long |>
#   filter(dataset == "RS23") |>
#   arrange(Trait)

# bl_df_long |>
#     group_by(dataset, Trait) |>
#     count()

# G matrix and scaled G matrix (From paper Gs = MM'/r, r = trace G/m where m = number of genotypes))
# M: individuals x markers
# Tim uses pedicure package here
# M <- markers[,-1] |> as.matrix()

# marker_var <- apply(M, 2, var, na.rm = TRUE)
# scaling <- mean(marker_var)

# G <- tcrossprod(M) / scaling

# Gscaled <- tcrossprod(M) / (diag(G) / nrow(G)) # Is this the same as pedicure??

# # gkeep and gdrop
# # Are gkeep the genotypes where their marker data was NOT imputed? (-1,0,1)
# # and are gdrop the genotypes where their marker data was imputed? (e.g. fractional values?)
# # Is this across all markers? If one SNP was imputed then you count it for gdrop?
# gdrop_ids <- setdiff(unique(bl_df$genotype), markers$genotype) # These are gdrop
# gkeep_ids <- intersect(unique(bl_df$genotype), markers$genotype)  # These are gkeeps

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

# levels(bl_df_valid$genotype) |>  length()
# levels(bl_df_valid$gkeep) |>  length()
# levels(bl_df_valid$gdrop) |>  length()

# bl_df
# str(bl_df_valid)

# # ------- Univariate baseline --------

# # Genetic model : diag
# # Residual model : a separate spatial model for each Trait
# #########################################################
# # gkeep/gdrop a trick to use to maximise data usage, rather than throwing away missing data
# #TODO: There is a paper, animals
# bl_ss.df <- bl_df_valid |> filter(dataset == "RS23") |> arrange(Trait)

# bl_uni.b.asr <- asreml(pheno ~ Trait + Trait:gdrop,
#                     random = ~ diag(Trait):gkeep +
#                                at(Trait):column +
#                                at(Trait):row,
#                     residual = ~ dsum( ~ ar1(column):ar1(row) | Trait),
#                     na.action=na.method(x='include', y='include'),
#                     data = bl_ss.df,
#                     maxit = 20,
#                     workspace = "1000mb")


# # H2_Oakey(bl_uni.b.asr, "gkeep")
# n_g <- bl_uni.b.asr$noeff[["Trait:gkeep"]]
# vc_g <- bl_uni.b.asr$vparameters[["Trait:gkeep!Trait_blcklgsur1"]] * bl_uni.b.asr$sigma2
# vcov_g <- predict(bl_uni.b.asr, # I think I'm getting it wrong here. What am I predicting??
#                   classify = "gkeep",
#                   only = "gkeep",
#                   vcov = TRUE,
#                   trace = FALSE
# )$vcov

# Gg_inv <- diag(1 / vc_g, nrow = n_g, ncol = n_g)

# M <- diag(n_g) - (Gg_inv %*% vcov_g)
# eM <- eigen(M)

# ##########################################################

# # ------- Bivariate baseline --------

# # Genetic model : corgh
# # Residual model : corgh three-way separable
# #########################################################

# bl_biv.b.asr <- asreml(pheno ~ Trait + Trait:gdrop,
#                     random = ~ corgh(Trait):gkeep +
#                       at(Trait):column +
#                       at(Trait):row,
#                     residual =  ~ corgh(Trait):ar1(column):ar1(row),
#                     na.action=na.method(x='include', y='include'),
#                     data = bl_ss.df,
#                     maxit = 20,
#                     workspace = "1400mb")

##### Looking at Tim's Rdata

# load(here("ignore/FONTI/timtam.RData"))
#
#
# glimpse(df)
# df |> count(Trait)
# df |> count(genotype)
# df |> count(gkeep)
# df |> count(gdrop)
#
# M_df <- M |> as_tibble()
# M_df$genotype <- rownames(M)
#
# M_df |> filter(str_detect(genotype, c("G585|G586"))) # They don't exist in M
#
# df |> filter(str_detect(genotype, c("G585|G586")))
#
# setdiff(levels(df$genotype), M_df$genotype) # These are gdrop
# intersect(levels(df$genotype), M_df$genotype) |> length() # These are gkeeps
# levels(df$gkeep)
# levels(df$gdrop)


