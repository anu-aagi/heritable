# test_that("Refactored asreml works",{
#   skip()
# lettuce_phenotypes <- subset(heritable::lettuce_phenotypes, !is.na(y))
# model <-  lme4::lmer(y ~  rep + (1 | gen * loc), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1 | gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1 | gen) + (1|gen:loc), data = lettuce_phenotypes)
#
# # when i fit a conterpart model of
# model <-  lme4::lmer(y ~  rep + (1 | gen * loc), data = lettuce_phenotypes)
# should it be
# model <-  lme4::lmer(y ~  rep + gen + gen:loc + (1 | loc), data = lettuce_phenotypes)
# # or
# model <-  lme4::lmer(y ~  rep + gen + gen:loc + loc, data = lettuce_phenotypes)
#
# # if i want base line heritability, should the conterpart model be
# model <-  lme4::lmer(y ~  rep + gen + gen:loc + (1 | loc), data = lettuce_phenotypes)
# # or
# model <-  lme4::lmer(y ~  rep + gen + (1 | loc:gen) + (1 | loc), data = lettuce_phenotypes)
#
#
# H2(model, "gen", marginal = FALSE)
# H2(model, "gen", marginal = TRUE)
# H2(model, "gen", stratification = data.frame(loc = "L1"))
#
# lettuce_subset <- subset(heritable::lettuce_phenotypes, loc == "L1")
# model <-  lme4::lmer(y ~  rep + (1 | gen), data = lettuce_subset)
# H2(model, "gen", marginal = FALSE)


#
# H2(model, "gen", marginal = TRUE, type = "BLUE")
# H2(model, "gen", marginal = TRUE, stratification = data.frame(loc = "L2"))
# h2(model, "gen", marginal = TRUE, type = "BLUE")
#
# vc <- var_comp(model, "gen", marginal = TRUE, calc_V = TRUE, calc_C22 = FALSE)
# V <- vc$V
# Z <- vc$Z
# G <- vc$G
# idx <- vc$idx
# X <- vc$X
# y <- vc$y
# m <- vc$m
#
# X_tilde <- cbind(X, Z[, idx, drop=FALSE])
# P <- ginv_sym_sparse(crossprod(X_tilde)) %*% t(X_tilde)
# P <- t(X_tilde) %*% ginv_sym_sparse(tcrossprod(X_tilde)) %*% X_tilde
# c <- c(0,0,0, m[,1] - m[,3])
# c - P %*% c
#
#
# lettuce_subset <- subset(heritable::lettuce_phenotypes, loc == "L2")
# model <-  lme4::lmer(y ~  rep + (1 | gen), data = lettuce_subset)
# H2(model, "gen", marginal = TRUE, type = "BLUE")
# h2(model, "gen", marginal = TRUE, type = "BLUE")

#
# model <-  lme4::lmer(y ~  rep + (0+ loc | gen ) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (loc | gen ) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1 | gen * loc) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1 | gen * loc) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1+ 1+ loc || gen ) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1| gen), data = lettuce_subset)

# require(devtools)
# load_all()
# model <-  lme4::lmer(y ~  rep + (1 | gen * loc), data = lettuce_phenotypes)
# system.time(
#   confint(H2(model, "gen"), B = 20)
# )
# system.time(
#   confint(H2(model, "gen"), B = 20, parallel = "multicore", ncpus = 10)
# )
# system.time(
#   confint(H2(model, "gen"), B = 20, parallel = "snow", ncpus = 2)
# )
# H2(model, "gen", method = "Delta")
#

#
# h2_Standard.lmerMod(model, "gen")
# h2_Oakey.lmerMod(model, "gen")
# h2_Delta_pairwise.lmerMod(model, "gen", type = "BLUE")
# h2_Delta_pairwise.lmerMod(model, "gen", type = "BLUP")
#

# model <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc + gen + gen:loc ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
# H2(model, "gen", marginal = TRUE)
# H2(model, "gen", marginal = TRUE, stratification = data.frame(loc = "L3"))
#



# system.time(
#   confint(H2(model, "gen"), B = 20)
# )
# system.time(
#   confint(H2(model, "gen"), B = 20, parallel = "multicore", ncpus = 10)
# )
# system.time(
#   confint(H2(model, "gen"), B = 100, parallel = "snow", ncpus = 3)
# )
# var_comp(model, "gen", marginal = TRUE, stratification = data.frame(loc = "L1"))
#
#
#
#
#
# M <- as.matrix(lettuce_markers[, -1] + 1)
# N <- nrow(M)
# pm <- colSums(M) / (2 * N) # allele freq per marker (diploid X)
# pm <- pmin(pmax(pm, 1e-6), 1 - 1e-6) # guard against 0 or 1
# W <- sweep(M, 2, 2 * pm, "-")
# W <- sweep(W, 2, sqrt(2 * pm * (1 - pm)), "/")
# G <- tcrossprod(W) / ncol(M)
# Ginv <- MASS::ginv(G)
# dimnames(G) <- dimnames(Ginv) <- list(lettuce_markers$gen, lettuce_markers$gen)
# attr(Ginv, "INVERSE") <- TRUE

# lettuce_G <- lettuce_GRM
# #
# require(devtools)
# load_all()
#
# model <- asreml::asreml(
#   fixed = y ~ rep,
#   random = ~ vm(gen, lettuce_G) * idv(loc) ,
#   data = lettuce_phenotypes)
#
# h2_Cullis(model, "gen" , marginal = FALSE, source = list(lettuce_G = lettuce_G), solver ="direct")
# .Machine$double.eps
#
# my_h2 <- heritable::h2(model, "gen", marginal = FALSE, source = list(lettuce_G = lettuce_G))
# my_ci <- confint(my_h2, B = 10, parallel  = "snow")
#
# my_h2 <- heritable::h2(model, "gen", marginal = FALSE, source = list(lettuce_G = lettuce_G))
# system.time(
#   {
#     confint(my_h2, B = 10)
#   }
# )
#
# my_h2 <- heritable::h2(model, "gen", marginal = FALSE, source = list(lettuce_G = lettuce_G), solver = "LMM")
# system.time(
#   {
#     confint(my_h2, B = 10)
#   }
# )


# model$G.param
#
# asreml::asreml.options(
# )
#
#
# model <- asreml::asreml(
#   fixed = y ~ rep,
#   random = ~ gen*loc ,
#   data = lettuce_phenotypes)
# my_h2 <- heritable::h2(model, "gen", marginal = FALSE)
# confint(my_h2, B = 10, parallel  = "snow", ncpus = 3)
# phrase_G
# phrase_vm
# ??H2
# model$G.param$`vm(gen, source = lettuce_GRM):loc`[[1]]
#
# check_GRM_exists(model, "gen", return = TRUE)
# lettuce_GRM <- lettuce_GRM
# rm(lettuce_GRM)
#
# pull_terms(model)
# pull_terms_without_specials(model)
#
# heritable::h2(model, "gen", marginal = FALSE, source = list(lettuce_GRM = lettuce_GRM))

# model <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc +vm(gen, lettuce_GRM) + gen:loc ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
# confint(h2(model, "gen", source = lettuce_GRM), B = 10)
#
# lettuce_GRM <- lettuce_GRM
# h2(model,"gen")
# rm(lettuce_GRM)
# h2(model,"gen")
# h2(model, "gen", source = lettuce_GRM)
#
#
# model <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc +vm(gen, Ginv) + gen:loc ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
# h2(model, "gen", source = G)
#
# model <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc +vm(gen, G, "PSD") + gen:loc ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
#
#
#
# })
