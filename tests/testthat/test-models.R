# test_that("Refactored asreml works",{
#   skip()
#
# lettuce_subset <- subset(lettuce_phenotypes, loc == "L2")
#
# model <-  lme4::lmer(y ~  rep + (0+ loc | gen ) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (loc | gen ) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1 | gen * loc) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1 | gen * loc) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1+ 1+ loc || gen ) + (1| gen), data = lettuce_phenotypes)
# model <-  lme4::lmer(y ~  rep + (1| gen), data = lettuce_subset)
# model <-  lme4::lmer(y ~  rep + (1 | gen * loc), data = lettuce_phenotypes)
# system.time(
#   confint(H2(model, "gen"), B = 20)
# )
# system.time(
#   confint(H2(model, "gen"), B = 20, parallel = "multicore", ncpus = 10)
# )
# system.time(
#   confint(H2(model, "gen"), B = 100, parallel = "snow", ncpus = 2)
# )


#
# h2_Standard.lmerMod(model, "gen")
# h2_Oakey.lmerMod(model, "gen")
# h2_Delta_pairwise.lmerMod(model, "gen", type = "BLUE")
# h2_Delta_pairwise.lmerMod(model, "gen", type = "BLUP")
#
#
# model <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc + gen + gen:loc ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
# system.time(
#   confint(H2(model, "gen"), B = 20)
# )
# system.time(
#   confint(H2(model, "gen"), B = 20, parallel = "multicore", ncpus = 10)
# )
# system.time(
#   confint(H2(model, "gen"), B = 100, parallel = "snow", ncpus = 3)
# )
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
#
#
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
