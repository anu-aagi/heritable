# test_that("Refactored asreml works",{
#   skip()
#
#   lettuce_subset <- subset(lettuce_phenotypes, loc == "L2")
#   lettuce_GRM <- lettuce_GRM
#
#   model <-  lme4::lmer(y ~  rep + (0+ loc | gen ) + (1| gen), data = lettuce_phenotypes)
#   model <-  lme4::lmer(y ~  rep + (loc | gen ) + (1| gen), data = lettuce_phenotypes)
#   model <-  lme4::lmer(y ~  rep + (1 | gen * loc) + (1| gen), data = lettuce_phenotypes)
#   model <-  lme4::lmer(y ~  rep + (1 | gen * loc) + (1| gen), data = lettuce_phenotypes)
#   model <-  lme4::lmer(y ~  rep + (1+ 1+ loc || gen ) + (1| gen), data = lettuce_phenotypes)
#   model <-  lme4::lmer(y ~  rep + (1| gen), data = lettuce_subset)
#
#   model <- asreml(
#     fixed = y ~ rep,
#     random =  ~ loc + gen + gen:loc ,
#     data = lettuce_phenotypes,
#     trace = FALSE,
#   )
#
#   model <- asreml(
#     fixed = y ~ rep,
#     random =  ~ loc + vm(gen, lettuce_GRM) + gen:loc ,
#     data = lettuce_phenotypes,
#     trace = FALSE,
#   )
#
#
# })
