# # ASREML
# # Genotype as random effect
# asreml_model_random <- asreml::asreml(
#   fixed = yield ~ rep,
#   random = ~ gen + rep:block,
#   data = agridat::john.alpha,
#   trace = FALSE
# )
#
# saveRDS(asreml_model_random, file = test_path("fixtures/asreml_model_random.rds"))
#
# # Genotype as fixed effect
# asreml_model_fixed <- asreml::asreml(
#   fixed = yield ~ rep + gen,
#   random = ~ rep:block,
#   data = agridat::john.alpha,
#   trace = FALSE
# )
#
# saveRDS(asreml_model_fixed, file = test_path("fixtures/asreml_model_fixed.rds"))
#
# # Genotype as both fixed and random effect
# asreml_model_both <- asreml::asreml(
#   fixed = yield ~ gen,
#   random = ~ gen:block,
#   data = agridat::john.alpha,
#   trace = FALSE
# )
#
# saveRDS(asreml_model_both, file = test_path("fixtures/asreml_model_both.rds"))
#
# # # Model that is not convegd
# # asreml_model_failed_converge <- asreml_model_random
# # model_failed_converge$converge <- FALSE
#
# # GxE models --------------------------------------------------------------
# asremL_model_g_by_e <- asreml::asreml(yield ~ year:loc:trial,
#                         random = ~ gen + gen:loc,
#                         data = agridat::adugna.sorghum,
#                         trace = FALSE
# )
#
# saveRDS(asremL_model_g_by_e, file = test_path("fixtures/asreml_model_g_by_e.rds"))
#
# # Narrow sense
asreml_model_grm <- asreml::asreml(y ~ rep,
                        random=~ vm(gen, lettuce_GRM),
                        data = lettuce_phenotypes |>
                          subset(loc == "L2"))
#
# saveRDS(asreml_model_grm, file = test_path("fixtures/asreml_model_grm.rds"))
#
# # Different residual structure
# asreml_model_R <- asreml::asreml(yield ~ 1,
#                         random = ~ gen + col + row,
#                         residual = ~ ar1(col):ar1(row),
#                         data = agridat::gilmour.serpentine |>
#                           transform(
#                             col = as.factor(col),
#                             row = as.factor(row)
#                           ),
#                         trace = FALSE
# )
#
# saveRDS(asreml_model_R, file = test_path("fixtures/asreml_model_R.rds"))
#
# # lme4
#
# # Genotype as random effect
# lmer_model_random <- lme4::lmer(data = agridat::john.alpha, formula = yield ~ rep + (1 | gen) + (1 | rep:block))
#
# saveRDS(lmer_model_random, file = test_path("fixtures/lmer_model_random.rds"))
#
# # Genotype as fixed effect
# lmer_model_fixed <- lme4::lmer(data = agridat::john.alpha, formula = yield ~ rep + gen + (1 | rep:block))
#
# saveRDS(lmer_model_fixed, file = test_path("fixtures/lmer_model_fixed.rds"))
#
# # Genotype as both fixed and random effect
# lmer_model_g_by_e <- lme4::lmer(data = agridat::adugna.sorghum, formula = yield ~ year:loc:trial + (1 | gen) + (1 | gen:loc))
#
# saveRDS(lmer_model_g_by_e, file = test_path("fixtures/lmer_model_g_by_e.rds"))
#
# #TODO: Narrow sense for lme4
#
