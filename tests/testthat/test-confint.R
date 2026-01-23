# Setup
# install.packages("C:/Users/yidid/Downloads/asreml_4.2.0.392.zip", repos = NULL, type = "win.binary")
# library(asreml)
# asreml.license.activate()

test_that("Confint works",{
  skip_if_not_installed("asreml")
  skip_on_cran()
  skip_on_ci()
  skip()

  require(asreml)

# lettuce_subset <- lettuce_phenotypes |>
#   dplyr::filter(loc == "L2")
#
# lme4
# lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
  lettuce_lme4 <- readRDS(test_path("fixtures/confint_lettuce_lme4.rds"))

# boot.obj <- lme4::bootMer(lettuce_lme4,
#                           FUN = function(fit) fit@beta,
#                           nsim = 2000, seed = 2,
#                           use.u = TRUE
# )
# boot::boot.ci(boot.obj, type = "basic")

H2_values_lme4 <- H2(lettuce_lme4, "gen", c("Standard")) # H2
ci_lme4 <- confint(H2_values_lme4, B = 1000, seed = 1)
ci_fix_lme4 <- confint(H2_values_lme4, B = 1000, seed = 1, random_effect = "conditional")

# attr(ci_lme4 , "boot_mod") |> boot::boot.ci(type = "basic")
# hist(attr(ci_lme4 , "boot_mod")$t)

# asreml
# N <- nrow(lettuce_subset)

# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ gen,
#   data = lettuce_subset,
#   trace = FALSE,
# )

lettuce_asreml <- readRDS(test_path("fixtures/confint_lettuce_asreml.rds"))

# get_fixed_fit_asreml(lettuce_asreml)
# boot.obj <- bootstrap_asreml(lettuce_asreml, function(fit) fit$sigma2, nsim = 10, use.u = TRUE)

H2_values_asreml <- H2(lettuce_asreml, "gen", c("Standard")) # H2
ci_asreml <- confint(H2_values_asreml, B = 1000, seed = 1)
ci_fix_asreml <- confint(H2_values_asreml, B = 1000, seed = 1, random_effect = "conditional")
# attr(ci_asreml, "boot_mod") |> boot::boot.ci(type = "basic")
# hist(attr(ci_asreml , "boot_mod")$t)

# Compare results
plot_df <- data.frame(lme4_resample =  attr(ci_lme4, "boot_mod")$t |> drop(),
                      asreml_resample = attr(ci_asreml, "boot_mod")$t |> drop(),
                      lme4_fix = attr(ci_fix_lme4, "boot_mod")$t |> drop(),
                      asreml_fix = attr(ci_fix_asreml, "boot_mod")$t |> drop())


expect_true(ks.test(plot_df$lme4_fix, plot_df$asreml_fix)$p.value > 0.05 & ks.test(plot_df$lme4_resample, plot_df$asreml_resample)$p.value > 0.05)

# plot_df <- reshape2::melt(plot_df)
# library(ggplot2)
# ggplot(plot_df) + geom_boxplot(aes(value, variable, fill = variable)) +
#   scale_fill_discrete("")+
#   theme_classic(base_size = 13)+
#   geom_vline(xintercept = H2_values_asreml, color = "red", linetype = "dashed")+
#   labs(x = "Bootstraped H2", y = "")
})



######################## bootstrap_asreml checks ###############################

#
# # Model 1
# pseudo_var1 <- sample(c("A","B"), size = N, replace = T) |> as.factor
# pseudo_var2 <- sample(c("A","B"), size = N, replace = T) |> as.factor
# lettuce_asreml <- asreml(
#   fixed = y ~ rep * pseudo_var1,
#   random =  ~ gen,
#   sparse = ~ pseudo_var2,
#   data = lettuce_subset,
#   trace = FALSE
# )
#
# get_fixed_fit_asreml(lettuce_asreml)
# asremlPlus::estimateV.asreml(lettuce_asreml)
# tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
# summary(lettuce_asreml)
#
# # Model 2
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ gen,
#   data = lettuce_subset,
#   trace = FALSE,
# )
#
# get_fixed_fit_asreml(lettuce_asreml)
# asremlPlus::estimateV.asreml(lettuce_asreml)
# tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
# summary(lettuce_asreml)
#
# # Model 3
# pseudo_var1 <- rnorm(N)
# lettuce_asreml <- asreml(
#   fixed = y ~ spl(pseudo_var1, 5),
#   random =  ~ gen,
#   data = lettuce_subset,
#   trace = FALSE,
# )
# get_fixed_fit_asreml(lettuce_asreml)
# asremlPlus::estimateV.asreml(lettuce_asreml)
# tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
# summary(lettuce_asreml)
#
# # Model 4
# lettuce_asreml_grm <- asreml(
#   fixed = y ~ loc,
#   random = ~ vm(gen, lettuce_GRM) + rep,
#   data = lettuce_phenotypes,
#   trace = FALSE
# )
# get_fixed_fit_asreml(lettuce_asreml_grm)
# asremlPlus::estimateV.asreml(lettuce_asreml_grm)
# tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
# summary(lettuce_asreml)
#
# # Try to break
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ ar1(gen),
#   data = lettuce_subset,
#   trace = FALSE,
# )
#
# bootstrap_asreml(lettuce_asreml, function(fit) fit$vparameters, nsim = 100, use.u = TRUE)
#
# # h2
# h2_values_asreml <- h2(lettuce_asreml, "gen")
#
# #https://asreml.kb.vsni.co.uk/knowledge-base/converting-variance-ratios-v-predict/
