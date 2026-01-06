# Setup
devtools::document()
devtools::load_all()
library(asreml)
library(dplyr)
library(lme4)
library(asremlPlus)
# install.packages("C:/Users/yidid/Downloads/asreml_4.2.0.392.zip", repos = NULL, type = "win.binary")
# library(asreml)
# asreml.license.activate()
lettuce_subset <- lettuce_phenotypes |>
  filter(loc == "L2")

# Test for lme4
lettuce_lme4 <- lmer(y ~ rep + (1 | gen), data = lettuce_subset)
H2_values_lme4 <- H2(lettuce_lme4, "gen", c("Cullis", "Standard"))
confint(H2_values_lme4)

# Test for asreml
N <- nrow(lettuce_subset)
lettuce_asreml <- asreml(
  fixed = y ~ rep,
  random =  ~ gen,
  data = lettuce_subset,
  trace = FALSE,
)
get_fixed_fit_asreml(lettuce_asreml)
bootstrap_asreml(lettuce_asreml, function(fit) fit$sigma2, nsim = 10, use.u = TRUE)
H2_values_asreml <- H2(lettuce_asreml, "gen", c("Cullis", "Standard"))
confint(H2_values_asreml)

# Model 1
pseudo_var1 <- sample(c("A","B"), size = N, replace = T) %>% as.factor
pseudo_var2 <- sample(c("A","B"), size = N, replace = T) %>% as.factor
lettuce_asreml <- asreml(
  fixed = y ~ rep * pseudo_var1,
  random =  ~ gen,
  sparse = ~ pseudo_var2,
  data = lettuce_subset,
  trace = FALSE
)

get_fixed_fit_asreml(lettuce_asreml)
asremlPlus::estimateV.asreml(lettuce_asreml)
tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
summary(lettuce_asreml)

# Model 2
lettuce_asreml <- asreml(
  fixed = y ~ rep,
  random =  ~ gen,
  data = lettuce_subset,
  trace = FALSE,
)

get_fixed_fit_asreml(lettuce_asreml)
asremlPlus::estimateV.asreml(lettuce_asreml)
tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
summary(lettuce_asreml)

# Model 3
pseudo_var1 <- rnorm(N)
lettuce_asreml <- asreml(
  fixed = y ~ spl(pseudo_var1, 5),
  random =  ~ gen,
  data = lettuce_subset,
  trace = FALSE,
)
get_fixed_fit_asreml(lettuce_asreml)
asremlPlus::estimateV.asreml(lettuce_asreml)
tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
summary(lettuce_asreml)

# Model 4
lettuce_asreml_grm <- asreml(
  fixed = y ~ loc,
  random = ~ vm(gen, lettuce_GRM) + rep,
  data = lettuce_phenotypes,
  trace = FALSE
)
get_fixed_fit_asreml(lettuce_asreml_grm)
asremlPlus::estimateV.asreml(lettuce_asreml_grm)
tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
summary(lettuce_asreml)

# Try to break
lettuce_asreml <- asreml(
  fixed = y ~ rep,
  random =  ~ ar1(gen),
  data = lettuce_subset,
  trace = FALSE,
)

bootstrap_asreml(lettuce_asreml, function(fit) fit$vparameters, nsim = 100, use.u = TRUE)

# h2
h2_values_asreml <- h2(lettuce_asreml, "gen")

#https://asreml.kb.vsni.co.uk/knowledge-base/converting-variance-ratios-v-predict/
