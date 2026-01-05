# Function
H2_confint <- function(heritable,
                       B = 100,
                       seed = NULL,
                       sample.random = TRUE,
                       level = 0.95,
                       type = c("basic", "norm", "perc"),
                       ...){
  UseMethod("H2_confint")
}


H2_confint.heritable <- function(heritable,
                                 B = 100,
                                 seed = NULL,
                                 sample.random = TRUE,
                                 level = 0.95,
                                 type = c("basic", "norm", "perc"),
                                 ...){

  # basic: bias corrected percentile interval
  # norm: bias corrected normal interval
  # perc: percentile interval
  type <- match.arg(type)

  method <- names(heritable)
  target <- attr(heritable, "target")
  model <- attr(heritable, "model") # Get the model

  H2_use <- function(x){
    heritable::H2(x, target, method, options = list(check = FALSE))
  }

  if(inherits(model, "lmerMod")){
    H2.boot <- lme4::bootMer(model, H2_use, nsim = B, seed = seed,
                             use.u = sample.random, ...)
  } else {

  }

  confint(H2.boot, level = level, type = type)
}



# Test
library(dplyr)
library(lme4)
library(asremlPlus)
# install.packages("C:/Users/yidid/Downloads/asreml_4.2.0.392.zip", repos = NULL, type = "win.binary")
# library(asreml)
# asreml.license.activate()

lettuce_subset <- lettuce_phenotypes |>
  filter(loc == "L2")
lettuce_lme4 <- lmer(y ~ rep + (1 | gen), data = lettuce_subset)
lettuce_asreml <- asreml(
  fixed = y ~ rep,
  random = ~ gen,
  data = lettuce_subset,
  trace = FALSE
)

H2_values_lme4 <- H2(lettuce_lme4, "gen", c("Cullis", "Standard"))
H2_confint(H2_values_lme4)

H2_values_asreml <- H2(lettuce_asreml, "gen", c("Cullis", "Standard"))
H2_confint(H2_values_asreml)


LT <- estimateV.asreml(lettuce_asreml, matrix = "V") %>%
  Matrix::Matrix() %>%
  Matrix::chol() %>%
  t
LT %*% rnorm(267)



(lettuce_subset$y - lettuce_asreml$residuals) %>% sort

asremlPlus::simulate.asreml

lettuce_asreml$linear.predictors %>% unique() %>% length
lettuce_asreml$coefficients$fixed %>% unique() %>% length
lettuce_asreml$coefficients$random %>% unique() %>% length


plot(lettuce_asreml$linear.predictors)
fitted(lettuce_asreml)

plot(fitted.values(lettuce_asreml) %>% sort)



lettuce_asreml <- asreml(
  fixed = y ~ rep*sample(c(1,2), 267, replace = T),
  random = ~ gen,
  data = lettuce_subset,
  trace = FALSE
)

unique(fitted.values(lettuce_asreml)) %>% length

lettuce_asreml$coefficients$random %>% unique

lettuce_asreml$coefficients %>% str

lettuce_asreml$coefficients
tcrossprod(model.matrix(~0+lettuce_asreml$mf$gen))
summary(lettuce_asreml)

model.matrix(formula(lettuce_asreml$call$fixed), data = as.data.frame(lettuce_asreml$mf))

estimateV.asreml


lettuce_asreml$call$random
lettuce_asreml$coefficients$random


estimateV.asreml(lettuce_asreml, )
estimateV.asreml


lettuce_asreml$formulae$random
formula(lettuce_asreml$call$random) %>% terms
simulate.asreml
