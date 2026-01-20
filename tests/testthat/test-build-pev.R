require(asreml)

asreml.options(design = TRUE)

lettuce_subset <- lettuce_subset <- lettuce_phenotypes |>
     dplyr::filter(loc == "L2")

lettuce_asreml <- asreml(
  fixed = y ~ rep,
  random =  ~ gen,
  data = lettuce_subset,
  trace = FALSE,
)
N <- nrow(lettuce_subset)
Ng <- nrow(lettuce_asreml$coefficients$random)
PEV <- predict(lettuce_asreml,
                  classify = "gen",
                  only = "gen",
                  vcov = TRUE,
                  trace = FALSE
)$vcov


y <- lettuce_subset$y
G <- diag(lettuce_asreml$sigma2 * lettuce_asreml$vparameters["gen"], Ng)
R <- diag(lettuce_asreml$sigma2 * lettuce_asreml$vparameters["units!R"], N)
X <- lettuce_asreml$design[,rownames(lettuce_asreml$coefficients$fixed)]
Z <- lettuce_asreml$design[,rownames(lettuce_asreml$coefficients$random)]
V <- R + Z %*% G %*% t(Z)
Vinv <- solve(V)
P <- Vinv - Vinv %*% X %*% solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv
C22 <- G - G %*% t(Z) %*% P %*% Z %*% G


G - (diag(Ng) - tcrossprod(rep(1,Ng))/Ng) %*% G %*% t(Z) %*% P %*% Z %*% G %*% (diag(Ng) - tcrossprod(rep(1,Ng))/Ng) - PEV
