test_that("asreml PEV definition",{
  skip()
  require(asreml)

  asreml.options(design = TRUE)

  lettuce_subset <- lettuce_subset <- lettuce_phenotypes |>
       dplyr::filter(loc == "L2")

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ vm(gen, lettuce_GRM),
    data = lettuce_subset,
    trace = FALSE,
  )
  N <- nrow(lettuce_subset)
  Ng <- nrow(lettuce_asreml$coefficients$random)
  PEV <- predict(lettuce_asreml,
                    classify = "vm(gen, lettuce_GRM)",
                    only = "vm(gen, lettuce_GRM)",
                    vcov = TRUE,
                    trace = FALSE
  )$vcov



  y <- lettuce_subset$y
  G <- lettuce_asreml$sigma2 * lettuce_asreml$vparameters["vm(gen, lettuce_GRM)"] * lettuce_GRM
  R <- diag(lettuce_asreml$sigma2 * lettuce_asreml$vparameters["units!R"], N)
  X <- lettuce_asreml$design[,rownames(lettuce_asreml$coefficients$fixed)]
  Z <- lettuce_asreml$design[,rownames(lettuce_asreml$coefficients$random)]
  V <- R + Z %*% G %*% t(Z)
  Vinv <- MASS::ginv(as.matrix(V))
  P <- Vinv - Vinv %*% X %*% MASS::ginv(as.matrix(t(X) %*% Vinv %*% X)) %*% t(X) %*% Vinv
  C22 <- G - G %*% t(Z) %*% P %*% Z %*% G


  expect_true(
    all(C22 - PEV < 1e-6)
  )
})
