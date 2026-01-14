test_that("Check GRM specification", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  skip_on_ci()

  require(asreml)

  M <- as.matrix(lettuce_markers[, -1] + 1)
  N <- nrow(M)
  pm <- colSums(M) / (2 * N) # allele freq per marker (diploid X)
  pm <- pmin(pmax(pm, 1e-6), 1 - 1e-6) # guard against 0 or 1
  W <- sweep(M, 2, 2 * pm, "-")
  W <- sweep(W, 2, sqrt(2 * pm * (1 - pm)), "/")
  G <- tcrossprod(W) / ncol(M)
  Ginv <- MASS::ginv(G)
  dimnames(G) <- dimnames(Ginv) <- list(lettuce_markers$gen, lettuce_markers$gen)
  attr(Ginv, "INVERSE") <- TRUE

  # Error in asreml::asreml(y ~ rep, random = ~vm(gen, G), data = subset(lettuce_phenotypes,  :
  # Error   : Iteration failed; ifault: 2020
  #  Warning :  Unexpected data when GRM diagonal element is ZERO
  model1 <- asreml::asreml(y ~ rep,
                                       random = ~ vm(gen, Ginv),
                                       data = lettuce_phenotypes |>
                                         subset(loc == "L2")
  )



  Ginv2 <- data.frame(Row = as.vector(row(Ginv)[!upper.tri(Ginv)]),
                      Column = as.vector(col(Ginv)[!upper.tri(Ginv)]),
                      value = as.vector(Ginv[!upper.tri(Ginv)])) |>
    dplyr::arrange(Row, Column)
  attr(Ginv2, "rowNames") <- lettuce_markers$gen
  attr(Ginv2, "INVERSE") <- TRUE

  model2 <- asreml::asreml(y ~ rep,
                           random = ~ vm(gen, Ginv2),
                           data = lettuce_phenotypes |>
                             subset(loc == "L2")
  )


  model3 <- asreml::asreml(y ~ rep,
                           random = ~ vm(gen, G, singG = "PSD"),
                           data = lettuce_phenotypes |>
                             subset(loc == "L2")
  )
  expect_equal(
    model1$sigma2,
    model2$sigma2,
    tolerance = 1e-5
  )
})
