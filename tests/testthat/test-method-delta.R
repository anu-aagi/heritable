test_that("delta method works", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  # Genotype as random effect
  model_random <- asreml::asreml(
    fixed = yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  # Genotype as fixed effect
  model_fixed <- asreml::asreml(
    fixed = yield ~ rep + gen,
    random = ~ rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  # Genotype as both fixed and random effect
  model_both <- asreml::asreml(
    fixed = yield ~ gen,
    random = ~ gen:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  target <- "gen"

  expect_equal(H2_Delta(model_random, target), 0.8090841)
  expect_error(H2_Delta(model_fixed, target, aggregate = "harmonic"))
  expect_lt(H2_Delta(model_random, target, type = "BLUP", aggregate = "harmonic"), H2_Delta(model_random, target, type = "BLUP"))
  expect_error(h2_Delta(model_both, target, aggregate = "harmonic"))

  res_ls <- H2_Delta_by_genotype(model_random, target = "gen", type = "BLUP")
  expect_named(H2_Delta_by_genotype(model_random, target, type = "BLUP"), levels(model_random$mf$gen))
  expect_true(length(res_ls) == length(levels(model_random$mf[["gen"]])))

  H2_mat <- H2_Delta_pairwise(model_random, target = "gen", type = "BLUP")
  expect_true(is.matrix(as.matrix(H2_mat)))
  expected <- rowMeans(as.matrix(H2_mat), na.rm = TRUE)

  # compare numeric values
  expect_equal(as.numeric(unlist(res_ls)), as.numeric(expected), tolerance = 1e-7)
})
