test_that("H2_Reliability_parameters computes r^2 = 1 - diag(C22)/diag(G)", {
  G_g <- diag(c(0.15, 0.20, 0.10), 3, 3)
  C22_g <- matrix(
    c(
      0.08, 0.01, 0.00,
      0.01, 0.07, 0.01,
      0.00, 0.01, 0.09
    ),
    nrow = 3, byrow = TRUE
  )

  expected <- 1 - diag(C22_g) / diag(G_g)
  expect_equal(H2_Reliability_parameters(G_g, C22_g), expected)

  # Off-diagonals are ignored (only diagonals matter)
  C22_no_offdiag <- diag(diag(C22_g))
  expect_equal(
    H2_Reliability_parameters(G_g, C22_g),
    H2_Reliability_parameters(G_g, C22_no_offdiag)
  )
})

test_that("overall reliability reproduces Schmidt et al. (2019), Fig. 1, Example 1", {
  # Independent ground truth (not a self-consistency check): Example 1 of
  # Schmidt et al. (2019) analyses the `agridat::john.alpha` oat trial *as if*
  # it were a CRD (Eq. 16), i.e. ignoring rep / incomplete-block structure:
  #   y = mu + g_i + e.
  # For that balanced setting the paper reports (Fig. 1) an overall reliability
  # of 0.556, constant across genotypes, alongside Standard = Cullis = 0.58.
  skip_if_not_installed("agridat")
  skip_if_not_installed("lme4")

  data("john.alpha", package = "agridat")
  model <- lme4::lmer(yield ~ (1 | gen), data = john.alpha)

  # r-bar^2 to the paper's printed precision (3 dp).
  expect_equal(round(H2_Reliability(model, "gen"), 3), 0.556)

  # r^2_i is constant across genotypes and equals r-bar^2 (paper: "constant ~0.556").
  r2i <- H2_Reliability_by_genotype(model, "gen")
  expect_equal(round(unname(r2i), 3), rep(0.556, length(r2i)))

  # Anchor that the model matches the paper's Example 1 setup (Standard = 0.58, 2 dp).
  expect_equal(round(H2_Standard(model, "gen"), 2), 0.58)
})

test_that("reliability (lme4) matches manual PEV/G computation", {
  lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
  model <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)

  vc <- var_comp(model, "gen", calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE)
  manual <- 1 - diag(as.matrix(vc$C22_g)) / diag(as.matrix(vc$G_g))

  by_geno <- H2_Reliability_by_genotype(model, "gen")

  expect_equal(unname(by_geno), unname(manual))
  expect_equal(names(by_geno), vc$gnames)
  expect_length(by_geno, vc$n_g)
})

test_that("overall reliability is the mean of the per-genotype values", {
  lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
  model <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)

  by_geno <- H2_Reliability_by_genotype(model, "gen")

  expect_equal(H2_Reliability(model, "gen"), mean(by_geno))
  # h2_ and H2_ share the backend and agree for lme4 models
  expect_equal(h2_Reliability(model, "gen"), H2_Reliability(model, "gen"))
})

test_that("reliability is bounded in [0, 1]", {
  lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
  model <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)

  r2 <- H2_Reliability(model, "gen")
  expect_gte(r2, 0)
  expect_lte(r2, 1)
})

test_that("reliability sits just below Cullis for this balanced dataset", {
  # Regression pin, NOT a general theorem: for iid genotypes the exact
  # relation is Cullis - mean(r^2) = mean(off-diagonal C22) / sigma^2_g, so
  # r^2 <= Cullis holds only when the mean off-diagonal prediction-error
  # covariance is non-negative -- true for this balanced RCBD, not guaranteed
  # in general (cf. Schmidt et al. 2019, Fig. 1, which shows the same ordering
  # for their balanced examples).
  lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
  model <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)

  expect_lte(H2_Reliability(model, "gen"), H2_Cullis(model, "gen"))
})

test_that("Reliability is reachable via the h2()/H2() wrappers", {
  lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
  model <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)

  H2_out <- H2(model, "gen", method = c("Cullis", "Reliability"))
  expect_named(H2_out, c("Cullis", "Reliability"))
  expect_equal(unname(H2_out["Reliability"]), H2_Reliability(model, "gen"))

  h2_out <- h2(model, "gen", method = "Reliability")
  expect_equal(unname(h2_out["Reliability"]), h2_Reliability(model, "gen"))
})

test_that("Reliability is not part of the default method set", {
  lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
  model <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)

  expect_false("Reliability" %in% names(H2(model, "gen")))
  expect_false("Reliability" %in% names(h2(model, "gen")))
})

test_that("reliability (asreml) matches manual PEV/G computation", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()

  lettuce_asreml <- readRDS(test_path("fixtures/lettuce_asreml.rds"))

  vc <- var_comp(lettuce_asreml, "gen", calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE)
  manual <- 1 - diag(as.matrix(vc$C22_g)) / diag(as.matrix(vc$G_g))

  by_geno <- H2_Reliability_by_genotype(lettuce_asreml, "gen")

  expect_equal(unname(by_geno), unname(manual))
  expect_equal(H2_Reliability(lettuce_asreml, "gen"), mean(by_geno))
})
