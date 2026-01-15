test_that("h2 heritability works", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()

  asreml_model_random <- readRDS(file = test_path("fixtures/asreml_model_random.rds"))
  asreml_model_grm <- readRDS(file = test_path("fixtures/asreml_model_grm.rds"))
  data("lettuce_GRM")

  # From Schmidt et al 2019 Fig 1.
  # Delta
  expect_equal(h2_Delta(asreml_model_grm, target = "gen", type = "BLUP"), 0.871, tolerance = 1e-3)
  # TODO: h2 BLUE values not matching
  expect_equal(h2_Delta(asreml_model_grm, target = "gen", type = "BLUE"), 0.818, tolerance = 1e-3)

  # Oakey
  # TODO: h2 Oakey values not matching
  expect_equal(h2_Oakey(asreml_model_grm, target = "gen"), 0.586, tolerance = 1e-3)

  # Structural checks
  expect_named(h2(asreml_model_grm, target = "gen"), c("Oakey", "Delta"))
  expect_s4_class(h2_Delta_pairwise(asreml_model_grm, target = "gen", type = "BLUP"), "dspMatrix")
  expect_type(h2_Delta_by_genotype(asreml_model_grm, target = "gen", type = "BLUP"), "list")
  expect_named(h2(asreml_model_grm, target = "gen"), c("Oakey", "Delta"))
  expect_length(h2(asreml_model_grm, target = "gen"), 2)
  expect_named(h2_Delta_by_genotype(asreml_model_grm, target = "gen", type = "BLUP"), levels(asreml_model_grm$mf$gen), ignore.order = TRUE)
  expect_equal(nrow(h2_Delta_pairwise(asreml_model_grm, target = "gen", type = "BLUP")), length(levels(asreml_model_grm$mf$gen)))
})

test_that("VanRadden GRM", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()

  # library(sommer)
  # lettuce_GRM_vanradden <- sommer::A.mat(as.matrix(lettuce_markers[, -1]))
  # dimnames(lettuce_GRM_vanradden) <- list(lettuce_markers$gen, lettuce_markers$gen)
  # saveRDS(lettuce_GRM_vanradden, file = test_path("fixtures/lettuce_GRM_vanradden.rds"))
  # lettuce_GRM_vanradden <- readRDS(file = test_path("fixtures/lettuce_GRM_vanradden.rds"))

  # Hand calculate vanRadden
  M <- as.matrix(lettuce_markers[, -1] + 1)
  N <- nrow(M)
  pm <- colSums(M) / (2 * N) # allele freq per marker (diploid X)
  pm <- pmin(pmax(pm, 1e-6), 1 - 1e-6) # guard against 0 or 1
  W <- sweep(M, 2, 2 * pm, "-")
  W <- sweep(W, 2, sqrt(2 * pm * (1 - pm)), "/")
  G <- tcrossprod(W) / ncol(M)
  Ginv <- MASS::ginv(G)
  dimnames(Ginv) <- list(lettuce_markers$gen, lettuce_markers$gen)
  attr(Ginv, "INVERSE") <- TRUE

  # # Crashes on Fonti's MBP
  # asreml_model_gr_vr <- asreml::asreml(y ~ rep,
  #                                      random = ~ vm(gen, Ginv),
  #                                      data = lettuce_phenotypes |>
  #                                        subset(loc == "L2")
  # )
  #
  # matrixcalc::is.singular.matrix(Ginv)
  # det(Ginv)
  # qr(Ginv)$rank < ncol(Ginv)
  #
  # saveRDS(asreml_model_gr_vr, test_path("fixtures/lettuce_asreml_vradden_grm.rds"))
  asreml_model_gr_vr <- readRDS(test_path("fixtures/lettuce_asreml_vradden_grm.rds"))

  # h2(asreml_model_gr_vr, "gen", source = Ginv)
})

test_that("Refactoring delta parameter functions works", {
  G_g <- matrix(c(
    0.5, 0.2, 0.2,
    0.2, 0.6, 0.3,
    0.2, 0.3, 0.7
  ), nrow = 3, byrow = TRUE)
  vd_matrix <- matrix(c(
    0.1, 0.15, 0.2,
    0.15, 0.12, 0.18,
    0.2, 0.18, 0.14
  ), nrow = 3, byrow = TRUE)

  expect_equal(
    h2_Delta_BLUP_parameters(G_g, vd_matrix),
    h2_Delta_parameters(G_g, vd_matrix, type = "BLUP")
  )

  expect_equal(
    h2_Delta_BLUE_parameters(G_g, vd_matrix),
    h2_Delta_parameters(G_g, vd_matrix, type = "BLUE")
  )

  vc_g <- 0.01
  vd_matrix <- matrix(c(NA, 0.2, 0.2, NA), 2, 2)

  expect_equal(
    H2_Delta_BLUE_parameters(vc_g, vd_matrix),
    H2_Delta_parameters(vc_g, vd_matrix, type = "BLUE")
  )

  expect_equal(
    H2_Delta_BLUP_parameters(vc_g, vd_matrix),
    H2_Delta_parameters(vc_g, vd_matrix, type = "BLUP")
  )
})

# test_that("Alternative way to get sigma a",{
#   model <- readRDS(file = test_path("fixtures/asreml_model_grm.rds"))
#   target <- "gen"
#
#   vm <- target_vm_term_asreml(model, target)
#   n_g <- model$noeff[[vm$target_vm]]
#   vc_g <- model$vparameters[[vm$target_vm]] * model$sigma2 * semivariance(vm$GRM)
#   # PEV
#   # Why is this giving me zero?
#   vdBLUP_mat <- predict(model,
#                     classify = vm$target_vm,
#                     only = vm$target_vm,
#                     vcov = TRUE,
#                     trace = FALSE
#   )$vcov
#
#   dim(vdBLUP_mat)
#
#   sigma_a_fernando_gonzales <- vc_g + (sum(diag(as.matrix(vdBLUP_mat))) / n_g)
#   model$vparameters
#
# }
# )

test_that("Try another package", {
  # install.packages("heritabilty")
  library(heritability)
  #
  data(LD)
  data(K_atwell)
  # ?K_atwell
  #
  # # Heritability estimation for all observations:
  # out <- marker_h2(data.vector=LD$LD,geno.vector=LD$genotype,
  #                 covariates=LD[,4:8],K=K_atwell)
  #
  # out$h2

  G_inv <- MASS::ginv(K_atwell)
  dimnames(G_inv) <- dimnames(K_atwell)

  LD_fit <- asreml(LD ~ rep1 + rep2 + rep3 + rep4 + rep5,
         random = ~vm(genotype, G_inv, singG = "PSD"),
         data = LD)


  h2(LD_fit, "genotype", source = G_inv) # error
  h2_Delta(LD_fit, "genotype")

})
