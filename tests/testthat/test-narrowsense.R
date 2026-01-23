test_that("h2 heritability works", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()
  skip()

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
  expect_equal(heritable:::h2_Oakey(asreml_model_grm, target = "gen"), 0.586, tolerance = 1e-3)

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
#                     classify = target,
#                     only = target,
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

test_that("Try GPT simulation", {
  # oakey_true_from_matrices_eigen <- function(X, Z, G_marker, sigma_g2, sigma_e2,
  #                                            tol_eig = 1e-8, tol_G = 1e-10) {
  #   n <- nrow(Z)
  #   m <- ncol(Z)
  #
  #   # Genetic covariance for g: G = sigma_g2 * G_marker
  #   G <- sigma_g2 * G_marker
  #   G <- (G + t(G))/2
  #
  #   # V = R + Z G Z'
  #   V <- diag(sigma_e2, n) + Z %*% G %*% t(Z)
  #   Vinv <- solve(V)
  #
  #   # Pv = Vinv - Vinv X (X' Vinv X)^-1 X' Vinv
  #   XtVinvX <- t(X) %*% Vinv %*% X
  #   XtVinvX_inv <- solve(XtVinvX)
  #   Pv <- Vinv - Vinv %*% X %*% XtVinvX_inv %*% t(X) %*% Vinv
  #   Pv <- (Pv + t(Pv))/2
  #
  #   # Build symmetric similar matrix S = G^{1/2} Z' Pv Z G^{1/2}
  #   eg <- eigen(G, symmetric = TRUE)
  #   d <- pmax(eg$values, 0)
  #   U <- eg$vectors
  #   Ghalf <- U %*% diag(sqrt(d)) %*% t(U)
  #
  #   S <- Ghalf %*% (t(Z) %*% Pv %*% Z) %*% Ghalf
  #   S <- (S + t(S))/2
  #
  #   lam <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  #
  #   # Oakey: drop (near) zero eigenvalues and average remaining
  #   lam_pos <- lam[lam > tol_eig]
  #   H2 <- mean(lam_pos)
  #
  #   list(H2 = H2, eigenvalues = lam, kept = lam_pos)
  # }
  #
  #
  # simulate_grm_data <- function(n_gen = 50, n_markers = 300, n_rep = 3,
  #                               sigma_g2 = 1.0, sigma_e2 = 1.0, seed = 1) {
  #   set.seed(seed)
  #
  #   # Marker-derived GRM (simple)
  #   M <- matrix(rbinom(n_gen * n_markers, 2, 0.5), n_gen, n_markers)
  #   M <- scale(M, center = TRUE, scale = TRUE)
  #   G <- tcrossprod(M) / n_markers
  #   G <- (G + t(G)) / 2
  #
  #   # Design: n_rep obs per genotype
  #   N <- n_gen * n_rep
  #   geno <- factor(rep(seq_len(n_gen), each = n_rep))
  #
  #   Z <- model.matrix(~ 0 + geno)   # N x n_gen
  #   X <- matrix(1, nrow = nrow(Z), ncol = 1) # intercept only
  #
  #   # Simulate u ~ N(0, sigma_g2 * G)
  #   eg <- eigen(G, symmetric = TRUE)
  #   u <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0))) %*% rnorm(n_gen)
  #   u <- as.numeric(u) * sqrt(sigma_g2)
  #
  #   e <- rnorm(N, sd = sqrt(sigma_e2))
  #   y <- drop(Z %*% u) + rnorm(nrow(Z), sd = sqrt(sigma_e2))
  #
  #
  #   dat <- data.frame(y = y, genotype = factor(paste0("geno",geno)))
  #   list(dat = dat, G = G, X = X, Z = Z,
  #        sigma_g2 = sigma_g2, sigma_e2 = sigma_e2)
  # }
  #
  #
  # sim <- simulate_grm_data()
  #
  # truth <- oakey_true_from_matrices_eigen(sim$X, sim$Z, sim$G, sim$sigma_g2, sim$sigma_e2)
  # truth$H2
  #
  # G_inv <- MASS::ginv(sim$G)
  # dimnames(G_inv) <- list(dimnames(sim$Z)[[2]], dimnames(sim$Z)[[2]])
  #
  # model <- asreml::asreml(y ~ 1,
  #                 random = ~ vm(genotype, G_inv, singG="PSD"),
  #                 data = sim$dat)
  #
  #
  # # Compare:
  # c(true = truth$H2, estimated = h2_Oakey(model, "genotype"))
  # #  true estimated
  # # 0.7248769 0.6755711
  #
  # Gg_inv_true <- (1/sim$sigma_g2 * sim$sigma_e2) * G_inv
  #
  # # Substitute truth from simulation and see if h2 can recover from model
  # model_plugged_in <- model
  # model_plugged_in$G.param$`vm(genotype, G_inv, singG = "PSD")`$variance$initial <- sim$sigma_g2
  # model_plugged_in$R.param$units$variance$initial  <- sim$sigma_e2
  #
  # # asreml.options(fixgammas = TRUE)
  # # model_plugged_in <- update(model_plugged_in)     # refit, but keep variance parameters fixed
  # # asreml.options(fixgammas = FALSE)  # turn it back off afterwards
  #
  # c(true = truth$H2,
  #   fixedVC = h2_Oakey(model_plugged_in, "genotype"))
  # # true   fixedVC
  # # 0.7248769 0.2702343

})
