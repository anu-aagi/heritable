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

test_that("Alternative way to get sigma a",{
  model <- readRDS(file = test_path("fixtures/asreml_model_grm.rds"))
  target <- "gen"

  vm <- target_vm_term_asreml(model, target)
  n_g <- model$noeff[[vm$target_vm]]
  vc_g <- model$vparameters[[vm$target_vm]] * model$sigma2 * semivariance(vm$GRM)
  # PEV
  # Why is this giving me zero?
  vdBLUP_mat <- predict(model,
                    classify = target,
                    only = target,
                    vcov = TRUE,
                    trace = FALSE
  )$vcov

  dim(vdBLUP_mat)

  sigma_a_fernando_gonzales <- vc_g + (sum(diag(as.matrix(vdBLUP_mat))) / n_g)
  model$vparameters

}
)

test_that("Try GPT simulation", {
  oakey_true_from_matrices <- function(X, Z, G, sigma_g2, sigma_e2, tol = 1e-10) {
    n <- ncol(Z)

    # R^{-1} for iid residuals
    Rinv <- diag(1 / sigma_e2, nrow(Z))

    # (sigma_g2 * G)^{-1} using eigen (handles PSD)
    eg <- eigen(G, symmetric = TRUE)
    d <- eg$values
    U <- eg$vectors
    dinv <- ifelse(d > tol, 1 / d, 0)
    Ginv <- U %*% diag(dinv) %*% t(U)
    Ginv <- (Ginv + t(Ginv)) / 2

    Ginvg <- (1 / sigma_g2) * Ginv

    XtRinvX <- t(X) %*% Rinv %*% X
    XtRinvZ <- t(X) %*% Rinv %*% Z
    ZtRinvX <- t(Z) %*% Rinv %*% X
    ZtRinvZ <- t(Z) %*% Rinv %*% Z

    K <- rbind(
      cbind(XtRinvX, XtRinvZ),
      cbind(ZtRinvX, ZtRinvZ + Ginvg)
    )

    Kinv <- MASS::ginv(K)
    p <- ncol(X)
    Cuu <- Kinv[(p + 1):(p + n), (p + 1):(p + n)]
    Cuu <- (Cuu + t(Cuu)) / 2  # symmetrize

    # average PEV of pairwise differences
    # mean over i<j of (Cii + Cjj - 2Cij)
    diagC <- diag(Cuu)
    # Efficient mean pairwise difference variance:
    # mean_{i<j}(Cii + Cjj - 2Cij)
    nG <- length(diagC)
    sumCiiCjj <- (nG - 1) * sum(diagC) * 2  # sum over ordered pairs i!=j of (Cii + Cjj)
    sumOff <- sum(Cuu) - sum(diagC)         # sum_{i!=j} Cij / 2? careful: Cuu includes both
    # Let's do it directly robustly:
    idx <- which(upper.tri(Cuu))
    pev_diff_bar <- mean(diagC[row(Cuu)[idx]] + diagC[col(Cuu)[idx]] - 2 * Cuu[idx])

    H2 <- 1 - pev_diff_bar / (2 * sigma_g2)

    list(H2 = H2, pev = Cuu, pev_diff_bar = pev_diff_bar)
  }


  simulate_grm_data <- function(n_gen = 50, n_markers = 300, n_rep = 3,
                                sigma_g2 = 1.0, sigma_e2 = 1.0, seed = 1) {
    set.seed(seed)

    # Marker-derived GRM (simple)
    M <- matrix(rbinom(n_gen * n_markers, 2, 0.5), n_gen, n_markers)
    M <- scale(M, center = TRUE, scale = TRUE)
    G <- tcrossprod(M) / n_markers
    G <- (G + t(G)) / 2

    # Design: n_rep obs per genotype
    N <- n_gen * n_rep
    geno <- factor(rep(seq_len(n_gen), each = n_rep))

    Z <- model.matrix(~ 0 + geno)   # N x n_gen
    X <- X <- matrix(1, nrow = nrow(Z), ncol = 1) # intercept only

    # Simulate u ~ N(0, sigma_g2 * G)
    eg <- eigen(G, symmetric = TRUE)
    u <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0))) %*% rnorm(n_gen)
    u <- as.numeric(u) * sqrt(sigma_g2)

    e <- rnorm(N, sd = sqrt(sigma_e2))
    y <- drop(Z %*% u) + rnorm(nrow(Z), sd = sqrt(sigma_e2))


    # browser()
    dat <- data.frame(y = y, genotype = factor(paste0("geno",geno)))
    list(dat = dat, G = G, X = X, Z = Z,
         sigma_g2 = sigma_g2, sigma_e2 = sigma_e2)
  }


  sim <- simulate_grm_data()

  truth <- oakey_true_from_matrices(sim$X, sim$Z, sim$G, sim$sigma_g2, sim$sigma_e2)
  truth$H2

  G_inv <- MASS::ginv(sim$G)
  dimnames(G_inv) <- list(dimnames(sim$Z)[[2]], dimnames(sim$Z)[[2]])

  model <- asreml(y ~ 1,
                  random = ~ vm(genotype, G_inv, singG="PSD"),
                  data = sim$dat,
                  ai.sing = TRUE)


  # Compare:
  c(true = truth$H2, estimated = h2_Oakey(model, "genotype"))

  G_inv <- MASS::ginv(sim$G)
  Gg_inv_true <- (1/sim$sigma_g2) * G_inv

  H2_eig_true <- H2_Oakey_parameters(Gg_inv_true, truth$pev)
})
