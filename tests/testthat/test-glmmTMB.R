# The glmmTMB backend must reproduce the lme4 backend exactly: fitting the same
# Gaussian model with lme4 (REML) and glmmTMB (REML) has to yield numerically
# identical heritability for every estimator. This is a self-contained ground
# truth that needs no external reference values.

fit_pair <- function(formula) {
  list(
    lme4 = lme4::lmer(formula, data = john.alpha),
    tmb  = glmmTMB::glmmTMB(formula, data = john.alpha, REML = TRUE)
  )
}

expect_backends_equal <- function(fn, models, ..., tolerance = 1e-4) {
  r_l <- suppressMessages(fn(models$lme4, "gen", ...))
  r_t <- suppressMessages(fn(models$tmb, "gen", ...))
  expect_equal(r_t, r_l, tolerance = tolerance, ignore_attr = TRUE)
}

# lme4 exposes the relative covariance factor Lambda, so the unconditional
# random-effect covariance G = sigma^2 * Lambda Lambda^T can be built WITHOUT any
# heritable code. This is an independent oracle for build_G_glmmTMB(), which must
# reconstruct the same G from VarCorr()$cond because glmmTMB exposes no Lambda
# (and getME(., "Gp") errors). The two backends order the columns of Z
# differently, so we compare the order-invariant spectrum (sorted eigenvalues).
G_from_lme4 <- function(m) {
  as.matrix(stats::sigma(m)^2 * Matrix::tcrossprod(lme4::getME(m, "Lambda")))
}
spectrum <- function(M) {
  sort(eigen(as.matrix(M), only.values = TRUE, symmetric = TRUE)$values)
}

john.alpha <- within(agridat::john.alpha, {
  gen   <- factor(gen)
  rep   <- factor(rep)
  block <- factor(block)
  x     <- as.numeric(scale(plot)) # deterministic numeric covariate for slopes
})

test_that("glmmTMB backend reproduces lme4 for a single random effect", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  m <- fit_pair(yield ~ rep + (1 | gen))

  expect_backends_equal(H2_Standard, m)
  expect_backends_equal(H2_Cullis,   m)
  expect_backends_equal(H2_Oakey,    m)
  expect_backends_equal(H2_Piepho,   m)
  expect_backends_equal(H2_Delta_pairwise, m, type = "BLUP")
  expect_backends_equal(H2_Delta_pairwise, m, type = "BLUE")
})

test_that("glmmTMB backend reproduces lme4 for multiple random effects", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  m <- fit_pair(yield ~ rep + (1 | gen) + (1 | block))

  expect_backends_equal(H2_Standard, m)
  expect_backends_equal(H2_Cullis,   m)
  expect_backends_equal(H2_Oakey,    m)
  expect_backends_equal(H2_Piepho,   m)
})

test_that("glmmTMB backend reproduces lme4 for an interaction random effect", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  m <- fit_pair(yield ~ rep + (1 | gen) + (1 | gen:block))

  expect_backends_equal(H2_Standard, m)
  expect_backends_equal(H2_Cullis,   m)
  expect_backends_equal(H2_Oakey,    m)
  expect_backends_equal(H2_Piepho,   m)
})

test_that("narrow-sense estimators also match across backends", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  m <- fit_pair(yield ~ rep + (1 | gen))

  expect_backends_equal(h2_Standard, m)
  expect_backends_equal(h2_Cullis,   m)
  expect_backends_equal(h2_Oakey,    m)
  expect_backends_equal(h2_Piepho,   m)
})

test_that("unsupported glmmTMB model structures are rejected informatively", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  d <- john.alpha

  # This first backend supports only plain random-effect intercepts (1 | group).
  # Non-Gaussian extras and any non-intercept random-effect structure must fail
  # early with an informative error rather than a cryptic downstream failure or a
  # silently wrong result. Random slopes and structured G-side covariance are the
  # deferred follow-up; here we only pin that they are rejected, not computed.
  zi    <- glmmTMB::glmmTMB(yield ~ rep + (1 | gen), data = d,
                            ziformula = ~1, REML = TRUE)
  disp  <- glmmTMB::glmmTMB(yield ~ rep + (1 | gen), data = d,
                            dispformula = ~rep, REML = TRUE)
  slope <- suppressWarnings(
    glmmTMB::glmmTMB(yield ~ rep + (1 + x | gen), data = d, REML = TRUE))
  dbar  <- suppressWarnings(
    glmmTMB::glmmTMB(yield ~ rep + (1 + x || gen), data = d, REML = TRUE))

  # A structured G-side covariance on a (here non-target) effect, representative
  # of cs()/ar1()/toep()/us(): all share the same `blockCode != "us"` rejection.
  set.seed(1)
  P <- 4L; B <- 20L
  ds     <- expand.grid(p = factor(seq_len(P)), b = factor(seq_len(B)))
  ds$gen <- factor(sample(rep(seq_len(10), length.out = nrow(ds))))
  ds$y   <- stats::rnorm(nrow(ds))
  cs_m   <- suppressWarnings(
    glmmTMB::glmmTMB(y ~ (1 | gen) + cs(p + 0 | b), data = ds, REML = TRUE))

  expect_error(H2(zi,    "gen"), "[Zz]ero-inflat")
  expect_error(H2(disp,  "gen"), "dispformula")
  expect_error(H2(slope, "gen"), "[Ii]ntercept")
  expect_error(H2(dbar,  "gen"), "[Ii]ntercept")
  expect_error(H2(cs_m,  "gen"), "[Ii]ntercept")
})

test_that("intercept-only guard closes the single-column slope and diag holes", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  d <- john.alpha

  # ACCEPT: for a single column, diag() is identical to (1 | gen), so it must be
  # accepted and give the same heritability, not rejected as "structured".
  m_diag <- suppressWarnings(
    glmmTMB::glmmTMB(yield ~ rep + diag(1 | gen), data = d, REML = TRUE))
  m_ref  <- glmmTMB::glmmTMB(yield ~ rep + (1 | gen), data = d, REML = TRUE)
  expect_equal(
    suppressMessages(as.numeric(H2_Standard(m_diag, "gen"))),
    suppressMessages(as.numeric(H2_Standard(m_ref,  "gen"))),
    tolerance = 1e-6
  )

  # REJECT: a lone single-column random slope on the target evades a naive
  # "single-column us" gate; it must be rejected with the intercept-only message,
  # not mis-diagnosed later as "target not found".
  slope_t <- suppressWarnings(
    glmmTMB::glmmTMB(yield ~ x + (0 + x | gen), data = d, REML = TRUE))
  expect_error(H2(slope_t, "gen"), "[Ii]ntercept")

  # REJECT: the pre-expanded uncorrelated form (1 | gen) + (0 + x | gen) is the
  # same model as (1 + x || gen) and previously slipped through as two
  # single-column blocks and crashed cryptically downstream.
  dbar_exp <- suppressWarnings(
    glmmTMB::glmmTMB(yield ~ x + (1 | gen) + (0 + x | gen), data = d, REML = TRUE))
  expect_error(H2(dbar_exp, "gen"), "[Ii]ntercept")

  # REJECT: a single-column slope on a NON-target effect is outside the
  # intercepts-only scope and must be rejected too (its estimate happened to
  # match lme4, but the backend does not promise this class yet).
  slope_side <- suppressWarnings(
    glmmTMB::glmmTMB(yield ~ x + (1 | gen) + (0 + x | rep), data = d, REML = TRUE))
  expect_error(H2(slope_side, "gen"), "[Ii]ntercept")
})

test_that("prior weights warn but do not block glmmTMB estimation", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  # Prior weights are not yet accounted for (residual is treated as sigma^2 * I);
  # consistent handling across backends is tracked upstream (#52). For now the
  # backend warns and still returns an estimate rather than hard-blocking.
  set.seed(1)
  d <- john.alpha
  d$w <- runif(nrow(d), 0.5, 2)
  wt  <- glmmTMB::glmmTMB(yield ~ rep + (1 | gen), data = d,
                          weights = w, REML = TRUE)

  expect_warning(res <- H2(wt, "gen"), "weights")
  expect_true(all(is.finite(as.numeric(res))))
})

test_that("non-Gaussian glmmTMB models are rejected with an informative error", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  pois <- glmmTMB::glmmTMB(round(yield) ~ rep + (1 | gen),
                           data = john.alpha, family = stats::poisson())

  expect_error(H2(pois, "gen"), "Gaussian")
})

test_that("build_G_glmmTMB reconstructs the same G as lme4's Lambda route", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("agridat")

  # Direct proof of the G reconstruction itself (not just the downstream h2):
  # build_G_glmmTMB() must equal the G that lme4 yields via its own Lambda factor
  # on the same data. Compared spectrally because the two backends order Z
  # differently (glmmTMB keeps formula order, lme4 sorts terms by n levels).
  scenarios <- list(
    yield ~ rep + (1 | gen),
    yield ~ rep + (1 | gen) + (1 | block),
    yield ~ rep + (1 | gen) + (1 | gen:block)
  )
  for (f in scenarios) {
    ml <- lme4::lmer(f, data = john.alpha)
    mt <- glmmTMB::glmmTMB(f, data = john.alpha, REML = TRUE)
    Gt <- build_G_glmmTMB(mt)
    Gl <- G_from_lme4(ml)
    expect_equal(dim(Gt), dim(Gl))
    expect_equal(spectrum(Gt), spectrum(Gl), tolerance = 1e-4)
  }
})

test_that("build_G_glmmTMB carries the off-diagonal covariance of a k=2 block", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")

  # A 1x1 intercept Sigma cannot reveal a wrong kronecker layout; a genuine 2x2
  # slope block with real covariance can. Verify both the per-term Sigma that
  # build_G_glmmTMB() consumes and the assembled G against the lme4 oracle.
  set.seed(2024)
  L <- 30L; n_per <- 25L
  gen <- factor(rep(seq_len(L), each = n_per))
  x   <- stats::rnorm(L * n_per)
  S   <- matrix(c(1, 0.36, 0.36, 0.36), 2) # cov = 0.36 => correlation 0.6
  b   <- matrix(stats::rnorm(L * 2), L, 2) %*% chol(S)
  y   <- 10 + b[gen, 1] + b[gen, 2] * x + stats::rnorm(L * n_per, 0, 0.8)
  sim <- data.frame(y, x, gen)

  ml <- lme4::lmer(y ~ 1 + x + (1 + x | gen), data = sim)
  mt <- glmmTMB::glmmTMB(y ~ 1 + x + (1 + x | gen), data = sim, REML = TRUE)

  Sigma_l <- matrix(as.numeric(lme4::VarCorr(ml)$gen), 2, 2)
  Sigma_t <- matrix(as.numeric(glmmTMB::VarCorr(mt)$cond$gen), 2, 2)
  expect_equal(Sigma_t, Sigma_l, tolerance = 1e-3)      # off-diagonal preserved

  Gt <- build_G_glmmTMB(mt)
  expect_equal(dim(Gt), c(2L * L, 2L * L))
  expect_equal(spectrum(Gt), spectrum(G_from_lme4(ml)), tolerance = 1e-4)
})

# NOTE: glmmTMB's real advantage over lme4 for a Gaussian model -- carrying a
# structured G-side covariance (cs/ar1/toep/us/diag) on the random effects, and
# random slopes -- is deferred to a follow-up. Those cases are pinned as
# informative rejections in "unsupported glmmTMB model structures are rejected
# informatively" above; the acceptance tests that check they yield a valid H2
# live with that follow-up.
