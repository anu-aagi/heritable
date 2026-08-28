test_that("Inner checks are triggered", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  lme4_lettuce <- readRDS(test_path("fixtures/lettuce_lme4.rds"))
  lmer_model_random <- readRDS(test_path("fixtures/lmer_model_random.rds"))
  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))
  asreml_model_fixed <- readRDS(test_path("fixtures/asreml_model_fixed.rds"))
  target <- "gen"

  # Won't converge
  lmer_model_failed_converge <- lmer_model_random
  lmer_model_failed_converge@optinfo$conv$opt <- 10
  asreml_model_failed_converge <- asreml_model_random
  asreml_model_failed_converge$converge <- FALSE

  expect_error(H2(
    model = c(asreml_model_random, asreml_model_fixed),
    target = target
  ))
  expect_warning(H2(
    lmer_model_failed_converge,
    target = target,
    method = "Oakey"
  ))
  expect_warning(H2(
    asreml_model_failed_converge,
    target = target,
    method = "Oakey"
  ))
  expect_error(H2(model = asreml_model_random, target = "foo"))
  expect_error(H2(model = asreml_model_fixed, target = "gen"))

  # Target level
  expect_false(check_target_random(asreml_model_fixed, target))
  expect_true(check_target_random(asreml_model_random, target))
  expect_error(H2(model = asreml_model_random, target = c("foo", "rep")))

  # Method level
  expect_error(H2(asreml_model_fixed, target, "Oakey"))
  expect_error(H2(asreml_model_fixed, target, "Cullis"))
})

test_that("Non-Gaussian lme4 models are rejected with an informative error", {
  skip_if_not_installed("lme4")

  set.seed(1)
  ng <- 8L
  nrep <- 4L
  gen_effect <- stats::rnorm(ng, sd = 2)
  df <- data.frame(gen = factor(rep(seq_len(ng), each = nrep)))
  # Give the genotypes a real signal so the Gaussian fit is not singular
  df$y_num <- gen_effect[df$gen] + stats::rnorm(nrow(df))
  df$y_bin <- stats::rbinom(nrow(df), 1, 0.5)
  df$y_cnt <- stats::rpois(nrow(df), 5)

  lmm <- lme4::lmer(y_num ~ 1 + (1 | gen), data = df)
  glmm_bin <- suppressWarnings(
    lme4::glmer(y_bin ~ 1 + (1 | gen), data = df, family = binomial)
  )
  glmm_pois <- suppressWarnings(
    lme4::glmer(y_cnt ~ 1 + (1 | gen), data = df, family = poisson)
  )

  # GLMMs abort with an informative error, not a cryptic dispatch failure
  expect_error(H2(glmm_bin, target = "gen", method = "Cullis"), "identity link")
  expect_error(h2(glmm_bin, target = "gen", method = "Cullis"), "identity link")
  expect_error(
    H2(glmm_pois, target = "gen", method = "Cullis"),
    "identity link"
  )

  # The check itself: Gaussian identity passes, non-Gaussian aborts
  expect_error(check_model_family(lmm), NA)
  expect_error(check_model_family(glmm_bin), "identity link")

  # A Gaussian lme4 model is unaffected and still returns a value
  expect_true(is.finite(H2(lmm, target = "gen", method = "Cullis")))
})

test_that("We can find GRM", {
  asreml_model_grm <- readRDS(test_path("fixtures/asreml_model_grm.rds"))

  expect_error(check_GRM_exists(asreml_model_grm, "gen"))

  data("lettuce_GRM")

  expect_true(
    do.call(
      c,
      check_GRM_exists(
        asreml_model_grm,
        "gen",
        source = list(lettuce_GRM = lettuce_GRM)
      )
    )
  )
  expect_true(
    do.call(c, check_GRM_exists(asreml_model_grm, "gen"))
  )
})
