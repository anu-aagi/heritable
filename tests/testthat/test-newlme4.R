
test_that("Confint works",{
  skip()

# require(stringr)
# require(asreml)
# devtools::document()
lettuce_subset <- lettuce_phenotypes |>
  dplyr::filter(loc == "L2")
N <- nrow(lettuce_subset)

############################## Test for lme4 ###################################
pseudo_var <- sample(c("A", "B"), size = N, replace =TRUE)
pseudo_var2 <- sample(c("A", "B"), size = N, replace =TRUE)

# Fails!
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen * pseudo_var), data = lettuce_subset)
pull_terms(lettuce_lme4)
pull_terms_without_specials(lettuce_lme4)

# Fails!
lettuce_lme4 <- lme4::lmer(y ~ rep + (pseudo_var | gen ), data = lettuce_subset)
pull_terms(lettuce_lme4)
pull_terms_without_specials(lettuce_lme4)

# Fails!
lettuce_lme4 <- lme4::lmer(y ~ rep + (0 + pseudo_var | gen ), data = lettuce_subset)
pull_terms(lettuce_lme4)
pull_terms_without_specials(lettuce_lme4)

# Fails!
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen ) + (1|pseudo_var), data = lettuce_subset)
pull_terms(lettuce_lme4)
pull_terms_without_specials(lettuce_lme4)

# Fails!
lettuce_lme4 <- lme4::lmer(y ~ (rep + pseudo_var) + (1 | gen ), data = lettuce_subset)
pull_terms(lettuce_lme4)
pull_terms_without_specials(lettuce_lme4)


# pull_terms.lmerMod was corrected so it works now.

# Check cov construction 1
lettuce_lme4 <-  lme4::lmer(y ~ rep + (pseudo_var | gen ), data = lettuce_subset)
H2(lettuce_lme4, target = "gen")
PEV_from_lme4(lettuce_lme4)


# Check cov construction 2: Fails ! Multiple terms
lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen/pseudo_var ), data = lettuce_subset)
H2(lettuce_lme4, target = "gen")
PEV_from_lme4(lettuce_lme4)


# Check cov construction 3: Fails ! Multiple terms, Singular G
lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen*pseudo_var ), data = lettuce_subset)
H2(lettuce_lme4, target = "gen", method = "Cullis")
PEV_from_lme4(lettuce_lme4) # Fails! Modified PEV_from_lme4 to make it robust to singular G
geno_components_from_lme4(lettuce_lme4, "gen", PEV_from_lme4(lettuce_lme4))$C22_g


# Check whether the modifed one matches with the result before
lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen ), data = lettuce_subset)

old_PEV_from_lme4 <- function(model) {
  vc <- lme4::VarCorr(model)
  ngrps <- lme4::ngrps(model)
  # Note the index and kronecker order needs to be followed careful downstream
  Glist <- lapply(names(vc), function(agrp) {
    Matrix::kronecker(vc[[agrp]], diag(ngrps[[agrp]]))
  })
  G <- do.call(Matrix::bdiag, Glist)

  n <- nrow(model@frame)
  R <- diag(n) * stats::sigma(model)^2

  X <- as.matrix(lme4::getME(model, "X"))
  Z <- as.matrix(lme4::getME(model, "Z"))

  Ginv <-   Cinv <- tryCatch(
    solve(G),
    error = function(e) {
      warning("Matrix inversion fails due to singularity, use pseudo-inverse.")
      MASS::ginv(G)
    }
  )

  C11 <- t(X) %*% solve(R) %*% X
  C12 <- t(X) %*% solve(R) %*% Z
  C21 <- t(Z) %*% solve(R) %*% X
  C22 <- t(Z) %*% solve(R) %*% Z + Ginv

  C <- rbind(
    cbind(C11, C12),
    cbind(C21, C22)
  )
  Cinv <- tryCatch(
    solve(C),
    error = function(e) {
      warning("Matrix inversion fails due to singularity, use pseudo-inverse.")
      MASS::ginv(C)
    }
  )
  Cinv
} # Original
old_PEV_from_lme4(lettuce_lme4)
geno_components_from_lme4(lettuce_lme4, "gen") # Correct
H2(lettuce_lme4, target = "gen")


# 01-13
# To Do: debug, error message when genotype factor has interaction effects. (Done)
# To Do: warning message for (f | gen) (Done).
# To Do: identify interaction terms with genotype factor.
# To Do: modify the standard method in presense of interaction terms.


# Check multiple terms
lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen * pseudo_var ), data = lettuce_subset)
model_terms <- pull_terms_without_specials(lettuce_lme4) # All good
form <- as.formula(paste("~", paste(
  c(
    model_terms$fixed,
    model_terms$random
  ),
  collapse = " + "
)))
form <- update(form, paste0("~ . - ", "gen"))
fcts <- rownames(attr(terms(form), "factors")) # The fcts will capture the interaction, so not good
# for modelling MET.

lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen) + (1 | gen), data = lettuce_subset)
model_terms <- pull_terms_without_specials(lettuce_lme4) # All good
form <- as.formula(paste("~", paste(
  c(
    model_terms$fixed,
    model_terms$random
  ),
  collapse = " + "
)))
form <- update(form, paste0("~ . - ", "gen"))
fcts <- rownames(attr(terms(form), "factors")) # No captrues. Not good.


lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | rep)  + (pseudo_var || rep), data = lettuce_subset)
H2(lettuce_lme4 , "gen") # Error in PEV_from_lme4 due to discrepancy between:
vc <- lme4::VarCorr(model)
ngrps <- lme4::ngrps(model)


# Solution, further modify the random effect pull function, align the output to asremel
model_formula <- y ~ x * y + (1||a*b) + (1|a*b)
term_labels <- attr(terms(model_formula), "term.labels")
ran_trms_formula <-
  stringr::str_extract_all(deparse1(model_formula), "(?<=\\()[^|()]+\\|+[^|()]+(?=\\))")[[1]]
fixed_trms <- setdiff(term_labels, ran_trms_formula)
ran_trms <- lapply(ran_trms_formula, function(frm){
  frm <- stringr::str_split(frm, "\\|")[[1]] |> tail(n=1)
  as.formula(
    paste0("~", frm)
  ) |> terms() |> attr("term.labels")
})
ran_trms <- do.call(c, ran_trms)

# Solution, further modify PEV_from_lme4, find different way of calculating ngrps
# Var(y)=ZVar(b)Z⊤+Var(ε)=σ2(ZΛΛ⊤Z⊤+I)
Lambda <- lme4::getME(lettuce_lme4, "Lambda")
G <- tcrossprod(Lambda) * sigma(lettuce_lme4)^2
lme4::VarCorr(lettuce_lme4)[[3]] # Correct
# Hence I modified the PEV_from_lme4 function.
lettuce_lme4@flist


# Bug in fit_counterpart_model.lmerMod when the model multiple terms
# with target as grouping factors.
lettuce_lme4 <- lme4::lmer(y ~ rep + ( pseudo_var || gen) + (1|pseudo_var), data = lettuce_subset)
ran_frms <- reformulas::findbars(formula(lettuce_lme4))
use <- reformulas::findbars(formula(lettuce_lme4)) %>%
  sapply(., function(frm){
    frm <- stringr::str_split(deparse1(frm), " \\| ")[[1]] |> tail(n=1)
  }) == "gen"
ran_frms <- ran_frms[use] %>%
  sapply(., function(frm){
    paste0("(", deparse1(frm), ")")
  })
fit_counterpart_model.lmerMod(lettuce_lme4, "gen")

# 01-19
# Check simple model fit for borad sense heritability.
lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen * pseudo_var ), data = lettuce_subset)
pull_terms(lettuce_lme4)
check_model_specification(lettuce_lme4, "gen", "broad_sense")

lettuce_lme4 <-  lme4::lmer(y ~ rep + ( pseudo_var | gen  ), data = lettuce_subset)
check_model_specification(lettuce_lme4, "gen", "broad_sense")

lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen  ) + (0+ pseudo_var | gen  ), data = lettuce_subset)
check_model_specification(lettuce_lme4, "gen", "broad_sense")
geno_components_from_lme4(lettuce_lme4, "gen")
H2(lettuce_lme4, target = "gen")
pull_terms_without_specials(lettuce_lme4)

# 01-19 Renaming, standardise
H2_Standard(lettuce_lme4, target = "gen")

H2_Cullis.old <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  g <- var_comp(model, target)
  s2_g <- mean(Matrix::diag(g$G_g))

  P_mu <- Matrix::Diagonal(n = g$n_g, x = g$n_g) - 1
  vdBLUP_sum <- sum(Matrix::diag(P_mu %*% g$C22_g))
  vdBLUP_avg <- vdBLUP_sum * (2 / (g$n_g * (g$n_g - 1)))

  return(H2_Cullis_parameters(vdBLUP_avg, s2_g))
}
H2_Cullis.new <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  g <- var_comp(model, target)
  s2_g <- mean(Matrix::diag(g$G_g))
  n <- g$n_g
  C22_g <- g$C22_g

  # This is equivalent to delta <- var_diff(C22_g); delta_avg = mean(delta[lower.tri(delta)])
  delta_avg <-  (2 / (n * (n - 1))) * (n * sum(diag(C22_g)) - sum(C22_g))

  return(H2_Cullis_parameters(delta_avg, s2_g))
}
H2_Cullis.old(lettuce_lme4, target = "gen")
H2_Cullis.new(lettuce_lme4, target = "gen")
H2_Cullis(lettuce_lme4, target = "gen")
initial_checks(lettuce_lme4, target = "gen", NULL)

H2_Oakey(lettuce_lme4, target = "gen")

H2_Piepho(lettuce_lme4, target = "gen")

H2_Delta(lettuce_lme4, target = "gen")

H2_Delta_BLUE_pairwise.old <- function(
    model,
    target = NULL,
    options = NULL
) {
  initial_checks(model, target, options)

  conterpart <- fit_counterpart_model(model, target)

  # Extract vc_g and vc_e
  g <- var_comp(model, target, calc_C22 = FALSE)
  s2_g <- mean(diag(g$G_g))

  # Calculate mean variance of a difference between genotypes
  frm <-  as.formula(paste("pairwise ~", target))
  EMM_fit <- emmeans::emmeans(conterpart, specs = frm)$contrasts
  EMM_fit <- data.frame(EMM_fit)  # Get variance
  EMM_fit$var <- EMM_fit$SE^2

  # Take pairwise differences and turn into variance-covariance matrix
  gnames <- g$gnames
  n_g <- g$n_g

  # Start with empty variance matrix for differences
  delta <- matrix(0, nrow = n_g, ncol = n_g)
  dimnames(delta) <- list(gnames, gnames)

  # Fill in the pairwise variances from deltas
  for (i in 1:nrow(EMM_fit)) {
    # Extract genotype names from contrast column
    pair <- strsplit(as.character(EMM_fit$contrast[i]), " - ")[[1]]
    g1 <- pair[1]
    g2 <- pair[2]

    # Variance of difference: Var(g1 - g2) = Var(g1) + Var(g2) - 2*Cov(g1, g2)
    # Get covariance between g1 and g2 (0 by default, but can be specified)
    delta[g1, g2] <- EMM_fit$var[i]
    delta[g2, g1] <- EMM_fit$var[i] # symmetric
  }

  # H2 Delta BLUE
  H2_Delta_BLUE <- H2_Delta_BLUE_parameters(s2_g, delta)

  return(H2_Delta_BLUE)
}

H2_Delta_BLUE_pairwise.old(lettuce_lme4, target = "gen")

lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen ), data = lettuce_subset)
H2(lettuce_lme4, "gen")

# Get GxE terms
lettuce_lme4 <-  lme4::lmer(y ~ rep + ( 1 | gen * pseudo_var ), data = lettuce_subset)
ran_trms <- pull_terms_without_specials(lettuce_lme4)$random
target <- "gen"
pattern <- paste0("(^|:)", target, "($|:)")
grepl(pattern, ran_trms)

sigma2 <- sigma(lettuce_lme4)^2
Z <- lme4::getME(lettuce_lme4, "Z")
Lambda <- lme4::getME(lettuce_lme4, "Lambda")
G <- tcrossprod(Lambda) * sigma2
dimnames(G) <- list(colnames(Z), colnames(Z))


# 01-20 Try get all terms associated with gen, and then extract each of the component
lettuce_lme4 <-  lme4::lmer(y ~  (pseudo_var | gen ) + (1|pseudo_var), data = lettuce_subset)
lettuce_lme4 <-  lme4::lmer(y ~  (pseudo_var | pseudo_var * gen ) + (1|pseudo_var), data = lettuce_subset)
lettuce_lme4 <-  lme4::lmer(y ~  rep* pseudo_var2 + (1 | gen) + (1|pseudo_var), data = lettuce_subset)
lettuce_lme4 <-  lme4::lmer(y ~  (pseudo_var | gen ) + (1| gen), data = lettuce_subset)
lettuce_lme4 <-  lme4::lmer(y ~  rep + ( 0 + pseudo_var | gen ) + (1| gen), data = lettuce_subset)
lettuce_lme4 <-  lme4::lmer(y ~  rep + (1| gen), data = lettuce_subset)
lettuce_lme4 <-  lme4::lmer(y ~  rep + (1| gen) + (1|pseudo_var), data = lettuce_subset)

# Get variance componenets
maps <- map_target_terms(lettuce_lme4, "gen", reconstruct = FALSE)
H2(lettuce_lme4, "gen")
colnames(lme4::getME(lettuce_lme4, "Z")[,maps$idx]) ==
  names(maps$w)#Equal


maps <- map_target_terms(lmer_model_random, "gen", reconstruct = FALSE)
H2(lmer_model_random, "gen")
var_comp(lmer_model_random, "gen")


# lme4::ranef(lettuce_lme4)
# lme4::VarCorr(lettuce_lme4)
# lme4::ngrps(lettuce_lme4)
# lme4::getME(lettuce_lme4, "Z") %>% colnames()
# terms(y ~ x * z * w + (1|x))
# formula <- y ~ x * z * w + (1|x) + (1|z)
# grep("x",deparse1(formula))
# stringr::str_extract_all(deparse1(formula), "(?<=\\()[^|()]+\\|[^|()]+(?=\\))")
})
# terms(formula(y ~ (1|x*z)))

# require(asreml)
# lettuce_asreml <- asreml(
#   fixed = y ~ rep * as.factor(pseudo_var),
#   random =  ~ gen + ar1(gen) ,
#   data = lettuce_subset,
#   trace = FALSE,
# )
# fixed_trms <- attr(lettuce_asreml$formulae$fixed, "term.labels")
# ran_trms <- attr(lettuce_asreml$formulae$random, "term.labels")
# model_terms <- pull_terms_without_specials(lettuce_asreml) # All good
# form <- as.formula(paste("~", paste(
#   c(
#     model_terms$fixed,
#     model_terms$random
#   ),
#   collapse = " + "
# )))
# form <- update(form, paste0("~ . - ", "gen"))
# fcts <- rownames(attr(terms(form), "factors")) # The fcts will capture the interaction, so not good


# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ gen + rnorm(nrow(lettuce_subset))*factor(pseudo_var) ,
#   data = lettuce_subset,
#   trace = FALSE,
# )
