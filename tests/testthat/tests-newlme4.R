
test_that("Confint works",{
  skip()

require(stringr)
lettuce_subset <- lettuce_phenotypes |>
  dplyr::filter(loc == "L2")
N <- nrow(lettuce_subset)

############################## Test for lme4 ###################################
pseudo_var <- sample(c("A", "B"), size = N, replace =TRUE)

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
# To Do: debug, error message when genotype factor has interaction effects.
# To Do: warning message for (f | gen).
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
