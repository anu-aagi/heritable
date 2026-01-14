
test_that("Confint works",{
  skip()

devtools::document()
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
H2(lettuce_lme4, target = "gen") # Correct

# To Do: debug, error message when genotype factor has interaction effects.
# To Do: warning message for (f | gen).
# To Do: identify interaction terms with genotype factor.
# To Do: modify the standard method in presense of interaction terms.






# lme4::ranef(lettuce_lme4)
# lme4::VarCorr(lettuce_lme4)
# lme4::ngrps(lettuce_lme4)
# lme4::getME(lettuce_lme4, "Z") %>% colnames()
# terms(y ~ x * z * w + (1|x))
# formula <- y ~ x * z * w + (1|x) + (1|z)
# grep("x",deparse1(formula))
# stringr::str_extract_all(deparse1(formula), "(?<=\\()[^|()]+\\|[^|()]+(?=\\))")
})
