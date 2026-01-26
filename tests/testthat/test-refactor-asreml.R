test_that("Refactored asreml works",{
  skip_if_not_installed("lme4")
  skip()

  pseudo_var <- rnorm(nrow(lettuce_phenotypes))

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ diag(loc):vm(gen, lettuce_GRM),
    data = lettuce_phenotypes,
    trace = FALSE,
  )

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ gen+ vm(gen, lettuce_GRM),
    data = lettuce_phenotypes,
    trace = FALSE,
  )


  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ gen  + ar1(loc)*ar1(gen),
    data = lettuce_phenotypes,
    trace = FALSE,
  )


  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ gen  + diag(loc):gen,
    data = lettuce_phenotypes,
    trace = FALSE,
  )

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ ar1(gen),
    data = lettuce_phenotypes,
    trace = FALSE,
  )

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ gen  + gen:pseudo_var,
    data = lettuce_phenotypes,
    trace = FALSE,
  )


target <- "gen"

sapply(lettuce_asreml$G.param, function(x){
  facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
  any(do.call(c,facnam))
})

sapply(lettuce_asreml$G.param, function(x){
  facnam <- lapply(x[-1], function(y) names(y[["model"]]))
  paste0(do.call(c,facnam), collapse = ":")
})

sapply(lettuce_asreml$G.param, function(x){
  facnam <- lapply(x[-1], function(y) y[["model"]])
  paste0(do.call(c,facnam), collapse = ":")
})

phrase_G(lettuce_asreml$G.param$gen)



lettuce_asreml$coefficients

lettuce_asreml <- asreml(
  fixed = y ~ 1,
  random =  ~ gen  + rep:pol(pseudo_var, 3):gen,
  data = lettuce_phenotypes,
  trace = FALSE,
)

lettuce_asreml$G.param$`rep:pol(pseudo_var, 3):gen`


# build_new_Z.lmermod lacks check for whether the required grouping variable is provided.



new_data <- c("loc" = "L1")
target <- "gen"

mf <- lettuce_asreml$mf

grp_names <- sapply(lettuce_asreml$G.param, function(x){
  facnam <- lapply(x[-1], function(y) names(y[["model"]]))
  paste0(do.call(c,facnam), collapse = ":")
})

matched_grp <- sapply(lettuce_asreml$G.param, function(x){
  facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
  any(do.call(c,facnam))
}) |> which()

var_used <- sapply(lettuce_asreml$G.param, function(x){
  facnam <- lapply(x[-1], function(y) names(y[["model"]]))
  unique(do.call(c,facnam))
})
var_used <- do.call(c,var_used[matched_grp]) |> unique()

u <- lettuce_asreml$coefficients$random
Gp <- cumsum(attr(u, "terms")[,"n"])
Gp <- c(0, Gp)
idx <- do.call(c,
               lapply(matched_grp, function(x) (Gp[x]+1):Gp[x+1])
)

trms <- names(new_data)

for(trm in trms){
  if (!trm %in% colnames(mf)){
    cli::cli_abort("{.code {trm}} in {.code new_data} was not found in the model.")
  }

  if (is.factor(mf[[trm]]) && !new_data[trm] %in% levels(mf[[trm]])) {
    cli::cli_abort("Unknow level in {.code {trm}} detected: {.code {new_data[trm]}}.")
  }
}

missing_trms <- var_used[!var_used %in% trms]
missing_trms <- missing_trms[missing_trms != target]
if(length(missing_trms) > 0){
  cli::cli_abort("Terms {.code {missing_trms}} interact with {.code {target}} but were not provided.")
}

# Target levels
g      <- mf[[target]]
gnames <- levels(g)
n_g    <- nlevels(g)

new_data <- matrix(rep(new_data, n_g), ncol = n_g) |> t() |> data.frame()
colnames(new_data) <- trms
new_data[[target]] <- gnames

# Add factor level and contrast
for (trm in trms) {
  if (is.factor(mf[[trm]])) {
    new_data[[trm]] <- factor(new_data[[trm]], levels = levels(mf[[trm]]))
    stats::contrasts(new_data[[trm]]) <- stats::contrasts(mf[[trm]])
  }
}

Z_list <- list()

for (i in seq_along(matched_grp)) {
  g_idx <- matched_grp[i]
  term  <- grp_names[g_idx]
  frm <- as.formula(paste0("~0+", term))
  Z_list[[i]] <- Matrix::sparse.model.matrix(frm, new_data)
}

Z <- do.call(cbind,Z_list)
colnames(Z) <- names(u[idx,])

})
