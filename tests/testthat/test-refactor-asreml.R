test_that("Refactored asreml works",{
  skip_if_not_installed("lme4")
  skip()

  pseudo_var <- rnorm(nrow(lettuce_phenotypes))
  pseudo_var2 <- rnorm(nrow(lettuce_phenotypes))

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

asreml.options(design = TRUE)

lettuce_asreml <- asreml(
  fixed = y ~ 1,
  random =  ~ gen  + pol(pseudo_var, 2):gen,
  data = lettuce_phenotypes,
  trace = FALSE,
)

#lettuce_asreml$G.param$`rep:pol(pseudo_var, 3):gen`
#lettuce_asreml$G.param$`spl(pseudo_var, 3):pol(pseudo_var, 3):gen`

# build_new_Z.lmermod lacks check for whether the required grouping variable is provided.
build_new_Z(lettuce_asreml, "gen", data.frame("loc" = "L1"))
build_new_Z(lettuce_asreml, "gen", data.frame("pseudo_var" = 1))


new_data <- data.frame("pseudo_var" = 1)
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

required_var <- sapply(lettuce_asreml$G.param, function(x){
  facnam <- lapply(x[-1], function(y) names(y[["model"]]))
  unique(do.call(c,facnam))
})
required_var <- do.call(c,required_var[matched_grp]) |> unique()

u <- lettuce_asreml$coefficients$random
Gp <- cumsum(attr(u, "terms")[,"n"])
Gp <- c(0, Gp)
idx <- lapply(matched_grp, function(x) (Gp[x]+1):Gp[x+1])

trms <- colnames(new_data)

for(trm in trms){
  if (!trm %in% colnames(mf)){
    cli::cli_abort("{.code {trm}} in {.code new_data} was not found in the model.")
  }

  if (!is.factor(mf[[trm]]) && !inherits(new_data[,trm], "numeric")) {
    cli::cli_abort("{.code {trm}} should be a numeric.")
  }

  if (is.factor(mf[[trm]]) && !new_data[, trm] %in% levels(mf[[trm]])) {
    cli::cli_abort("Unknow level in {.code {trm}} detected: {.code {new_data[trm]}}.")
  }
}

missing_trms <- required_var[!required_var %in% trms]
missing_trms <- missing_trms[missing_trms != target]
if(length(missing_trms) > 0){
  cli::cli_abort("Terms {.code {missing_trms}} interact with {.code {target}} but were not provided.")
}

# Target levels
g      <- mf[[target]]
gnames <- levels(g)
n_g    <- nlevels(g)

if(length(trms) == 1){
  new_data <- data.frame(rep(new_data[,1], n_g))
  colnames(new_data) <- trms
} else {
  new_data <- rep(new_data, n_g)
}
new_data[[target]] <- factor(gnames, levels = gnames)


# Add factor level and contrast
var_type <- c()
for (trm in trms) {
  var_type[trm] <- class(mf[[trm]])
  if (is.factor(mf[[trm]])) {
    new_data[[trm]] <- factor(new_data[[trm]], levels = levels(mf[[trm]]))
    stats::contrasts(new_data[[trm]]) <- stats::contrasts(mf[[trm]])
  }
}
numeric_var <- trms[var_type == "numeric"]

Z_list <- list()

for (i in seq_along(matched_grp)) {
  g_idx <- matched_grp[i]
  term  <- grp_names[g_idx]
  vars <-  sapply(lettuce_asreml$G.param[[i]][-1], function(x) names(x[["model"]]))
  levs <-  sapply(lettuce_asreml$G.param[[i]][-1], function(x) x[["levels"]])
  specs <- sapply(lettuce_asreml$G.param[[i]][-1], function(x) x[["model"]])
  nlevs <- sapply(levs, length)

  numeric_var_idx <- which(vars %in% numeric_var)
  factor_idx <- which(!vars %in% numeric_var)
  stra_idx <- which(!vars %in% numeric_var & vars != target)
  target_idx <- which(vars == target)


  if (length(numeric_var_idx) > 0 && any(specs[numeric_var_idx] != "id")) {
    # Build a group level matrix
    combo_df <- do.call(
      expand.grid,
      c(rev(levs), list(stringsAsFactors = FALSE))
    )
    combo_df <- combo_df[, rev(seq_len(ncol(combo_df))), drop = FALSE]

    # Helper function, build factor matrix
    build_f_mat <- function(x, level){
      i <- seq_along(x)
      j <- match(x, level)
      keep <- !is.na(j)
      mm_grp <- Matrix::sparseMatrix(
        i    = i[keep],
        j    = j[keep],
        x    = rep(1, sum(keep)),
        dims = c(length(x), length(level))
      )
      mm_grp
    }

    # Helper function, build basis
    build_basis <- function(Z, stra){
      sapply(stra, function(s) {
        s <- intersect(colnames(Z), s)
        rowSums(Z[, s, drop = FALSE])
      })
    }

    # Group by numeric terms
    numeric_var_key <- apply(combo_df[, numeric_var_idx,drop = FALSE], 1, function(x) paste0(x, collapse = ":"))
    numeric_var_id <- match(numeric_var_key, unique(numeric_var_key))
    numeric_var_colname <- lapply(
      unique(numeric_var_id),
      function(id) rownames(u)[idx[[i]]][numeric_var_id == id]
    )

    # Choose best approximation row in mf, then get weights per group
    best_apprx <- sweep(
      data.frame(mf)[, vars[numeric_var_idx], drop = FALSE],
      2,
      new_data[1, vars[numeric_var_idx]],
      "-"
    )^2 |> rowSums() |> which.min()

    z <- build_basis(lettuce_asreml$design, numeric_var_colname)[best_apprx, numeric_var_id]

    if(length(stra_idx) > 1){
      new_strata <- new_data[1, vars[stra_idx]]
      stra_key <- apply(combo_df[, stra_idx, drop = FALSE], 1, function(x) paste0(x, collapse = ":"))
      z <- z * (stra_key == new_strata)
    }
    Z_list[[i]] <- t(z * build_f_mat(combo_df[, target_idx], gnames))

  } else {

    frm <- as.formula(paste0("~0+", term))
    Z_list[[i]] <- Matrix::sparse.model.matrix(frm, new_data)

  }

}

Z <- do.call(cbind,Z_list)
colnames(Z) <- names(u[do.call(c,idx),])

})


pheatmap::pheatmap(Z, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
colnames(Z)
