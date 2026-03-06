test_that("Refactored asreml works",{
  skip_if_not_installed("lme4")
  skip()

  lettuce_phenotypes_na_impute <- lettuce_phenotypes
  lettuce_phenotypes_na_impute$y[is.na(lettuce_phenotypes_na_impute$y)] <- 3
  lettuce_subset <- subset(lettuce_phenotypes, loc == "L2")
  table(lettuce_phenotypes_na_impute$rep, lettuce_phenotypes_na_impute$gen)

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ loc + gen + gen:loc ,
    data = lettuce_phenotypes_na_impute,
    trace = FALSE,
  )

  lettuce_lme4 <- lmer(
    y~  rep + (1|loc) +  (1|gen) +  (1|gen:loc),
    data = lettuce_phenotypes_na_impute
  )

  H2_Cullis(lettuce_asreml, "gen")
  H2_Cullis(lettuce_lme4, "gen")
  H2_Cullis(lettuce_asreml, "gen", stratification = data.frame(loc = "L1"))
  H2_Cullis(lettuce_lme4, "gen", stratification = data.frame(loc = "L1"))
  H2_Cullis(lettuce_asreml, "gen", marginal = FALSE)
  H2_Cullis(lettuce_lme4, "gen", marginal = FALSE)


  H2_Oakey(lettuce_asreml, "gen")
  H2_Oakey(lettuce_lme4, "gen")
  H2_Oakey(lettuce_asreml, "gen", stratification = data.frame(loc = "L1"))
  H2_Oakey(lettuce_lme4, "gen", stratification = data.frame(loc = "L1"))
  H2_Oakey(lettuce_asreml, "gen", marginal = FALSE)
  H2_Oakey(lettuce_lme4, "gen", marginal = FALSE)

  H2_Standard(lettuce_asreml, "gen")
  H2_Standard(lettuce_lme4, "gen")
  H2_Standard(lettuce_asreml, "gen", stratification = data.frame(loc = "L1"))
  H2_Standard(lettuce_lme4, "gen", stratification = data.frame(loc = "L1"))
  H2_Standard(lettuce_asreml, "gen", marginal = FALSE)
  H2_Standard(lettuce_lme4, "gen", marginal = FALSE)

  H2_Delta(lettuce_asreml, "gen")
  H2_Delta(lettuce_lme4, "gen")
  H2_Delta(lettuce_asreml, "gen", stratification = data.frame(loc = "L1"))
  H2_Delta(lettuce_lme4, "gen", stratification = data.frame(loc = "L1"))
  H2_Delta(lettuce_asreml, "gen", marginal = FALSE)
  H2_Delta(lettuce_lme4, "gen", marginal = FALSE)

  # Test counterpart
  H2_Delta(lettuce_asreml, "gen", marginal = FALSE, type = "BLUE")
  H2_Delta(lettuce_asreml, "gen", marginal = FALSE, type = "BLUP")
  H2_Delta(lettuce_lme4, "gen", marginal = FALSE, type = "BLUE")
  H2_Delta(lettuce_lme4, "gen", marginal = FALSE, type = "BLUP")


  H2(lettuce_asreml, "gen")
  H2(lettuce_asreml, "gen", stratification = data.frame(loc = "L1"))
  H2(lettuce_asreml, "gen", marginal = FALSE)
  H2(lettuce_lme4, "gen")


  # Test simple
  lettuce_subset <- subset(lettuce_phenotypes, loc == "L2")
  table(lettuce_subset$gen)
  table(lettuce_phenotypes_na_impute$gen)

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ gen ,
    data = lettuce_phenotypes_na_impute,
    trace = FALSE,
  )
  H2_Delta(lettuce_asreml, "gen", type = "BLUP")
  H2_Piepho(lettuce_asreml, "gen")
  H2_Standard(lettuce_asreml, "gen")
  H2_Cullis(lettuce_asreml, "gen")
  H2_Oakey(lettuce_asreml, "gen")
  H2(lettuce_asreml, "gen")

  lettuce_lme4 <- lmer(
    y~  rep +  (1|gen),
    data = lettuce_phenotypes_na_impute
  )
  H2_Delta(lettuce_lme4, "gen", type = "BLUP")
  H2_Piepho(lettuce_lme4, "gen")
  H2_Standard(lettuce_lme4, "gen")
  H2_Cullis(lettuce_lme4, "gen")
  H2_Oakey(lettuce_lme4, "gen")
  H2(lettuce_lme4, "gen")

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ loc + ar1(gen) + ar1(gen):loc ,
    data = lettuce_phenotypes_na_impute,
    trace = FALSE,
  )

  lettuce_asreml$G.param
  phrase_G(lettuce_asreml$G.param$`gen:loc`)
  lettuce_asreml$G.param$`gen:loc`$gen$initial
  lettuce_asreml$G.param$`gen:loc`$variance
  phrase_id(lettuce_asreml$G.param$`gen:loc`$loc)

  # Test vm
  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ loc + vm(gen, lettuce_GRM) ,
    data = lettuce_phenotypes_na_impute,
    trace = FALSE,
  )
  lettuce_GRM <- lettuce_GRM
  var_comp(lettuce_asreml, "gen")

#   pseudo_var <- rnorm(nrow(lettuce_phenotypes))
#   pseudo_var2 <- rnorm(nrow(lettuce_phenotypes))
#
#   lettuce_asreml <- asreml(
#     fixed = y ~ rep,
#     random =  ~ diag(loc):vm(gen, lettuce_GRM),
#     data = lettuce_phenotypes,
#     trace = FALSE,
#   )
#
#   lettuce_asreml <- asreml(
#     fixed = y ~ rep,
#     random =  ~ gen+ vm(gen, lettuce_GRM),
#     data = lettuce_phenotypes,
#     trace = FALSE,
#   )
#
#
#   lettuce_asreml <- asreml(
#     fixed = y ~ rep,
#     random =  ~ gen  + ar1(loc)*ar1(gen),
#     data = lettuce_phenotypes,
#     trace = FALSE,
#   )
#
#
#   lettuce_asreml <- asreml(
#     fixed = y ~ rep,
#     random =  ~ gen  + diag(loc):gen,
#     data = lettuce_phenotypes,
#     trace = FALSE,
#   )
#
#   lettuce_asreml <- asreml(
#     fixed = y ~ rep,
#     random =  ~ ar1(gen),
#     data = lettuce_phenotypes,
#     trace = FALSE,
#   )
#
#   lettuce_asreml <- asreml(
#     fixed = y ~ rep,
#     random =  ~ gen  + gen:pseudo_var,
#     data = lettuce_phenotypes,
#     trace = FALSE,
#   )
#
#
# target <- "gen"
#
# sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
#   any(do.call(c,facnam))
# })
#
# sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) names(y[["model"]]))
#   paste0(do.call(c,facnam), collapse = ":")
# })
#
# sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) y[["model"]])
#   paste0(do.call(c,facnam), collapse = ":")
# })
#
# phrase_G(lettuce_asreml$G.param$gen)
#
#
#
# lettuce_asreml$coefficients
#
# asreml.options(design = TRUE)
#
# lettuce_asreml <- asreml(
#   fixed = y ~ 1,
#   random =  ~ gen  + pol(pseudo_var, 2):gen,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
# #lettuce_asreml$G.param$`rep:pol(pseudo_var, 3):gen`
# #lettuce_asreml$G.param$`spl(pseudo_var, 3):pol(pseudo_var, 3):gen`
#
# # build_new_Z.lmermod lacks check for whether the required grouping variable is provided.
# build_new_Z(lettuce_asreml, "gen", data.frame("loc" = "L1"))
# build_new_Z(lettuce_asreml, "gen", data.frame("pseudo_var" = 1))
#
#
# new_data <- data.frame("pseudo_var" = 1)
# target <- "gen"
#
# mf <- lettuce_asreml$mf
#
# grp_names <- sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) names(y[["model"]]))
#   paste0(do.call(c,facnam), collapse = ":")
# })
#
# matched_grp <- sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
#   any(do.call(c,facnam))
# }) |> which()
#
# required_var <- sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) names(y[["model"]]))
#   unique(do.call(c,facnam))
# })
# required_var <- do.call(c,required_var[matched_grp]) |> unique()
#
# u <- lettuce_asreml$coefficients$random
# Gp <- cumsum(attr(u, "terms")[,"n"])
# Gp <- c(0, Gp)
# idx <- lapply(matched_grp, function(x) (Gp[x]+1):Gp[x+1])
#
# trms <- colnames(new_data)
#
# for(trm in trms){
#   if (!trm %in% colnames(mf)){
#     cli::cli_abort("{.code {trm}} in {.code new_data} was not found in the model.")
#   }
#
#   if (!is.factor(mf[[trm]]) && !inherits(new_data[,trm], "numeric")) {
#     cli::cli_abort("{.code {trm}} should be a numeric.")
#   }
#
#   if (is.factor(mf[[trm]]) && !new_data[, trm] %in% levels(mf[[trm]])) {
#     cli::cli_abort("Unknow level in {.code {trm}} detected: {.code {new_data[trm]}}.")
#   }
# }
#
# missing_trms <- required_var[!required_var %in% trms]
# missing_trms <- missing_trms[missing_trms != target]
# if(length(missing_trms) > 0){
#   cli::cli_abort("Terms {.code {missing_trms}} interact with {.code {target}} but were not provided.")
# }
#
# # Target levels
# g      <- mf[[target]]
# if(!is.factor(g)){
#   cli::cli_abort("The target {.code {target}} provided is not a factor.")
# }
# gnames <- levels(g)
# n_g    <- nlevels(g)
#
#
# if(length(trms) == 1){
#   new_data <- data.frame(rep(new_data[,1], n_g))
#   colnames(new_data) <- trms
# } else {
#   new_data <- rep(new_data, n_g)
# }
# new_data[[target]] <- factor(gnames, levels = gnames)
#
#
# # Add factor level and contrast
# var_type <- c()
# for (trm in trms) {
#   var_type[trm] <- class(mf[[trm]])
#   if (is.factor(mf[[trm]])) {
#     new_data[[trm]] <- factor(new_data[[trm]], levels = levels(mf[[trm]]))
#     stats::contrasts(new_data[[trm]]) <- stats::contrasts(mf[[trm]])
#   }
# }
# numeric_var <- trms[var_type == "numeric"]
#
# Z_list <- list()
#
# for (i in seq_along(matched_grp)) {
#   g_idx <- matched_grp[i]
#   term  <- grp_names[g_idx]
#   vars <-  sapply(lettuce_asreml$G.param[[i]][-1], function(x) names(x[["model"]]))
#   levs <-  sapply(lettuce_asreml$G.param[[i]][-1], function(x) x[["levels"]])
#   specs <- sapply(lettuce_asreml$G.param[[i]][-1], function(x) x[["model"]])
#   nlevs <- sapply(levs, length)
#
#   numeric_var_idx <- which(vars %in% numeric_var)
#   factor_idx <- which(!vars %in% numeric_var)
#   stra_idx <- which(!vars %in% numeric_var & vars != target)
#   target_idx <- which(vars == target)
#
#
#   if (length(numeric_var_idx) > 0 && any(specs[numeric_var_idx] != "id")) {
#     # Build a group level matrix
#     combo_df <- do.call(
#       expand.grid,
#       c(rev(levs), list(stringsAsFactors = FALSE))
#     )
#     combo_df <- combo_df[, rev(seq_len(ncol(combo_df))), drop = FALSE]
#
#     # Helper function, build factor matrix
#     build_f_mat <- function(x, level){
#       i <- seq_along(x)
#       j <- match(x, level)
#       keep <- !is.na(j)
#       mm_grp <- Matrix::sparseMatrix(
#         i    = i[keep],
#         j    = j[keep],
#         x    = rep(1, sum(keep)),
#         dims = c(length(x), length(level))
#       )
#       mm_grp
#     }
#
#     # Helper function, build basis
#     build_basis <- function(Z, stra){
#       sapply(stra, function(s) {
#         s <- intersect(colnames(Z), s)
#         rowSums(Z[, s, drop = FALSE])
#       })
#     }
#
#     # Group by numeric terms
#     numeric_var_key <- apply(combo_df[, numeric_var_idx,drop = FALSE], 1, function(x) paste0(x, collapse = ":"))
#     numeric_var_id <- match(numeric_var_key, unique(numeric_var_key))
#     numeric_var_colname <- lapply(
#       unique(numeric_var_id),
#       function(id) rownames(u)[idx[[i]]][numeric_var_id == id]
#     )
#
#     # Choose best approximation row in mf, then get weights per group
#     best_apprx <- sweep(
#       data.frame(mf)[, vars[numeric_var_idx], drop = FALSE],
#       2,
#       new_data[1, vars[numeric_var_idx]],
#       "-"
#     )^2 |> rowSums() |> which.min()
#
#     z <- build_basis(lettuce_asreml$design, numeric_var_colname)[best_apprx, numeric_var_id]
#
#     if(length(stra_idx) > 1){
#       new_strata <- new_data[1, vars[stra_idx]]
#       stra_key <- apply(combo_df[, stra_idx, drop = FALSE], 1, function(x) paste0(x, collapse = ":"))
#       z <- z * (stra_key == new_strata)
#     }
#     Z_list[[i]] <- t(z * build_f_mat(combo_df[, target_idx], gnames))
#
#   } else {
#
#     frm <- as.formula(paste0("~0+", term))
#     Z_list[[i]] <- Matrix::sparse.model.matrix(frm, new_data)
#
#   }
#
# }
#
# Z <- do.call(cbind,Z_list)
# colnames(Z) <- names(u[do.call(c,idx),])
#
# pheatmap::pheatmap(Z, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
# colnames(Z)
#
# ### mapper
# pseudo_var <- rnorm(nrow(lettuce_phenotypes))
#
# lettuce_asreml <- asreml(
#   fixed = y ~ 1,
#   random =  ~ loc:gen + gen,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
# target <- "gen"
# marginal <- TRUE
# map_target_terms(lettuce_asreml, "gen")
# map_target_terms(model2, "gen")
#
#
# mf <- lettuce_asreml$mf
# grp_names <- sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) names(y[["model"]]))
#   paste0(do.call(c,facnam), collapse = ":")
# })
#
# matched_grp <- sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
#   any(do.call(c,facnam))
# }) |> which()
#
# u <- lettuce_asreml$coefficients$random
# Gp <- cumsum(attr(u, "terms")[,"n"])
# Gp <- c(0, Gp)
# idx <- lapply(matched_grp, function(x) (Gp[x]+1):Gp[x+1])
# idx_all <- do.call(c, idx)
# terms <- rownames(u)[idx_all]
#
# g      <- mf[[target]]
# if(!is.factor(g)){
#   cli::cli_abort("The target {.code {target}} provided is not a factor.")
# }
# target_grp <- levels(g)
# n_tg  <- nlevels(g)
# n <- nrow(mf)
#
# w_list <- list()
# m_list <- list()
# intercept_idx <- c()
#
#
# for(i in seq_along(matched_grp)){
#   g_idx <- matched_grp[i]
#   term  <- grp_names[g_idx]
#   vars <-  sapply(lettuce_asreml$G.param[[i]][-1], function(x) names(x[["model"]]))
#   levs <-  sapply(lettuce_asreml$G.param[[i]][-1], function(x) x[["levels"]])
#   nlevs <- sapply(levs, length)
#   target_idx <- which(vars == target)
#   stra_idx <- which(vars != target)
#   grp <- rownames(u)[idx[[i]]]
#   n_g <- length(grp)
#
#   if(grp_names[g_idx] != target){
#
#     combo_df <- do.call(
#       expand.grid,
#       c(rev(levs), list(stringsAsFactors = FALSE))
#     )
#     combo_df <- combo_df[, rev(seq_len(ncol(combo_df))), drop = FALSE]
#     target_key <- combo_df[[target]]
#
#     if (marginal) {
#       stra_key <- apply(combo_df[, stra_idx,drop = FALSE], 1, function(x) paste0(x, collapse = ":"))
#       stra_id <- match(stra_key, unique(stra_key))
#       stra_colname <- lapply(
#         unique(stra_id),
#         function(id) grp[stra_id == id]
#       )
#
#       design <- lettuce_asreml$design
#
#       stra_w <-sapply(stra_colname, function(s) {
#         s <- intersect(colnames(design), s)
#         sum(design[, s, drop = FALSE])/n
#       })
#       w <- stra_w[stra_id]
#     } else {
#       w <- rep(1, prod(nlevs))
#     }
#
#     # Get intercept terms
#     intercept_idx <- c(intercept_idx, rep(0, prod(nlevs)))
#   } else {
#     target_key <- target_grp
#     w <- rep(1, n_g)
#     intercept_idx <- c(intercept_idx, rep(1, n_g))
#   }
#
#   # Helper function, build factor matrix
#   build_f_mat <- function(x, level){
#     i <- seq_along(x)
#     j <- match(x, level)
#     keep <- !is.na(j)
#     mm_grp <- Matrix::sparseMatrix(
#       i    = i[keep],
#       j    = j[keep],
#       x    = rep(1, sum(keep)),
#       dims = c(length(x), length(level))
#     )
#     mm_grp
#   }
#
#   m <- build_f_mat(target_key, target_grp)
#   m <- m * w
#   dimnames(m) <- list(grp, target_grp)
#   m_list[[i]] <- m
#   w_list[[i]] <- w
# }
#
# m <- do.call(rbind, m_list)
# w <- do.call(c, w_list)
#
# intercept <- intercept_idx==1
# list(m = m,
#      w = w,
#      idx = setNames(idx_all, terms),
#      intercept = setNames(intercept, terms))
#
#
#
# ### Build variance component
# lettuce_asreml <- asreml(
#   fixed = y ~ 1,
#   random =  ~ loc + gen,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
# target <- "gen"
#
#
# grp_names <- sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) y[["facnam"]])
#   paste0(do.call(c,facnam), collapse = ":")
# })
#
# matched_grp <- sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
#   any(do.call(c,facnam))
# }) |> which()
#
# grp_names
# lapply(matched_grp, function(x){
#   predict(lettuce_asreml,
#           classify = grp_names[1],
#           vcov = TRUE,
#           trace = FALSE
#   )$vcov
# })
#
#
# lettuce_asreml <- asreml(
#   fixed = y ~ 1,
#   random =  ~ loc + gen,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
#
# predict(lettuce_asreml,
#         classify = "loc:gen",
#         vcov = TRUE,
#         trace = FALSE
# )$vcov
#
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ gen + gen:loc,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
# predict(lettuce_asreml,
#                classify = "gen:loc",
#                only = "gen:loc",
#                vcov = TRUE,
#                trace = FALSE
# )$vcov
#
#  predict(lettuce_asreml,
#          classify = "ar1(gen)",
#          only = "ar1(gen)",
#          vcov = TRUE,
#          trace = FALSE
#  )$vcov
#
#
#
#  M <- as.matrix(lettuce_markers[, -1] + 1)
#  N <- nrow(M)
#  pm <- colSums(M) / (2 * N) # allele freq per marker (diploid X)
#  pm <- pmin(pmax(pm, 1e-6), 1 - 1e-6) # guard against 0 or 1
#  W <- sweep(M, 2, 2 * pm, "-")
#  W <- sweep(W, 2, sqrt(2 * pm * (1 - pm)), "/")
#  G <- tcrossprod(W) / ncol(M)
#  Ginv <- MASS::ginv(G)
#  dimnames(G) <- dimnames(Ginv) <- list(lettuce_markers$gen, lettuce_markers$gen)
#  attr(Ginv, "INVERSE") <- TRUE
#
#
#  lettuce_asreml <- asreml(
#    fixed = y ~ rep,
#    random =  ~ vm(gen, G, singG = "PSD") ,
#    data = lettuce_phenotypes,
#    trace = FALSE,
#  )
#  map_target_terms(lettuce_asreml, "gen")
#  lettuce_asreml$coefficients$random
#  lettuce_asreml$G.param[[1]]
#
#  grp_names <- sapply(lettuce_asreml$G.param, function(x){
#    facnam <- lapply(x[-1], function(y) y[["facnam"]])
#    paste0(do.call(c,facnam), collapse = ":")
#  })
#  G_list <- sapply(seq_along(grp_names), function(x){
#    phrase_G(lettuce_asreml$G.param[[x]])
#  })
#  G <- Matrix::bdiag(G_list)
#  design <- lettuce_asreml$design
#  u <- lettuce_asreml$coefficients$random
#  u_names <- rownames(u)
#
#
#  lettuce_asreml <- asreml(
#    fixed = y ~ rep,
#    random =  ~ gen ,
#    data = lettuce_phenotypes,
#    trace = FALSE,
#  )
#
#
# lapply(lettuce_asreml$G.param, function(x){
#    facnam <- lapply(x[-1], function(y) paste0(y[["facnam"]], "_", y[["levels"]]))
#    facnam <- do.call(
#      expand.grid,
#      c(rev(facnam), list(stringsAsFactors = FALSE))
#    )
#    facnam <- facnam[, rev(seq_len(ncol(facnam))), drop = FALSE]
#    facnam
# })
#
#
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc:gen + gen ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ vm(gen, G, singG = "PSD") ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ gen + pol(pseudo_var,3):gen ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
#
#
# map_target_terms(lettuce_asreml, "gen")
# build_new_Z(lettuce_asreml, "gen", data.frame("loc" = "L1"))
#
# lettuce_asreml$coefficients$random
#
#
# df <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = sample(c("A","B"), 100, TRUE), x4 = sample(c("A","B"), 100, TRUE), x5 = sample(c("A","B"), 100, TRUE))
# model.matrix(~0+ x4:x3:x5, df)
#
#
#
#
#
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc + gen ,
#   data = lettuce_phenotypes,
#   trace = FALSE,
# )
# target <- "gen"
#
#
# mapper <- map_target_terms(lettuce_asreml, target)
# g <- mapper$idx
#
# grp_names <- sapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) y[["facnam"]])
#   paste0(do.call(c,facnam), collapse = ":")
# })
# ran_terms <- lapply(lettuce_asreml$G.param, function(x){
#   facnam <- lapply(x[-1], function(y) {
#     y <- paste0(y[["facnam"]], "_", y[["levels"]])
#     factor(y, levels = y)
#   }
#   )
#   facnam <- do.call(
#     expand.grid,
#     c(rev(facnam), list(stringsAsFactors = FALSE))
#   )
#   facnam <- facnam[, rev(seq_len(ncol(facnam))), drop = FALSE]
#   facnam <- apply(facnam, 1, function(x) paste0(x, collapse = ":"))
# })
# ran_terms <- do.call(c, ran_terms)
#
#
# G_list <- sapply(seq_along(grp_names), function(x){
#   phrase_G(lettuce_asreml$G.param[[x]])
# })
# G <- Matrix::bdiag(G_list)
# G <- G * lettuce_asreml$sigma2
# dimnames(G) <- list(ran_terms, ran_terms)
# Q <- ncol(G)
#
# design <- lettuce_asreml$design
# N <- nrow(design)
# Z <- Matrix::Matrix(0, N, Q)
# colnames(Z) <- ran_terms
# common_col <- intersect(colnames(design), ran_terms)
# Z[,common_col] <- design[, common_col]
#
# X <- design[,!colnames(design) %in% ran_terms]
#
#
# # Matches
# lettuce_phenotypes_na_impute <- lettuce_phenotypes
# lettuce_phenotypes_na_impute$y[is.na(lettuce_phenotypes_na_impute$y)] <- 3
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc + gen:loc ,
#   data = lettuce_phenotypes_na_impute,
#   trace = FALSE,
# )
# var_comp(lettuce_asreml, "gen")
#
#
# lettuce_lme4 <- lmer(
#   y~ rep + (1|loc) + (1|gen:loc),
#   data = lettuce_phenotypes_na_impute
# )
# var_comp(lettuce_lme4, "gen")
#
#
#
# lettuce_asreml <- asreml(
#   fixed = y ~ rep,
#   random =  ~ loc + ar1(gen) + gen:loc ,
#   data = lettuce_phenotypes_na_impute,
#   trace = FALSE,
# )
# check_model_specification.asreml(lettuce_asreml, "gen", "broad_sense")
#
# lettuce_lme4 <- lmer(
#   y~ rep + (1|loc) + (rep|gen) + (1|gen:loc),
#   data = lettuce_phenotypes_na_impute
# )
# check_model_specification(lettuce_lme4, "gen", "broad_sense")

})


test_that("VM extraction works", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()
  skip()

  lettuce_phenotypes_na_impute <- lettuce_phenotypes
  lettuce_phenotypes_na_impute$y[is.na(lettuce_phenotypes_na_impute$y)] <- 3

  lettuce_asreml <- asreml(
    fixed = y ~ rep,
    random =  ~ loc + vm(gen, lettuce_GRM) ,
    data = lettuce_phenotypes_na_impute,
    trace = FALSE,
  )

  var_comps <- var_comp(lettuce_asreml, "gen", calc_C22 = TRUE)

  expect_named(var_comps)
  expect_length(var_comps, 5)
  expect_equal(dim(var_comps$G_g), dim(lettuce_GRM))
})

