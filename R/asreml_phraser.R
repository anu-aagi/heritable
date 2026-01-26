phrase_G <- function(G_params){
  terms <- names(G_params[-1])
  specs <- sapply(G_params[-1], function(x) x[["model"]])
  known_specs <- c("id","diag", "ar1")
  unknown_specs <- unique(specs[!specs %in% known_specs])
  if(length(unknown_specs) > 0){
    cli::cli_abort("Don't know how to phrase {.code {unknown_specs}} yet.")
  }

  G_list <- list()
  v <- G_params[["variance"]][["initial"]]
  if(is.na(v)) v <- 1

  for(trm in terms){
    spec <- G_params[[trm]][["model"]]
    if (spec == "id") G_list[[trm]] <- phrase_id(G_params[[trm]]) * v
    if (spec == "diag") G_list[[trm]] <- phrase_diag(G_params[[trm]]) * v
    if (spec == "ar1") G_list[[trm]] <- phrase_ar1(G_params[[trm]]) * v
  }

  G_list
}

phrase_id <- function(G_params){
  spec <- G_params[["model"]]
  if(spec != "id"){
    cli::cli_abort("Wrong prase used.")
  }

  grp <- G_params[["levels"]]
  n_g <- length(G_params[["levels"]])

  G <- Matrix::Diagonal(n = n_g)
  dimnames(G) <- list(grp, grp)
  G
}

phrase_diag <- function(G_params){
  spec <- G_params[["model"]]
  if(spec != "diag"){
    cli::cli_abort("Wrong prase used.")
  }

  grp <- G_params[["levels"]]
  n_g <- length(G_params[["levels"]])

  G <- Matrix::Diagonal(x = G_params[["initial"]])
  dimnames(G) <- list(grp, grp)
  G
}

phrase_ar1<- function(G_params){
  spec <- G_params[["model"]]
  if(spec != "ar1"){
    cli::cli_abort("Wrong prase used.")
  }

  grp <- G_params[["levels"]]
  n_g <- length(G_params[["levels"]])

  idx <- 0:(n_g-1)
  rho <- G_params[["initial"]]
  G <- outer(idx, idx, function(i, j) rho^abs(i - j))

  dimnames(G) <- list(grp, grp)
  G
}

