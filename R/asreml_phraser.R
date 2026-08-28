phrase_G <- function(G_params, ...) {
  terms <- names(G_params[-1])
  specs <- sapply(G_params[-1], function(x) x[["model"]])
  known_specs <- c("id", "idv", "diag", "ar1", "vm")
  unknown_specs <- unique(specs[!specs %in% known_specs])
  if (length(unknown_specs) > 0) {
    cli::cli_abort("Don't know how to phrase {.code {unknown_specs}} yet.")
  }

  G_list <- list()
  v <- G_params[["variance"]][["initial"]]
  if (is.na(v)) {
    v <- 1
  }

  for (trm in terms) {
    spec <- G_params[[trm]][["model"]]
    if (spec == "id") {
      G_list[[trm]] <- phrase_id(G_params[[trm]])
    }
    if (spec == "idv") {
      G_list[[trm]] <- phrase_idv(G_params[[trm]])
    }
    if (spec == "diag") {
      G_list[[trm]] <- phrase_diag(G_params[[trm]])
    }
    if (spec == "ar1") {
      G_list[[trm]] <- phrase_ar1(G_params[[trm]])
    }
    if (spec == "vm") G_list[[trm]] <- phrase_vm(G_params[[trm]], ...)
  }

  Reduce(function(X, Y) Matrix::kronecker(X, Y), G_list) * v
}

phrase_id <- function(G_params, ...) {
  spec <- G_params[["model"]]
  if (spec != "id") {
    cli::cli_abort("Wrong phrase used.")
  }

  grp <- G_params[["levels"]]
  n_g <- length(G_params[["levels"]])

  G <- Matrix::Diagonal(n = n_g)
  dimnames(G) <- list(grp, grp)
  G
}

phrase_idv <- function(G_params, ...) {
  spec <- G_params[["model"]]
  if (spec != "idv") {
    cli::cli_abort("Wrong phrase used.")
  }

  grp <- G_params[["levels"]]
  n_g <- length(G_params[["levels"]])

  G <- Matrix::Diagonal(n = n_g)
  G <- G * G_params[["initial"]]
  dimnames(G) <- list(grp, grp)
  G
}

phrase_diag <- function(G_params, ...) {
  spec <- G_params[["model"]]
  if (spec != "diag") {
    cli::cli_abort("Wrong phrase used.")
  }

  grp <- G_params[["levels"]]
  n_g <- length(G_params[["levels"]])

  G <- Matrix::Diagonal(x = G_params[["initial"]])
  dimnames(G) <- list(grp, grp)
  G
}

phrase_ar1 <- function(G_params, ...) {
  spec <- G_params[["model"]]
  if (spec != "ar1") {
    cli::cli_abort("Wrong phrase used.")
  }

  grp <- G_params[["levels"]]
  n_g <- length(G_params[["levels"]])

  idx <- 0:(n_g - 1)
  rho <- G_params[["initial"]]
  G <- outer(idx, idx, function(i, j) rho^abs(i - j))

  dimnames(G) <- list(grp, grp)
  G
}

phrase_vm <- function(G_params, source = list(), ...) {
  spec <- G_params[["model"]]
  if (spec != "vm") {
    cli::cli_abort("Wrong phrase used.")
  }
  # Identify vm parameter
  vm_par <- G_params[["facnam"]]

  # Parameter estimation

  # TODO: capture when there is SingG argument.
  # Capture the name of the GRM object
  name_GRM <- stringr::str_match(
    vm_par,
    "source\\s*=\\s*([^\\s,\\)]+)"
  )[, 2]

  if (is.na(name_GRM)) {
    name_GRM <- stringr::str_match(
      vm_par,
      "vm\\([^,]+,\\s*([^\\s,\\)]+)"
    )[, 2]
  }

  if (!is.null(source[[name_GRM]])) {
    GRM_source <- source[[name_GRM]]
  } else if (exists(name_GRM, envir = .GlobalEnv, inherits = FALSE)) {
    GRM_source <- get(name_GRM, envir = .GlobalEnv, inherits = FALSE)
  } else {
    cli::cli_abort(
      "GRM used in model was neither provided nor found in the global environment."
    )
  }

  # TODO: inverted GRM source will not be inverted back, debug in the future.
  # Convert to appropriate matrix format
  if (is.data.frame(GRM_source) && ncol(GRM_source) == 3) {
    # TODO: Make be singular
    # Sparse triplet format - convert to matrix
    GRM <- sp2Matrix(GRM_source)
    if (inherits(GRM_source, "ginv") || isTRUE(attr(GRM_source, "INVERSE"))) {
      GRM <- ginv_sym_sparse(GRM)
    }
  } else if (is.matrix(GRM_source) || inherits(GRM_source, "Matrix")) {
    # Regular matrix or Matrix object
    GRM <- GRM_source
    if (inherits(GRM_source, "ginv") || isTRUE(attr(GRM_source, "INVERSE"))) {
      GRM <- ginv_sym_sparse(GRM)
    }
  } else {
    cli::cli_abort("Cannot get the source {.code {name_GRM}} for vm().")
  }

  dimnames(GRM) <- c(list(G_params[["levels"]]), list(G_params[["levels"]]))

  GRM
}
