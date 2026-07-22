#' Bayesian linear regression from sufficient statistics
#'
#' Fits the same prior models as [blm()] using cross-products instead of the
#' original response and predictor matrix.
#'
#' @param n Number of observations used to form the sufficient statistics.
#' @param XtX A finite, symmetric predictor cross-product matrix.
#' @param Xty A finite predictor-response cross-product vector.
#' @param ETA A prior specification or named list of prior blocks. Each block
#'   must contain `model` and may contain `indices`, an integer or character
#'   vector selecting columns of `XtX`. A single block may omit `indices` and
#'   then uses every predictor. Multiple blocks must partition all predictors.
#'   Prior parameters and `standardize` are the same as in [blm()].
#' @param yty Optional finite, nonnegative response sum of squares. It is
#'   required when `residual_var = NULL`.
#' @param X_means,y_mean Optional predictor means and response mean. They must
#'   be supplied together. When omitted, the model is fitted without an
#'   intercept.
#' @inheritParams blm
#'
#' @return A fitted object with the same block-specific posterior summaries as
#'   [blm()]. Intercept components are present only when both `X_means` and
#'   `y_mean` are supplied.
#'
#' @details If means are supplied, `XtX`, `Xty`, and `yty` are interpreted as
#'   uncentered cross-products and are centered internally. Without means,
#'   they are used as supplied. If any block requests standardization without
#'   means, a warning reminds the user that the supplied cross-products should
#'   already represent centered or standardized variables.
#'
#' @export
#'
#' @examples
#' X <- cbind(x1 = 1:20, x2 = rep(c(0, 1), 10))
#' y <- 1 + 2 * X[, "x1"] - X[, "x2"]
#' fit <- blm_ss(
#'   n = nrow(X),
#'   XtX = crossprod(X),
#'   Xty = drop(crossprod(X, y)),
#'   ETA = list(model = "Normal"),
#'   yty = sum(y^2),
#'   X_means = colMeans(X),
#'   y_mean = mean(y),
#'   residual_var = 1,
#'   iterations = 100,
#'   burnin = 50,
#'   seed = 123
#' )
blm_ss <- function(n, XtX, Xty, ETA, yty = NULL, X_means = NULL,
                   y_mean = NULL, residual_var = NULL,
                   residual_shape = NULL, residual_scale = NULL,
                   iterations = 4000L, burnin = 1000L, thin = 1L,
                   seed = NULL, version = c("Rcpp", "R"), verbose = FALSE,
                   nchains = 1L, store_samples = TRUE,
                   store_coefficient_cov = TRUE) {
  version <- match.arg(version)
  nchains <- .validate_nchains(nchains)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(store_samples) || length(store_samples) != 1L ||
      is.na(store_samples)) {
    stop("`store_samples` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(store_coefficient_cov) ||
      length(store_coefficient_cov) != 1L || is.na(store_coefficient_cov)) {
    stop("`store_coefficient_cov` must be TRUE or FALSE.", call. = FALSE)
  }

  statistics <- .validate_sufficient_statistics(
    n, XtX, Xty, yty, X_means, y_mean
  )
  n <- statistics$n
  XtX <- statistics$XtX
  Xty <- statistics$Xty
  yty <- statistics$yty
  X_means <- statistics$X_means
  y_mean <- statistics$y_mean
  predictor_names <- statistics$predictor_names
  fit_intercept <- !is.null(X_means)

  if (is.null(yty) && is.null(residual_var)) {
    stop(
      "`residual_var` must be fixed when `yty` is not supplied.",
      call. = FALSE
    )
  }
  .validate_residual_specification(
    residual_var, residual_shape, residual_scale
  )

  normalized <- .normalize_ss_eta(ETA, predictor_names, residual_var)
  blocks <- normalized$blocks
  source_indices <- normalized$source_indices
  if (!fit_intercept) {
    warning(
      "`X_means` and `y_mean` were not supplied; fitting without an intercept.",
      call. = FALSE
    )
    if (any(vapply(blocks, `[[`, logical(1), "standardize"))) {
      warning(
        paste0(
          "`standardize = TRUE` without `X_means`; the supplied ",
          "cross-products should be centered or standardized."
        ),
        call. = FALSE
      )
    }
  }

  centered_XtX <- if (fit_intercept) {
    XtX - n * tcrossprod(X_means)
  } else {
    XtX
  }
  centered_Xty <- if (fit_intercept) {
    Xty - n * X_means * y_mean
  } else {
    Xty
  }
  centered_yty <- if (is.null(yty)) {
    NULL
  } else if (fit_intercept) {
    yty - n * y_mean^2
  } else {
    yty
  }

  variance_reference <- if (fit_intercept) {
    pmax(1, abs(diag(XtX)), abs(n * X_means^2))
  } else {
    pmax(1, abs(diag(XtX)))
  }
  variance_tolerance <- 100 * .Machine$double.eps * variance_reference
  constant_predictors <- diag(centered_XtX) <= variance_tolerance
  if (any(constant_predictors)) {
    stop(
      sprintf(
        "The sufficient statistics contain constant predictor(s): %s.",
        paste(predictor_names[constant_predictors], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  predictor_scales <- lapply(seq_along(blocks), function(block_index) {
    indices <- source_indices[[block_index]]
    if (!blocks[[block_index]]$standardize) return(rep(1, length(indices)))
    variances <- diag(centered_XtX)[indices] / (n - 1)
    if (any(!is.finite(variances)) || any(variances <= 0)) {
      stop(
        sprintf(
          "ETA block `%s` has a predictor with nonpositive variance.",
          names(blocks)[block_index]
        ),
        call. = FALSE
      )
    }
    sqrt(variances)
  })
  for (block_index in seq_along(blocks)) {
    blocks[[block_index]]$predictor_scale <- predictor_scales[[block_index]]
  }

  source_order <- unlist(source_indices, use.names = FALSE)
  scale_order <- unlist(predictor_scales, use.names = FALSE)
  working_XtX <- centered_XtX[source_order, source_order, drop = FALSE] /
    outer(scale_order, scale_order)
  working_Xty <- centered_Xty[source_order] / scale_order
  pseudo <- .crossproducts_to_pseudo(working_XtX, working_Xty, centered_yty)

  block_sizes <- lengths(source_indices)
  block_ends <- cumsum(block_sizes)
  block_starts <- block_ends - block_sizes + 1L
  block_indices <- Map(seq.int, block_starts, block_ends)
  block_model <- vapply(blocks, `[[`, integer(1), "model_code")
  block_id <- rep.int(seq_along(blocks), block_sizes)
  internal_names <- unlist(lapply(seq_along(blocks), function(block_index) {
    paste0(names(blocks)[block_index], "::", blocks[[block_index]]$predictor_names)
  }))
  colnames(pseudo$x) <- internal_names

  normal_shape <- vapply(blocks, `[[`, numeric(1), "normal_shape")
  normal_scale <- vapply(blocks, `[[`, numeric(1), "normal_scale")
  pi_alpha <- vapply(blocks, `[[`, numeric(1), "pi_alpha")
  pi_beta <- vapply(blocks, `[[`, numeric(1), "pi_beta")
  slab_shape <- vapply(blocks, `[[`, numeric(1), "slab_shape")
  slab_scale <- vapply(blocks, `[[`, numeric(1), "slab_scale")
  global_scale <- vapply(blocks, `[[`, numeric(1), "global_scale")
  local_a <- vapply(blocks, function(block) block$local_shape[1L], numeric(1))
  local_b <- vapply(blocks, function(block) block$local_shape[2L], numeric(1))
  sampler_arguments <- list(
    y = pseudo$y,
    x = pseudo$x,
    residual_shape = if (is.null(residual_shape)) 1 else residual_shape,
    residual_scale = if (is.null(residual_scale)) 1 else residual_scale,
    residual_var = residual_var,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    block_id = block_id,
    block_model = block_model,
    normal_shape = normal_shape,
    normal_scale = normal_scale,
    pi_alpha = pi_alpha,
    pi_beta = pi_beta,
    slab_shape = slab_shape,
    slab_scale = slab_scale,
    global_scale = global_scale,
    local_a = local_a,
    local_b = local_b,
    store_samples = store_samples,
    store_coefficient_cov = store_coefficient_cov,
    effective_n = n,
    fit_intercept = fit_intercept,
    intercept_x_mean = if (fit_intercept) {
      X_means[source_order] / scale_order
    } else {
      rep(0, length(source_order))
    },
    intercept_y_mean = if (fit_intercept) y_mean else 0
  )
  run_chains <- function(progressor = NULL) {
    .run_blm_chains(
      sampler_arguments, version, nchains, seed, block_model, progressor
    )
  }
  samples <- if (verbose) {
    progressr::with_progress({
      progress <- progressr::progressor(steps = nchains * iterations)
      run_chains(progress)
    }, enable = TRUE)
  } else {
    run_chains()
  }

  .assemble_blm_result(
    blocks, block_indices, samples, nchains, store_samples,
    store_coefficient_cov, fit_intercept
  )
}

.validate_residual_specification <- function(residual_var, residual_shape,
                                             residual_scale) {
  if (!is.null(residual_var)) {
    if (!is.null(residual_shape) || !is.null(residual_scale)) {
      stop(
        "Supply either `residual_var` or the inverse-gamma prior, not both.",
        call. = FALSE
      )
    }
    .validate_variance(residual_var, "residual_var")
  } else {
    if (is.null(residual_shape) || is.null(residual_scale)) {
      stop(
        paste0(
          "`residual_shape` and `residual_scale` are required when ",
          "`residual_var` is NULL."
        ),
        call. = FALSE
      )
    }
    .validate_variance(residual_shape, "residual_shape")
    .validate_variance(residual_scale, "residual_scale")
  }
}

.validate_sufficient_statistics <- function(n, XtX, Xty, yty, X_means,
                                            y_mean) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n != floor(n) || n < 2) {
    stop("`n` must be an integer of at least two.", call. = FALSE)
  }
  if (!is.matrix(XtX) || !is.numeric(XtX) || nrow(XtX) < 1L ||
      nrow(XtX) != ncol(XtX) || anyNA(XtX) || any(!is.finite(XtX))) {
    stop("`XtX` must be a finite numeric square matrix.", call. = FALSE)
  }
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(XtX)))
  if (max(abs(XtX - t(XtX))) > tolerance) {
    stop("`XtX` must be symmetric.", call. = FALSE)
  }
  XtX <- (XtX + t(XtX)) / 2
  eigenvalues <- eigen(XtX, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < -tolerance) {
    stop("`XtX` must be positive semidefinite.", call. = FALSE)
  }
  p <- ncol(XtX)
  if (is.matrix(Xty) && is.numeric(Xty) &&
      identical(dim(Xty), c(p, 1L))) {
    Xty <- drop(Xty)
  }
  if (!is.numeric(Xty) || !is.atomic(Xty) || is.object(Xty) ||
      !is.null(dim(Xty)) || length(Xty) != p || anyNA(Xty) ||
      any(!is.finite(Xty))) {
    stop("`Xty` must be a finite numeric vector matching `XtX`.",
         call. = FALSE)
  }
  if (!is.null(yty) && (!is.numeric(yty) || length(yty) != 1L ||
      is.na(yty) || !is.finite(yty) || yty < 0)) {
    stop("`yty` must be a finite, nonnegative numeric scalar.", call. = FALSE)
  }
  if (xor(is.null(X_means), is.null(y_mean))) {
    stop("`X_means` and `y_mean` must be supplied together.", call. = FALSE)
  }
  if (!is.null(X_means)) {
    if (!is.numeric(X_means) || !is.atomic(X_means) || is.object(X_means) ||
        !is.null(dim(X_means)) || length(X_means) != p || anyNA(X_means) ||
        any(!is.finite(X_means))) {
      stop("`X_means` must be a finite numeric vector matching `XtX`.",
           call. = FALSE)
    }
    if (!is.numeric(y_mean) || length(y_mean) != 1L || is.na(y_mean) ||
        !is.finite(y_mean)) {
      stop("`y_mean` must be a finite numeric scalar.", call. = FALSE)
    }
  }
  predictor_names <- colnames(XtX)
  if (!is.null(predictor_names) && !is.null(rownames(XtX)) &&
      !identical(rownames(XtX), predictor_names)) {
    stop("The row and column names of `XtX` must match.", call. = FALSE)
  }
  if (is.null(predictor_names)) predictor_names <- names(Xty)
  if (is.null(predictor_names)) predictor_names <- paste0("x", seq_len(p))
  if (length(predictor_names) != p || anyNA(predictor_names) ||
      any(predictor_names == "") || anyDuplicated(predictor_names)) {
    stop("Predictor names must be nonempty and unique.", call. = FALSE)
  }
  list(
    n = as.integer(n), XtX = XtX, Xty = as.numeric(Xty), yty = yty,
    X_means = if (is.null(X_means)) NULL else as.numeric(X_means),
    y_mean = y_mean, predictor_names = predictor_names
  )
}

.normalize_ss_eta <- function(ETA, predictor_names, residual_var) {
  if (!is.list(ETA) || length(ETA) < 1L) {
    stop("`ETA` must be a non-empty list.", call. = FALSE)
  }
  if ("model" %in% names(ETA)) ETA <- list(ETA1 = ETA)
  if (!all(vapply(ETA, is.list, logical(1)))) {
    stop("`ETA` must contain prior specifications.", call. = FALSE)
  }
  block_names <- names(ETA)
  if (is.null(block_names)) {
    block_names <- paste0("ETA", seq_along(ETA))
  } else {
    missing_name <- is.na(block_names) | block_names == ""
    block_names[missing_name] <- paste0("ETA", which(missing_name))
    block_names <- make.unique(block_names)
  }
  names(ETA) <- block_names
  p <- length(predictor_names)
  source_indices <- lapply(seq_along(ETA), function(block_index) {
    indices <- ETA[[block_index]]$indices
    if (is.null(indices)) {
      if (length(ETA) > 1L) {
        stop("Every ETA block must supply `indices` when using multiple blocks.",
             call. = FALSE)
      }
      return(seq_len(p))
    }
    if (is.character(indices)) {
      if (anyNA(indices) || any(!indices %in% predictor_names)) {
        stop("Character `indices` must match columns of `XtX`.", call. = FALSE)
      }
      indices <- match(indices, predictor_names)
    }
    if (!is.numeric(indices) || !is.atomic(indices) || is.object(indices) ||
        !is.null(dim(indices)) || length(indices) < 1L || anyNA(indices) ||
        any(!is.finite(indices)) || any(indices != floor(indices)) ||
        any(indices < 1L | indices > p) || anyDuplicated(indices)) {
      stop("`indices` must select unique columns of `XtX`.", call. = FALSE)
    }
    as.integer(indices)
  })
  all_indices <- unlist(source_indices, use.names = FALSE)
  if (length(all_indices) != p || anyDuplicated(all_indices) ||
      !setequal(all_indices, seq_len(p))) {
    stop("ETA block `indices` must partition all predictors exactly once.",
         call. = FALSE)
  }

  standardize <- vapply(ETA, function(specification) {
    value <- specification$standardize
    if (is.null(value)) value <- TRUE
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("Each `standardize` value must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }, logical(1))
  parser_eta <- lapply(seq_along(ETA), function(block_index) {
    specification <- ETA[[block_index]]
    specification$indices <- NULL
    k <- length(source_indices[[block_index]])
    specification$X <- matrix(rep(c(0, 1), k), nrow = 2L)
    colnames(specification$X) <- predictor_names[source_indices[[block_index]]]
    specification$standardize <- FALSE
    specification
  })
  names(parser_eta) <- block_names
  blocks <- .normalize_eta(parser_eta, 2L, residual_var)
  for (block_index in seq_along(blocks)) {
    blocks[[block_index]]$standardize <- unname(standardize[block_index])
    blocks[[block_index]]$predictor_names <-
      predictor_names[source_indices[[block_index]]]
  }
  list(blocks = blocks, source_indices = source_indices)
}

.crossproducts_to_pseudo <- function(XtX, Xty, yty = NULL) {
  decomposition <- eigen(XtX, symmetric = TRUE)
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(decomposition$values)))
  positive <- decomposition$values > tolerance
  coordinates <- drop(crossprod(decomposition$vectors, Xty))
  if (any(abs(coordinates[!positive]) >
          sqrt(tolerance) * max(1, sqrt(sum(Xty^2))))) {
    stop("`Xty` is incompatible with `XtX`.", call. = FALSE)
  }
  minimum_yty <- if (any(positive)) {
    sum(coordinates[positive]^2 / decomposition$values[positive])
  } else {
    0
  }
  joint_tolerance <- sqrt(.Machine$double.eps) * max(1, minimum_yty)
  if (is.null(yty)) {
    yty <- minimum_yty
  } else if (yty < minimum_yty - joint_tolerance) {
    stop("`yty` is incompatible with `XtX` and `Xty`.", call. = FALSE)
  } else {
    yty <- max(yty, minimum_yty)
  }
  joint <- rbind(cbind(XtX, Xty), c(Xty, yty))
  joint <- (joint + t(joint)) / 2
  joint_decomposition <- eigen(joint, symmetric = TRUE)
  joint_tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(joint_decomposition$values)))
  if (min(joint_decomposition$values) < -joint_tolerance) {
    stop("The supplied sufficient statistics are not jointly valid.",
         call. = FALSE)
  }
  retained <- joint_decomposition$values > joint_tolerance
  rank <- sum(retained)
  if (rank == 0L) {
    values <- matrix(0, nrow = 2L, ncol = ncol(joint))
  } else {
    factor <- sqrt(joint_decomposition$values[retained]) *
      t(joint_decomposition$vectors[, retained, drop = FALSE])
    basis_input <- cbind(
      rep(1, rank + 1L),
      diag(rank + 1L)[, seq_len(rank), drop = FALSE]
    )
    zero_mean_basis <- qr.Q(qr(basis_input))[, -1L, drop = FALSE]
    values <- zero_mean_basis %*% factor
  }
  list(
    x = values[, seq_len(ncol(XtX)), drop = FALSE],
    y = values[, ncol(values)]
  )
}
