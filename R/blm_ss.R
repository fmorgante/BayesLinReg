#' Bayesian linear regression from sufficient statistics
#'
#' Fits the same prior models as [blm()] directly from cross-products instead
#' of the original response and predictor matrix.
#'
#' @param n Number of observations used to form the sufficient statistics.
#' @param XtX A finite, symmetric predictor cross-product matrix, or a nonempty
#'   list of such matrices representing an exactly block-diagonal
#'   cross-product. Dense base-R matrices and compressed sparse `Matrix`
#'   objects of class `dgCMatrix` or `dsCMatrix` are supported. List and sparse
#'   inputs require `version = "Rcpp"`.
#' @param XtX_storage Storage policy for symmetric sparse matrices, supplied
#'   either directly or within list input. `"speed"` expands them to general
#'   sparse matrices, `"memory"` chooses the smaller exact representation, and
#'   `"auto"` expands them only while the estimated total representation
#'   remains within `XtX_memory_limit`.
#' @param XtX_memory_limit Positive approximate byte limit used by
#'   `XtX_storage = "auto"`. It controls internal Gram storage, not posterior
#'   draws or other fit allocations.
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
#' @param reference_response_var Optional positive reference response variance
#'   used to calibrate blocks that specify `expected_pve`. When omitted, it is
#'   recovered from `yty` and `y_mean`; if that is not possible, it must be
#'   supplied explicitly.
#' @param check_psd If `TRUE`, use a full eigendecomposition to verify that the
#'   centered `XtX` is positive semidefinite and that `Xty` and `yty` are
#'   jointly compatible with it. The default, `FALSE`, avoids this
#'   \eqn{O(p^3)} validation cost. Symmetry and basic input checks are always
#'   performed. Sparse or list `XtX` is materialized as one dense matrix for
#'   this optional validation.
#' @param nthreads Number of threads used within one Rcpp chain for list `XtX`.
#'   Values greater than one require `nchains = 1` and zero working predictor
#'   means. The default preserves the serial sampler and its RNG sequence.
#' @param likelihood_df Optional positive integer no greater than `n`, giving
#'   the likelihood dimension used in residual-variance and PVE calculations.
#'   The default is `n - 1` when `X_means` and `y_mean` identify an intercept,
#'   and `n` otherwise. For centered sufficient statistics supplied without
#'   means, set `likelihood_df = n - 1` explicitly.
#' @inheritParams blm
#'
#' @return A fitted object with the same block-specific posterior summaries as
#'   [blm()]. Intercept components are present only when both `X_means` and
#'   `y_mean` are supplied. For sparse `XtX`, the result records the selected
#'   storage representation. For list `XtX`, it also records the block-diagonal
#'   assumption and block sizes.
#'
#' @details If means are supplied, `XtX`, `Xty`, and `yty` are interpreted as
#'   uncentered cross-products and are centered internally. Without means,
#'   they are used as supplied. If any block requests standardization without
#'   means, a warning reminds the user that the supplied cross-products should
#'   already represent centered or standardized variables.
#'   For list `XtX`, predictor order is the concatenation of the matrices in
#'   list order and every omitted cross-block product is assumed to be exactly
#'   zero. Gram blocks and `ETA` prior blocks are independent: either may split
#'   or span the other. When predictor names are used, every Gram block must
#'   provide matching row and column names that are unique across the list.
#'   Centering is still exact because its global dense rank-one correction is
#'   maintained separately.
#'   When `expected_pve` is used and `reference_response_var` is omitted, the
#'   response variance is recovered from `yty` and `y_mean`; otherwise an
#'   explicit positive `reference_response_var` is required.
#'   When every block supplies `expected_pve`, omitting `residual_scale`
#'   calibrates it from their total expected PVE exactly as in [blm()].
#'   A `"Fixed"` block uses a flat prior and is jointly rank-checked with every
#'   other fixed block using the centered and standardized fixed-predictor
#'   submatrix of `XtX`. This check is always performed and is separate from
#'   the optional full-matrix `check_psd` validation.
#'   Optional posterior PVE calculations use the centered cross-products and
#'   the definitions in [blm()]. Their cost is quadratic in block size for a
#'   dense `XtX` and proportional to the relevant stored entries for sparse
#'   input, at each retained draw.
#'
#'   Coefficients are sampled with a right-hand-side update that maintains
#'   \eqn{X'y-X'X\beta}; individual-level pseudo-observations are not formed.
#'   For sparse `XtX`, the Rcpp sampler traverses only stored entries.
#'   `XtX_storage` can expand a symmetric `dsCMatrix` once to a general
#'   `dgCMatrix`, or retain a lower triangle and stream each ascending Gibbs
#'   sweep without a reverse adjacency index. This applies to both direct and
#'   list input. The full right-hand-side state is reconstructed once after
#'   each streaming sweep. The dense rank-one centering correction is
#'   maintained separately rather than materialized.
#'   With list input and zero working predictor means, `nthreads > 1` updates
#'   separate Gram blocks concurrently through `RcppParallel`. Updates remain
#'   sequential within each Gram block. Shared prior hyperparameters and the
#'   residual variance are updated after the parallel coefficient sweep. For
#'   GlobalLocal blocks, the conditionally independent local-variance GIG
#'   updates also run concurrently across Gram blocks.
#'   After [set.seed()], threaded runs are reproducible and use a different RNG
#'   stream from the serial sampler.
#'   Set `check_psd = TRUE` when the supplied statistics are not known to come
#'   from a valid common data matrix and response vector.
#'   When residual variance is learned, the sampler also checks reconstructed
#'   residual SSE values during its existing periodic state refresh and whenever
#'   an incremental value becomes suspiciously negative. Values within a
#'   scale-aware floating-point tolerance are clamped to zero; materially
#'   negative values produce an incompatibility error. This runtime guard is
#'   inexpensive but does not prove that `XtX` is positive semidefinite, so it
#'   does not replace `check_psd = TRUE`.
#'
#' @export
#'
#' @examples
#' X <- cbind(x1 = 1:20, x2 = rep(c(0, 1), 10))
#' y <- 1 + 2 * X[, "x1"] - X[, "x2"]
#' set.seed(123)
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
#'   burnin = 50
#' )
blm_ss <- function(n, XtX, Xty, ETA, yty = NULL, X_means = NULL,
                   y_mean = NULL, residual_var = NULL,
                   residual_shape = NULL, residual_scale = NULL,
                   reference_response_var = NULL,
                   iterations = 4000L, burnin = 1000L, thin = 1L,
                   version = c("Rcpp", "R"), verbose = FALSE,
                   nchains = 1L, nthreads = 1L, store_samples = FALSE,
                   store_coefficient_cov = FALSE, check_psd = FALSE,
                   XtX_storage = c("auto", "speed", "memory"),
                   XtX_memory_limit = 1024^3,
                   compute_pve = FALSE,
                   pve_type = c("standalone", "allocated"),
                   likelihood_df = NULL) {
  version <- match.arg(version)
  XtX_storage <- match.arg(XtX_storage)
  pve_controls <- .validate_pve_controls(compute_pve, pve_type)
  compute_pve <- pve_controls$compute_pve
  pve_type <- pve_controls$pve_type
  list_XtX <- is.list(XtX) && !is.matrix(XtX)
  sparse_XtX <- inherits(XtX, "sparseMatrix") ||
    (list_XtX && any(vapply(XtX, inherits, logical(1), "sparseMatrix")))
  if ((sparse_XtX || list_XtX) && version != "Rcpp") {
    stop("Sparse or list `XtX` requires `version = \"Rcpp\"`.", call. = FALSE)
  }
  nchains <- .validate_nchains(nchains)
  nthreads <- .validate_nthreads(nthreads)
  if (nthreads > 1L && nchains != 1L) {
    stop("`nthreads > 1` requires `nchains = 1`.", call. = FALSE)
  }
  if (nthreads > 1L && (!list_XtX || version != "Rcpp")) {
    stop(
      "`nthreads > 1` requires list `XtX` and `version = \"Rcpp\"`.",
      call. = FALSE
    )
  }
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
  if (!is.logical(check_psd) || length(check_psd) != 1L || is.na(check_psd)) {
    stop("`check_psd` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(XtX_memory_limit) || length(XtX_memory_limit) != 1L ||
      is.na(XtX_memory_limit) || !is.finite(XtX_memory_limit) ||
      XtX_memory_limit <= 0) {
    stop("`XtX_memory_limit` must be a positive finite number.", call. = FALSE)
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
  list_XtX <- statistics$list_input
  gram_block_names <- statistics$block_names
  gram_block_sizes <- statistics$block_sizes
  fit_intercept <- !is.null(X_means)
  likelihood_df <- .resolve_likelihood_df(n, fit_intercept, likelihood_df)

  if (is.null(yty) && is.null(residual_var)) {
    stop(
      "`residual_var` must be fixed when `yty` is not supplied.",
      call. = FALSE
    )
  }
  normalized <- .normalize_ss_eta(ETA, predictor_names, residual_var, n)
  blocks <- normalized$blocks
  source_indices <- normalized$source_indices
  has_expected_pve <- any(vapply(
    blocks, function(block) !is.null(block$expected_pve), logical(1)
  ))
  if (!is.null(reference_response_var)) {
    .validate_variance(reference_response_var, "reference_response_var")
    if (!has_expected_pve) {
      stop(
        "`reference_response_var` requires at least one `expected_pve` block.",
        call. = FALSE
      )
    }
  }
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

  XtX_diagonal <- if (list_XtX) {
    unlist(lapply(XtX, function(block) {
      if (inherits(block, "sparseMatrix")) Matrix::diag(block) else diag(block)
    }), use.names = FALSE)
  } else if (sparse_XtX) {
    Matrix::diag(XtX)
  } else {
    diag(XtX)
  }
  centered_diagonal <- if (fit_intercept) {
    XtX_diagonal - n * X_means^2
  } else {
    XtX_diagonal
  }
  centered_XtX <- if (list_XtX || sparse_XtX) {
    NULL
  } else if (fit_intercept) {
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
    pmax(1, abs(XtX_diagonal), abs(n * X_means^2))
  } else {
    pmax(1, abs(XtX_diagonal))
  }
  variance_tolerance <- 100 * .Machine$double.eps * variance_reference
  constant_predictors <- centered_diagonal <= variance_tolerance
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
    variances <- centered_diagonal[indices] / (n - 1)
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
  if (has_expected_pve) {
    if (is.null(reference_response_var)) {
      if (!fit_intercept || is.null(centered_yty)) {
        stop(
          paste0(
            "`reference_response_var` is required for `expected_pve` ",
            "calibration when centered response variance cannot be recovered."
          ),
          call. = FALSE
        )
      }
      reference_response_var <- centered_yty / (n - 1)
      .validate_variance(
        reference_response_var,
        "(yty - n * y_mean^2) / (n - 1)"
      )
    }
    predictor_variance_sums <- vapply(
      seq_along(blocks),
      function(block_index) {
        indices <- source_indices[[block_index]]
        sum(
          centered_diagonal[indices] /
            predictor_scales[[block_index]]^2 / (n - 1)
        )
      },
      numeric(1)
    )
    blocks <- .calibrate_eta_priors(
      blocks, predictor_variance_sums, reference_response_var, n
    )
  }
  residual_prior <- .prepare_residual_prior(
    residual_var, residual_shape, residual_scale, blocks,
    reference_response_var
  )
  residual_var <- residual_prior$residual_var
  residual_shape <- residual_prior$residual_shape
  residual_scale <- residual_prior$residual_scale

  layout <- .prepare_block_layout(blocks, source_indices, predictor_scales)
  source_order <- layout$source_order
  scale_order <- layout$scale_order
  block_indices <- layout$block_indices
  block_model <- layout$block_model
  internal_names <- layout$internal_names

  gram_plan <- NULL
  gram_indices <- NULL
  gram_types <- NULL
  direct_streaming_XtX <- FALSE
  if (list_XtX) {
    p <- length(predictor_names)
    original_scale <- numeric(p)
    original_scale[source_order] <- scale_order
    sampler_position <- integer(p)
    sampler_position[source_order] <- seq_len(p)
    gram_ends <- cumsum(gram_block_sizes)
    gram_starts <- gram_ends - gram_block_sizes + 1L
    original_gram_indices <- Map(seq.int, gram_starts, gram_ends)
    gram_indices <- lapply(original_gram_indices, function(indices) {
      sampler_position[indices]
    })
    gram_plan <- .plan_gram_storage(
      XtX, gram_block_names, XtX_storage, XtX_memory_limit
    )
    working_XtX <- Map(function(matrix, indices) {
      .scale_gram_block(matrix, 1 / original_scale[indices])
    }, gram_plan$blocks, original_gram_indices)
    names(working_XtX) <- gram_block_names
    gram_plan$blocks <- working_XtX
    gram_types <- gram_plan$types
    working_center <- if (fit_intercept) {
      sqrt(n) * X_means[source_order] / scale_order
    } else {
      numeric(p)
    }
    if (nthreads > 1L && any(working_center != 0)) {
      stop(
        paste0(
          "`nthreads > 1` requires zero working predictor means; use zero ",
          "`X_means` or fit without an intercept."
        ),
        call. = FALSE
      )
    }
  } else if (sparse_XtX) {
    working_XtX <- XtX[source_order, source_order, drop = FALSE]
    gram_plan <- .plan_gram_storage(
      list(XtX = working_XtX), "XtX", XtX_storage, XtX_memory_limit
    )
    working_XtX <- .scale_gram_block(
      gram_plan$blocks[[1L]], 1 / scale_order
    )
    gram_plan$blocks[[1L]] <- working_XtX
    gram_types <- gram_plan$types
    direct_streaming_XtX <- identical(gram_types[[1L]], 2L)
    working_center <- if (fit_intercept) {
      sqrt(n) * X_means[source_order] / scale_order
    } else {
      numeric(length(source_order))
    }
  } else {
    working_XtX <- centered_XtX[source_order, source_order, drop = FALSE] /
      outer(scale_order, scale_order)
    working_center <- numeric(length(source_order))
  }
  working_Xty <- centered_Xty[source_order] / scale_order
  if (!list_XtX) {
    dimnames(working_XtX) <- list(internal_names, internal_names)
  }
  names(working_Xty) <- internal_names
  fixed_indices <- .fixed_predictor_indices(blocks, block_indices)
  fixed_gram <- if (list_XtX && length(fixed_indices)) {
    .materialize_gram_blocks(
      working_XtX, gram_indices, working_center, selected = fixed_indices
    )
  } else if (list_XtX) {
    matrix(numeric(), 0L, 0L)
  } else if (sparse_XtX) {
    working_XtX - Matrix::tcrossprod(working_center)
  } else {
    working_XtX
  }
  .validate_fixed_gram(
    fixed_gram,
    if (list_XtX) seq_along(fixed_indices) else fixed_indices,
    internal_names[fixed_indices]
  )
  if (check_psd) {
    validation_XtX <- if (list_XtX) {
      .materialize_gram_blocks(
        working_XtX, gram_indices, working_center
      )
    } else if (sparse_XtX) {
      as.matrix(working_XtX) - tcrossprod(working_center)
    } else {
      working_XtX
    }
    .validate_working_crossproducts(
      validation_XtX, working_Xty, centered_yty
    )
  }

  sampler_arguments <- .prepare_sampler_arguments(
    blocks = blocks,
    layout = layout,
    y = numeric(),
    x = matrix(numeric(), nrow = 0L, ncol = length(working_Xty)),
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    residual_var = residual_var,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    store_samples = store_samples,
    store_coefficient_cov = store_coefficient_cov,
    compute_pve = compute_pve,
    pve_type = pve_type,
    effective_n = n,
    likelihood_df = likelihood_df,
    fit_intercept = fit_intercept,
    intercept_x_mean = if (fit_intercept) {
      X_means[source_order] / scale_order
    } else {
      rep(0, length(source_order))
    },
    intercept_y_mean = if (fit_intercept) y_mean else 0
  )
  sampler_arguments$XtX <- if (direct_streaming_XtX) {
    list(XtX = working_XtX)
  } else {
    working_XtX
  }
  sampler_arguments$XtX_center <- working_center
  sampler_arguments$Xty <- working_Xty
  sampler_arguments$yty <- centered_yty
  if (list_XtX || direct_streaming_XtX) {
    sampler_arguments$XtX_indices <- if (list_XtX) {
      gram_indices
    } else {
      list(seq_along(working_Xty))
    }
    sampler_arguments$XtX_types <- gram_types
    sampler_arguments$nthreads <- if (list_XtX) nthreads else 1L
  }
  samples <- .run_prepared_sampler(
    sampler_arguments, version, nchains, block_model, verbose, iterations
  )

  result <- .assemble_blm_result(
    blocks, block_indices, samples, nchains, store_samples,
    store_coefficient_cov, fit_intercept,
    compute_pve = compute_pve, pve_type = pve_type,
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    residual_scale_calibrated =
      residual_prior$residual_scale_calibrated,
    expected_pve_total = residual_prior$expected_pve_total,
    reference_response_var = reference_response_var,
    reference_residual_var = residual_prior$reference_residual_var,
    likelihood_df = likelihood_df
  )
  if (list_XtX) {
    result$XtX_representation <- "block_diagonal"
    result$XtX_number_of_blocks <- length(working_XtX)
    result$XtX_block_sizes <- stats::setNames(
      gram_block_sizes, gram_block_names
    )
    result$XtX_storage <- gram_plan$metadata
    result$XtX_cross_block_assumption <- "zero"
    result$nthreads <- nthreads
  } else if (sparse_XtX) {
    result$XtX_representation <- gram_plan$metadata$representation[[1L]]
    result$XtX_storage <- gram_plan$metadata
  }
  result
}

.validate_sufficient_statistics <- function(n, XtX, Xty, yty, X_means,
                                            y_mean) {
  n <- .validate_bounded_integer(
    n, "n", minimum = 2,
    message = paste0(
      "`n` must be an integer of at least two and no larger than ",
      "`.Machine$integer.max`."
    )
  )
  list_input <- is.list(XtX) && !is.matrix(XtX)
  if (list_input) {
    if (length(XtX) < 1L) {
      stop("List `XtX` must contain at least one Gram matrix.", call. = FALSE)
    }
    block_names <- names(XtX)
    if (is.null(block_names)) block_names <- paste0("block", seq_along(XtX))
    missing_names <- is.na(block_names) | block_names == ""
    block_names[missing_names] <- paste0("block", which(missing_names))
    block_names <- make.unique(block_names)
    blocks <- Map(function(block, name) {
      .validate_gram_block(block, sprintf("`XtX[[\"%s\"]]`", name), TRUE)
    }, XtX, block_names)
    names(blocks) <- block_names
  } else {
    block_names <- "XtX"
    blocks <- list(XtX = .validate_gram_block(XtX, "`XtX`", FALSE))
  }
  block_sizes <- vapply(blocks, ncol, integer(1))
  p <- sum(block_sizes)
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
  supplied_block_names <- vapply(
    blocks, function(block) !is.null(colnames(block)), logical(1)
  )
  if (list_input && any(supplied_block_names) && !all(supplied_block_names)) {
    stop(
      "List `XtX` must name predictors in every matrix or in none of them.",
      call. = FALSE
    )
  }
  predictor_names <- if (all(supplied_block_names)) {
    unlist(lapply(blocks, colnames), use.names = FALSE)
  } else {
    names(Xty)
  }
  if (is.null(predictor_names)) predictor_names <- paste0("x", seq_len(p))
  if (length(predictor_names) != p || anyNA(predictor_names) ||
      any(predictor_names == "") || anyDuplicated(predictor_names)) {
    stop("Predictor names must be nonempty and unique.", call. = FALSE)
  }
  if (!is.null(names(Xty)) &&
      !identical(names(Xty), predictor_names)) {
    stop("Names of `Xty` must match the `XtX` predictor names and order.",
         call. = FALSE)
  }
  if (!is.null(X_means) && !is.null(names(X_means)) &&
      !identical(names(X_means), predictor_names)) {
    stop("Names of `X_means` must match the `XtX` predictor names and order.",
         call. = FALSE)
  }
  list(
    n = as.integer(n), XtX = if (list_input) blocks else blocks[[1L]],
    Xty = as.numeric(Xty), yty = yty,
    X_means = if (is.null(X_means)) NULL else as.numeric(X_means),
    y_mean = y_mean, predictor_names = predictor_names,
    list_input = list_input, block_names = block_names,
    block_sizes = unname(block_sizes)
  )
}

.validate_gram_block <- function(matrix, label, strict_names) {
  sparse <- inherits(matrix, "sparseMatrix")
  supported_sparse <- inherits(matrix, c("dgCMatrix", "dsCMatrix"))
  valid_dense <- is.matrix(matrix) && is.numeric(matrix)
  finite <- if (sparse) {
    supported_sparse && !anyNA(matrix@x) && all(is.finite(matrix@x))
  } else {
    valid_dense && !anyNA(matrix) && all(is.finite(matrix))
  }
  if ((!valid_dense && !supported_sparse) || nrow(matrix) < 1L ||
      nrow(matrix) != ncol(matrix) || !finite) {
    stop(sprintf("%s must be a finite numeric square matrix.", label),
         call. = FALSE)
  }
  row_names <- rownames(matrix)
  column_names <- colnames(matrix)
  if (!is.null(row_names)) row_names <- as.character(row_names)
  if (!is.null(column_names)) column_names <- as.character(column_names)
  if ((strict_names && xor(is.null(row_names), is.null(column_names))) ||
      (!is.null(row_names) && !is.null(column_names) &&
       !identical(row_names, column_names))) {
    stop(sprintf("The row and column names of %s must match.", label),
         call. = FALSE)
  }
  maximum <- if (sparse && length(matrix@x) == 0L) {
    0
  } else if (sparse) {
    max(abs(matrix@x))
  } else {
    max(abs(matrix))
  }
  tolerance <- sqrt(.Machine$double.eps) * max(1, maximum)
  asymmetry <- if (inherits(matrix, "dsCMatrix")) {
    0
  } else {
    difference <- if (sparse) matrix - Matrix::t(matrix) else matrix - t(matrix)
    if (sparse && length(difference@x) == 0L) 0 else max(abs(difference))
  }
  if (asymmetry > tolerance) {
    stop(sprintf("%s must be symmetric.", label), call. = FALSE)
  }
  if (!inherits(matrix, "dsCMatrix")) {
    matrix <- if (sparse) {
      (matrix + Matrix::t(matrix)) / 2
    } else {
      (matrix + t(matrix)) / 2
    }
  }
  matrix
}

.scale_gram_block <- function(matrix, inverse_scale) {
  if (is.matrix(matrix)) {
    return(matrix * tcrossprod(inverse_scale))
  }
  column <- rep.int(seq_len(ncol(matrix)), diff(matrix@p))
  row <- matrix@i + 1L
  matrix@x <- matrix@x * inverse_scale[row] * inverse_scale[column]
  matrix
}

.gram_storage_estimates <- function(matrix) {
  p <- ncol(matrix)
  if (is.matrix(matrix)) {
    bytes <- 8 * as.double(p) * p
    return(c(dense = bytes, general = bytes, symmetric = bytes))
  }
  nonzeros <- length(matrix@x)
  if (inherits(matrix, "dgCMatrix")) {
    bytes <- 12 * as.double(nonzeros) + 4 * (p + 1)
    return(c(dense = Inf, general = bytes, symmetric = bytes))
  }
  column <- rep.int(seq_len(p), diff(matrix@p))
  diagonal <- sum(matrix@i + 1L == column)
  edges <- nonzeros - diagonal
  general <- 12 * (diagonal + 2 * as.double(edges)) + 4 * (p + 1)
  symmetric <- 12 * as.double(nonzeros) + 4 * (p + 1)
  c(dense = Inf, general = general, symmetric = symmetric)
}

.plan_gram_storage <- function(blocks, block_names, strategy, memory_limit) {
  input_class <- vapply(blocks, function(block) class(block)[1L], character(1))
  estimates <- lapply(blocks, .gram_storage_estimates)
  symmetric <- vapply(blocks, inherits, logical(1), "dsCMatrix")
  representation <- ifelse(
    vapply(blocks, is.matrix, logical(1)), "dense", "general_sparse"
  )
  if (strategy == "memory") {
    use_symmetric <- vapply(seq_along(blocks), function(index) {
      symmetric[index] &&
        estimates[[index]]["symmetric"] < estimates[[index]]["general"]
    }, logical(1))
    representation[use_symmetric] <- "symmetric_streaming"
  } else if (strategy == "auto") {
    speed_bytes <- vapply(seq_along(blocks), function(index) {
      estimates[[index]]["general"]
    }, numeric(1))
    memory_bytes <- vapply(seq_along(blocks), function(index) {
      estimates[[index]]["symmetric"]
    }, numeric(1))
    total <- sum(speed_bytes)
    candidates <- which(symmetric & memory_bytes < speed_bytes)
    if (total > memory_limit && length(candidates)) {
      savings <- speed_bytes[candidates] - memory_bytes[candidates]
      for (index in candidates[order(savings, decreasing = TRUE)]) {
        representation[index] <- "symmetric_streaming"
        total <- total - (speed_bytes[index] - memory_bytes[index])
        if (total <= memory_limit) break
      }
    }
  }

  converted <- Map(function(block, selected) {
    if (selected == "symmetric_streaming") {
      # A lower triangle lets the ascending Gibbs scan update only the current
      # and not-yet-visited coordinates without constructing reverse adjacency.
      return(Matrix::forceSymmetric(block, uplo = "L"))
    }
    if (selected != "general_sparse" || !inherits(block, "dsCMatrix")) {
      return(block)
    }
    block <- methods::as(block, "generalMatrix")
    if (!inherits(block, "dgCMatrix")) {
      block <- methods::as(block, "dgCMatrix")
    }
    block
  }, blocks, representation)
  names(converted) <- block_names
  types <- ifelse(
    representation == "dense", 0L,
    ifelse(representation == "general_sparse", 1L, 2L)
  )
  estimated_bytes <- vapply(seq_along(blocks), function(index) {
    if (representation[index] == "symmetric_streaming") {
      estimates[[index]]["symmetric"]
    } else if (representation[index] == "general_sparse") {
      estimates[[index]]["general"]
    } else {
      estimates[[index]]["dense"]
    }
  }, numeric(1))
  if (sum(estimated_bytes) > memory_limit && strategy == "auto") {
    warning(
      sprintf(
        paste0(
          "Estimated internal `XtX` storage is %.1f MiB, above the ",
          "%.1f MiB `XtX_memory_limit`; no more compact exact representation ",
          "is available for the supplied blocks."
        ),
        sum(estimated_bytes) / 1024^2, memory_limit / 1024^2
      ),
      call. = FALSE
    )
  }
  metadata <- data.frame(
    block = block_names,
    predictors = vapply(blocks, ncol, integer(1)),
    input_class = input_class,
    representation = representation,
    estimated_bytes = estimated_bytes,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  list(blocks = converted, types = as.integer(types), metadata = metadata)
}

.materialize_gram_blocks <- function(blocks, indices, center,
                                     selected = seq_along(center)) {
  selected <- as.integer(selected)
  result <- matrix(0, length(selected), length(selected))
  for (block_index in seq_along(blocks)) {
    positions <- match(indices[[block_index]], selected, nomatch = 0L)
    keep <- which(positions > 0L)
    if (!length(keep)) next
    local <- as.matrix(blocks[[block_index]][keep, keep, drop = FALSE])
    result[positions[keep], positions[keep]] <- local
  }
  result - tcrossprod(center[selected])
}

.normalize_ss_eta <- function(ETA, predictor_names, residual_var, n) {
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
  blocks <- .normalize_eta(
    parser_eta, 2L, residual_var, calibration_n = n
  )
  for (block_index in seq_along(blocks)) {
    blocks[[block_index]]$standardize <- unname(standardize[block_index])
    blocks[[block_index]]$predictor_names <-
      predictor_names[source_indices[[block_index]]]
  }
  list(blocks = blocks, source_indices = source_indices)
}

.validate_working_crossproducts <- function(XtX, Xty, yty = NULL) {
  decomposition <- eigen(XtX, symmetric = TRUE)
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(decomposition$values)))
  if (min(decomposition$values) < -tolerance) {
    stop("The centered `XtX` must be positive semidefinite.", call. = FALSE)
  }
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
  if (!is.null(yty) && yty < minimum_yty - joint_tolerance) {
    stop("`yty` is incompatible with `XtX` and `Xty`.", call. = FALSE)
  }
  invisible(NULL)
}
