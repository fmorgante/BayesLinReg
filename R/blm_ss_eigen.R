#' Bayesian linear regression from an eigen representation
#'
#' Fits the same models as [blm()] using eigenpairs of the centered predictor
#' cross-product matrix. The supplied representation may be exact or truncated.
#'
#' @param n Number of observations used to form the sufficient statistics.
#' @param XtX_eigenvectors A finite numeric matrix with one row per predictor
#'   and one column per retained eigenvector. The eigenvectors must describe
#'   the centered predictor cross-product matrix. Alternatively, a list of
#'   matrices representing an exactly block-diagonal centered cross-product.
#' @param XtX_eigenvalues A positive finite numeric vector corresponding to
#'   the columns of `XtX_eigenvectors`, or a matching list of vectors.
#' @param XtX_prop_var A number, or for list input one value per block, in
#'   `(0, 1]` declaring the proportion of the corresponding centered-`XtX`
#'   variance represented by the supplied eigenpairs.
#'   Use `1` only when all positive eigenpairs have been retained.
#' @param Xty A finite predictor-response cross-product vector on the original
#'   predictor scale.
#' @param ETA A prior specification or named list of prior blocks, with the
#'   same `indices`, prior parameters, and standardization options as [blm_ss()].
#' @param yty Optional finite, nonnegative response sum of squares. It is
#'   required when `residual_var = NULL`.
#' @param X_means,y_mean Optional predictor and response means. If supplied,
#'   `Xty` and `yty` are centered internally and an intercept is returned.
#'   The eigenpairs themselves must already represent centered `XtX`.
#' @param reference_response_var Optional positive response variance used for
#'   `expected_pve` calibration. It is recovered from `yty` and `y_mean` when
#'   possible.
#' @param check_eigenvectors If `TRUE`, check that the supplied eigenvectors
#'   are mutually orthonormal. The default avoids an
#'   \eqn{O(pq^2)} cross-product.
#' @param nthreads Number of threads used within one Rcpp chain for list eigen
#'   input. Values above one require `nchains = 1`.
#' @inheritParams blm
#'
#' @return A `blm_fit` object. In addition to the usual fields, it records
#'   `XtX_eigen_rank`, `XtX_prop_var`, `XtX_approximate`, and
#'   `residual_sse_offset`. List input also records block sizes, block-specific
#'   ranks and retained fractions, the cross-block assumption, and `nthreads`.
#'
#' @details Let \eqn{G=X_c'X_c}, \eqn{s=X_c'y_c}, and suppose the supplied
#'   eigenpairs give \eqn{G_q=U_q\Lambda_qU_q'}. The sampler uses
#'   \deqn{Q=\Lambda_q^{1/2}U_q',\qquad
#'     w=\Lambda_q^{-1/2}U_q's}
#'   and maintains the transformed residual \eqn{w-Q\beta}. With
#'   `XtX_prop_var = 1`, this is an exact re-expression of the coefficient
#'   likelihood, subject to numerical precision. When `yty` is supplied, the
#'   response variation orthogonal to the predictor space is retained as a
#'   beta-independent SSE offset, so learning the residual variance also
#'   matches [blm_ss()].
#'
#'   Values of `XtX_prop_var` below one explicitly request an approximate
#'   posterior based on the supplied eigenspace. Components of `Xty` outside
#'   that space do not affect coefficient sampling and are included in the
#'   beta-independent SSE offset. Standardization uses the diagonal of
#'   \eqn{G_q}; aggressive truncation can therefore make a predictor appear
#'   constant. All `"Fixed"` predictors must be full rank within the retained
#'   eigenspace; this condition is checked before sampling.
#'
#'   This entry point is implemented only with the Rcpp sampler. It does not
#'   compute an eigendecomposition internally, allowing a precomputed
#'   representation to be reused without an \eqn{O(p^3)} initialization.
#'   Optional posterior PVE calculations use the retained representation, at
#'   a cost of approximately \eqn{O(qp)} per retained draw across all blocks.
#'
#'   With list input, each matrix and corresponding eigenvalue vector defines
#'   one block of the centered cross-product matrix; all cross-block products
#'   are assumed to be exactly zero. Predictor order is obtained by
#'   concatenating the matrix rows. The implementation stores each
#'   block-specific transformed design and residual separately, requiring
#'   approximately \eqn{\sum_b p_b q_b} values instead of a global
#'   \eqn{p\sum_b q_b} pseudo-design. ETA prior blocks may cross eigen blocks.
#'   When `nthreads > 1`, coefficient sweeps are sequential within each eigen
#'   block and concurrent across blocks. Shared prior and residual parameters
#'   are updated after the workers join. Nonzero `X_means` are supported
#'   because the supplied eigenpairs must already represent centered `XtX`;
#'   independently decomposed uncentered blocks are not valid substitutes.
#'
#' @export
#'
#' @examples
#' X <- cbind(x1 = 1:20, x2 = rep(c(0, 1), 10))
#' y <- 1 + 2 * X[, "x1"] - X[, "x2"]
#' X_centered <- sweep(X, 2, colMeans(X), FUN = "-")
#' decomposition <- eigen(crossprod(X_centered), symmetric = TRUE)
#' keep <- decomposition$values > 1e-10
#' fit <- blm_ss_eigen(
#'   n = nrow(X),
#'   XtX_eigenvectors = decomposition$vectors[, keep, drop = FALSE],
#'   XtX_eigenvalues = decomposition$values[keep],
#'   XtX_prop_var = 1,
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
blm_ss_eigen <- function(
    n, XtX_eigenvectors, XtX_eigenvalues, XtX_prop_var, Xty, ETA,
    yty = NULL, X_means = NULL, y_mean = NULL, residual_var = NULL,
    residual_shape = NULL, residual_scale = NULL,
    reference_response_var = NULL,
    iterations = 4000L, burnin = 1000L, thin = 1L, seed = NULL,
    verbose = FALSE, nchains = 1L, nthreads = 1L, store_samples = TRUE,
    store_coefficient_cov = TRUE, check_eigenvectors = FALSE,
    compute_pve = FALSE,
    pve_type = c("standalone", "allocated")) {
  pve_controls <- .validate_pve_controls(compute_pve, pve_type)
  compute_pve <- pve_controls$compute_pve
  pve_type <- pve_controls$pve_type
  controls <- list(
    verbose = verbose,
    store_samples = store_samples,
    store_coefficient_cov = store_coefficient_cov,
    check_eigenvectors = check_eigenvectors
  )
  for (control_name in names(controls)) {
    value <- controls[[control_name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", control_name), call. = FALSE)
    }
  }
  nchains <- .validate_nchains(nchains)
  nthreads <- .validate_nthreads(nthreads)
  if (nthreads > 1L && nchains != 1L) {
    stop("`nthreads > 1` requires `nchains = 1`.", call. = FALSE)
  }
  list_eigen <- is.list(XtX_eigenvectors) && !is.matrix(XtX_eigenvectors)
  if (nthreads > 1L && !list_eigen) {
    stop(
      "`nthreads > 1` requires list `XtX_eigenvectors`.",
      call. = FALSE
    )
  }
  statistics <- .validate_eigen_sufficient_statistics(
    n, XtX_eigenvectors, XtX_eigenvalues, XtX_prop_var, Xty, yty,
    X_means, y_mean, check_eigenvectors
  )
  n <- statistics$n
  eigenvectors <- statistics$eigenvectors
  eigenvalues <- statistics$eigenvalues
  XtX_prop_var <- statistics$prop_var
  Xty <- statistics$Xty
  yty <- statistics$yty
  X_means <- statistics$X_means
  y_mean <- statistics$y_mean
  predictor_names <- statistics$predictor_names
  list_eigen <- statistics$block_input
  fit_intercept <- !is.null(X_means)

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
          "eigen representation should be centered or standardized."
        ),
        call. = FALSE
      )
    }
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
  yty_tolerance <- sqrt(.Machine$double.eps) *
    max(1, if (is.null(yty)) 0 else abs(yty))
  if (!is.null(centered_yty) && centered_yty < -yty_tolerance) {
    stop("The centered `yty` must be nonnegative.", call. = FALSE)
  }
  if (!is.null(centered_yty)) centered_yty <- max(0, centered_yty)

  eigenvector_blocks <- if (list_eigen) eigenvectors else list(eigenvectors)
  eigenvalue_blocks <- if (list_eigen) eigenvalues else list(eigenvalues)
  prop_var_blocks <- if (list_eigen) XtX_prop_var else c(XtX_prop_var)
  source_eigen_indices <- statistics$eigen_indices
  projected_blocks <- lapply(seq_along(eigenvector_blocks), function(block) {
    drop(crossprod(
      eigenvector_blocks[[block]], centered_Xty[source_eigen_indices[[block]]]
    ))
  })
  transformed_y_blocks <- Map(
    function(projected, values) projected / sqrt(values),
    projected_blocks, eigenvalue_blocks
  )
  projected_centered_Xty <- numeric(length(centered_Xty))
  for (block in seq_along(eigenvector_blocks)) {
    indices <- source_eigen_indices[[block]]
    projected_centered_Xty[indices] <- drop(
      eigenvector_blocks[[block]] %*% projected_blocks[[block]]
    )
    if (prop_var_blocks[block] == 1) {
      projection_error <- centered_Xty[indices] -
        projected_centered_Xty[indices]
      projection_tolerance <- sqrt(.Machine$double.eps) *
        max(1, sqrt(sum(centered_Xty[indices]^2)))
      if (sqrt(sum(projection_error^2)) > projection_tolerance) {
        stop(
          sprintf(
            paste0(
              "`Xty` has a component outside the supplied exact eigenspace ",
              "for eigen block `%s`; the eigenpairs and sufficient ",
              "statistics are incompatible."
            ),
            names(eigenvector_blocks)[block]
          ),
          call. = FALSE
        )
      }
    }
  }
  transformed_y_norm <- sum(vapply(
    transformed_y_blocks, function(value) sum(value^2), numeric(1)
  ))
  residual_sse_offset <- if (is.null(centered_yty)) {
    0
  } else {
    centered_yty - transformed_y_norm
  }
  offset_tolerance <- sqrt(.Machine$double.eps) *
    max(1, if (is.null(centered_yty)) 0 else centered_yty)
  if (residual_sse_offset < -offset_tolerance) {
    stop(
      paste0(
        "`yty` is incompatible with `Xty` and the supplied ",
        "eigen representation."
      ),
      call. = FALSE
    )
  }
  residual_sse_offset <- max(0, residual_sse_offset)

  approximate_diagonal <- numeric(length(centered_Xty))
  for (block in seq_along(eigenvector_blocks)) {
    approximate_diagonal[source_eigen_indices[[block]]] <- drop(
      eigenvector_blocks[[block]]^2 %*% eigenvalue_blocks[[block]]
    )
  }
  diagonal_tolerance <- 100 * .Machine$double.eps *
    pmax(1, approximate_diagonal)
  constant_predictors <- approximate_diagonal <= diagonal_tolerance
  if (any(constant_predictors)) {
    stop(
      sprintf(
        "The eigen representation contains constant predictor(s): %s.",
        paste(predictor_names[constant_predictors], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  predictor_scales <- lapply(seq_along(blocks), function(block_index) {
    indices <- source_indices[[block_index]]
    if (!blocks[[block_index]]$standardize) return(rep(1, length(indices)))
    sqrt(approximate_diagonal[indices] / (n - 1))
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
          approximate_diagonal[indices] /
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

  source_order <- unlist(source_indices, use.names = FALSE)
  scale_order <- unlist(predictor_scales, use.names = FALSE)
  inverse_source_order <- match(seq_along(source_order), source_order)
  if (list_eigen) {
    transformed_X_blocks <- lapply(seq_along(eigenvector_blocks), function(b) {
      indices <- source_eigen_indices[[b]]
      design <- sqrt(eigenvalue_blocks[[b]]) * t(eigenvector_blocks[[b]])
      sweep(
        design, 2L, scale_order[inverse_source_order[indices]], FUN = "/"
      )
    })
    eigen_internal_indices <- lapply(
      source_eigen_indices,
      function(indices) inverse_source_order[indices]
    )
    transformed_X <- NULL
  } else {
    transformed_X <- sqrt(eigenvalues) * t(eigenvectors)
    transformed_X <- sweep(
      transformed_X[, source_order, drop = FALSE],
      2L, scale_order, FUN = "/"
    )
  }
  block_sizes <- lengths(source_indices)
  block_ends <- cumsum(block_sizes)
  block_starts <- block_ends - block_sizes + 1L
  block_indices <- Map(seq.int, block_starts, block_ends)
  block_model <- vapply(blocks, `[[`, integer(1), "model_code")
  block_id <- rep.int(seq_along(blocks), block_sizes)
  internal_names <- unlist(lapply(seq_along(blocks), function(block_index) {
    paste0(names(blocks)[block_index], "::", blocks[[block_index]]$predictor_names)
  }))
  fixed_indices <- .fixed_predictor_indices(blocks, block_indices)
  fixed_design <- if (list_eigen && length(fixed_indices) > 0L) {
    do.call(rbind, lapply(seq_along(transformed_X_blocks), function(b) {
      answer <- matrix(
        0, nrow = nrow(transformed_X_blocks[[b]]),
        ncol = length(fixed_indices)
      )
      selected <- match(eigen_internal_indices[[b]], fixed_indices, nomatch = 0L)
      keep <- selected > 0L
      answer[, selected[keep]] <- transformed_X_blocks[[b]][, keep, drop = FALSE]
      answer
    }))
  } else {
    transformed_X
  }
  if (!is.null(transformed_X)) colnames(transformed_X) <- internal_names
  .validate_fixed_design(
    fixed_design,
    if (list_eigen) seq_along(fixed_indices) else fixed_indices,
    internal_names[fixed_indices]
  )

  sampler_arguments <- list(
    y = if (list_eigen) numeric() else transformed_y_blocks[[1L]],
    x = if (list_eigen) {
      matrix(numeric(), nrow = 0L, ncol = length(source_order))
    } else {
      transformed_X
    },
    residual_shape = if (is.null(residual_shape)) 1 else residual_shape,
    residual_scale = if (is.null(residual_scale)) 1 else residual_scale,
    residual_var = residual_var,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    block_id = block_id,
    block_model = block_model,
    normal_shape = vapply(blocks, `[[`, numeric(1), "normal_shape"),
    normal_scale = vapply(blocks, `[[`, numeric(1), "normal_scale"),
    pi_alpha = vapply(blocks, `[[`, numeric(1), "pi_alpha"),
    pi_beta = vapply(blocks, `[[`, numeric(1), "pi_beta"),
    spike_var_shape = vapply(
      blocks, `[[`, numeric(1), "spike_var_shape"
    ),
    spike_var_scale = vapply(
      blocks, `[[`, numeric(1), "spike_var_scale"
    ),
    global_scale = vapply(blocks, `[[`, numeric(1), "global_scale"),
    local_a = vapply(
      blocks, function(block) block$local_shape[1L], numeric(1)
    ),
    local_b = vapply(
      blocks, function(block) block$local_shape[2L], numeric(1)
    ),
    multi_gamma = lapply(blocks, `[[`, "multi_gamma"),
    multi_pi_alpha = lapply(blocks, `[[`, "multi_pi_alpha"),
    multi_var_shape = vapply(
      blocks, `[[`, numeric(1), "multi_var_shape"
    ),
    multi_var_scale = vapply(
      blocks, `[[`, numeric(1), "multi_var_scale"
    ),
    store_samples = store_samples,
    store_coefficient_cov = store_coefficient_cov,
    compute_pve = compute_pve,
    pve_type = pve_type,
    effective_n = n,
    fit_intercept = fit_intercept,
    intercept_x_mean = if (fit_intercept) {
      X_means[source_order] / scale_order
    } else {
      rep(0, length(source_order))
    },
    intercept_y_mean = if (fit_intercept) y_mean else 0,
    center_observations = FALSE,
    residual_sse_offset = residual_sse_offset
  )
  if (list_eigen) {
    sampler_arguments$Xty <- projected_centered_Xty[source_order] / scale_order
    sampler_arguments$yty <- if (is.null(centered_yty)) 0 else centered_yty
    sampler_arguments$eigen_X <- transformed_X_blocks
    sampler_arguments$eigen_y <- transformed_y_blocks
    sampler_arguments$eigen_indices <- eigen_internal_indices
    sampler_arguments$nthreads <- nthreads
    sampler_arguments$residual_sse_offset <- 0
  }
  run_chains <- function(progressor = NULL) {
    .run_blm_chains(
      sampler_arguments, "Rcpp", nchains, seed, block_model, progressor
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
    reference_residual_var = residual_prior$reference_residual_var
  )
  eigen_rank_by_block <- vapply(eigenvalue_blocks, length, integer(1))
  retained_trace <- vapply(eigenvalue_blocks, sum, numeric(1))
  aggregate_prop_var <- if (list_eigen) {
    sum(retained_trace) / sum(retained_trace / prop_var_blocks)
  } else {
    prop_var_blocks[[1L]]
  }
  result$XtX_eigen_rank <- sum(eigen_rank_by_block)
  result$XtX_prop_var <- aggregate_prop_var
  result$XtX_approximate <- any(prop_var_blocks < 1)
  result$residual_sse_offset <- residual_sse_offset
  if (list_eigen) {
    result$XtX_representation <- "block_diagonal_eigen"
    result$XtX_number_of_blocks <- length(eigenvector_blocks)
    result$XtX_block_sizes <- stats::setNames(
      vapply(eigenvector_blocks, nrow, integer(1)),
      names(eigenvector_blocks)
    )
    result$XtX_eigen_rank_by_block <- stats::setNames(
      eigen_rank_by_block, names(eigenvector_blocks)
    )
    result$XtX_prop_var_by_block <- stats::setNames(
      prop_var_blocks, names(eigenvector_blocks)
    )
    result$XtX_cross_block_assumption <- "zero"
    result$nthreads <- nthreads
  }
  result
}

.validate_eigen_sufficient_statistics <- function(
    n, eigenvectors, eigenvalues, prop_var, Xty, yty, X_means, y_mean,
    check_eigenvectors) {
  block_input <- is.list(eigenvectors) && !is.matrix(eigenvectors)
  if (block_input) {
    if (length(eigenvectors) < 1L) {
      stop(
        "List `XtX_eigenvectors` must contain at least one matrix.",
        call. = FALSE
      )
    }
    if (!is.list(eigenvalues) || length(eigenvalues) != length(eigenvectors)) {
      stop(
        paste0(
          "List `XtX_eigenvalues` must have the same length as ",
          "`XtX_eigenvectors`."
        ),
        call. = FALSE
      )
    }
    block_names <- names(eigenvectors)
    if (is.null(block_names)) block_names <- paste0("block", seq_along(eigenvectors))
    missing_name <- is.na(block_names) | block_names == ""
    block_names[missing_name] <- paste0("block", which(missing_name))
    block_names <- make.unique(block_names)
    names(eigenvectors) <- block_names
    eigenvalue_names <- names(eigenvalues)
    if (!is.null(eigenvalue_names)) {
      if (anyDuplicated(eigenvalue_names) ||
          !setequal(eigenvalue_names, block_names)) {
        stop(
          "Named eigenvalue lists must match the eigenvector block names.",
          call. = FALSE
        )
      }
      eigenvalues <- eigenvalues[block_names]
    }
    names(eigenvalues) <- block_names
    if (is.list(prop_var)) prop_var <- unlist(prop_var, use.names = TRUE)
    if (!is.numeric(prop_var) || anyNA(prop_var) ||
        any(!is.finite(prop_var))) {
      stop("`XtX_prop_var` must contain finite numeric values.", call. = FALSE)
    }
    if (length(prop_var) == 1L) {
      prop_var <- rep(unname(prop_var), length(eigenvectors))
    }
    if (length(prop_var) != length(eigenvectors)) {
      stop(
        "List eigen input requires one `XtX_prop_var` value per block.",
        call. = FALSE
      )
    }
    if (!is.null(names(prop_var)) && any(names(prop_var) != "")) {
      if (anyDuplicated(names(prop_var)) ||
          !setequal(names(prop_var), block_names)) {
        stop(
          "Named `XtX_prop_var` values must match the eigen block names.",
          call. = FALSE
        )
      }
      prop_var <- prop_var[block_names]
    }
    if (any(prop_var <= 0 | prop_var > 1)) {
      stop(
        "Every `XtX_prop_var` value must be between zero (exclusive) and one.",
        call. = FALSE
      )
    }
    valid_matrix <- vapply(eigenvectors, function(value) {
      is.matrix(value) && is.numeric(value) && nrow(value) > 0L
    }, logical(1))
    if (!all(valid_matrix)) {
      stop(
        "Every `XtX_eigenvectors` block must be a nonempty numeric matrix.",
        call. = FALSE
      )
    }
    block_sizes <- vapply(eigenvectors, nrow, integer(1))
    p <- sum(block_sizes)
    if (is.matrix(Xty) && is.numeric(Xty) && identical(dim(Xty), c(p, 1L))) {
      Xty <- drop(Xty)
    }
    if (!is.numeric(Xty) || !is.atomic(Xty) || is.object(Xty) ||
        !is.null(dim(Xty)) || length(Xty) != p || anyNA(Xty) ||
        any(!is.finite(Xty))) {
      stop("`Xty` must be a finite numeric vector with one value per predictor.",
           call. = FALSE)
    }
    if (xor(is.null(X_means), is.null(y_mean))) {
      stop("Supply `X_means` and `y_mean` together.", call. = FALSE)
    }
    if (!is.null(X_means) && length(X_means) != p) {
      stop("`X_means` must have one finite value per predictor.", call. = FALSE)
    }
    ends <- cumsum(block_sizes)
    starts <- ends - block_sizes + 1L
    eigen_indices <- Map(seq.int, starts, ends)
    supplied_predictor_names <- !is.null(names(Xty)) ||
      any(vapply(eigenvectors, function(value) !is.null(rownames(value)), logical(1)))
    validated <- lapply(seq_along(eigenvectors), function(block) {
      indices <- eigen_indices[[block]]
      .validate_eigen_sufficient_statistics(
        n, eigenvectors[[block]], eigenvalues[[block]], prop_var[[block]],
        Xty[indices], yty,
        if (is.null(X_means)) NULL else X_means[indices], y_mean,
        check_eigenvectors
      )
    })
    predictor_names <- unlist(lapply(validated, `[[`, "predictor_names"))
    if (!is.null(names(Xty))) {
      missing_predictor_name <- is.na(predictor_names) | predictor_names == ""
      predictor_names[missing_predictor_name] <- names(Xty)[missing_predictor_name]
    }
    predictor_names <- if (supplied_predictor_names) {
      make.unique(predictor_names)
    } else {
      paste0("x", seq_len(p))
    }
    eigenvectors <- lapply(validated, `[[`, "eigenvectors")
    eigenvalues <- lapply(validated, `[[`, "eigenvalues")
    names(eigenvectors) <- names(eigenvalues) <- block_names
    return(list(
      n = validated[[1L]]$n,
      eigenvectors = eigenvectors,
      eigenvalues = eigenvalues,
      prop_var = stats::setNames(as.numeric(prop_var), block_names),
      Xty = as.numeric(Xty),
      yty = validated[[1L]]$yty,
      X_means = if (is.null(X_means)) NULL else as.numeric(X_means),
      y_mean = if (is.null(y_mean)) NULL else as.numeric(y_mean),
      predictor_names = predictor_names,
      eigen_indices = eigen_indices,
      block_input = TRUE
    ))
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n != floor(n) || n < 2) {
    stop("`n` must be an integer of at least two.", call. = FALSE)
  }
  if (!is.matrix(eigenvectors) || !is.numeric(eigenvectors) ||
      nrow(eigenvectors) < 1L || ncol(eigenvectors) < 1L ||
      ncol(eigenvectors) > nrow(eigenvectors) || anyNA(eigenvectors) ||
      any(!is.finite(eigenvectors))) {
    stop(
      "`XtX_eigenvectors` must be a finite numeric p-by-q matrix with q <= p.",
      call. = FALSE
    )
  }
  q <- ncol(eigenvectors)
  if (!is.numeric(eigenvalues) || !is.atomic(eigenvalues) ||
      is.object(eigenvalues) || !is.null(dim(eigenvalues)) ||
      length(eigenvalues) != q || anyNA(eigenvalues) ||
      any(!is.finite(eigenvalues)) || any(eigenvalues <= 0)) {
    stop(
      paste0(
        "`XtX_eigenvalues` must contain one positive finite value per ",
        "eigenvector."
      ),
      call. = FALSE
    )
  }
  if (!is.numeric(prop_var) || length(prop_var) != 1L ||
      is.na(prop_var) || !is.finite(prop_var) ||
      prop_var <= 0 || prop_var > 1) {
    stop("`XtX_prop_var` must be between zero (exclusive) and one.",
         call. = FALSE)
  }
  if (check_eigenvectors) {
    orthogonality_error <- max(abs(
      crossprod(eigenvectors) - diag(q)
    ))
    tolerance <- sqrt(.Machine$double.eps) * max(1, nrow(eigenvectors))
    if (orthogonality_error > tolerance) {
      stop("`XtX_eigenvectors` must have orthonormal columns.",
           call. = FALSE)
    }
  }

  p <- nrow(eigenvectors)
  if (is.matrix(Xty) && is.numeric(Xty) &&
      identical(dim(Xty), c(p, 1L))) {
    Xty <- drop(Xty)
  }
  if (!is.numeric(Xty) || !is.atomic(Xty) || is.object(Xty) ||
      !is.null(dim(Xty)) || length(Xty) != p || anyNA(Xty) ||
      any(!is.finite(Xty))) {
    stop("`Xty` must be a finite numeric vector with one value per predictor.",
         call. = FALSE)
  }
  if (!is.null(yty) &&
      (!is.numeric(yty) || length(yty) != 1L || is.na(yty) ||
       !is.finite(yty) || yty < 0)) {
    stop("`yty` must be NULL or a finite nonnegative number.",
         call. = FALSE)
  }
  if (xor(is.null(X_means), is.null(y_mean))) {
    stop("Supply `X_means` and `y_mean` together.", call. = FALSE)
  }
  if (!is.null(X_means)) {
    if (!is.numeric(X_means) || !is.atomic(X_means) ||
        is.object(X_means) || !is.null(dim(X_means)) ||
        length(X_means) != p || anyNA(X_means) ||
        any(!is.finite(X_means))) {
      stop("`X_means` must have one finite value per predictor.",
           call. = FALSE)
    }
    if (!is.numeric(y_mean) || length(y_mean) != 1L || is.na(y_mean) ||
        !is.finite(y_mean)) {
      stop("`y_mean` must be a finite number.", call. = FALSE)
    }
  }
  predictor_names <- rownames(eigenvectors)
  if (is.null(predictor_names)) predictor_names <- names(Xty)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("x", seq_len(p))
  } else {
    missing_name <- is.na(predictor_names) | predictor_names == ""
    predictor_names[missing_name] <- paste0("x", which(missing_name))
    predictor_names <- make.unique(predictor_names)
  }
  storage.mode(eigenvectors) <- "double"
  list(
    n = as.integer(n),
    eigenvectors = eigenvectors,
    eigenvalues = as.numeric(eigenvalues),
    prop_var = as.numeric(prop_var),
    Xty = as.numeric(Xty),
    yty = if (is.null(yty)) NULL else as.numeric(yty),
    X_means = if (is.null(X_means)) NULL else as.numeric(X_means),
    y_mean = if (is.null(y_mean)) NULL else as.numeric(y_mean),
    predictor_names = predictor_names,
    eigen_indices = list(seq_len(p)),
    block_input = FALSE
  )
}
