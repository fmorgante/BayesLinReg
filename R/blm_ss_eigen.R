#' Bayesian linear regression from an eigen representation
#'
#' Fits the same models as [blm()] using eigenpairs of the centered predictor
#' cross-product matrix. The supplied representation may be exact or truncated.
#'
#' @param n Number of observations used to form the sufficient statistics.
#' @param XtX_eigenvectors A finite numeric matrix with one row per predictor
#'   and one column per retained eigenvector. The eigenvectors must describe
#'   the centered predictor cross-product matrix.
#' @param XtX_eigenvalues A positive finite numeric vector corresponding to
#'   the columns of `XtX_eigenvectors`.
#' @param XtX_prop_var A number in `(0, 1]` declaring the proportion of the
#'   total centered-`XtX` variance represented by the supplied eigenpairs.
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
#' @inheritParams blm
#'
#' @return A `blm_fit` object. In addition to the usual fields, it records
#'   `XtX_eigen_rank`, `XtX_prop_var`, `XtX_approximate`, and
#'   `residual_sse_offset`.
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
#'   constant.
#'
#'   This entry point is implemented only with the Rcpp sampler. It does not
#'   compute an eigendecomposition internally, allowing a precomputed
#'   representation to be reused without an \eqn{O(p^3)} initialization.
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
    verbose = FALSE, nchains = 1L, store_samples = TRUE,
    store_coefficient_cov = TRUE, check_eigenvectors = FALSE) {
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

  projected_Xty <- drop(crossprod(eigenvectors, centered_Xty))
  transformed_y <- projected_Xty / sqrt(eigenvalues)
  if (XtX_prop_var == 1) {
    projection_error <- centered_Xty -
      drop(eigenvectors %*% projected_Xty)
    projection_tolerance <- sqrt(.Machine$double.eps) *
      max(1, sqrt(sum(centered_Xty^2)))
    if (sqrt(sum(projection_error^2)) > projection_tolerance) {
      stop(
        paste0(
          "`Xty` has a component outside the supplied exact eigenspace; ",
          "the eigenpairs and sufficient statistics are incompatible."
        ),
        call. = FALSE
      )
    }
  }
  residual_sse_offset <- if (is.null(centered_yty)) {
    0
  } else {
    centered_yty - sum(transformed_y^2)
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

  approximate_diagonal <- drop(
    eigenvectors^2 %*% eigenvalues
  )
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
  transformed_X <- sqrt(eigenvalues) * t(eigenvectors)
  transformed_X <- sweep(
    transformed_X[, source_order, drop = FALSE],
    2L, scale_order, FUN = "/"
  )
  block_sizes <- lengths(source_indices)
  block_ends <- cumsum(block_sizes)
  block_starts <- block_ends - block_sizes + 1L
  block_indices <- Map(seq.int, block_starts, block_ends)
  block_model <- vapply(blocks, `[[`, integer(1), "model_code")
  block_id <- rep.int(seq_along(blocks), block_sizes)
  internal_names <- unlist(lapply(seq_along(blocks), function(block_index) {
    paste0(names(blocks)[block_index], "::", blocks[[block_index]]$predictor_names)
  }))
  colnames(transformed_X) <- internal_names

  sampler_arguments <- list(
    y = transformed_y,
    x = transformed_X,
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
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    residual_scale_calibrated =
      residual_prior$residual_scale_calibrated,
    expected_pve_total = residual_prior$expected_pve_total,
    reference_response_var = reference_response_var,
    reference_residual_var = residual_prior$reference_residual_var
  )
  result$XtX_eigen_rank <- length(eigenvalues)
  result$XtX_prop_var <- XtX_prop_var
  result$XtX_approximate <- XtX_prop_var < 1
  result$residual_sse_offset <- residual_sse_offset
  result
}

.validate_eigen_sufficient_statistics <- function(
    n, eigenvectors, eigenvalues, prop_var, Xty, yty, X_means, y_mean,
    check_eigenvectors) {
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
    predictor_names = predictor_names
  )
}
