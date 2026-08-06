#' Compute regression sufficient statistics from GWAS results
#'
#' Reconstructs centered linear-regression cross-products from marginal GWAS
#' effect estimates, their standard errors, and signed LD correlations using
#' the finite-sample transformation employed by SuSiE-RSS.
#'
#' @param beta A finite numeric vector of marginal ordinary least-squares
#'   effect estimates.
#' @param se A positive finite numeric vector of standard errors corresponding
#'   to `beta`.
#' @param LD A finite, symmetric LD correlation matrix, or a nonempty list of
#'   such matrices representing exactly zero cross-block LD. Predictors follow
#'   matrix order or the concatenated list order and must have the same
#'   effect-allele orientation as `beta`. Dense base-R matrices and compressed
#'   sparse `Matrix` objects of class `dgCMatrix` or `dsCMatrix` are supported.
#'   Supply signed correlations, not squared correlations.
#' @param n Number of observations used to calculate the GWAS statistics.
#' @param response_var Optional positive sample variance of the response. It is
#'   required for `scale = "original"`.
#' @param scale Scale on which to construct the sufficient statistics.
#'   `"auto"` uses the original coefficient scale when `response_var` is
#'   supplied and otherwise standardizes both predictors and response.
#' @param residual_df Positive residual degrees of freedom used by the marginal
#'   regressions. The SuSiE-RSS ordinary least-squares conversion uses the
#'   default `n - 2`.
#' @param output Whether to return cross-products for [blm_ss()] or a complete
#'   raw eigendecomposition for preparing input to [blm_ss_eigen()]. For list
#'   `LD`, each scaled block is decomposed independently.
#'
#' @return For `output = "sufficient"`, a list containing `n`, `XtX`, `Xty`,
#'   `yty`, zero `X_means` and `y_mean`, and `reference_response_var`. `XtX`
#'   is a matrix for matrix `LD` and a same-structure list for list `LD`. Pass
#'   `reference_response_var` to [blm_ss()] only when at least one `ETA` block
#'   specifies `expected_pve`. For `output = "eigen"`, `XtX` is replaced by
#'   `XtX_eigenvectors_raw`, `XtX_eigenvalues_raw`, and
#'   `XtX_eigenvalue_tolerance`. These are matrices and a scalar for matrix
#'   `LD`, or matching named lists and a named vector for list `LD`. Raw
#'   eigenpairs are returned without filtering or modification and therefore
#'   are not necessarily valid inputs to [blm_ss_eigen()], which requires
#'   strictly positive eigenvalues.
#'
#' @details Let \eqn{z_j=\hat b_j/\hat s_j} and
#'   \deqn{a_j=(n-1)/(z_j^2+\nu),}
#'   where \eqn{\nu} is `residual_df`. On the standardized scale, the function
#'   returns
#'   \deqn{X'X=(n-1)R,\quad
#'         X'y=\sqrt{n-1}\,\sqrt{a}\,z,\quad y'y=n-1.}
#'   When `response_var` is supplied and the original scale is requested, it
#'   sets
#'   \deqn{d_j=response\_var\,a_j/\hat s_j^2,}
#'   \deqn{X'X=D^{1/2}RD^{1/2},\quad X'y=D\hat b,\quad
#'         y'y=(n-1)response\_var.}
#'
#'   These are exact centered cross-products when the marginal statistics are
#'   ordinary least-squares results from a common sample and `LD` is the
#'   corresponding in-sample predictor correlation matrix. With reference LD,
#'   meta-analysis statistics, variable per-variant sample sizes, mixed-model
#'   statistics, or other non-OLS inputs, they are approximate working
#'   cross-products. In particular, learning the residual variance is generally
#'   safer with in-sample LD; fixing it is generally safer with reference LD.
#'   For list `LD`, every omitted cross-block correlation is assumed to be
#'   exactly zero and the returned `XtX` has the same list structure. Scaling
#'   is performed within each block, and symmetric sparse blocks retain their
#'   triangular representation for [blm_ss()] storage planning.
#'
#'   The function performs only basic input validation. It does not diagnose
#'   allele flips or broader incompatibility between the GWAS statistics and
#'   LD. Eigen output uses a complete dense eigendecomposition. Matrix input
#'   requires \eqn{O(p^2)} memory and \eqn{O(p^3)} time; list input processes
#'   one scaled block at a time and requires peak eigendecomposition storage
#'   for the largest block. Every eigenpair is returned, including zero and
#'   negative eigenvalues. Block-specific warnings identify negative
#'   eigenvalues. Users must choose the eigenpairs passed to [blm_ss_eigen()]
#'   and calculate the corresponding `XtX_prop_var` values.
#'
#' @references
#' Zou Y, Carbonetto P, Wang G, Stephens M (2022). Fine-mapping from summary
#' data with the "Sum of Single Effects" model. PLoS Genetics 18:e1010299.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' X <- scale(matrix(rnorm(200), 40, 5), center = TRUE, scale = FALSE)
#' y <- drop(scale(rnorm(40), center = TRUE, scale = FALSE))
#' marginal_beta <- drop(crossprod(X, y)) / colSums(X^2)
#' marginal_se <- sqrt(
#'   (sum(y^2) - marginal_beta^2 * colSums(X^2)) /
#'     ((nrow(X) - 2) * colSums(X^2))
#' )
#' ss <- compute_ss_from_gwas(
#'   marginal_beta, marginal_se, stats::cor(X), nrow(X),
#'   response_var = stats::var(y)
#' )
compute_ss_from_gwas <- function(
    beta, se, LD, n, response_var = NULL,
    scale = c("auto", "standardized", "original"), residual_df = n - 2,
    output = c("sufficient", "eigen")) {
  scale <- match.arg(scale)
  output <- match.arg(output)

  if (!is.numeric(beta) || !is.atomic(beta) || is.object(beta) ||
      !is.null(dim(beta)) || length(beta) < 1L || anyNA(beta) ||
      any(!is.finite(beta))) {
    stop("`beta` must be a nonempty finite numeric vector.", call. = FALSE)
  }
  if (!is.numeric(se) || !is.atomic(se) || is.object(se) ||
      !is.null(dim(se)) || length(se) != length(beta) || anyNA(se) ||
      any(!is.finite(se)) || any(se <= 0)) {
    stop("`se` must be a positive finite numeric vector matching `beta`.",
         call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n != floor(n) || n < 2) {
    stop("`n` must be an integer of at least two.", call. = FALSE)
  }
  if (!is.numeric(residual_df) || length(residual_df) != 1L ||
      is.na(residual_df) || !is.finite(residual_df) || residual_df <= 0 ||
      residual_df > n - 1) {
    stop("`residual_df` must be positive and no greater than `n - 1`.",
         call. = FALSE)
  }
  if (!is.null(response_var) &&
      (!is.numeric(response_var) || length(response_var) != 1L ||
       is.na(response_var) || !is.finite(response_var) ||
       response_var <= 0)) {
    stop("`response_var` must be NULL or a positive finite number.",
         call. = FALSE)
  }

  if (scale == "auto") {
    scale <- if (is.null(response_var)) "standardized" else "original"
  }
  if (scale == "original" && is.null(response_var)) {
    stop("`response_var` is required for `scale = \"original\"`.",
         call. = FALSE)
  }
  if (scale == "standardized" && !is.null(response_var)) {
    stop("`response_var` must be NULL for `scale = \"standardized\"`.",
         call. = FALSE)
  }

  p <- length(beta)
  beta_names <- names(beta)
  se_names <- names(se)
  validated_LD <- .validate_gwas_ld(LD, p, beta_names, se_names)
  LD <- validated_LD$LD
  list_LD <- validated_LD$list_input
  predictor_names <- validated_LD$predictor_names
  block_sizes <- validated_LD$block_sizes
  block_names <- validated_LD$block_names
  beta <- as.numeric(beta)
  se <- as.numeric(se)
  z <- beta / se
  if (any(!is.finite(z)) || any(!is.finite(z^2))) {
    stop("The `beta / se` z-scores are too large to represent safely.",
         call. = FALSE)
  }
  adjustment <- (n - 1) / (z^2 + residual_df)
  block_ends <- cumsum(block_sizes)
  block_starts <- c(1L, block_ends[-length(block_ends)] + 1L)
  block_indices <- Map(seq.int, block_starts, block_ends)

  if (scale == "standardized") {
    scale_block <- function(block, indices) (n - 1) * block
    Xty <- sqrt(n - 1) * sqrt(adjustment) * z
    yty <- n - 1
    reference_response_var <- 1
  } else {
    XtX_diagonal <- response_var * adjustment / se^2
    root_diagonal <- sqrt(XtX_diagonal)
    scale_block <- function(block, indices) {
      .scale_gram_block(block, root_diagonal[indices])
    }
    Xty <- XtX_diagonal * beta
    yty <- (n - 1) * response_var
    reference_response_var <- response_var
  }
  if (anyNA(Xty) || any(!is.finite(Xty)) || !is.finite(yty)) {
    stop("The reconstructed sufficient statistics are not finite.",
         call. = FALSE)
  }

  names(Xty) <- predictor_names
  X_means <- stats::setNames(numeric(p), predictor_names)

  common <- list(
    n = as.integer(n), Xty = Xty, yty = as.numeric(yty),
    X_means = X_means, y_mean = 0,
    reference_response_var = as.numeric(reference_response_var)
  )
  if (output == "sufficient") {
    XtX <- if (list_LD) {
      Map(scale_block, LD, block_indices)
    } else {
      scale_block(LD, block_indices[[1L]])
    }
    finite_XtX <- if (list_LD) {
      all(vapply(XtX, .finite_gwas_matrix, logical(1)))
    } else {
      .finite_gwas_matrix(XtX)
    }
    if (!finite_XtX) {
      stop("The reconstructed sufficient statistics are not finite.",
           call. = FALSE)
    }
    if (list_LD) {
      XtX <- Map(function(block, indices) {
        names <- predictor_names[indices]
        dimnames(block) <- list(names, names)
        block
      }, XtX, block_indices)
      names(XtX) <- block_names
    } else {
      dimnames(XtX) <- list(predictor_names, predictor_names)
    }
    return(c(list(XtX = XtX), common))
  }

  decompositions <- Map(function(block, indices, block_name) {
    scaled_block <- scale_block(block, indices)
    if (!.finite_gwas_matrix(scaled_block)) {
      stop("The reconstructed sufficient statistics are not finite.",
           call. = FALSE)
    }
    .decompose_gwas_gram(
      scaled_block, predictor_names[indices],
      if (list_LD) sprintf("eigen block `%s`", block_name) else "`XtX`"
    )
  }, if (list_LD) LD else list(LD), block_indices, block_names)

  if (list_LD) {
    eigenvectors <- lapply(decompositions, `[[`, "vectors")
    eigenvalues <- lapply(decompositions, `[[`, "values")
    eigenvalue_tolerance <- vapply(
      decompositions, `[[`, numeric(1), "tolerance"
    )
    names(eigenvectors) <- names(eigenvalues) <-
      names(eigenvalue_tolerance) <- block_names
  } else {
    eigenvectors <- decompositions[[1L]]$vectors
    eigenvalues <- decompositions[[1L]]$values
    eigenvalue_tolerance <- decompositions[[1L]]$tolerance
  }
  c(
    list(
      XtX_eigenvectors_raw = eigenvectors,
      XtX_eigenvalues_raw = eigenvalues,
      XtX_eigenvalue_tolerance = eigenvalue_tolerance
    ),
    common
  )
}

.validate_gwas_ld <- function(LD, p, beta_names, se_names) {
  list_input <- is.list(LD) && !is.matrix(LD)
  if (list_input) {
    if (length(LD) < 1L) {
      stop("List `LD` must contain at least one correlation matrix.",
           call. = FALSE)
    }
    block_names <- names(LD)
    if (is.null(block_names)) block_names <- paste0("block", seq_along(LD))
    missing_names <- is.na(block_names) | block_names == ""
    block_names[missing_names] <- paste0("block", which(missing_names))
    block_names <- make.unique(block_names)
    blocks <- Map(function(block, name) {
      .validate_gwas_ld_block(
        block, sprintf("`LD[[\"%s\"]]`", name), strict_names = TRUE
      )
    }, LD, block_names)
    names(blocks) <- block_names
  } else {
    block_names <- "LD"
    blocks <- list(LD = .validate_gwas_ld_block(
      LD, "`LD`", strict_names = TRUE
    ))
  }

  block_sizes <- vapply(blocks, ncol, integer(1))
  if (sum(block_sizes) != p) {
    stop("The total dimension of `LD` must match `beta` and `se`.",
         call. = FALSE)
  }
  supplied_block_names <- vapply(
    blocks, function(block) !is.null(colnames(block)), logical(1)
  )
  if (list_input && any(supplied_block_names) && !all(supplied_block_names)) {
    stop("List `LD` must name predictors in every matrix or in none of them.",
         call. = FALSE)
  }
  LD_names <- if (all(supplied_block_names)) {
    unlist(lapply(blocks, colnames), use.names = FALSE)
  } else {
    NULL
  }
  supplied_names <- Filter(Negate(is.null), list(
    beta = beta_names, se = se_names, LD = LD_names
  ))
  if (length(supplied_names) > 1L &&
      !all(vapply(supplied_names[-1L], identical, logical(1),
                  supplied_names[[1L]]))) {
    stop("Predictor names in `beta`, `se`, and `LD` must match in order.",
         call. = FALSE)
  }
  predictor_names <- if (length(supplied_names)) {
    supplied_names[[1L]]
  } else {
    paste0("x", seq_len(p))
  }
  if (length(predictor_names) != p || anyNA(predictor_names) ||
      any(predictor_names == "") || anyDuplicated(predictor_names)) {
    stop("Predictor names must be nonempty and unique.", call. = FALSE)
  }

  block_ends <- cumsum(block_sizes)
  block_starts <- c(1L, block_ends[-length(block_ends)] + 1L)
  blocks <- Map(function(block, start, end) {
    names <- predictor_names[seq.int(start, end)]
    dimnames(block) <- list(names, names)
    block
  }, blocks, block_starts, block_ends)
  names(blocks) <- block_names
  list(
    LD = if (list_input) blocks else blocks[[1L]],
    list_input = list_input,
    predictor_names = predictor_names,
    block_sizes = unname(block_sizes), block_names = block_names
  )
}

.validate_gwas_ld_block <- function(matrix, label, strict_names) {
  matrix <- .validate_gram_block(matrix, label, strict_names)
  diagonal <- if (inherits(matrix, "sparseMatrix")) {
    Matrix::diag(matrix)
  } else {
    diag(matrix)
  }
  if (any(abs(diagonal - 1) > 1e-6)) {
    stop(sprintf("The diagonal of %s must equal one.", label), call. = FALSE)
  }
  maximum <- if (inherits(matrix, "sparseMatrix")) {
    if (length(matrix@x)) max(abs(matrix@x)) else 0
  } else {
    max(abs(matrix))
  }
  if (maximum > 1 + 1e-6) {
    stop(sprintf("Every correlation in %s must be between -1 and 1.", label),
         call. = FALSE)
  }
  matrix
}

.finite_gwas_matrix <- function(matrix) {
  if (inherits(matrix, "sparseMatrix")) {
    !anyNA(matrix@x) && all(is.finite(matrix@x))
  } else {
    !anyNA(matrix) && all(is.finite(matrix))
  }
}

.decompose_gwas_gram <- function(matrix, predictor_names, label) {
  decomposition <- eigen(as.matrix(matrix), symmetric = TRUE)
  rownames(decomposition$vectors) <- predictor_names
  colnames(decomposition$vectors) <- paste0(
    "eigen", seq_along(decomposition$values)
  )
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(decomposition$values)))
  negative <- decomposition$values < 0
  if (any(negative)) {
    materially_negative <- decomposition$values < -tolerance
    warning(
      sprintf(
        paste0(
          "The reconstructed %s has %d negative eigenvalue(s) ",
          "(minimum %.6g); %d are below -%.6g. All eigenpairs were ",
          "returned unchanged."
        ),
        label, sum(negative), min(decomposition$values),
        sum(materially_negative), tolerance
      ),
      call. = FALSE
    )
  }
  list(
    vectors = decomposition$vectors,
    values = decomposition$values,
    tolerance = tolerance
  )
}
