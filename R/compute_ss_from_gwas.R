#' Compute regression sufficient statistics from GWAS results
#'
#' Reconstructs centered linear-regression cross-products from marginal GWAS
#' effect estimates, their standard errors, and a signed LD correlation matrix
#' using the finite-sample transformation employed by SuSiE-RSS.
#'
#' @param beta A finite numeric vector of marginal ordinary least-squares
#'   effect estimates.
#' @param se A positive finite numeric vector of standard errors corresponding
#'   to `beta`.
#' @param LD A finite, symmetric LD correlation matrix whose rows and columns
#'   have the same order and effect-allele orientation as `beta`. Dense base-R
#'   matrices and compressed sparse `Matrix` objects of class `dgCMatrix` or
#'   `dsCMatrix` are supported. Supply signed correlations, not squared
#'   correlations.
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
#'   raw eigendecomposition for preparing input to [blm_ss_eigen()].
#'
#' @return For `output = "sufficient"`, a list containing `n`, `XtX`, `Xty`,
#'   `yty`, zero `X_means` and `y_mean`, and `reference_response_var`. Pass
#'   `reference_response_var` to [blm_ss()] only when at least one `ETA` block
#'   specifies `expected_pve`. For `output = "eigen"`, `XtX` is replaced by
#'   `XtX_eigenvectors_raw`, `XtX_eigenvalues_raw`, and
#'   `XtX_eigenvalue_tolerance`. Raw eigenpairs are returned without filtering
#'   or modification and therefore are not necessarily valid inputs to
#'   [blm_ss_eigen()], which requires strictly positive eigenvalues.
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
#'
#'   The function performs only basic input validation. It does not diagnose
#'   allele flips or broader incompatibility between the GWAS statistics and
#'   LD. Eigen output uses a complete dense eigendecomposition, requiring
#'   \eqn{O(p^2)} memory and \eqn{O(p^3)} time. Every eigenpair is returned,
#'   including zero and negative eigenvalues. A warning is produced if any
#'   negative eigenvalues are present. Users must choose the eigenpairs passed
#'   to [blm_ss_eigen()] and calculate the corresponding `XtX_prop_var`.
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

  sparse_LD <- inherits(LD, "sparseMatrix")
  supported_sparse <- inherits(LD, c("dgCMatrix", "dsCMatrix"))
  valid_dense <- is.matrix(LD) && is.numeric(LD)
  finite_LD <- if (sparse_LD) {
    supported_sparse && !anyNA(LD@x) && all(is.finite(LD@x))
  } else {
    valid_dense && !anyNA(LD) && all(is.finite(LD))
  }
  p <- length(beta)
  if ((!valid_dense && !supported_sparse) || nrow(LD) != p ||
      ncol(LD) != p || !finite_LD) {
    stop("`LD` must be a finite numeric square matrix matching `beta`.",
         call. = FALSE)
  }

  maximum_LD <- if (sparse_LD && length(LD@x) == 0L) {
    0
  } else if (sparse_LD) {
    max(abs(LD@x))
  } else {
    max(abs(LD))
  }
  symmetry_tolerance <- sqrt(.Machine$double.eps) * max(1, maximum_LD)
  asymmetry <- if (inherits(LD, "dsCMatrix")) {
    0
  } else {
    difference <- if (sparse_LD) LD - Matrix::t(LD) else LD - t(LD)
    if (sparse_LD && length(difference@x) == 0L) {
      0
    } else {
      max(abs(difference))
    }
  }
  if (asymmetry > symmetry_tolerance) {
    stop("`LD` must be symmetric.", call. = FALSE)
  }
  if (!inherits(LD, "dsCMatrix")) {
    LD <- if (sparse_LD) {
      (LD + Matrix::t(LD)) / 2
    } else {
      (LD + t(LD)) / 2
    }
  }

  correlation_tolerance <- 1e-6
  LD_diagonal <- if (sparse_LD) Matrix::diag(LD) else diag(LD)
  if (any(abs(LD_diagonal - 1) > correlation_tolerance)) {
    stop("The diagonal of `LD` must equal one.", call. = FALSE)
  }
  if (maximum_LD > 1 + correlation_tolerance) {
    stop("Entries of `LD` must be correlations between -1 and 1.",
         call. = FALSE)
  }

  beta_names <- names(beta)
  se_names <- names(se)
  LD_row_names <- rownames(LD)
  LD_column_names <- colnames(LD)
  if (xor(is.null(LD_row_names), is.null(LD_column_names)) ||
      (!is.null(LD_row_names) && !identical(LD_row_names, LD_column_names))) {
    stop("The row and column names of `LD` must match.", call. = FALSE)
  }
  supplied_names <- Filter(Negate(is.null), list(
    beta = beta_names, se = se_names, LD = LD_column_names
  ))
  if (length(supplied_names) > 1L &&
      !all(vapply(supplied_names[-1L], identical, logical(1),
                  supplied_names[[1L]]))) {
    stop("Names of `beta`, `se`, and `LD` must match when supplied.",
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

  beta <- as.numeric(beta)
  se <- as.numeric(se)
  z <- beta / se
  if (any(!is.finite(z)) || any(!is.finite(z^2))) {
    stop("The `beta / se` z-scores are too large to represent safely.",
         call. = FALSE)
  }
  adjustment <- (n - 1) / (z^2 + residual_df)

  if (scale == "standardized") {
    XtX <- (n - 1) * LD
    Xty <- sqrt(n - 1) * sqrt(adjustment) * z
    yty <- n - 1
    reference_response_var <- 1
  } else {
    XtX_diagonal <- response_var * adjustment / se^2
    root_diagonal <- sqrt(XtX_diagonal)
    if (sparse_LD) {
      scaling <- Matrix::Diagonal(x = root_diagonal)
      XtX <- scaling %*% LD %*% scaling
    } else {
      XtX <- LD * tcrossprod(root_diagonal)
    }
    Xty <- XtX_diagonal * beta
    yty <- (n - 1) * response_var
    reference_response_var <- response_var
  }
  finite_XtX <- if (sparse_LD) {
    !anyNA(XtX@x) && all(is.finite(XtX@x))
  } else {
    !anyNA(XtX) && all(is.finite(XtX))
  }
  if (!finite_XtX || anyNA(Xty) || any(!is.finite(Xty)) ||
      !is.finite(yty)) {
    stop("The reconstructed sufficient statistics are not finite.",
         call. = FALSE)
  }

  dimnames(XtX) <- list(predictor_names, predictor_names)
  names(Xty) <- predictor_names
  X_means <- stats::setNames(numeric(p), predictor_names)

  common <- list(
    n = as.integer(n), Xty = Xty, yty = as.numeric(yty),
    X_means = X_means, y_mean = 0,
    reference_response_var = as.numeric(reference_response_var)
  )
  if (output == "sufficient") {
    return(c(list(XtX = XtX), common))
  }

  decomposition <- eigen(as.matrix(XtX), symmetric = TRUE)
  rownames(decomposition$vectors) <- predictor_names
  colnames(decomposition$vectors) <- paste0(
    "eigen", seq_along(decomposition$values)
  )
  eigenvalue_tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(decomposition$values)))
  negative <- decomposition$values < 0
  if (any(negative)) {
    materially_negative <- decomposition$values < -eigenvalue_tolerance
    warning(
      sprintf(
        paste0(
          "The reconstructed `XtX` has %d negative eigenvalue(s) ",
          "(minimum %.6g); %d are below -%.6g. All eigenpairs were ",
          "returned unchanged."
        ),
        sum(negative), min(decomposition$values), sum(materially_negative),
        eigenvalue_tolerance
      ),
      call. = FALSE
    )
  }
  c(
    list(
      XtX_eigenvectors_raw = decomposition$vectors,
      XtX_eigenvalues_raw = decomposition$values,
      XtX_eigenvalue_tolerance = eigenvalue_tolerance
    ),
    common
  )
}
