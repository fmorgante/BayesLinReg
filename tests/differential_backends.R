library(BayesLinReg)

# Randomized, reproducible differential checks protect the reference R and
# optimized Rcpp implementations against behavioral drift across refactors.
models <- c("Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab", "Fixed")
exact_rng_models <- c("Normal", "SpikeSlab", "Fixed")

for (case in seq_len(3L)) {
  set.seed(900 + case)
  n <- 25L
  p <- 4L
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- drop(X %*% c(0.8, -0.4, 0, 0) + rnorm(n))

  for (model in models) {
    eta <- list(X = X, model = model)
    set.seed(1200 + case)
    reference <- blm(
      y, eta, residual_var = 1, iterations = 120, burnin = 40,
      version = "R"
    )
    set.seed(1200 + case)
    optimized <- blm(
      y, eta, residual_var = 1, iterations = 120, burnin = 40,
      version = "Rcpp"
    )

    difference <- max(abs(
      reference$ETA[[1L]]$coefficient_mean -
        optimized$ETA[[1L]]$coefficient_mean
    ))
    tolerance <- if (model %in% exact_rng_models) 1e-10 else 0.5
    stopifnot(
      is.finite(difference),
      difference < tolerance,
      identical(reference$store_samples, FALSE),
      identical(reference$store_coefficient_cov, FALSE),
      identical(optimized$store_samples, FALSE),
      identical(optimized$store_coefficient_cov, FALSE)
    )
  }
}

# Random exact sufficient statistics agree across dense, general sparse,
# block-streaming, full-eigen, and LD-native representations. These cases also
# exercise posterior PVE and within-block coefficient covariance.
summary_models <- c("Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab")
for (case in seq_len(4L)) {
  set.seed(1700 + case)
  n <- 40L
  p <- 6L
  X <- scale(matrix(rnorm(n * p), nrow = n), center = TRUE, scale = TRUE)
  y <- drop(scale(rnorm(n), center = TRUE, scale = TRUE))
  predictor_names <- paste0("v", seq_len(p))
  colnames(X) <- predictor_names
  XtX <- crossprod(X)
  Xty <- drop(crossprod(X, y))
  names(Xty) <- predictor_names
  yty <- sum(y^2)
  decomposition <- eigen(XtX, symmetric = TRUE)
  model <- summary_models[[case]]
  ETA <- list(model = model)
  arguments <- list(
    n = n, Xty = Xty, ETA = ETA, yty = yty, residual_var = 1,
    iterations = 70L, burnin = 20L, store_samples = TRUE,
    store_coefficient_cov = TRUE, compute_pve = TRUE
  )

  set.seed(1800 + case)
  dense <- suppressWarnings(do.call(blm_ss, c(list(XtX = XtX), arguments)))
  set.seed(1800 + case)
  sparse <- suppressWarnings(do.call(
    blm_ss,
    c(list(XtX = Matrix::Matrix(XtX, sparse = TRUE)), arguments)
  ))
  set.seed(1800 + case)
  blocked <- suppressWarnings(do.call(
    blm_ss, c(list(XtX = list(all = XtX)), arguments)
  ))
  set.seed(1800 + case)
  eigen_fit <- suppressWarnings(do.call(
    blm_ss_eigen,
    c(list(
      XtX_eigenvectors = decomposition$vectors,
      XtX_eigenvalues = decomposition$values,
      XtX_prop_var = 1
    ), arguments)
  ))

  marginal_beta <- Xty / diag(XtX)
  marginal_se <- sqrt(
    (yty - marginal_beta^2 * diag(XtX)) /
      ((n - 2) * diag(XtX))
  )
  variants <- data.frame(
    CHR = 1, ID = predictor_names, POS = seq_len(p),
    A1 = rep(c("A", "C"), length.out = p),
    A0 = rep(c("C", "A"), length.out = p)
  )
  gwas <- transform(
    variants, N = n, BETA = marginal_beta, SE = marginal_se
  )
  ld <- as_blm_ld(stats::cor(X), variants)
  set.seed(1800 + case)
  ld_fit <- blm_gwas(
    gwas, ld, ETA, residual_var = 1,
    iterations = 70L, burnin = 20L, store_samples = TRUE,
    store_coefficient_cov = TRUE, compute_pve = TRUE
  )

  fits <- list(
    sparse = sparse, blocked = blocked, eigen = eigen_fit, ld = ld_fit
  )
  comparison_messages <- lapply(fits, function(fit) {
    comparisons <- list(
      ETA = all.equal(dense$ETA, fit$ETA, tolerance = 1e-8),
      total_pve = all.equal(
        dense$total_pve_samples, fit$total_pve_samples, tolerance = 1e-8
      ),
      cross_block_pve = all.equal(
        dense$cross_block_pve_samples, fit$cross_block_pve_samples,
        tolerance = 1e-8
      )
    )
    failed <- !vapply(comparisons, isTRUE, logical(1))
    if (!any(failed)) return(character())
    paste0(
      names(comparisons)[failed], ": ",
      vapply(comparisons[failed], paste, character(1), collapse = " ")
    )
  })
  representations_agree <- !lengths(comparison_messages)
  if (!all(representations_agree)) {
    stop(sprintf(
      "Random differential case %d (%s) disagreed: %s.",
      case, model,
      paste(vapply(
        names(fits)[!representations_agree], function(name) {
          detail <- paste(comparison_messages[[name]], collapse = "; ")
          paste0(name, " [", detail, "]")
        }, character(1)
      ), collapse = "; ")
    ))
  }
}
