library(BayesLinReg)

blm_public <- BayesLinReg::blm
blm <- function(..., store_samples = TRUE, store_coefficient_cov = TRUE) {
  blm_public(..., store_samples = store_samples,
             store_coefficient_cov = store_coefficient_cov)
}
blm_ss_public <- BayesLinReg::blm_ss
blm_ss <- function(..., store_samples = TRUE, store_coefficient_cov = TRUE) {
  blm_ss_public(..., store_samples = store_samples,
                store_coefficient_cov = store_coefficient_cov)
}
blm_ss_eigen_public <- BayesLinReg::blm_ss_eigen
blm_ss_eigen <- function(..., store_samples = TRUE,
                         store_coefficient_cov = TRUE) {
  blm_ss_eigen_public(..., store_samples = store_samples,
                      store_coefficient_cov = store_coefficient_cov)
}

set.seed(801)
n <- 50
X <- matrix(rnorm(n * 6), nrow = n)
colnames(X) <- paste0("x", seq_len(ncol(X)))
y <- drop(1 + X %*% c(1.5, -0.5, 0, -1, 0, 0.25) + rnorm(n))
eta <- list(
  first = list(X = X[, 1:2, drop = FALSE], model = "Normal"),
  second = list(X = X[, 3:6, drop = FALSE], model = "SpikeSlab")
)

# PVE is evaluated at retained draws and agrees with direct fitted-value
# calculations and online summaries in both sampling implementations.
for (sampler_version in c("Rcpp", "R")) {
  set.seed(802)
  baseline <- blm(
    y, eta, residual_shape = 2, residual_scale = 1,
    iterations = 90, burnin = 30, thin = 2,
    version = sampler_version
  )
  for (pve_type in c("standalone", "allocated")) {
    set.seed(802)
    stored <- blm(
      y, eta, residual_shape = 2, residual_scale = 1,
      iterations = 90, burnin = 30, thin = 2,
      version = sampler_version, compute_pve = TRUE,
      pve_type = pve_type
    )
    set.seed(802)
    online <- blm(
      y, eta, residual_shape = 2, residual_scale = 1,
      iterations = 90, burnin = 30, thin = 2,
      version = sampler_version, compute_pve = TRUE,
      pve_type = pve_type, store_samples = FALSE
    )
    stopifnot(
      identical(stored$pve_type, pve_type),
      length(stored$total_pve_samples) == 30L,
      isTRUE(all.equal(
        baseline$ETA$first$coefficient_samples,
        stored$ETA$first$coefficient_samples
      )),
      isTRUE(all.equal(stored$total_pve_mean, online$total_pve_mean)),
      isTRUE(all.equal(stored$total_pve_var, online$total_pve_var)),
      isTRUE(all.equal(
        stored$cross_block_pve_mean, online$cross_block_pve_mean
      )),
      isTRUE(all.equal(
        stored$ETA$first$pve_mean, online$ETA$first$pve_mean
      )),
      isTRUE(all.equal(
        stored$ETA$second$pve_var, online$ETA$second$pve_var
      ))
    )

    centered_X <- sweep(X, 2L, colMeans(X), FUN = "-")
    first_fitted <- centered_X[, 1:2, drop = FALSE] %*%
      t(stored$ETA$first$coefficient_samples)
    second_fitted <- centered_X[, 3:6, drop = FALSE] %*%
      t(stored$ETA$second$coefficient_samples)
    total_fitted <- first_fitted + second_fitted
    total_variance <- colSums(total_fitted^2) / (n - 1)
    denominator <- total_variance + stored$residual_var_samples
    direct_total <- total_variance / denominator
    direct_standalone <- rbind(
      colSums(first_fitted^2) / (n - 1) / denominator,
      colSums(second_fitted^2) / (n - 1) / denominator
    )
    direct_allocated <- rbind(
      colSums(first_fitted * total_fitted) / (n - 1) / denominator,
      colSums(second_fitted * total_fitted) / (n - 1) / denominator
    )
    direct_cross <- direct_total - colSums(direct_standalone)
    reported_blocks <- rbind(
      stored$ETA$first$pve_samples,
      stored$ETA$second$pve_samples
    )
    accounted_total <- if (pve_type == "standalone") {
      colSums(reported_blocks) + stored$cross_block_pve_samples
    } else {
      colSums(reported_blocks)
    }
    stopifnot(
      max(abs(stored$total_pve_samples - direct_total)) < 1e-12,
      max(abs(stored$cross_block_pve_samples - direct_cross)) < 1e-12,
      max(abs(reported_blocks - if (pve_type == "standalone") {
        direct_standalone
      } else {
        direct_allocated
      })) < 1e-12,
      max(abs(accounted_total - stored$total_pve_samples)) < 1e-12
    )
  }
  stopifnot(
    is.null(baseline$total_pve_mean),
    is.null(baseline$ETA$first$pve_mean)
  )
}

pve_diagnostics <- assess_convergence(stored, plot = FALSE)
stopifnot(all(c(
  "pve_first", "pve_second", "total_pve", "cross_block_pve"
) %in% names(pve_diagnostics$effective_sample_size)))

# Dense, sparse, and exact eigen sufficient statistics return the same draws
# and PVE summaries.
XtX <- crossprod(X)
Xty <- crossprod(X, y)
yty <- sum(y^2)
ss_eta <- list(
  first = list(indices = 1:2, model = "Normal"),
  second = list(indices = 3:6, model = "SpikeSlab")
)
ss_arguments <- list(
  n = n, Xty = Xty, ETA = ss_eta, yty = yty,
  X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
  iterations = 70, burnin = 30, compute_pve = TRUE
)
set.seed(803)
dense <- do.call(blm_ss, c(list(XtX = XtX), ss_arguments))
set.seed(803)
sparse <- do.call(
  blm_ss,
  c(list(XtX = Matrix::Matrix(XtX, sparse = TRUE)), ss_arguments)
)
sparse_without_storage <- sparse
sparse_without_storage$XtX_representation <- NULL
sparse_without_storage$XtX_storage <- NULL
centered_X <- sweep(X, 2L, colMeans(X), FUN = "-")
decomposition <- eigen(crossprod(centered_X), symmetric = TRUE)
keep <- decomposition$values >
  sqrt(.Machine$double.eps) * max(decomposition$values)
set.seed(803)
eigen_fit <- blm_ss_eigen(
  n, decomposition$vectors[, keep, drop = FALSE],
  decomposition$values[keep], 1, Xty, ETA = ss_eta, yty = yty,
  X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
  iterations = 70, burnin = 30, compute_pve = TRUE
)
stopifnot(
  isTRUE(all.equal(dense, sparse_without_storage, tolerance = 1e-8)),
  max(abs(dense$total_pve_samples - eigen_fit$total_pve_samples)) < 1e-8,
  max(abs(
    dense$ETA$first$pve_samples - eigen_fit$ETA$first$pve_samples
  )) < 1e-8,
  max(abs(
    dense$cross_block_pve_samples - eigen_fit$cross_block_pve_samples
  )) < 1e-8
)

# Invalid PVE controls fail before sampling.
invalid_pve <- list(
  function() blm(y, eta, residual_var = 1, compute_pve = NA),
  function() blm(y, eta, residual_var = 1, pve_type = "unknown")
)
stopifnot(all(vapply(
  invalid_pve,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))
