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
