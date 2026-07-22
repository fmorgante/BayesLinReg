library(BayesLinReg)

set.seed(501)
n <- 70
X <- matrix(rnorm(n * 5), nrow = n)
colnames(X) <- paste0("x", seq_len(ncol(X)))
y <- drop(1 + X %*% c(1.5, -1, 0, 0.5, 0) + rnorm(n))
XtX <- crossprod(X)
Xty <- crossprod(X, y)
yty <- sum(y^2)

# Sufficient-statistics fits reproduce raw-data fits for every prior and engine.
for (sampler_version in c("Rcpp", "R")) {
  for (model in c(
    "Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab"
  )) {
    raw_fit <- blm(
      y,
      ETA = list(X = X, model = model),
      residual_shape = 2,
      residual_scale = 1,
      iterations = 300,
      burnin = 100,
      seed = 502,
      version = sampler_version
    )
    ss_fit <- blm_ss(
      n = n,
      XtX = XtX,
      Xty = Xty,
      yty = yty,
      X_means = colMeans(X),
      y_mean = mean(y),
      ETA = list(model = model),
      residual_shape = 2,
      residual_scale = 1,
      iterations = 300,
      burnin = 100,
      seed = 502,
      version = sampler_version
    )
    stopifnot(
      isTRUE(all.equal(raw_fit, ss_fit, tolerance = 1e-8))
    )
  }
}

# Multiple blocks select and reorder columns using integer or character indices.
raw_eta <- list(
  fixed = list(X = X[, c("x1", "x3")], model = "Normal"),
  selection = list(X = X[, "x2", drop = FALSE], model = "SpikeSlab"),
  shrinkage = list(X = X[, c("x4", "x5")], model = "GlobalLocal")
)
ss_eta <- list(
  fixed = list(indices = c("x1", "x3"), model = "Normal"),
  selection = list(indices = 2, model = "SpikeSlab"),
  shrinkage = list(indices = c(4, 5), model = "GlobalLocal")
)
raw_multi <- blm(
  y, ETA = raw_eta, residual_var = 1,
  iterations = 250, burnin = 100, seed = 503, version = "Rcpp"
)
ss_multi <- blm_ss(
  n, XtX, Xty, ETA = ss_eta, yty = yty,
  X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
  iterations = 250, burnin = 100, seed = 503, version = "Rcpp"
)
stopifnot(isTRUE(all.equal(raw_multi, ss_multi, tolerance = 1e-8)))

# Compressed sparse cross-products use the separate RcppEigen path. Centering
# is represented implicitly even when it makes the logical Gram matrix dense.
set.seed(511)
sparse_n <- 80
sparse_p <- 5
sparse_group <- sample(c(seq_len(sparse_p), NA_integer_), sparse_n, TRUE)
sparse_X <- Matrix::sparseMatrix(
  i = which(!is.na(sparse_group)),
  j = sparse_group[!is.na(sparse_group)],
  x = runif(sum(!is.na(sparse_group)), 0.5, 1.5),
  dims = c(sparse_n, sparse_p),
  dimnames = list(NULL, paste0("s", seq_len(sparse_p)))
)
sparse_y <- drop(
  0.75 + as.matrix(sparse_X) %*% c(1, -0.5, 0, 0.25, -0.75) +
    rnorm(sparse_n, sd = 0.5)
)
sparse_XtX <- Matrix::crossprod(sparse_X)
sparse_Xty <- as.numeric(Matrix::crossprod(sparse_X, sparse_y))
sparse_yty <- sum(sparse_y^2)
sparse_args <- list(
  n = sparse_n,
  Xty = sparse_Xty,
  yty = sparse_yty,
  X_means = Matrix::colMeans(sparse_X),
  y_mean = mean(sparse_y),
  ETA = list(
    first = list(indices = c("s3", "s1"), model = "Normal"),
    second = list(
      indices = c("s2", "s4", "s5"), model = "SpikeMultiSlab"
    )
  ),
  residual_shape = 2,
  residual_scale = 1,
  iterations = 180,
  burnin = 60,
  seed = 512,
  version = "Rcpp"
)
sparse_dense_fit <- do.call(
  blm_ss,
  c(list(XtX = as.matrix(sparse_XtX)), sparse_args)
)
sparse_symmetric_fit <- do.call(
  blm_ss,
  c(list(XtX = sparse_XtX), sparse_args)
)
sparse_general_fit <- do.call(
  blm_ss,
  c(list(XtX = methods::as(sparse_XtX, "generalMatrix")), sparse_args)
)
stopifnot(
  inherits(sparse_XtX, "dsCMatrix"),
  isTRUE(all.equal(
    sparse_dense_fit, sparse_symmetric_fit, tolerance = 1e-8
  )),
  isTRUE(all.equal(
    sparse_dense_fit, sparse_general_fit, tolerance = 1e-8
  ))
)

for (sparse_model in c(
  "Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab"
)) {
  model_args <- utils::modifyList(sparse_args, list(
    residual_var = 1,
    residual_shape = NULL,
    residual_scale = NULL,
    iterations = 120,
    burnin = 40,
    seed = 513
  ))
  model_args$ETA <- list(model = sparse_model)
  dense_model_fit <- do.call(
    blm_ss,
    c(list(XtX = as.matrix(sparse_XtX)), model_args)
  )
  sparse_model_fit <- do.call(
    blm_ss,
    c(list(XtX = sparse_XtX), model_args)
  )
  stopifnot(isTRUE(all.equal(
    dense_model_fit, sparse_model_fit, tolerance = 1e-8
  )))
}

# Optional dense PSD validation and summary-only storage also work through the
# sparse entry point.
sparse_checked <- do.call(
  blm_ss,
  c(
    list(XtX = sparse_XtX),
    utils::modifyList(sparse_args, list(
      check_psd = TRUE,
      store_samples = FALSE,
      store_coefficient_cov = FALSE,
      iterations = 100,
      burnin = 40
    ))
  )
)
stopifnot(
  identical(sparse_checked$store_samples, FALSE),
  is.null(sparse_checked$ETA$first$coefficient_cov)
)

sparse_r_error <- try(do.call(
  blm_ss,
  c(
    list(XtX = sparse_XtX),
    utils::modifyList(sparse_args, list(version = "R"))
  )
), silent = TRUE)
stopifnot(
  inherits(sparse_r_error, "try-error"),
  grepl("requires `version = \\\"Rcpp\\\"`", sparse_r_error)
)

asymmetric_sparse <- methods::as(sparse_XtX, "generalMatrix")
asymmetric_sparse[1, 2] <- 0.25
sparse_symmetry_error <- try(do.call(
  blm_ss,
  c(list(XtX = asymmetric_sparse), sparse_args)
), silent = TRUE)
stopifnot(
  inherits(sparse_symmetry_error, "try-error"),
  grepl("must be symmetric", sparse_symmetry_error)
)

# yty is unnecessary when the residual variance is fixed, including when
# means are used to fit an intercept.
raw_fixed <- blm(
  y, ETA = list(X = X, model = "SpikeSlab"), residual_var = 1,
  iterations = 200, burnin = 80, seed = 507, version = "Rcpp"
)
ss_fixed <- blm_ss(
  n, XtX, Xty, ETA = list(model = "SpikeSlab"),
  X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
  iterations = 200, burnin = 80, seed = 507, version = "Rcpp"
)
stopifnot(isTRUE(all.equal(raw_fixed, ss_fixed, tolerance = 1e-8)))

# Without yty, a fixed residual variance supports every prior, including
# SpikeSlab, and summary-only storage remains available.
fixed_spike <- suppressWarnings(blm_ss(
  n, XtX, Xty,
  ETA = list(model = "SpikeSlab"),
  residual_var = 1,
  iterations = 150,
  burnin = 50,
  seed = 504,
  store_samples = FALSE,
  store_coefficient_cov = FALSE
))
stopifnot(
  identical(fixed_spike$residual_var_mean, 1),
  identical(fixed_spike$residual_var_var, 0),
  is.null(fixed_spike$intercept_mean),
  is.null(fixed_spike$ETA$ETA1$coefficient_cov),
  length(fixed_spike$ETA$ETA1$coefficient_var) == ncol(X)
)

# Omitting means produces the no-intercept and standardization warnings.
warnings <- character()
no_intercept_fit <- withCallingHandlers(
  blm_ss(
    n, XtX, Xty, ETA = list(model = "Normal"), residual_var = 1,
    iterations = 100, burnin = 40, seed = 505
  ),
  warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  length(warnings) == 2L,
  any(grepl("without an intercept", warnings, fixed = TRUE)),
  any(grepl("should be centered or standardized", warnings, fixed = TRUE)),
  is.null(no_intercept_fit$intercept_samples),
  !"intercept" %in% names(assess_convergence(
    no_intercept_fit, plot = FALSE
  )$rhat)
)

# yty permits learning the residual variance in a no-intercept model.
no_intercept_r <- suppressWarnings(blm_ss(
  n, XtX, Xty, ETA = list(model = "Normal"), yty = yty,
  residual_shape = 2, residual_scale = 1,
  iterations = 150, burnin = 50, seed = 508, version = "R"
))
no_intercept_rcpp <- suppressWarnings(blm_ss(
  n, XtX, Xty, ETA = list(model = "Normal"), yty = yty,
  residual_shape = 2, residual_scale = 1,
  iterations = 150, burnin = 50, seed = 508, version = "Rcpp"
))
stopifnot(
  is.null(no_intercept_r$intercept_mean),
  all(no_intercept_r$residual_var_samples > 0),
  max(abs(
    no_intercept_r$ETA$ETA1$coefficient_mean -
      no_intercept_rcpp$ETA$ETA1$coefficient_mean
  )) < 0.2,
  abs(
    no_intercept_r$residual_var_mean - no_intercept_rcpp$residual_var_mean
  ) < 0.2
)

# Rank-deficient cross-products do not require a pseudo-design factorization.
rank_X <- cbind(
  x1 = seq_len(40),
  x2 = rep(c(-1, 1), 20),
  x3 = seq_len(40) + rep(c(-1, 1), 20)
)
rank_y <- drop(2 + rank_X %*% c(0.5, -0.25, 0.1) + rnorm(40, sd = 0.2))
rank_fits <- lapply(c("R", "Rcpp"), function(sampler_version) {
  blm_ss(
    nrow(rank_X), crossprod(rank_X), crossprod(rank_X, rank_y),
    ETA = list(model = "Normal"), yty = sum(rank_y^2),
    X_means = colMeans(rank_X), y_mean = mean(rank_y),
    residual_shape = 2, residual_scale = 1,
    iterations = 150, burnin = 50, seed = 509,
    version = sampler_version
  )
})
stopifnot(
  all(vapply(rank_fits, function(fit) {
    all(is.finite(fit$ETA$ETA1$coefficient_mean)) &&
      all(fit$residual_var_samples > 0)
  }, logical(1))),
  max(abs(
    rank_fits[[1]]$ETA$ETA1$coefficient_mean -
      rank_fits[[2]]$ETA$ETA1$coefficient_mean
  )) < 0.2
)

warnings <- character()
invisible(withCallingHandlers(
  blm_ss(
    n, XtX, Xty,
    ETA = list(model = "Normal", standardize = FALSE),
    residual_var = 1, iterations = 50, burnin = 20, seed = 506
  ),
  warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
))
stopifnot(
  length(warnings) == 1L,
  grepl("without an intercept", warnings, fixed = TRUE)
)

# Full PSD and joint-compatibility validation is opt-in.
unchecked_fit <- blm_ss(
  n, XtX, Xty, yty = 0, ETA = list(model = "Normal"),
  X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
  iterations = 50, burnin = 20, seed = 510
)
stopifnot(inherits(unchecked_fit, "blm_fit"))

# Invalid or incomplete sufficient statistics are rejected.
invalid_calls <- list(
  function() blm_ss(
    n, XtX, Xty, ETA = list(model = "Normal"),
    residual_shape = 2, residual_scale = 1
  ),
  function() blm_ss(
    n, XtX, Xty, ETA = list(model = "Normal"),
    X_means = colMeans(X), residual_var = 1
  ),
  function() blm_ss(
    n, XtX, Xty, ETA = list(model = "Normal"),
    y_mean = mean(y), residual_var = 1
  ),
  function() blm_ss(
    n, XtX, Xty[-1], ETA = list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss(
    n, XtX + upper.tri(XtX), Xty,
    ETA = list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss(
    n, XtX, Xty,
    ETA = list(
      first = list(indices = 1:3, model = "Normal"),
      second = list(indices = 3:5, model = "Normal")
    ),
    residual_var = 1
  ),
  function() blm_ss(
    n, XtX, Xty, yty = 0,
    ETA = list(model = "Normal"),
    residual_shape = 2, residual_scale = 1, check_psd = TRUE
  ),
  function() blm_ss(
    10, matrix(c(1, 2, 2, 1), 2), c(0, 0),
    ETA = list(model = "Normal", standardize = FALSE),
    residual_var = 1, check_psd = TRUE
  ),
  function() blm_ss(
    n, XtX, Xty, ETA = list(model = "Normal"), residual_var = 1,
    check_psd = NA
  )
)
stopifnot(all(vapply(
  invalid_calls,
  function(call) inherits(suppressWarnings(try(call(), silent = TRUE)),
                           "try-error"),
  logical(1)
)))
