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
    set.seed(502)
    raw_fit <- blm(
      y,
      ETA = list(X = X, model = model),
      residual_shape = 2,
      residual_scale = 1,
      iterations = 300,
      burnin = 100,

      version = sampler_version
    )
    set.seed(502)
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
set.seed(503)
raw_multi <- blm(
  y, ETA = raw_eta, residual_var = 1,
  iterations = 250, burnin = 100, version = "Rcpp"
)
set.seed(503)
ss_multi <- blm_ss(
  n, XtX, Xty, ETA = ss_eta, yty = yty,
  X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
  iterations = 250, burnin = 100, version = "Rcpp"
)
stopifnot(isTRUE(all.equal(raw_multi, ss_multi, tolerance = 1e-8)))

# A list of Gram matrices represents exact block-diagonal cross-products.
# Gram blocks and prior blocks are independent: each prior block below spans
# both Gram blocks while also splitting predictors within each one.
gram_a <- matrix(c(
  12, 2,
  2, 10
), 2L, dimnames = list(c("g1", "g2"), c("g1", "g2")))
gram_b <- matrix(c(
  9, -1,
  -1, 8
), 2L, dimnames = list(c("g3", "g4"), c("g3", "g4")))
block_gram <- matrix(0, 4L, 4L)
block_gram[1:2, 1:2] <- gram_a
block_gram[3:4, 3:4] <- gram_b
dimnames(block_gram) <- list(paste0("g", 1:4), paste0("g", 1:4))
block_Xty <- stats::setNames(c(1, 2, -1, 0.5), paste0("g", 1:4))
block_means <- stats::setNames(c(0.1, -0.1, 0.05, -0.05), paste0("g", 1:4))
crossing_eta <- list(
  normal = list(indices = c("g1", "g3"), model = "Normal"),
  selection = list(indices = c("g2", "g4"), model = "SpikeMultiSlab")
)
block_arguments <- list(
  n = 20L, Xty = block_Xty, ETA = crossing_eta, yty = 20,
  X_means = block_means, y_mean = 0.02,
  residual_shape = 2, residual_scale = 1,
  iterations = 120, burnin = 40,
  compute_pve = TRUE, check_psd = TRUE
)
set.seed(512)
single_gram_fit <- do.call(
  blm_ss, c(list(XtX = block_gram), block_arguments)
)
set.seed(512)
list_gram_fit <- do.call(
  blm_ss,
  c(list(XtX = list(region_a = gram_a, region_b = gram_b)), block_arguments)
)
stopifnot(
  isTRUE(all.equal(single_gram_fit$ETA, list_gram_fit$ETA,
                   tolerance = 1e-10)),
  isTRUE(all.equal(single_gram_fit$intercept_mean,
                   list_gram_fit$intercept_mean, tolerance = 1e-10)),
  isTRUE(all.equal(single_gram_fit$residual_var_mean,
                   list_gram_fit$residual_var_mean, tolerance = 1e-10)),
  isTRUE(all.equal(single_gram_fit$total_pve_mean,
                   list_gram_fit$total_pve_mean, tolerance = 1e-10)),
  identical(list_gram_fit$XtX_representation, "block_diagonal"),
  identical(list_gram_fit$XtX_number_of_blocks, 2L),
  identical(unname(list_gram_fit$XtX_block_sizes), c(2L, 2L)),
  identical(list_gram_fit$XtX_cross_block_assumption, "zero"),
  identical(list_gram_fit$XtX_storage$representation, c("dense", "dense"))
)

# Symmetric sparse blocks can stream one lower triangle or be expanded for
# speed. Both representations produce the same chain up to numerical drift.
set.seed(513)
symmetric_dense <- lapply(c("left", "right"), function(name) {
  values <- crossprod(matrix(rnorm(24), 6L, 4L)) + diag(4)
  predictor_names <- paste0(substr(name, 1L, 1L), seq_len(4L))
  dimnames(values) <- list(predictor_names, predictor_names)
  values
})
names(symmetric_dense) <- c("left", "right")
symmetric_sparse <- Map(function(matrix, uplo) {
  methods::as(
    Matrix::forceSymmetric(Matrix::Matrix(matrix, sparse = TRUE), uplo = uplo),
    "dsCMatrix"
  )
}, symmetric_dense, c("U", "L"))
symmetric_Xty <- stats::setNames(
  seq(-0.7, 0.7, length.out = 8L),
  unlist(lapply(symmetric_dense, colnames), use.names = FALSE)
)
symmetric_arguments <- list(
  n = 30L, XtX = symmetric_sparse, Xty = symmetric_Xty,
  ETA = list(model = "SpikeSlab"), residual_var = 1,
  X_means = stats::setNames(numeric(8L), names(symmetric_Xty)), y_mean = 0,
  iterations = 100, burnin = 40, compute_pve = TRUE
)
set.seed(514)
memory_fit <- do.call(
  blm_ss, c(symmetric_arguments, list(XtX_storage = "memory"))
)
set.seed(514)
speed_fit <- do.call(
  blm_ss, c(symmetric_arguments, list(XtX_storage = "speed"))
)
set.seed(514)
auto_fit <- do.call(
  blm_ss,
  c(symmetric_arguments, list(XtX_storage = "auto", XtX_memory_limit = 280))
)
stopifnot(
  isTRUE(all.equal(memory_fit$ETA, speed_fit$ETA, tolerance = 1e-10)),
  isTRUE(all.equal(memory_fit$ETA, auto_fit$ETA, tolerance = 1e-10)),
  isTRUE(all.equal(memory_fit$total_pve_mean, speed_fit$total_pve_mean,
                   tolerance = 1e-10)),
  identical(memory_fit$XtX_storage$representation,
            rep("symmetric_streaming", 2L)),
  identical(speed_fit$XtX_storage$representation,
            rep("general_sparse", 2L)),
  identical(auto_fit$XtX_storage$representation,
            rep("symmetric_streaming", 2L)),
  sum(auto_fit$XtX_storage$estimated_bytes) <= 280
)

# The separate dense centering correction remains exact with streaming blocks.
centered_arguments <- utils::modifyList(symmetric_arguments, list(
  X_means = stats::setNames(seq(-0.2, 0.2, length.out = 8L),
                            names(symmetric_Xty)),
  y_mean = 0.15, yty = 35, iterations = 80, burnin = 30
))
set.seed(5141)
centered_memory <- do.call(
  blm_ss, c(centered_arguments, list(XtX_storage = "memory"))
)
set.seed(5141)
centered_speed <- do.call(
  blm_ss, c(centered_arguments, list(XtX_storage = "speed"))
)
stopifnot(
  isTRUE(all.equal(centered_memory$ETA, centered_speed$ETA,
                   tolerance = 1e-10)),
  isTRUE(all.equal(centered_memory$total_pve_mean,
                   centered_speed$total_pve_mean, tolerance = 1e-10))
)

# With zero working means, Gram blocks are conditionally independent during
# the coefficient sweep. Threaded results are reproducible across repeated
# runs and thread counts, including prior blocks that cross Gram blocks.
threaded_eta <- list(
  normal = list(indices = c("l1", "r1"), model = "Normal"),
  selection = list(indices = c("l2", "r2"), model = "SpikeSlab"),
  shrinkage = list(indices = c("l3", "r3"), model = "GlobalLocal"),
  multi = list(indices = c("l4", "r4"), model = "SpikeMultiSlab")
)
threaded_arguments <- utils::modifyList(symmetric_arguments, list(
  residual_var = NULL, residual_shape = 2,
  residual_scale = 1, yty = 30, iterations = 120, burnin = 40
))
threaded_arguments$ETA <- threaded_eta
set.seed(515)
threaded_memory <- do.call(
  blm_ss,
  c(threaded_arguments, list(XtX_storage = "memory", nthreads = 2))
)
set.seed(515)
threaded_repeat <- do.call(
  blm_ss,
  c(threaded_arguments, list(XtX_storage = "memory", nthreads = 2))
)
set.seed(515)
threaded_three <- do.call(
  blm_ss,
  c(threaded_arguments, list(XtX_storage = "memory", nthreads = 3))
)
set.seed(515)
threaded_speed <- do.call(
  blm_ss,
  c(threaded_arguments, list(XtX_storage = "speed", nthreads = 2))
)
stopifnot(
  identical(threaded_memory, threaded_repeat),
  identical(threaded_memory$ETA, threaded_three$ETA),
  identical(threaded_memory$residual_var_samples,
            threaded_three$residual_var_samples),
  isTRUE(all.equal(threaded_memory$ETA, threaded_speed$ETA,
                   tolerance = 1e-10)),
  isTRUE(all.equal(threaded_memory$residual_var_samples,
                   threaded_speed$residual_var_samples, tolerance = 1e-10)),
  identical(threaded_memory$nthreads, 2L),
  identical(threaded_three$nthreads, 3L)
)

# List-specific validation prevents ambiguous predictor ordering and unsupported
# use of the reference R sampler.
unnamed_gram <- unname(gram_b)
dimnames(unnamed_gram) <- NULL
invalid_gram_list_calls <- list(
  function() blm_ss(
    20, list(), numeric(), list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss(
    20, list(gram_a, unnamed_gram), block_Xty,
    list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss(
    20, list(gram_a, gram_b), block_Xty,
    list(model = "Normal"), residual_var = 1, version = "R"
  ),
  function() blm_ss(
    20, list(gram_a, gram_b), block_Xty,
    list(model = "Normal"), residual_var = 1, XtX_memory_limit = 0
  ),
  function() blm_ss(
    20, list(gram_a, gram_b), block_Xty,
    list(model = "Normal"), residual_var = 1, nthreads = 0
  ),
  function() blm_ss(
    20, block_gram, block_Xty,
    list(model = "Normal"), residual_var = 1, nthreads = 2
  ),
  function() blm_ss(
    20, list(gram_a, gram_b), block_Xty,
    list(model = "Normal"), residual_var = 1,
    nthreads = 2, nchains = 2
  ),
  function() blm_ss(
    20, list(gram_a, gram_b), block_Xty,
    list(model = "Normal"), residual_var = 1,
    X_means = block_means, y_mean = 0, nthreads = 2
  )
)
stopifnot(all(vapply(
  invalid_gram_list_calls,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))

# Expected-sparsity calibration uses the original `n`, not the parser's
# two-row placeholder, and is identical for raw and sufficient-statistic fits.
calibration_eta_raw <- list(
  X = X, model = "GlobalLocal", expected_nonzero = 1.5,
  reference_residual_var = 0.64
)
calibration_eta_ss <- list(
  model = "GlobalLocal", expected_nonzero = 1.5,
  reference_residual_var = 0.64
)
calibration_scale <- 1.5 / (ncol(X) - 1.5) * 0.8 / sqrt(n)
set.seed(503)
calibrated_raw <- blm(
  y, ETA = calibration_eta_raw, residual_var = 1,
  iterations = 80, burnin = 30
)
set.seed(503)
calibrated_ss <- blm_ss(
  n, XtX, Xty, ETA = calibration_eta_ss, yty = yty,
  X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
  iterations = 80, burnin = 30
)
stopifnot(
  identical(calibrated_raw$ETA$ETA1$global_scale, calibration_scale),
  identical(calibrated_ss$ETA$ETA1$global_scale, calibration_scale),
  isTRUE(all.equal(calibrated_raw, calibrated_ss, tolerance = 1e-8))
)

# Raw and sufficient-statistic expected-PVE calibration agree for every model,
# including calibration of the residual inverse-gamma scale from centered yty.
for (model in c(
  "Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab"
)) {
  raw_specification <- list(
    X = X, model = model, expected_pve = 0.25
  )
  ss_specification <- list(model = model, expected_pve = 0.25)
  if (model == "GlobalLocal") {
    raw_specification$expected_nonzero <- 1
    ss_specification$expected_nonzero <- 1
  }
  set.seed(504)
  pve_raw <- blm(
    y, ETA = raw_specification, residual_shape = 3,
    iterations = 80, burnin = 30
  )
  set.seed(504)
  pve_ss <- blm_ss(
    n, XtX, Xty, ETA = ss_specification, yty = yty,
    X_means = colMeans(X), y_mean = mean(y), residual_shape = 3,
    iterations = 80, burnin = 30
  )
  expected_residual_var <- 0.75 * stats::var(y)
  stopifnot(
    isTRUE(all.equal(pve_raw, pve_ss, tolerance = 1e-8)),
    isTRUE(all.equal(
      pve_raw$residual_scale, (3 - 1) * expected_residual_var
    )),
    identical(pve_raw$residual_scale_calibrated, TRUE),
    identical(pve_ss$residual_scale_calibrated, TRUE),
    identical(pve_raw$expected_pve_total, 0.25),
    isTRUE(all.equal(
      pve_ss$reference_residual_var, expected_residual_var
    ))
  )
}

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

  version = "Rcpp"
)
set.seed(512)
sparse_dense_fit <- do.call(
  blm_ss,
  c(list(XtX = as.matrix(sparse_XtX)), sparse_args)
)
set.seed(512)
sparse_symmetric_fit <- do.call(
  blm_ss,
  c(list(XtX = sparse_XtX), sparse_args)
)
set.seed(512)
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
    burnin = 40
  ))
  model_args$ETA <- list(model = sparse_model)
  set.seed(513)
  dense_model_fit <- do.call(
    blm_ss,
    c(list(XtX = as.matrix(sparse_XtX)), model_args)
  )
  set.seed(513)
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
set.seed(507)
raw_fixed <- blm(
  y, ETA = list(X = X, model = "SpikeSlab"), residual_var = 1,
  iterations = 200, burnin = 80, version = "Rcpp"
)
set.seed(507)
ss_fixed <- blm_ss(
  n, XtX, Xty, ETA = list(model = "SpikeSlab"),
  X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
  iterations = 200, burnin = 80, version = "Rcpp"
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
    iterations = 100, burnin = 40
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
set.seed(508)
no_intercept_r <- suppressWarnings(blm_ss(
  n, XtX, Xty, ETA = list(model = "Normal"), yty = yty,
  residual_shape = 2, residual_scale = 1,
  iterations = 150, burnin = 50, version = "R"
))
set.seed(508)
no_intercept_rcpp <- suppressWarnings(blm_ss(
  n, XtX, Xty, ETA = list(model = "Normal"), yty = yty,
  residual_shape = 2, residual_scale = 1,
  iterations = 150, burnin = 50, version = "Rcpp"
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
  set.seed(509)
  blm_ss(
    nrow(rank_X), crossprod(rank_X), crossprod(rank_X, rank_y),
    ETA = list(model = "Normal"), yty = sum(rank_y^2),
    X_means = colMeans(rank_X), y_mean = mean(rank_y),
    residual_shape = 2, residual_scale = 1,
    iterations = 150, burnin = 50,
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
    residual_var = 1, iterations = 50, burnin = 20
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
  iterations = 50, burnin = 20
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
    n, XtX, stats::setNames(as.numeric(Xty), rev(colnames(XtX))),
    ETA = list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss(
    n, XtX, Xty, ETA = list(model = "Normal"),
    X_means = stats::setNames(colMeans(X), rev(colnames(XtX))),
    y_mean = mean(y), residual_var = 1
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
  ),
  function() blm_ss(
    n, XtX, Xty,
    ETA = list(model = "Normal", expected_pve = 0.2),
    residual_var = 1
  ),
  function() blm_ss(
    n, XtX, Xty, yty = yty,
    X_means = colMeans(X), y_mean = mean(y),
    ETA = list(model = "Normal"), residual_shape = 2
  ),
  function() blm_ss(
    n, XtX, Xty, yty = yty,
    X_means = colMeans(X), y_mean = mean(y),
    ETA = list(model = "Normal", expected_pve = 0.2),
    residual_shape = 1
  )
)
stopifnot(all(vapply(
  invalid_calls,
  function(call) inherits(suppressWarnings(try(call(), silent = TRUE)),
                           "try-error"),
  logical(1)
)))
