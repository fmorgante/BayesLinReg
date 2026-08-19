library(BayesLinReg)

# With a list dsCMatrix, a reordered ETA can put the serial sampler order at
# odds with the block's lower-triangular storage order.  The memory path must
# instead follow the block-local order and agree with the corresponding full
# sparse transition after coefficients are aligned by name.
predictor_names <- c("x1", "x2")
gram <- matrix(
  c(4, 1.5,
    1.5, 3),
  nrow = 2L,
  dimnames = list(predictor_names, predictor_names)
)
symmetric_gram <- methods::as(
  Matrix::forceSymmetric(Matrix::Matrix(gram, sparse = TRUE), uplo = "L"),
  "dsCMatrix"
)
common <- list(
  n = 20L,
  XtX = list(region = symmetric_gram),
  Xty = stats::setNames(c(0.6, -0.4), predictor_names),
  residual_var = 1,
  iterations = 90L,
  burnin = 30L,
  store_samples = TRUE,
  nthreads = 1L
)

set.seed(2601)
memory_fit <- suppressWarnings(do.call(blm_ss, c(
  common,
  list(
    ETA = list(normal = list(
      indices = c("x2", "x1"), model = "Normal", var = 1,
      standardize = FALSE
    )),
    XtX_storage = "memory"
  )
)))
set.seed(2601)
speed_fit <- suppressWarnings(do.call(blm_ss, c(
  common,
  list(
    ETA = list(normal = list(
      indices = c("x1", "x2"), model = "Normal", var = 1,
      standardize = FALSE
    )),
    XtX_storage = "speed"
  )
)))

stopifnot(
  identical(memory_fit$XtX_storage$representation, "symmetric_streaming"),
  identical(speed_fit$XtX_storage$representation, "general_sparse"),
  isTRUE(all.equal(
    memory_fit$ETA$normal$coefficient_samples[, predictor_names],
    speed_fit$ETA$normal$coefficient_samples,
    tolerance = 1e-10
  ))
)
