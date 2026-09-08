library(BayesLinReg)

ids <- c("rs1", "rs2")
R <- matrix(c(1, 0.95, 0.95, 1), 2L)
dimnames(R) <- list(ids, ids)
variants <- data.frame(
  CHR = 1, ID = ids, POS = 1:2,
  A1 = "A", A0 = "C"
)
ld <- as_blm_ld(R, variants)
gwas <- transform(
  variants, N = 1000, BETA = c(0.3, -0.3), SE = 0.05
)
original_ld <- ld
original_gwas <- gwas

set.seed(1501)
diagnostics <- diagnose_gwas_ld(
  gwas, ld, window_variants = 2L, overlap_variants = 0L,
  store_variant_report = "all"
)
set.seed(1501)
diagnostics_repeat <- diagnose_gwas_ld(
  gwas, ld, window_variants = 2L, overlap_variants = 0L,
  store_variant_report = "all"
)
stopifnot(
  identical(diagnostics, diagnostics_repeat),
  inherits(diagnostics, "blm_gwas_ld_diagnostics"),
  identical(ld, original_ld),
  identical(gwas, original_gwas),
  diagnostics$summary[["retained_variants"]] == 2,
  diagnostics$summary[["assessed_variants"]] == 2,
  diagnostics$summary[["flagged_variants"]] == 2,
  all(diagnostics$variant_report$conditional_outlier),
  all(diagnostics$variant_report$possible_allele_flip),
  all(diagnostics$variant_report$status == "possible_allele_flip"),
  all(diagnostics$variant_report$tagging > 0.8),
  all(diagnostics$variant_report$p_value < 5e-8),
  nrow(diagnostics$variant_failures) == 2L,
  diagnostics$block_report$failed_windows == 0L,
  isFALSE(diagnostics$filtering_applied),
  identical(
    diagnostics$score_correlation_assumption,
    "gwas_score_correlation_equals_reference_ld"
  )
)

# Dense packed and symmetric sparse native LD use the same bounded-window
# extraction semantics.
sparse_R <- Matrix::forceSymmetric(Matrix::Matrix(R, sparse = TRUE), "L")
sparse_ld <- as_blm_ld(sparse_R, variants)
set.seed(1501)
sparse_diagnostics <- diagnose_gwas_ld(
  gwas, sparse_ld, window_variants = 2L, overlap_variants = 0L,
  store_variant_report = "all"
)
stopifnot(isTRUE(all.equal(
  diagnostics$variant_report,
  sparse_diagnostics$variant_report,
  tolerance = 0
)))

# Default storage retains failures without duplicating the full genome-wide
# report, and explicit allele metadata orientation is honored.
reversed <- gwas
reversed$A1[1] <- "C"
reversed$A0[1] <- "A"
reversed$BETA[1] <- -reversed$BETA[1]
set.seed(1502)
failure_only <- diagnose_gwas_ld(
  reversed, ld, window_variants = 2L, overlap_variants = 0L
)
stopifnot(
  is.null(failure_only$variant_report),
  all(failure_only$variant_failures$possible_allele_flip),
  identical(
    failure_only$variant_failures$ld_aligned_z,
    diagnostics$variant_failures$ld_aligned_z
  ),
  failure_only$variant_failures$effect_orientation[1] == -1,
  failure_only$harmonization_report[["flipped"]] == 1
)

# Failure-only storage preserves a typed zero-row table when no variant is
# flagged.
independent_ld <- as_blm_ld(diag(2), variants)
set.seed(1504)
no_failures <- diagnose_gwas_ld(
  gwas, independent_ld, window_variants = 2L, overlap_variants = 0L
)
stopifnot(
  no_failures$summary[["flagged_variants"]] == 0,
  nrow(no_failures$variant_failures) == 0L,
  identical(names(no_failures$variant_failures),
            names(diagnostics$variant_failures))
)

# Strong sample-size differences are reported rather than silently interpreted
# as calibrated conditional tests.
heterogeneous <- gwas
heterogeneous$N <- c(1000, 700)
set.seed(1503)
heterogeneous_diagnostics <- suppressWarnings(diagnose_gwas_ld(
  heterogeneous, ld, window_variants = 2L, overlap_variants = 0L,
  store_variant_report = "all"
))
stopifnot(
  heterogeneous_diagnostics$summary[["dissimilar_n_blocks"]] == 1,
  isFALSE(heterogeneous_diagnostics$block_report$similar_sample_sizes),
  all(!heterogeneous_diagnostics$variant_report$similar_sample_sizes)
)

# Invalid diagnostic controls and eigen LD input fail explicitly.
eigen_ld <- as_blm_ld_eigen(ld, prop_var = 1)
errors <- list(
  try(diagnose_gwas_ld(gwas, eigen_ld), silent = TRUE),
  try(diagnose_gwas_ld(gwas, ld, window_variants = 1), silent = TRUE),
  try(diagnose_gwas_ld(gwas, ld, overlap_variants = -1), silent = TRUE),
  try(diagnose_gwas_ld(gwas, ld, ld_shrink = 1), silent = TRUE),
  try(diagnose_gwas_ld(gwas, ld, p_threshold = 0), silent = TRUE),
  try(diagnose_gwas_ld(gwas, ld, min_tagging = 1), silent = TRUE)
)
stopifnot(all(vapply(errors, inherits, logical(1), "try-error")))

# BLAS thread requests are interpreted conservatively by vendor.
requested_blas_threads <- getFromNamespace(
  ".requested_blas_threads", "BayesLinReg"
)
stopifnot(
  requested_blas_threads("Intel oneMKL", c(MKL_NUM_THREADS = "1")) == 1L,
  requested_blas_threads("Intel oneMKL", c(MKL_NUM_THREADS = "8")) == 8L,
  requested_blas_threads(
    "Intel oneMKL",
    c(MKL_NUM_THREADS = "8", MKL_DOMAIN_NUM_THREADS = "MKL_BLAS=1")
  ) == 1L,
  requested_blas_threads(
    "OpenBLAS", c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "4")
  ) == 4L,
  requested_blas_threads("reference BLAS", character()) == 1L,
  is.na(requested_blas_threads("unknown", character()))
)

# Compare the fused kernel with direct R conditional calculations over several
# overlapping windows.
parallel_R <- matrix(0.2, 6L, 6L)
diag(parallel_R) <- 1
parallel_variants <- data.frame(
  CHR = 1, ID = paste0("parallel", seq_len(6L)), POS = seq_len(6L),
  A1 = "A", A0 = "C"
)
parallel_ld <- as_blm_ld(parallel_R, parallel_variants)
parallel_gwas <- transform(
  parallel_variants, N = 1000, BETA = seq(-0.2, 0.3, length.out = 6L),
  SE = 0.05
)
reference_diagnostic <- function(R, z, window_variants, overlap_variants,
                                 ld_shrink) {
  p <- length(z)
  answer <- list(
    predicted_z = rep(NA_real_, p), conditional_variance = rep(NA_real_, p),
    tagging = rep(NA_real_, p), statistic = rep(NA_real_, p),
    flip_log_likelihood_ratio = rep(NA_real_, p)
  )
  for (start in seq.int(1L, p, by = window_variants)) {
    core <- seq.int(start, min(p, start + window_variants - 1L))
    expanded <- seq.int(
      max(1L, core[[1L]] - overlap_variants),
      min(p, core[[length(core)]] + overlap_variants)
    )
    permutation <- sample.int(length(expanded))
    group <- integer(length(expanded))
    group[permutation] <- rep_len(1:2, length(expanded))
    for (target_group in 1:2) {
      target <- which(expanded %in% core & group == target_group)
      predictor <- which(group != target_group)
      if (!length(target) || !length(predictor)) next
      covariance <- (1 - ld_shrink) *
        R[expanded[predictor], expanded[predictor], drop = FALSE]
      diag(covariance) <- 1
      cross <- (1 - ld_shrink) *
        R[expanded[target], expanded[predictor], drop = FALSE]
      fitted <- as.numeric(cross %*% solve(covariance, z[expanded[predictor]]))
      tagging <- rowSums(cross * t(solve(covariance, t(cross))))
      variance <- 1 - tagging
      global <- expanded[target]
      standardized <- (z[global] - fitted) / sqrt(variance)
      answer$predicted_z[global] <- fitted
      answer$conditional_variance[global] <- variance
      answer$tagging[global] <- tagging
      answer$statistic[global] <- standardized^2
      answer$flip_log_likelihood_ratio[global] <-
        -2 * z[global] * fitted / variance
    }
  }
  answer
}
set.seed(1505)
reference_windows <- reference_diagnostic(
  parallel_R, parallel_gwas$BETA / parallel_gwas$SE, 2L, 1L, 0.05
)
set.seed(1505)
serial_windows <- diagnose_gwas_ld(
  parallel_gwas, parallel_ld, window_variants = 2L,
  overlap_variants = 1L, nthreads = 1L, store_variant_report = "all"
)
for (name in names(reference_windows)) {
  stopifnot(isTRUE(all.equal(
    serial_windows$variant_report[[name]], reference_windows[[name]],
    tolerance = 1e-12
  )))
}

# Independent windows are reproducible whenever nested BLAS parallelism is not
# requested.
build_information <- blm_build_info()
parallel_permitted <- !build_information$eigen_blas || identical(
  requested_blas_threads(build_information$blas, Sys.getenv()), 1L
)

# A single native call queues windows across LD blocks, so threading remains
# available when every block contributes only one window.
second_ids <- c("rs3", "rs4")
second_variants <- transform(variants, CHR = 2, ID = second_ids)
second_R <- R
dimnames(second_R) <- list(second_ids, second_ids)
multi_block_ld <- as_blm_ld(
  list(first = R, second = second_R),
  list(first = variants, second = second_variants)
)
multi_block_gwas <- rbind(
  gwas,
  transform(gwas, CHR = 2, ID = second_ids, BETA = rev(BETA))
)
set.seed(1506)
multi_block_serial <- diagnose_gwas_ld(
  multi_block_gwas, multi_block_ld,
  window_variants = 2000L, overlap_variants = 0L,
  nthreads = 1L, store_variant_report = "all"
)
stopifnot(
  identical(multi_block_serial$block_report$windows, c(1L, 1L)),
  nrow(multi_block_serial$variant_report) == 4L
)
if (parallel_permitted) {
  set.seed(1505)
  parallel_windows <- diagnose_gwas_ld(
    parallel_gwas, parallel_ld, window_variants = 2L,
    overlap_variants = 1L, nthreads = 2L, store_variant_report = "all"
  )
  stopifnot(
    identical(serial_windows$variant_report, parallel_windows$variant_report),
    identical(serial_windows$block_report, parallel_windows$block_report)
  )

  set.seed(1506)
  multi_block_parallel <- diagnose_gwas_ld(
    multi_block_gwas, multi_block_ld,
    window_variants = 2000L, overlap_variants = 0L,
    nthreads = 2L, store_variant_report = "all"
  )
  stopifnot(
    identical(
      multi_block_serial$variant_report,
      multi_block_parallel$variant_report
    ),
    identical(
      multi_block_serial$block_report,
      multi_block_parallel$block_report
    )
  )
}
