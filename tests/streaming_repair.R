library(BayesLinReg)

# The incremental one-sided repair agrees with the diagnostic full-refresh
# fallback across both compressed LD layouts, serial and parallel sweeps, and
# the 100-iteration numerical refresh boundary.
make_variants <- function(ids, chromosome) {
  data.frame(
    CHR = chromosome, ID = ids, POS = seq_along(ids),
    A1 = rep(c("A", "C"), length.out = length(ids)),
    A0 = rep(c("C", "A"), length.out = length(ids))
  )
}

ids_a <- paste0("a", seq_len(6L))
ids_b <- paste0("b", seq_len(6L))
interval_R <- outer(
  seq_len(6L), seq_len(6L), function(i, j) 0.35^abs(i - j)
)
irregular_R <- diag(6L)
irregular_R[4L, 1L] <- irregular_R[1L, 4L] <- 0.2
irregular_R[6L, 2L] <- irregular_R[2L, 6L] <- -0.15
irregular_R[5L, 3L] <- irregular_R[3L, 5L] <- 0.1
dimnames(interval_R) <- list(ids_a, ids_a)
dimnames(irregular_R) <- list(ids_b, ids_b)
variants_a <- make_variants(ids_a, 1)
variants_b <- make_variants(ids_b, 2)
ld <- as_blm_ld(
  list(interval = interval_R, indexed = irregular_R),
  list(interval = variants_a, indexed = variants_b)
)
stopifnot(identical(
  ld$block_table$storage,
  c("interval_triangular", "indexed_triangular")
))

gwas <- transform(
  rbind(variants_a, variants_b),
  N = 500L,
  BETA = seq(-0.08, 0.09, length.out = 12L),
  SE = 0.05
)
ETA <- list(model = "SpikeMultiSlab")
common <- list(
  gwas = gwas, ld = ld, ETA = ETA,
  residual_shape = 2, residual_scale = 1,
  iterations = 130L, burnin = 30L,
  store_samples = TRUE, compute_pve = TRUE
)

fit_with_repair <- function(seed, nthreads, full_refresh) {
  if (full_refresh) {
    Sys.setenv(BAYESLINREG_FULL_STREAMING_REFRESH = "1")
  } else {
    Sys.unsetenv("BAYESLINREG_FULL_STREAMING_REFRESH")
  }
  set.seed(seed)
  do.call(blm_gwas, c(common, list(nthreads = nthreads)))
}

incremental_serial <- fit_with_repair(3101, 1L, FALSE)
full_serial <- fit_with_repair(3101, 1L, TRUE)
incremental_parallel <- fit_with_repair(3102, 2L, FALSE)
full_parallel <- fit_with_repair(3102, 2L, TRUE)
Sys.unsetenv("BAYESLINREG_FULL_STREAMING_REFRESH")

compare_fit <- function(incremental, full) {
  isTRUE(all.equal(
    incremental$ETA, full$ETA, tolerance = 1e-9
  )) &&
    isTRUE(all.equal(
      incremental$residual_var_samples,
      full$residual_var_samples,
      tolerance = 1e-9
    )) &&
    isTRUE(all.equal(
      incremental$total_pve_samples,
      full$total_pve_samples,
      tolerance = 1e-9
    ))
}

stopifnot(
  compare_fit(incremental_serial, full_serial),
  compare_fit(incremental_parallel, full_parallel)
)
