# Benchmark the incremental one-sided LD repair against the internal
# full-reconstruction fallback. Run from an installed BayesLinReg source tree.
library(BayesLinReg)

block_size <- 2000L
block_count <- 8L
iterations <- 120L
nthreads <- 4L
repetitions <- 3L

rho <- 0.95
correlation <- rho^abs(outer(
  seq_len(block_size), seq_len(block_size), `-`
))
ids <- paste0("v", seq_len(block_size * block_count))
R_list <- stats::setNames(
  rep(list(correlation), block_count),
  paste0("block", seq_len(block_count))
)
variants <- lapply(seq_len(block_count), function(block) {
  indices <- ((block - 1L) * block_size + 1L):(block * block_size)
  data.frame(
    CHR = block, ID = ids[indices], POS = indices,
    A1 = "A", A0 = "C"
  )
})
names(variants) <- names(R_list)

ld <- as_blm_ld(R_list, variants)
gwas <- transform(
  do.call(rbind, variants), N = 100000L, BETA = 0, SE = 0.01
)
rownames(gwas) <- NULL
ETA <- list(model = "SpikeMultiSlab")

run_fit <- function(full_refresh, seed) {
  if (full_refresh) {
    Sys.setenv(BAYESLINREG_FULL_STREAMING_REFRESH = "1")
  } else {
    Sys.unsetenv("BAYESLINREG_FULL_STREAMING_REFRESH")
  }
  set.seed(seed)
  system.time(blm_gwas(
    gwas, ld, ETA, residual_var = 1,
    iterations = iterations, burnin = 20L,
    nthreads = nthreads
  ))[["elapsed"]]
}

full <- incremental <- numeric(repetitions)
for (repetition in seq_len(repetitions)) {
  full[repetition] <- run_fit(TRUE, 900L + repetition)
  incremental[repetition] <- run_fit(FALSE, 900L + repetition)
}
Sys.unsetenv("BAYESLINREG_FULL_STREAMING_REFRESH")

result <- data.frame(
  kernel = c("full_refresh", "incremental_repair"),
  median_seconds = c(median(full), median(incremental))
)
result$speedup_relative_to_full <- result$median_seconds[[1L]] /
  result$median_seconds
print(result, row.names = FALSE)
