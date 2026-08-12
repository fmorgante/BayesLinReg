# Benchmark LD construction, compressed filtering, and ETA-order-independent
# fitting. This script is intentionally not run during package checks.
#
# Usage after installing BayesLinReg:
#   source(system.file(
#     "benchmarks", "benchmark_ld_preprocessing.R",
#     package = "BayesLinReg"
#   ))
#
# Environment variables optionally control the workload:
#   BLM_BENCH_P=20000 BLM_BENCH_BAND=50 BLM_BENCH_BLOCKS=20 Rscript ...

library(BayesLinReg)
library(Matrix)

integer_environment <- function(name, default) {
  value <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < 1L) {
    stop(sprintf("Environment variable `%s` must be a positive integer.", name))
  }
  value
}

p <- integer_environment("BLM_BENCH_P", 20000L)
if (p < 2L) stop("`BLM_BENCH_P` must be at least two.")
bandwidth <- min(integer_environment("BLM_BENCH_BAND", 50L), p - 1L)
number_of_blocks <- min(integer_environment("BLM_BENCH_BLOCKS", 20L), p)

make_variants <- function(ids, chromosome) {
  data.frame(
    CHR = chromosome,
    ID = ids,
    POS = seq_along(ids),
    A1 = rep(c("A", "C"), length.out = length(ids)),
    A0 = rep(c("C", "A"), length.out = length(ids))
  )
}

band_offsets <- seq.int(0L, bandwidth)
band_values <- lapply(band_offsets, function(offset) {
  rep(0.8^abs(offset), p - abs(offset))
})
interval_R <- bandSparse(
  p, p, k = band_offsets, diagonals = band_values, symmetric = TRUE
)
interval_ids <- paste0("interval", seq_len(p))
interval_variants <- make_variants(interval_ids, 1L)

block_sizes <- rep.int(p %/% number_of_blocks, number_of_blocks)
block_sizes[seq_len(p %% number_of_blocks)] <-
  block_sizes[seq_len(p %% number_of_blocks)] + 1L
block_R <- lapply(block_sizes, function(size) {
  local_bandwidth <- min(bandwidth, size - 1L)
  offsets <- seq.int(0L, local_bandwidth)
  values <- lapply(offsets, function(offset) {
    rep(0.8^abs(offset), size - abs(offset))
  })
  bandSparse(size, size, k = offsets, diagonals = values, symmetric = TRUE)
})
block_ids <- Map(
  function(size, block) paste0("block", block, "_", seq_len(size)),
  block_sizes, seq_along(block_sizes)
)
block_variants <- Map(make_variants, block_ids, seq_along(block_ids))
names(block_R) <- names(block_variants) <- paste0("block", seq_along(block_R))

irregular_R <- interval_R
if (length(irregular_R@x)) {
  retain <- seq_along(irregular_R@x) %% 5L == 0L |
    irregular_R@i == rep.int(seq_len(ncol(irregular_R)) - 1L,
                             diff(irregular_R@p))
  irregular_R@x[!retain] <- 0
  irregular_R <- drop0(irregular_R)
  diag(irregular_R) <- 1
}

timed <- function(expression) {
  gc()
  timing <- system.time(value <- force(expression))
  list(value = value, elapsed = unname(timing[["elapsed"]]))
}

interval_build <- timed(as_blm_ld(interval_R, interval_variants))
block_build <- timed(as_blm_ld(block_R, block_variants))
irregular_build <- timed(as_blm_ld(irregular_R, interval_variants))

retained <- seq_len(p) %% 20L != 0L
filter_time <- timed(BayesLinReg:::.subset_blm_ld(
  interval_build$value, which(retained)
))

benchmark_gwas <- transform(
  interval_variants, N = 1000L, BETA = 0, SE = 0.05
)
ordered_eta <- list(
  first = list(indices = interval_ids[seq_len(p %/% 2L)], model = "Normal"),
  second = list(
    indices = interval_ids[seq.int(p %/% 2L + 1L, p)], model = "SpikeSlab"
  )
)
reversed_eta <- lapply(ordered_eta, function(specification) {
  specification$indices <- rev(specification$indices)
  specification
})
set.seed(1)
ordered_fit <- timed(blm_gwas(
  benchmark_gwas, interval_build$value, ordered_eta,
  residual_var = 1, iterations = 3L, burnin = 1L
))
set.seed(1)
reversed_fit <- timed(blm_gwas(
  benchmark_gwas, interval_build$value, reversed_eta,
  residual_var = 1, iterations = 3L, burnin = 1L
))

results <- data.frame(
  operation = c(
    "construct_interval", "construct_many_blocks", "construct_irregular",
    "filter_five_percent", "fit_ordered_eta", "fit_reversed_eta"
  ),
  predictors = c(p, p, p, sum(retained), p, p),
  elapsed_seconds = c(
    interval_build$elapsed, block_build$elapsed, irregular_build$elapsed,
    filter_time$elapsed, ordered_fit$elapsed, reversed_fit$elapsed
  ),
  object_megabytes = c(
    as.numeric(object.size(interval_build$value)),
    as.numeric(object.size(block_build$value)),
    as.numeric(object.size(irregular_build$value)),
    as.numeric(object.size(filter_time$value)),
    as.numeric(object.size(ordered_fit$value)),
    as.numeric(object.size(reversed_fit$value))
  ) / 1024^2,
  computational_blocks = c(
    length(interval_build$value$blocks), length(block_build$value$blocks),
    length(irregular_build$value$blocks), length(filter_time$value$blocks),
    length(interval_build$value$blocks), length(interval_build$value$blocks)
  )
)

print(results, row.names = FALSE, digits = 4)
