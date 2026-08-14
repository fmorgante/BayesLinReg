# Benchmark the independent native Hoermann--Leydold GIG sampler against the
# current GIGrvg backend. Run after installing BayesLinReg.

sample_size <- as.integer(Sys.getenv("BAYESLINREG_GIG_BENCH_N", "200000"))
repetitions <- as.integer(Sys.getenv("BAYESLINREG_GIG_BENCH_REPS", "7"))

cases <- data.frame(
  case = c("small parameter", "unshifted ratio", "shifted ratio", "negative lambda"),
  lambda = c(0.2, 0.8, 5, -0.5),
  chi = c(0.01, 0.25, 4, 2),
  psi = c(0.01, 1, 4, 3),
  stringsAsFactors = FALSE
)

time_sampler <- function(sampler, lambda, chi, psi) {
  times <- numeric(repetitions)
  for (iteration in seq_len(repetitions)) {
    set.seed(8100 + iteration)
    times[iteration] <- system.time(
      sampler(sample_size, lambda, chi, psi)
    )[["elapsed"]]
  }
  median(times)
}

results <- lapply(seq_len(nrow(cases)), function(index) {
  parameters <- cases[index, ]
  native_time <- time_sampler(
    BayesLinReg:::draw_gig_native_rcpp_cpp,
    parameters$lambda, parameters$chi, parameters$psi
  )
  reference_time <- time_sampler(
    BayesLinReg:::draw_gig_rcpp_cpp,
    parameters$lambda, parameters$chi, parameters$psi
  )
  data.frame(
    case = parameters$case,
    sample_size = sample_size,
    native_seconds = native_time,
    GIGrvg_seconds = reference_time,
    native_speedup = reference_time / native_time
  )
})

print(do.call(rbind, results), row.names = FALSE)
