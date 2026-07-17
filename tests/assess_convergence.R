library(blm)

X <- cbind(
  first = 1:20,
  second = rep(c(0, 1), 10)
)
y <- 1 + 2 * X[, "first"] - X[, "second"]

fit_one <- multiple_blm(
  y, X, 10,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 100,
  burnin = 40,
  seed = 101,
  version = "R"
)
fit_two <- multiple_blm(
  y, X, 10,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 100,
  burnin = 40,
  seed = 102,
  version = "R"
)

combined_fit <- fit_one
combined_fit$coefficient_samples <- rbind(
  fit_one$coefficient_samples,
  fit_two$coefficient_samples
)
combined_fit$intercept_samples <- c(
  fit_one$intercept_samples,
  fit_two$intercept_samples
)
combined_fit$residual_var_samples <- c(
  fit_one$residual_var_samples,
  fit_two$residual_var_samples
)
combined_fit$chain_id <- rep.int(1:2, c(60L, 60L))

# Include spike-and-slab pi draws.
set.seed(103)
combined_fit$pi_samples <- stats::rbeta(120, 2, 2)

diagnostics <- assess_convergence(combined_fit, plot = FALSE)
expected_parameters <- c(
  "intercept", "residual_var", "pi"
)
stopifnot(
  identical(names(diagnostics), c(
    "rhat", "geweke", "effective_sample_size", "nchains",
    "draws_per_chain"
  )),
  identical(names(diagnostics$rhat), expected_parameters),
  identical(dim(diagnostics$geweke), c(3L, 2L)),
  identical(rownames(diagnostics$geweke), expected_parameters),
  identical(
    names(diagnostics$effective_sample_size),
    expected_parameters
  ),
  diagnostics$nchains == 2L,
  diagnostics$draws_per_chain == 60L,
  all(is.finite(diagnostics$rhat)),
  all(diagnostics$effective_sample_size > 0)
)

# Trace plotting covers every assessed parameter.
trace_file <- tempfile(fileext = ".pdf")
grDevices::pdf(trace_file)
invisible(assess_convergence(combined_fit, plot = TRUE))
grDevices::dev.off()
stopifnot(file.info(trace_file)$size > 0)

single_chain_diagnostics <- assess_convergence(fit_one, plot = FALSE)
stopifnot(
  single_chain_diagnostics$nchains == 1L,
  all(is.na(single_chain_diagnostics$rhat))
)

known_fit <- multiple_blm(y, X, 10, residual_var = 1)
stopifnot(
  inherits(try(assess_convergence(known_fit), silent = TRUE), "try-error"),
  inherits(
    try(assess_convergence(fit_one, plot = NA), silent = TRUE),
    "try-error"
  )
)
