library(BayesLinReg)

X <- cbind(
  first = 1:20,
  second = rep(c(0, 1), 10),
  third = sin(seq_len(20)),
  fourth = cos(seq_len(20))
)
y <- 1 + 2 * X[, "first"] - X[, "second"]

diagnostic_eta <- list(
  normal = list(
    X = X[, "first", drop = FALSE], model = "Normal",
    var_shape = 2, var_scale = 10
  ),
  selection = list(
    X = X[, "second", drop = FALSE], model = "SpikeSlab",
    slab_shape = 2, slab_scale = 10
  ),
  shrinkage = list(
    X = X[, c("third", "fourth")], model = "GlobalLocal"
  )
)

fit_one <- blm(
  y,
  ETA = diagnostic_eta,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 100,
  burnin = 40,
  seed = 101,
  version = "R"
)
fit_two <- blm(
  y,
  ETA = diagnostic_eta,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 100,
  burnin = 40,
  seed = 102,
  version = "R"
)

combined_fit <- fit_one
for (block_name in names(combined_fit$ETA)) {
  combined_fit$ETA[[block_name]]$coefficient_samples <- rbind(
    fit_one$ETA[[block_name]]$coefficient_samples,
    fit_two$ETA[[block_name]]$coefficient_samples
  )
}
combined_fit$ETA$selection$pi_samples <- c(
  fit_one$ETA$selection$pi_samples,
  fit_two$ETA$selection$pi_samples
)
combined_fit$ETA$normal$normal_var_samples <- c(
  fit_one$ETA$normal$normal_var_samples,
  fit_two$ETA$normal$normal_var_samples
)
combined_fit$ETA$selection$slab_var_samples <- c(
  fit_one$ETA$selection$slab_var_samples,
  fit_two$ETA$selection$slab_var_samples
)
combined_fit$ETA$shrinkage$tau_sq_samples <- c(
  fit_one$ETA$shrinkage$tau_sq_samples,
  fit_two$ETA$shrinkage$tau_sq_samples
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

diagnostics <- assess_convergence(combined_fit, plot = FALSE)
expected_parameters <- c(
  "intercept",
  "residual_var",
  "normal_var_normal",
  "pi_selection",
  "slab_var_selection",
  "tau_sq_shrinkage"
)
stopifnot(
  identical(names(diagnostics), c(
    "rhat", "geweke", "effective_sample_size", "nchains",
    "draws_per_chain"
  )),
  identical(names(diagnostics$rhat), expected_parameters),
  identical(dim(diagnostics$geweke), c(6L, 2L)),
  identical(rownames(diagnostics$geweke), expected_parameters),
  identical(
    names(diagnostics$effective_sample_size),
    expected_parameters
  ),
  diagnostics$nchains == 2L,
  diagnostics$draws_per_chain == 60L,
  all(is.finite(diagnostics$rhat)),
  all(diagnostics$effective_sample_size > 0),
  !any(grepl("coefficient|inclusion|local_var|::", expected_parameters))
)

# Trace plotting covers every assessed parameter.
trace_file <- tempfile(fileext = ".pdf")
grDevices::pdf(trace_file)
trace_titles <- BayesLinReg:::.plot_blm_traces(
  BayesLinReg:::.as_blm_mcmc_list(combined_fit)
)
invisible(assess_convergence(combined_fit, plot = TRUE))
grDevices::dev.off()
stopifnot(
  file.info(trace_file)$size > 0,
  identical(trace_titles, paste("Trace of", expected_parameters)),
  "Trace of normal_var_normal" %in% trace_titles,
  "Trace of pi_selection" %in% trace_titles,
  "Trace of slab_var_selection" %in% trace_titles
)

single_chain_diagnostics <- assess_convergence(fit_one, plot = FALSE)
stopifnot(
  single_chain_diagnostics$nchains == 1L,
  all(is.na(single_chain_diagnostics$rhat))
)

known_fit <- blm(
  y,
  ETA = list(X = X, model = "Normal", var_shape = 2, var_scale = 10),
  residual_var = 1,
  iterations = 100,
  burnin = 40,
  seed = 103
)
summary_only_fit <- blm(
  y,
  ETA = list(X = X, model = "Normal"),
  residual_var = 1,
  iterations = 100,
  burnin = 40,
  seed = 104,
  store_samples = FALSE
)
stopifnot(
  assess_convergence(known_fit, plot = FALSE)$nchains == 1L,
  grepl(
    "store_samples = TRUE",
    as.character(try(assess_convergence(summary_only_fit), silent = TRUE)),
    fixed = TRUE
  ),
  inherits(
    try(assess_convergence(fit_one, plot = NA), silent = TRUE),
    "try-error"
  )
)
