# Direct distributional checks for the independent native GIG sampler. The
# cases exercise the small-parameter, unshifted ratio-of-uniforms, shifted
# ratio-of-uniforms, reciprocal negative-lambda, and unequal-scale paths.
gig_cases <- data.frame(
  case = c(
    "zero lambda small parameter",
    "positive small parameter",
    "unshifted ratio",
    "shifted ratio",
    "negative lambda",
    "unequal scale small parameter",
    "unequal scale shifted ratio"
  ),
  lambda = c(0, 0.2, 0.8, 5, -0.5, 0.5, -3),
  chi = c(0.01, 0.01, 0.25, 4, 2, 1e-3, 16),
  psi = c(0.01, 0.01, 1, 4, 3, 10, 0.25),
  stringsAsFactors = FALSE
)
gig_sample_size <- 30000L

# E[X^r] = (chi / psi)^(r / 2) K_(lambda + r)(sqrt(chi psi)) /
# K_lambda(sqrt(chi psi)). Exponentially scaled Bessel functions share the
# same exponential factor, which cancels in their ratio.
gig_moment <- function(order, lambda, chi, psi) {
  beta <- sqrt(chi * psi)
  exp(0.5 * order * (log(chi) - log(psi))) *
    besselK(beta, lambda + order, expon.scaled = TRUE) /
    besselK(beta, lambda, expon.scaled = TRUE)
}

native_gig <- BayesLinReg:::draw_gig_native_rcpp_cpp
for (case_index in seq_len(nrow(gig_cases))) {
  parameters <- gig_cases[case_index, ]
  set.seed(9100 + case_index)
  native_draws <- native_gig(
    gig_sample_size, parameters$lambda, parameters$chi, parameters$psi
  )
  stopifnot(
    all(is.finite(native_draws)),
    all(native_draws > 0)
  )

  # Check the first two raw moments using their exact Monte Carlo standard
  # errors, derived from moments through order four. Six standard errors is
  # deliberately conservative for stable checks across platforms and RNGs.
  moment_z <- vapply(1:2, function(order) {
    expected <- gig_moment(
      order, parameters$lambda, parameters$chi, parameters$psi
    )
    expected_square <- gig_moment(
      2 * order, parameters$lambda, parameters$chi, parameters$psi
    )
    standard_error <- sqrt(
      (expected_square - expected^2) / gig_sample_size
    )
    (mean(native_draws^order) - expected) / standard_error
  }, numeric(1))

  # Compare the complete empirical distribution with the established GIGrvg
  # implementation. The fixed D threshold is much wider than ordinary
  # two-sample fluctuation at this sample size, avoiding p-value flakiness.
  set.seed(10100 + case_index)
  reference_draws <- GIGrvg::rgig(
    gig_sample_size, parameters$lambda, parameters$chi, parameters$psi
  )
  reference_distance <- unname(suppressWarnings(stats::ks.test(
    native_draws, reference_draws, exact = FALSE
  )$statistic))

  # If Y ~ GIG(-lambda, psi, chi), then 1/Y has the target distribution. This
  # independently exercises the native negative-lambda reciprocal transform.
  set.seed(11100 + case_index)
  reciprocal_draws <- 1 / native_gig(
    gig_sample_size, -parameters$lambda, parameters$psi, parameters$chi
  )
  reciprocal_distance <- unname(suppressWarnings(stats::ks.test(
    native_draws, reciprocal_draws, exact = FALSE
  )$statistic))

  stopifnot(
    max(abs(moment_z)) < 6,
    reference_distance < 0.025,
    reciprocal_distance < 0.025
  )
}
