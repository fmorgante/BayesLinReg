#ifndef BAYESLINREG_GIG_SAMPLER_H
#define BAYESLINREG_GIG_SAMPLER_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace bayeslinreg {

// Independent implementation of the rejection algorithms described by
// Hoermann and Leydold (2014), Statistics and Computing 24:547--557.
// The RNG must provide uniform() draws strictly between zero and one.

inline double gig_log_kernel(
    const double x,
    const double lambda,
    const double beta) {
  return (lambda - 1.0) * std::log(x) -
    0.5 * beta * (x + 1.0 / x);
}

inline double gig_mode(
    const double lambda,
    const double beta) {
  const double order = lambda - 1.0;
  const double absolute_order = std::abs(order);
  const double maximum = std::max(absolute_order, beta);
  const double minimum = std::min(absolute_order, beta);
  const double ratio = minimum / maximum;
  const double radius = maximum * std::sqrt(1.0 + ratio * ratio);
  if (order >= 0.0) {
    return (radius + order) / beta;
  }
  return beta / (radius - order);
}

inline double gig_scaled_kernel(
    const double x,
    const double lambda,
    const double beta,
    const double log_kernel_at_mode) {
  const double log_value =
    gig_log_kernel(x, lambda, beta) - log_kernel_at_mode;
  return std::exp(std::min(0.0, log_value));
}

template <typename UniformRng>
double sample_gig_small_parameter(
    const double lambda,
    const double beta,
    UniformRng& rng) {
  const double one_minus_lambda = 1.0 - lambda;
  const double mode = gig_mode(lambda, beta);
  const double log_mode = gig_log_kernel(mode, lambda, beta);
  const double x0 = beta / one_minus_lambda;
  const double upper = 2.0 / beta;
  const double x_star = std::max(x0, upper);
  const double log_x0 = std::log(x0);
  const double log_upper = std::log(upper);

  const double log_area1 = log_x0;
  const double log_k2 = -beta - log_mode;
  double log_area2 = -std::numeric_limits<double>::infinity();
  if (x0 < upper) {
    const double log_ratio = log_upper - log_x0;
    if (std::abs(lambda) < 1e-8) {
      log_area2 = log_k2 + std::log(log_ratio);
    } else {
      const double log_power_difference = lambda * log_upper +
        std::log1p(-std::exp(lambda * (log_x0 - log_upper)));
      log_area2 = log_k2 + log_power_difference - std::log(lambda);
    }
  }

  const double log_k3 = (lambda - 1.0) * std::log(x_star) - log_mode;
  const double log_area3 = std::log(2.0 / beta) + log_k3 -
    0.5 * beta * x_star;
  const double maximum_log_area = std::max({
    log_area1, log_area2, log_area3
  });
  const double area1 = std::exp(log_area1 - maximum_log_area);
  const double area2 = std::exp(log_area2 - maximum_log_area);
  const double area3 = std::exp(log_area3 - maximum_log_area);
  const double total_area = area1 + area2 + area3;

  for (int attempt = 0; attempt < 1000000; ++attempt) {
    const double log_accept_uniform = std::log(rng.uniform());
    const double selector = total_area * rng.uniform();
    double x = 0.0;
    double log_hat = 0.0;

    if (selector <= area1) {
      x = x0 * rng.uniform();
    } else if (selector <= area1 + area2) {
      const double mixture_uniform = rng.uniform();
      if (std::abs(lambda) < 1e-8) {
        x = std::exp(log_x0 + mixture_uniform * (log_upper - log_x0));
      } else {
        const double lower_term = std::log1p(-mixture_uniform) +
          lambda * log_x0;
        const double upper_term = std::log(mixture_uniform) +
          lambda * log_upper;
        const double maximum_term = std::max(lower_term, upper_term);
        const double log_power = maximum_term + std::log(
          std::exp(lower_term - maximum_term) +
            std::exp(upper_term - maximum_term)
        );
        x = std::exp(log_power / lambda);
      }
      log_hat = log_k2 + (lambda - 1.0) * std::log(x);
    } else {
      x = x_star - 2.0 / beta * std::log(rng.uniform());
      log_hat = log_k3 - 0.5 * beta * x;
    }

    const double target_log = gig_log_kernel(x, lambda, beta) - log_mode;
    if (std::isfinite(x) && x > 0.0 && std::isfinite(target_log) &&
        log_accept_uniform + log_hat <= std::min(0.0, target_log)) {
      return x;
    }
  }
  throw std::runtime_error("The native GIG rejection sampler did not converge.");
}

template <typename UniformRng>
double sample_gig_ratio_unshifted(
    const double lambda,
    const double beta,
    UniformRng& rng) {
  const double mode = gig_mode(lambda, beta);
  const double log_mode = gig_log_kernel(mode, lambda, beta);
  const double x_plus = (
    (1.0 + lambda) +
      std::sqrt((1.0 + lambda) * (1.0 + lambda) + beta * beta)
  ) / beta;
  const double v_plus = 1.0;
  const double u_plus = x_plus * std::sqrt(
    gig_scaled_kernel(x_plus, lambda, beta, log_mode)
  );

  for (int attempt = 0; attempt < 1000000; ++attempt) {
    const double u = u_plus * rng.uniform();
    const double v = v_plus * rng.uniform();
    const double x = u / v;
    if (std::isfinite(x) && x > 0.0 &&
        v * v <= gig_scaled_kernel(x, lambda, beta, log_mode)) {
      return x;
    }
  }
  throw std::runtime_error("The native GIG ratio-of-uniforms sampler did not converge.");
}

template <typename UniformRng>
double sample_gig_ratio_shifted(
    const double lambda,
    const double beta,
    UniformRng& rng) {
  const double mode = gig_mode(lambda, beta);
  const double log_mode = gig_log_kernel(mode, lambda, beta);
  const double a = -2.0 * (lambda + 1.0) / beta - mode;
  const double b = 2.0 * (lambda - 1.0) * mode / beta - 1.0;
  const double c = mode;
  const double cubic_p = b - a * a / 3.0;
  const double cubic_q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
  const double acos_argument = std::clamp(
    -0.5 * cubic_q * std::sqrt(-27.0 / (cubic_p * cubic_p * cubic_p)),
    -1.0,
    1.0
  );
  const double phi = std::acos(acos_argument);
  const double root_scale = std::sqrt(-4.0 * cubic_p / 3.0);
  const double pi = std::acos(-1.0);
  const double x_minus = root_scale * std::cos(phi / 3.0 + 4.0 * pi / 3.0) -
    a / 3.0;
  const double x_plus = root_scale * std::cos(phi / 3.0) - a / 3.0;
  const double u_minus = (x_minus - mode) * std::sqrt(
    gig_scaled_kernel(x_minus, lambda, beta, log_mode)
  );
  const double u_plus = (x_plus - mode) * std::sqrt(
    gig_scaled_kernel(x_plus, lambda, beta, log_mode)
  );

  if (!std::isfinite(u_minus) || !std::isfinite(u_plus) ||
      !(u_minus < u_plus)) {
    // The unshifted construction remains exact and avoids the cubic setup's
    // loss of range when beta is extremely close to zero.
    return sample_gig_ratio_unshifted(lambda, beta, rng);
  }

  for (int attempt = 0; attempt < 1000000; ++attempt) {
    const double u = u_minus + (u_plus - u_minus) * rng.uniform();
    const double v = rng.uniform();
    const double x = u / v + mode;
    if (std::isfinite(x) && x > 0.0 &&
        v * v <= gig_scaled_kernel(x, lambda, beta, log_mode)) {
      return x;
    }
  }
  throw std::runtime_error(
    "The native GIG shifted ratio-of-uniforms sampler did not converge."
  );
}

template <typename UniformRng>
double sample_gig_symmetric(
    const double lambda,
    const double beta,
    UniformRng& rng) {
  if (lambda < 1.0 &&
      beta < std::min(0.5, (2.0 / 3.0) * std::sqrt(1.0 - lambda))) {
    return sample_gig_small_parameter(lambda, beta, rng);
  }
  // The unshifted construction has much cheaper setup and remains efficient
  // throughout this central region. The shifted construction is reserved for
  // larger parameters, where centering materially improves acceptance.
  if (lambda <= 2.0 && beta <= 3.0) {
    return sample_gig_ratio_unshifted(lambda, beta, rng);
  }
  return sample_gig_ratio_shifted(lambda, beta, rng);
}

template <typename UniformRng>
double sample_gig(
    const double lambda,
    const double chi,
    const double psi,
    UniformRng& rng) {
  if (!std::isfinite(lambda) || !std::isfinite(chi) ||
      !std::isfinite(psi) || chi <= 0.0 || psi <= 0.0) {
    throw std::domain_error(
      "GIG parameters require finite lambda and positive finite chi and psi."
    );
  }

  const double sqrt_chi = std::sqrt(chi);
  const double sqrt_psi = std::sqrt(psi);
  const double beta = sqrt_chi * sqrt_psi;
  if (!(beta > 0.0) || !std::isfinite(beta)) {
    throw std::overflow_error(
      "The transformed native GIG scale is outside the finite double range."
    );
  }
  const bool reciprocal = lambda < 0.0;
  const double symmetric_draw = sample_gig_symmetric(
    std::abs(lambda), beta, rng
  );
  const double scale = sqrt_chi / sqrt_psi;
  double draw = reciprocal
    ? scale / symmetric_draw
    : scale * symmetric_draw;
  if (!(draw > 0.0) || !std::isfinite(draw)) {
    // Preserve range when a large scale and a compensating symmetric draw
    // would overflow or underflow if combined directly.
    const double log_scale = 0.5 * (std::log(chi) - std::log(psi));
    const double log_draw = log_scale +
      (reciprocal ? -std::log(symmetric_draw) : std::log(symmetric_draw));
    draw = std::exp(log_draw);
  }
  if (!(draw > 0.0) || !std::isfinite(draw)) {
    throw std::overflow_error(
      "The native GIG draw is outside the finite positive double range."
    );
  }
  return draw;
}

}  // namespace bayeslinreg

#endif
