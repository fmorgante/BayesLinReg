#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <GIGrvg.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

typedef SEXP (*gig_sampler_type)(int, double, double, double);

double draw_gig(const double lambda, const double chi, const double psi) {
  if (!std::isfinite(lambda) || !std::isfinite(chi) ||
      !std::isfinite(psi) || chi <= 0.0 || psi <= 0.0) {
    Rcpp::stop(
      "GIG parameters require finite lambda and positive finite chi and psi."
    );
  }

  static gig_sampler_type sampler = NULL;
  if (sampler == NULL) {
    sampler = reinterpret_cast<gig_sampler_type>(
      R_GetCCallable("GIGrvg", "do_rgig")
    );
  }

  SEXP result = PROTECT(sampler(1, lambda, chi, psi));
  const double draw = REAL(result)[0];
  UNPROTECT(1);
  if (!std::isfinite(draw) || draw <= 0.0) {
    Rcpp::stop("The GIG sampler returned a non-positive or non-finite draw.");
  }
  return draw;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector draw_gig_rcpp_cpp(
    const int n,
    const double lambda,
    const double chi,
    const double psi) {
  Rcpp::RNGScope scope;
  if (n < 1) {
    Rcpp::stop("GIG sample size must be a positive integer.");
  }
  Rcpp::NumericVector draws(n);
  for (int index = 0; index < n; ++index) {
    draws[index] = draw_gig(lambda, chi, psi);
  }
  return draws;
}

// [[Rcpp::export]]
Rcpp::List blm_gibbs_rcpp_cpp(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& X,
    const double residual_shape,
    const double residual_scale,
    const int iterations,
    const int burnin,
    const int thin,
    const Rcpp::Function& progress_callback,
    const Rcpp::IntegerVector& block_id,
    const Rcpp::IntegerVector& block_model,
    const Rcpp::NumericVector& normal_shape,
    const Rcpp::NumericVector& normal_scale,
    const Rcpp::NumericVector& pi_alpha,
    const Rcpp::NumericVector& pi_beta,
    const Rcpp::NumericVector& slab_shape,
    const Rcpp::NumericVector& slab_scale,
    const Rcpp::NumericVector& global_scale,
    const Rcpp::NumericVector& local_a,
    const Rcpp::NumericVector& local_b,
    const bool learn_residual_var,
    const double fixed_residual_var) {
  Rcpp::RNGScope scope;

  const int n = y.size();
  const int p = X.ncol();
  const int number_of_blocks = block_model.size();
  const int number_of_draws = (iterations - burnin - 1) / thin + 1;

  std::vector<double> x_mean(p, 0.0);
  double y_mean = 0.0;
  for (int i = 0; i < n; ++i) {
    y_mean += y[i];
    for (int j = 0; j < p; ++j) {
      x_mean[j] += X(i, j);
    }
  }
  y_mean /= n;
  for (int j = 0; j < p; ++j) {
    x_mean[j] /= n;
  }

  Rcpp::NumericMatrix x_centered(n, p);
  Rcpp::NumericVector y_centered(n);
  for (int i = 0; i < n; ++i) {
    y_centered[i] = y[i] - y_mean;
    for (int j = 0; j < p; ++j) {
      x_centered(i, j) = X(i, j) - x_mean[j];
    }
  }

  std::vector<double> x_squared(p, 0.0);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) {
      x_squared[j] += x_centered(i, j) * x_centered(i, j);
    }
  }

  Rcpp::NumericMatrix coefficient_samples(number_of_draws, p);
  Rcpp::NumericVector intercept_samples(number_of_draws);
  Rcpp::NumericVector residual_var_samples(number_of_draws);
  Rcpp::NumericMatrix normal_var_samples(number_of_draws, number_of_blocks);
  Rcpp::IntegerMatrix inclusion_samples(number_of_draws, p);
  Rcpp::NumericMatrix pi_samples(number_of_draws, number_of_blocks);
  Rcpp::NumericMatrix slab_var_samples(number_of_draws, number_of_blocks);
  Rcpp::NumericMatrix local_var_samples(number_of_draws, p);
  Rcpp::NumericMatrix tau_sq_samples(number_of_draws, number_of_blocks);
  std::vector<double> coefficient(p, 0.0);
  std::vector<int> inclusion(p, 1);
  std::vector<double> local_var(p, 1.0);
  std::vector<double> local_aux(p, 1.0);
  std::vector<double> residuals(n);
  for (int i = 0; i < n; ++i) {
    residuals[i] = y_centered[i];
  }

  double residual_var = learn_residual_var
    ? residual_scale / (residual_shape + 1.0)
    : fixed_residual_var;
  std::vector<double> pi(number_of_blocks, 0.5);
  std::vector<double> normal_var(number_of_blocks, 1.0);
  std::vector<double> slab_var(number_of_blocks, 1.0);
  std::vector<double> tau_sq(number_of_blocks, 1.0);
  std::vector<double> global_aux(number_of_blocks, 1.0);
  bool has_normal = false;
  bool has_spike_slab = false;
  bool has_global_local = false;
  for (int block = 0; block < number_of_blocks; ++block) {
    if (block_model[block] == 0) {
      has_normal = true;
      normal_var[block] =
        normal_scale[block] / (normal_shape[block] + 1.0);
    }
    if (block_model[block] == 1) {
      has_spike_slab = true;
      pi[block] = pi_alpha[block] / (pi_alpha[block] + pi_beta[block]);
      slab_var[block] = slab_scale[block] / (slab_shape[block] + 1.0);
    }
    if (block_model[block] == 2) {
      has_global_local = true;
      tau_sq[block] = global_scale[block] * global_scale[block];
    }
  }
  const double posterior_shape =
    residual_shape + static_cast<double>(n - 1) / 2.0;
  int retained_index = 0;
  int next_progress_percent = 10;
  int last_reported_iteration = 0;

  for (int iteration = 1; iteration <= iterations; ++iteration) {
    // Update each coefficient from its univariate conditional normal.
    for (int j = 0; j < p; ++j) {
      const int block = block_id[j] - 1;
      const int model = block_model[block];
      double conditional_numerator = 0.0;
      for (int i = 0; i < n; ++i) {
        residuals[i] += x_centered(i, j) * coefficient[j];
        conditional_numerator += x_centered(i, j) * residuals[i];
      }
      const double prior_precision = model == 2
        ? 1.0 / tau_sq[block] / local_var[j]
        : (model == 1 ? 1.0 / slab_var[block] : 1.0 / normal_var[block]);
      const double conditional_var = 1.0 / (
        x_squared[j] / residual_var + prior_precision
      );
      const double conditional_mean =
        conditional_var * conditional_numerator / residual_var;
      if (model == 1) {
        const double epsilon = std::numeric_limits<double>::epsilon();
        const double bounded_pi = std::min(
          std::max(pi[block], epsilon),
          1.0 - epsilon
        );
        const double log_inclusion_odds =
          std::log(bounded_pi) - std::log1p(-bounded_pi) +
          0.5 * std::log(conditional_var / slab_var[block]) +
          conditional_mean * conditional_mean / (2.0 * conditional_var);
        const double inclusion_probability =
          log_inclusion_odds >= 0.0
            ? 1.0 / (1.0 + std::exp(-log_inclusion_odds))
            : std::exp(log_inclusion_odds) /
                (1.0 + std::exp(log_inclusion_odds));
        inclusion[j] = static_cast<int>(
          R::rbinom(1.0, inclusion_probability)
        );
      }
      if (model != 1 || inclusion[j] == 1) {
        coefficient[j] = R::rnorm(
          conditional_mean,
          std::sqrt(conditional_var)
        );
      } else {
        coefficient[j] = 0.0;
      }
      for (int i = 0; i < n; ++i) {
        residuals[i] -= x_centered(i, j) * coefficient[j];
      }
    }

    if (has_normal) {
      for (int block = 0; block < number_of_blocks; ++block) {
        if (block_model[block] != 0) {
          continue;
        }
        int block_size = 0;
        double coefficient_sum_of_squares = 0.0;
        for (int j = 0; j < p; ++j) {
          if (block_id[j] - 1 == block) {
            coefficient_sum_of_squares += coefficient[j] * coefficient[j];
            ++block_size;
          }
        }
        const double normal_posterior_scale =
          normal_scale[block] + 0.5 * coefficient_sum_of_squares;
        normal_var[block] = 1.0 / R::rgamma(
          normal_shape[block] + 0.5 * block_size,
          1.0 / normal_posterior_scale
        );
      }
    }

    if (has_spike_slab) {
      for (int block = 0; block < number_of_blocks; ++block) {
        if (block_model[block] != 1) {
          continue;
        }
        int number_included = 0;
        int block_size = 0;
        double included_sum_of_squares = 0.0;
        for (int j = 0; j < p; ++j) {
          if (block_id[j] - 1 == block) {
            number_included += inclusion[j];
            if (inclusion[j] == 1) {
              included_sum_of_squares += coefficient[j] * coefficient[j];
            }
            ++block_size;
          }
        }
        pi[block] = R::rbeta(
          pi_alpha[block] + number_included,
          pi_beta[block] + block_size - number_included
        );
        const double slab_posterior_scale =
          slab_scale[block] + 0.5 * included_sum_of_squares;
        slab_var[block] = 1.0 / R::rgamma(
          slab_shape[block] + 0.5 * number_included,
          1.0 / slab_posterior_scale
        );
      }
    }

    if (has_global_local) {
      for (int block = 0; block < number_of_blocks; ++block) {
        if (block_model[block] != 2) {
          continue;
        }
        int block_size = 0;
        for (int j = 0; j < p; ++j) {
          if (block_id[j] - 1 != block) {
            continue;
          }
          ++block_size;
          const double raw_chi =
            coefficient[j] * coefficient[j] / tau_sq[block];
          const double chi = std::max(
            raw_chi,
            std::numeric_limits<double>::min()
          );
          local_var[j] = draw_gig(
            local_a[block] - 0.5,
            chi,
            2.0 * local_aux[j]
          );
          local_aux[j] = R::rgamma(
            local_a[block] + local_b[block],
            1.0 / (1.0 + local_var[j])
          );
        }

        double tau_rate = 1.0 / global_aux[block];
        for (int j = 0; j < p; ++j) {
          if (block_id[j] - 1 == block) {
            tau_rate += coefficient[j] * coefficient[j] /
              (2.0 * local_var[j]);
          }
        }
        tau_sq[block] = 1.0 / R::rgamma(
          (static_cast<double>(block_size) + 1.0) / 2.0,
          1.0 / tau_rate
        );
        const double global_aux_rate =
          1.0 / (global_scale[block] * global_scale[block]) +
          1.0 / tau_sq[block];
        global_aux[block] =
          1.0 / R::rgamma(1.0, 1.0 / global_aux_rate);
      }
    }

    if (learn_residual_var) {
      double sum_squared_residuals = 0.0;
      for (int i = 0; i < n; ++i) {
        sum_squared_residuals += residuals[i] * residuals[i];
      }
      const double posterior_scale =
        residual_scale + 0.5 * sum_squared_residuals;
      residual_var = 1.0 / R::rgamma(
        posterior_shape,
        1.0 / posterior_scale
      );
    }

    if (iteration > burnin &&
        (iteration - burnin - 1) % thin == 0) {
      double intercept_mean = y_mean;
      for (int j = 0; j < p; ++j) {
        const int block = block_id[j] - 1;
        const int model = block_model[block];
        coefficient_samples(retained_index, j) = coefficient[j];
        if (model == 1) {
          inclusion_samples(retained_index, j) = inclusion[j];
        }
        if (model == 2) {
          local_var_samples(retained_index, j) = local_var[j];
        }
        intercept_mean -= x_mean[j] * coefficient[j];
      }
      intercept_samples[retained_index] = R::rnorm(
        intercept_mean,
        std::sqrt(residual_var / n)
      );
      residual_var_samples[retained_index] = residual_var;
      if (has_normal) {
        for (int block = 0; block < number_of_blocks; ++block) {
          if (block_model[block] == 0) {
            normal_var_samples(retained_index, block) = normal_var[block];
          }
        }
      }
      if (has_spike_slab) {
        for (int block = 0; block < number_of_blocks; ++block) {
          if (block_model[block] == 1) {
            pi_samples(retained_index, block) = pi[block];
            slab_var_samples(retained_index, block) = slab_var[block];
          }
        }
      }
      if (has_global_local) {
        for (int block = 0; block < number_of_blocks; ++block) {
          if (block_model[block] == 2) {
            tau_sq_samples(retained_index, block) = tau_sq[block];
          }
        }
      }
      ++retained_index;
    }

    if (next_progress_percent <= 100) {
      int threshold = static_cast<int>(std::ceil(
        iterations * next_progress_percent / 100.0
      ));
      if (iteration >= threshold) {
        do {
          next_progress_percent += 10;
          if (next_progress_percent > 100) {
            break;
          }
          threshold = static_cast<int>(std::ceil(
            iterations * next_progress_percent / 100.0
          ));
        } while (iteration >= threshold);
        progress_callback(
          iteration - last_reported_iteration,
          iteration
        );
        last_reported_iteration = iteration;
      }
    }

    if (iteration % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("coefficient_samples") = coefficient_samples,
    Rcpp::Named("intercept_samples") = intercept_samples,
    Rcpp::Named("residual_var_samples") = residual_var_samples,
    Rcpp::Named("normal_var_samples") = normal_var_samples,
    Rcpp::Named("inclusion_samples") = inclusion_samples,
    Rcpp::Named("pi_samples") = pi_samples,
    Rcpp::Named("slab_var_samples") = slab_var_samples,
    Rcpp::Named("local_var_samples") = local_var_samples,
    Rcpp::Named("tau_sq_samples") = tau_sq_samples
  );
}
