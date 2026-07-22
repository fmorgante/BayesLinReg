#include <RcppEigen.h>
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

class DenseSummaryMatrix {
 public:
  explicit DenseSummaryMatrix(const Rcpp::NumericMatrix& matrix)
    : matrix_(matrix) {}

  int cols() const { return matrix_.ncol(); }

  double diagonal(const int j) const { return matrix_(j, j); }

  double corrected_value(
      const std::vector<double>& rhs,
      const int j,
      const double) const {
    return rhs[j];
  }

  void update(
      std::vector<double>& rhs,
      const int j,
      const double change) const {
    for (int k = 0; k < matrix_.nrow(); ++k) {
      rhs[k] -= matrix_(k, j) * change;
    }
  }

  void multiply(
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const {
    const int p = matrix_.ncol();
    std::fill(fitted.begin(), fitted.end(), 0.0);
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < p; ++k) {
        fitted[j] += matrix_(j, k) * coefficient[k];
      }
    }
  }

  double center_dot(const std::vector<double>&) const { return 0.0; }

  double centered_fitted(
      const double fitted,
      const int,
      const double) const {
    return fitted;
  }

  void update_center_dot(double&, const int, const double) const {}

 private:
  const Rcpp::NumericMatrix& matrix_;
};

class SparseSummaryMatrix {
 public:
  SparseSummaryMatrix(
      const Eigen::MappedSparseMatrix<double>& matrix,
      const Rcpp::NumericVector& center)
    : matrix_(matrix), center_(center) {
    if (center_.size() != matrix_.cols()) {
      Rcpp::stop("Sparse centering vector does not match `XtX`.");
    }
  }

  int cols() const { return matrix_.cols(); }

  double diagonal(const int j) const {
    return matrix_.coeff(j, j) - center_[j] * center_[j];
  }

  double corrected_value(
      const std::vector<double>& rhs,
      const int j,
      const double center_dot) const {
    return rhs[j] + center_[j] * center_dot;
  }

  void update(
      std::vector<double>& rhs,
      const int j,
      const double change) const {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator entry(matrix_, j);
         entry; ++entry) {
      rhs[entry.row()] -= entry.value() * change;
    }
  }

  void multiply(
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const {
    Eigen::Map<const Eigen::VectorXd> coefficient_map(
      coefficient.data(), coefficient.size()
    );
    Eigen::Map<Eigen::VectorXd> fitted_map(fitted.data(), fitted.size());
    fitted_map.noalias() = matrix_ * coefficient_map;
  }

  double center_dot(const std::vector<double>& coefficient) const {
    double result = 0.0;
    for (int j = 0; j < center_.size(); ++j) {
      result += center_[j] * coefficient[j];
    }
    return result;
  }

  double centered_fitted(
      const double fitted,
      const int j,
      const double center_dot) const {
    return fitted - center_[j] * center_dot;
  }

  void update_center_dot(
      double& center_dot,
      const int j,
      const double change) const {
    center_dot += center_[j] * change;
  }

 private:
  const Eigen::MappedSparseMatrix<double>& matrix_;
  const Rcpp::NumericVector& center_;
};

}  // namespace

// [[Rcpp::depends(RcppEigen)]]

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

template <typename SummaryMatrix>
Rcpp::List blm_gibbs_core(
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
    const Rcpp::NumericVector& spike_var_shape,
    const Rcpp::NumericVector& spike_var_scale,
    const Rcpp::NumericVector& global_scale,
    const Rcpp::NumericVector& local_a,
    const Rcpp::NumericVector& local_b,
    const Rcpp::List& multi_gamma_list,
    const Rcpp::List& multi_pi_alpha_list,
    const Rcpp::NumericVector& multi_var_shape,
    const Rcpp::NumericVector& multi_var_scale,
    const bool learn_residual_var,
    const double fixed_residual_var,
    const bool store_samples,
    const bool store_coefficient_cov,
    const int effective_n,
    const bool fit_intercept,
    const Rcpp::NumericVector& intercept_x_mean,
    const double intercept_y_mean,
    const bool use_sufficient_statistics,
    const SummaryMatrix& summary_XtX,
    const Rcpp::NumericVector& summary_Xty,
    const double summary_yty) {
  Rcpp::RNGScope scope;

  const int n = y.size();
  const int p = use_sufficient_statistics ? summary_XtX.cols() : X.ncol();
  const int number_of_blocks = block_model.size();
  const int number_of_draws = (iterations - burnin - 1) / thin + 1;

  std::vector<double> x_mean(p, 0.0);
  double y_mean = 0.0;
  if (!use_sufficient_statistics) {
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
  }

  Rcpp::NumericMatrix x_centered(n, p);
  Rcpp::NumericVector y_centered(n);
  if (!use_sufficient_statistics) {
    for (int i = 0; i < n; ++i) {
      y_centered[i] = y[i] - y_mean;
      for (int j = 0; j < p; ++j) {
        x_centered(i, j) = X(i, j) - x_mean[j];
      }
    }
  }

  std::vector<double> x_squared(p, 0.0);
  for (int j = 0; j < p; ++j) {
    if (use_sufficient_statistics) {
      x_squared[j] = summary_XtX.diagonal(j);
    } else {
      for (int i = 0; i < n; ++i) {
        x_squared[j] += x_centered(i, j) * x_centered(i, j);
      }
    }
  }

  const int stored_rows = store_samples ? number_of_draws : 0;
  Rcpp::NumericMatrix coefficient_samples(stored_rows, p);
  Rcpp::NumericVector intercept_samples(stored_rows);
  Rcpp::NumericVector residual_var_samples(stored_rows);
  Rcpp::NumericMatrix normal_var_samples(stored_rows, number_of_blocks);
  Rcpp::IntegerMatrix inclusion_samples(stored_rows, p);
  Rcpp::NumericMatrix pi_samples(stored_rows, number_of_blocks);
  Rcpp::NumericMatrix slab_var_samples(stored_rows, number_of_blocks);
  Rcpp::NumericMatrix local_var_samples(stored_rows, p);
  Rcpp::NumericMatrix tau_sq_samples(stored_rows, number_of_blocks);
  Rcpp::IntegerMatrix multi_component_samples(stored_rows, p);
  Rcpp::List multi_pi_samples(number_of_blocks);
  Rcpp::NumericMatrix multi_var_samples(stored_rows, number_of_blocks);
  Rcpp::NumericVector coefficient_sum(p);
  Rcpp::NumericVector coefficient_sum_sq(p);
  const int covariance_dimension = store_coefficient_cov ? p : 0;
  Rcpp::NumericMatrix coefficient_crossprod(
    covariance_dimension, covariance_dimension
  );
  double intercept_sum = 0.0;
  double intercept_sum_sq = 0.0;
  double residual_var_sum = 0.0;
  double residual_var_sum_sq = 0.0;
  Rcpp::NumericVector normal_var_sum(number_of_blocks);
  Rcpp::NumericVector normal_var_sum_sq(number_of_blocks);
  Rcpp::NumericVector inclusion_sum(p);
  Rcpp::NumericVector pi_sum(number_of_blocks);
  Rcpp::NumericVector pi_sum_sq(number_of_blocks);
  Rcpp::NumericVector slab_var_sum(number_of_blocks);
  Rcpp::NumericVector slab_var_sum_sq(number_of_blocks);
  Rcpp::NumericVector local_var_sum(p);
  Rcpp::NumericVector local_var_sum_sq(p);
  Rcpp::NumericVector tau_sq_sum(number_of_blocks);
  Rcpp::NumericVector tau_sq_sum_sq(number_of_blocks);
  Rcpp::List multi_component_sum(number_of_blocks);
  Rcpp::List multi_pi_sum(number_of_blocks);
  Rcpp::List multi_pi_sum_sq(number_of_blocks);
  Rcpp::NumericVector multi_var_sum(number_of_blocks);
  Rcpp::NumericVector multi_var_sum_sq(number_of_blocks);
  std::vector<double> coefficient(p, 0.0);
  std::vector<int> inclusion(p, 1);
  std::vector<int> multi_component(p, 0);
  std::vector<int> multi_local_index(p, -1);
  std::vector<double> local_var(p, 1.0);
  std::vector<double> local_aux(p, 1.0);
  std::vector<double> residuals(n);
  std::vector<double> corrected_rhs(p, 0.0);
  std::vector<double> fitted_crossproduct(p, 0.0);
  double center_dot = 0.0;
  double residual_sse = summary_yty;
  if (use_sufficient_statistics) {
    for (int j = 0; j < p; ++j) {
      corrected_rhs[j] = summary_Xty[j];
    }
  } else {
    for (int i = 0; i < n; ++i) {
      residuals[i] = y_centered[i];
    }
  }

  double residual_var = learn_residual_var
    ? residual_scale / (residual_shape + 1.0)
    : fixed_residual_var;
  std::vector<double> pi(number_of_blocks, 0.5);
  std::vector<double> normal_var(number_of_blocks, 1.0);
  std::vector<double> slab_var(number_of_blocks, 1.0);
  std::vector<double> tau_sq(number_of_blocks, 1.0);
  std::vector<double> global_aux(number_of_blocks, 1.0);
  std::vector<double> multi_var(number_of_blocks, 1.0);
  std::vector< std::vector<double> > multi_gamma(number_of_blocks);
  std::vector< std::vector<double> > multi_pi_alpha(number_of_blocks);
  std::vector< std::vector<double> > multi_pi(number_of_blocks);
  bool has_normal = false;
  bool has_spike_slab = false;
  bool has_global_local = false;
  bool has_spike_multi_slab = false;
  for (int block = 0; block < number_of_blocks; ++block) {
    if (block_model[block] == 0) {
      has_normal = true;
      normal_var[block] =
        normal_scale[block] / (normal_shape[block] + 1.0);
    }
    if (block_model[block] == 1) {
      has_spike_slab = true;
      pi[block] = pi_alpha[block] / (pi_alpha[block] + pi_beta[block]);
      slab_var[block] =
        spike_var_scale[block] / (spike_var_shape[block] + 1.0);
    }
    if (block_model[block] == 2) {
      has_global_local = true;
      tau_sq[block] = global_scale[block] * global_scale[block];
    }
    if (block_model[block] == 3) {
      has_spike_multi_slab = true;
      const Rcpp::NumericVector gamma_values = multi_gamma_list[block];
      const Rcpp::NumericVector alpha_values = multi_pi_alpha_list[block];
      multi_gamma[block] = Rcpp::as< std::vector<double> >(gamma_values);
      multi_pi_alpha[block] =
        Rcpp::as< std::vector<double> >(alpha_values);
      multi_pi[block].resize(alpha_values.size());
      double alpha_total = 0.0;
      for (int component = 0; component < alpha_values.size(); ++component) {
        alpha_total += alpha_values[component];
      }
      for (int component = 0; component < alpha_values.size(); ++component) {
        multi_pi[block][component] = alpha_values[component] / alpha_total;
      }
      multi_var[block] =
        multi_var_scale[block] / (multi_var_shape[block] + 1.0);
      if (store_samples) {
        multi_pi_samples[block] = Rcpp::NumericMatrix(
          stored_rows, alpha_values.size()
        );
      } else {
        int block_size = 0;
        for (int j = 0; j < p; ++j) {
          if (block_id[j] - 1 == block) {
            multi_local_index[j] = block_size;
            ++block_size;
          }
        }
        multi_component_sum[block] = Rcpp::NumericMatrix(
          block_size, alpha_values.size()
        );
        multi_pi_sum[block] = Rcpp::NumericVector(alpha_values.size());
        multi_pi_sum_sq[block] = Rcpp::NumericVector(alpha_values.size());
      }
    } else {
      multi_pi_samples[block] = R_NilValue;
      multi_component_sum[block] = R_NilValue;
      multi_pi_sum[block] = R_NilValue;
      multi_pi_sum_sq[block] = R_NilValue;
    }
  }
  const double posterior_shape =
    residual_shape +
      static_cast<double>(effective_n - (fit_intercept ? 1 : 0)) / 2.0;
  int retained_index = 0;
  int next_progress_percent = 10;
  int last_reported_iteration = 0;

  for (int iteration = 1; iteration <= iterations; ++iteration) {
    // Update each coefficient from its univariate conditional normal.
    for (int j = 0; j < p; ++j) {
      const int block = block_id[j] - 1;
      const int model = block_model[block];
      const double old_coefficient = coefficient[j];
      double conditional_numerator = 0.0;
      double partial_rhs = 0.0;
      if (use_sufficient_statistics) {
        partial_rhs = summary_XtX.corrected_value(
          corrected_rhs, j, center_dot
        ) + x_squared[j] * old_coefficient;
        conditional_numerator = partial_rhs;
      } else {
        for (int i = 0; i < n; ++i) {
          residuals[i] += x_centered(i, j) * old_coefficient;
          conditional_numerator += x_centered(i, j) * residuals[i];
        }
      }
      if (model == 3) {
        const int component_count = multi_gamma[block].size();
        std::vector<double> log_weights(component_count, 0.0);
        std::vector<double> conditional_vars(component_count, 0.0);
        std::vector<double> conditional_means(component_count, 0.0);
        double maximum_log_weight = -std::numeric_limits<double>::infinity();
        for (int component = 0; component < component_count; ++component) {
          log_weights[component] = std::log(std::max(
            multi_pi[block][component], std::numeric_limits<double>::min()
          ));
          if (component > 0) {
            const double prior_var =
              multi_gamma[block][component] * multi_var[block];
            conditional_vars[component] = 1.0 / (
              x_squared[j] / residual_var + 1.0 / prior_var
            );
            conditional_means[component] = conditional_vars[component] *
              conditional_numerator / residual_var;
            log_weights[component] +=
              0.5 * std::log(conditional_vars[component] / prior_var) +
              conditional_means[component] * conditional_means[component] /
                (2.0 * conditional_vars[component]);
          }
          maximum_log_weight = std::max(
            maximum_log_weight, log_weights[component]
          );
        }
        double weight_total = 0.0;
        for (int component = 0; component < component_count; ++component) {
          log_weights[component] = std::exp(
            log_weights[component] - maximum_log_weight
          );
          weight_total += log_weights[component];
        }
        const double threshold = R::runif(0.0, weight_total);
        double cumulative_weight = 0.0;
        int selected_component = component_count - 1;
        for (int component = 0; component < component_count; ++component) {
          cumulative_weight += log_weights[component];
          if (threshold <= cumulative_weight) {
            selected_component = component;
            break;
          }
        }
        multi_component[j] = selected_component;
        coefficient[j] = selected_component == 0
          ? 0.0
          : R::rnorm(
              conditional_means[selected_component],
              std::sqrt(conditional_vars[selected_component])
            );
      } else {
        const double prior_precision = model == 2
          ? 1.0 / tau_sq[block] / local_var[j]
          : (model == 1
              ? 1.0 / slab_var[block]
              : 1.0 / normal_var[block]);
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
      }
      if (use_sufficient_statistics) {
        const double coefficient_change = coefficient[j] - old_coefficient;
        if (coefficient_change != 0.0) {
          summary_XtX.update(corrected_rhs, j, coefficient_change);
          summary_XtX.update_center_dot(center_dot, j, coefficient_change);
          if (learn_residual_var) {
            residual_sse += -2.0 * coefficient_change * partial_rhs +
              (coefficient[j] * coefficient[j] -
               old_coefficient * old_coefficient) * x_squared[j];
          }
        }
      } else {
        for (int i = 0; i < n; ++i) {
          residuals[i] -= x_centered(i, j) * coefficient[j];
        }
      }
    }

    if (use_sufficient_statistics && iteration % 100 == 0) {
      double quadratic = 0.0;
      double linear = 0.0;
      summary_XtX.multiply(coefficient, fitted_crossproduct);
      center_dot = summary_XtX.center_dot(coefficient);
      for (int j = 0; j < p; ++j) {
        corrected_rhs[j] = summary_Xty[j] - fitted_crossproduct[j];
        linear += coefficient[j] * summary_Xty[j];
        quadratic += coefficient[j] * summary_XtX.centered_fitted(
          fitted_crossproduct[j], j, center_dot
        );
      }
      if (learn_residual_var) {
        residual_sse = summary_yty - 2.0 * linear + quadratic;
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
          spike_var_scale[block] + 0.5 * included_sum_of_squares;
        slab_var[block] = 1.0 / R::rgamma(
          spike_var_shape[block] + 0.5 * number_included,
          1.0 / slab_posterior_scale
        );
      }
    }

    if (has_spike_multi_slab) {
      for (int block = 0; block < number_of_blocks; ++block) {
        if (block_model[block] != 3) {
          continue;
        }
        const int component_count = multi_gamma[block].size();
        std::vector<int> counts(component_count, 0);
        int number_nonzero = 0;
        double scaled_sum_of_squares = 0.0;
        for (int j = 0; j < p; ++j) {
          if (block_id[j] - 1 != block) {
            continue;
          }
          const int component = multi_component[j];
          ++counts[component];
          if (component > 0) {
            ++number_nonzero;
            scaled_sum_of_squares += coefficient[j] * coefficient[j] /
              multi_gamma[block][component];
          }
        }
        double pi_total = 0.0;
        for (int component = 0; component < component_count; ++component) {
          multi_pi[block][component] = R::rgamma(
            multi_pi_alpha[block][component] + counts[component], 1.0
          );
          pi_total += multi_pi[block][component];
        }
        for (int component = 0; component < component_count; ++component) {
          multi_pi[block][component] /= pi_total;
        }
        const double posterior_scale = multi_var_scale[block] +
          0.5 * scaled_sum_of_squares;
        multi_var[block] = 1.0 / R::rgamma(
          multi_var_shape[block] + 0.5 * number_nonzero,
          1.0 / posterior_scale
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
      double sum_squared_residuals = residual_sse;
      if (!use_sufficient_statistics) {
        sum_squared_residuals = 0.0;
        for (int i = 0; i < n; ++i) {
          sum_squared_residuals += residuals[i] * residuals[i];
        }
      }
      const double posterior_scale =
        residual_scale + 0.5 * std::max(0.0, sum_squared_residuals);
      residual_var = 1.0 / R::rgamma(
        posterior_shape,
        1.0 / posterior_scale
      );
    }

    if (iteration > burnin &&
        (iteration - burnin - 1) % thin == 0) {
      double intercept_mean = intercept_y_mean;
      for (int j = 0; j < p; ++j) {
        const int block = block_id[j] - 1;
        const int model = block_model[block];
        if (store_samples) {
          coefficient_samples(retained_index, j) = coefficient[j];
          if (model == 1) {
            inclusion_samples(retained_index, j) = inclusion[j];
          }
          if (model == 2) {
            local_var_samples(retained_index, j) = local_var[j];
          }
          if (model == 3) {
            multi_component_samples(retained_index, j) =
              multi_component[j] + 1;
          }
        }
        intercept_mean -= intercept_x_mean[j] * coefficient[j];
      }
      const double intercept_draw = fit_intercept
        ? R::rnorm(
            intercept_mean,
            std::sqrt(residual_var / effective_n)
          )
        : 0.0;
      if (store_samples) {
        intercept_samples[retained_index] = intercept_draw;
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
        if (has_spike_multi_slab) {
          for (int block = 0; block < number_of_blocks; ++block) {
            if (block_model[block] != 3) {
              continue;
            }
            Rcpp::NumericMatrix block_pi_samples = multi_pi_samples[block];
            for (int component = 0;
                 component < static_cast<int>(multi_pi[block].size());
                 ++component) {
              block_pi_samples(retained_index, component) =
                multi_pi[block][component];
            }
            multi_var_samples(retained_index, block) = multi_var[block];
          }
        }
      } else {
        for (int j = 0; j < p; ++j) {
          coefficient_sum[j] += coefficient[j];
          coefficient_sum_sq[j] += coefficient[j] * coefficient[j];
          inclusion_sum[j] += inclusion[j];
          local_var_sum[j] += local_var[j];
          local_var_sum_sq[j] += local_var[j] * local_var[j];
          if (block_model[block_id[j] - 1] == 3) {
            const int block = block_id[j] - 1;
            Rcpp::NumericMatrix block_component_sum =
              multi_component_sum[block];
            block_component_sum(
              multi_local_index[j], multi_component[j]
            ) += 1.0;
          }
          if (store_coefficient_cov) {
            for (int k = 0; k < p; ++k) {
              coefficient_crossprod(j, k) += coefficient[j] * coefficient[k];
            }
          }
        }
        intercept_sum += intercept_draw;
        intercept_sum_sq += intercept_draw * intercept_draw;
        residual_var_sum += residual_var;
        residual_var_sum_sq += residual_var * residual_var;
        for (int block = 0; block < number_of_blocks; ++block) {
          normal_var_sum[block] += normal_var[block];
          normal_var_sum_sq[block] += normal_var[block] * normal_var[block];
          pi_sum[block] += pi[block];
          pi_sum_sq[block] += pi[block] * pi[block];
          slab_var_sum[block] += slab_var[block];
          slab_var_sum_sq[block] += slab_var[block] * slab_var[block];
          tau_sq_sum[block] += tau_sq[block];
          tau_sq_sum_sq[block] += tau_sq[block] * tau_sq[block];
          if (block_model[block] == 3) {
            Rcpp::NumericVector block_pi_sum = multi_pi_sum[block];
            Rcpp::NumericVector block_pi_sum_sq = multi_pi_sum_sq[block];
            for (int component = 0;
                 component < static_cast<int>(multi_pi[block].size());
                 ++component) {
              block_pi_sum[component] += multi_pi[block][component];
              block_pi_sum_sq[component] +=
                multi_pi[block][component] * multi_pi[block][component];
            }
            multi_var_sum[block] += multi_var[block];
            multi_var_sum_sq[block] += multi_var[block] * multi_var[block];
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

  if (store_samples) {
    return Rcpp::List::create(
      Rcpp::Named("coefficient_samples") = coefficient_samples,
      Rcpp::Named("intercept_samples") = intercept_samples,
      Rcpp::Named("residual_var_samples") = residual_var_samples,
      Rcpp::Named("normal_var_samples") = normal_var_samples,
      Rcpp::Named("inclusion_samples") = inclusion_samples,
      Rcpp::Named("pi_samples") = pi_samples,
      Rcpp::Named("slab_var_samples") = slab_var_samples,
      Rcpp::Named("local_var_samples") = local_var_samples,
      Rcpp::Named("tau_sq_samples") = tau_sq_samples,
      Rcpp::Named("multi_component_samples") = multi_component_samples,
      Rcpp::Named("multi_pi_samples") = multi_pi_samples,
      Rcpp::Named("multi_var_samples") = multi_var_samples
    );
  }
  Rcpp::List summaries = Rcpp::List::create(
    Rcpp::Named("number_of_draws") = number_of_draws,
    Rcpp::Named("coefficient_sum") = coefficient_sum,
    Rcpp::Named("coefficient_sum_sq") = coefficient_sum_sq,
    Rcpp::Named("intercept_sum") = intercept_sum,
    Rcpp::Named("intercept_sum_sq") = intercept_sum_sq,
    Rcpp::Named("residual_var_sum") = residual_var_sum,
    Rcpp::Named("residual_var_sum_sq") = residual_var_sum_sq,
    Rcpp::Named("normal_var_sum") = normal_var_sum,
    Rcpp::Named("normal_var_sum_sq") = normal_var_sum_sq,
    Rcpp::Named("inclusion_sum") = inclusion_sum,
    Rcpp::Named("pi_sum") = pi_sum,
    Rcpp::Named("pi_sum_sq") = pi_sum_sq,
    Rcpp::Named("slab_var_sum") = slab_var_sum,
    Rcpp::Named("slab_var_sum_sq") = slab_var_sum_sq,
    Rcpp::Named("local_var_sum") = local_var_sum,
    Rcpp::Named("local_var_sum_sq") = local_var_sum_sq,
    Rcpp::Named("tau_sq_sum") = tau_sq_sum,
    Rcpp::Named("tau_sq_sum_sq") = tau_sq_sum_sq,
    Rcpp::Named("multi_component_sum") = multi_component_sum,
    Rcpp::Named("multi_pi_sum") = multi_pi_sum,
    Rcpp::Named("multi_pi_sum_sq") = multi_pi_sum_sq,
    Rcpp::Named("multi_var_sum") = multi_var_sum,
    Rcpp::Named("multi_var_sum_sq") = multi_var_sum_sq
  );
  if (store_coefficient_cov) {
    summaries["coefficient_crossprod"] = coefficient_crossprod;
  }
  return summaries;
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
    const Rcpp::NumericVector& spike_var_shape,
    const Rcpp::NumericVector& spike_var_scale,
    const Rcpp::NumericVector& global_scale,
    const Rcpp::NumericVector& local_a,
    const Rcpp::NumericVector& local_b,
    const Rcpp::List& multi_gamma_list,
    const Rcpp::List& multi_pi_alpha_list,
    const Rcpp::NumericVector& multi_var_shape,
    const Rcpp::NumericVector& multi_var_scale,
    const bool learn_residual_var,
    const double fixed_residual_var,
    const bool store_samples,
    const bool store_coefficient_cov,
    const int effective_n,
    const bool fit_intercept,
    const Rcpp::NumericVector& intercept_x_mean,
    const double intercept_y_mean,
    const bool use_sufficient_statistics,
    const Rcpp::NumericMatrix& summary_XtX,
    const Rcpp::NumericVector& summary_Xty,
    const double summary_yty) {
  const DenseSummaryMatrix summary_matrix(summary_XtX);
  return blm_gibbs_core(
    y, X, residual_shape, residual_scale, iterations, burnin, thin,
    progress_callback, block_id, block_model, normal_shape, normal_scale,
    pi_alpha, pi_beta, spike_var_shape, spike_var_scale, global_scale,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, fit_intercept, intercept_x_mean,
    intercept_y_mean, use_sufficient_statistics, summary_matrix, summary_Xty,
    summary_yty
  );
}

// [[Rcpp::export]]
Rcpp::List blm_gibbs_sparse_rcpp_cpp(
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
    const Rcpp::NumericVector& spike_var_shape,
    const Rcpp::NumericVector& spike_var_scale,
    const Rcpp::NumericVector& global_scale,
    const Rcpp::NumericVector& local_a,
    const Rcpp::NumericVector& local_b,
    const Rcpp::List& multi_gamma_list,
    const Rcpp::List& multi_pi_alpha_list,
    const Rcpp::NumericVector& multi_var_shape,
    const Rcpp::NumericVector& multi_var_scale,
    const bool learn_residual_var,
    const double fixed_residual_var,
    const bool store_samples,
    const bool store_coefficient_cov,
    const int effective_n,
    const bool fit_intercept,
    const Rcpp::NumericVector& intercept_x_mean,
    const double intercept_y_mean,
    const Eigen::MappedSparseMatrix<double>& summary_XtX,
    const Rcpp::NumericVector& summary_center,
    const Rcpp::NumericVector& summary_Xty,
    const double summary_yty) {
  const Rcpp::NumericVector empty_y;
  const Rcpp::NumericMatrix empty_X(0, 0);
  const SparseSummaryMatrix summary_matrix(summary_XtX, summary_center);
  return blm_gibbs_core(
    empty_y, empty_X, residual_shape, residual_scale, iterations, burnin, thin,
    progress_callback, block_id, block_model, normal_shape, normal_scale,
    pi_alpha, pi_beta, spike_var_shape, spike_var_scale, global_scale,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, fit_intercept, intercept_x_mean,
    intercept_y_mean, true, summary_matrix, summary_Xty, summary_yty
  );
}
