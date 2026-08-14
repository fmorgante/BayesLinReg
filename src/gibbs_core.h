#ifndef BAYESLINREG_GIBBS_CORE_H
#define BAYESLINREG_GIBBS_CORE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include "parallel_blocks.h"

namespace bayeslinreg {

double draw_gig(double lambda, double chi, double psi);

inline void update_online_moment(
    const double value,
    const int count,
    double& mean,
    double& m2) {
  const double delta = value - mean;
  mean += delta / static_cast<double>(count);
  m2 += delta * (value - mean);
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
    const Rcpp::NumericVector& fixed_global_var,
    const Rcpp::NumericVector& local_a,
    const Rcpp::NumericVector& local_b,
    const Rcpp::List& multi_gamma_list,
    const Rcpp::List& multi_pi_alpha_list,
    const Rcpp::NumericVector& multi_var_shape,
    const Rcpp::NumericVector& multi_var_scale,
    const Rcpp::NumericVector& fixed_var,
    const bool learn_residual_var,
    const double fixed_residual_var,
    const bool store_samples,
    const bool store_coefficient_cov,
    const int effective_n,
    const int likelihood_df,
    const bool fit_intercept,
    const Rcpp::NumericVector& intercept_x_mean,
    const double intercept_y_mean,
    const bool use_sufficient_statistics,
    const SummaryMatrix& summary_XtX,
    const Rcpp::NumericVector& summary_Xty,
    const double summary_yty,
    const bool center_observations,
    const double residual_sse_offset,
    const bool compute_pve,
    const int pve_type_code,
    const int nthreads) {
  Rcpp::RNGScope scope;

  const int n = y.size();
  const int p = use_sufficient_statistics ? summary_XtX.cols() : X.ncol();
  const int number_of_blocks = block_model.size();
  const int number_of_draws = (iterations - burnin - 1) / thin + 1;
  const PveType pve_type = bayeslinreg::pve_type_from_code(pve_type_code);
  std::vector< std::vector<int> > block_predictors(number_of_blocks);
  std::vector<PriorModel> prior_models(number_of_blocks);
  std::vector<int> model_local_index(p, -1);
  std::vector<int> block_local_index(p, -1);
  std::array<int, bayeslinreg::prior_model_count> model_size = {};
  bool has_normal = false;
  bool has_spike_slab = false;
  bool has_global_local = false;
  bool has_spike_multi_slab = false;
  for (int block = 0; block < number_of_blocks; ++block) {
    const PriorModel model =
      bayeslinreg::prior_model_from_code(block_model[block]);
    prior_models[block] = model;
    has_normal = has_normal || model == PriorModel::Normal;
    has_spike_slab = has_spike_slab || model == PriorModel::SpikeSlab;
    has_global_local = has_global_local || model == PriorModel::GlobalLocal;
    has_spike_multi_slab =
      has_spike_multi_slab || model == PriorModel::SpikeMultiSlab;
  }
  for (int j = 0; j < p; ++j) {
    const int block = block_id[j] - 1;
    if (block < 0 || block >= number_of_blocks) {
      Rcpp::stop("Coefficient block identifiers are out of range.");
    }
    const PriorModel model = prior_models[block];
    block_local_index[j] =
      static_cast<int>(block_predictors[block].size());
    block_predictors[block].push_back(j);
    const int model_index = bayeslinreg::prior_model_index(model);
    model_local_index[j] = model_size[model_index];
    ++model_size[model_index];
  }

  std::vector<double> x_mean(p, 0.0);
  double y_mean = 0.0;
  if (!use_sufficient_statistics && center_observations) {
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

  Rcpp::NumericMatrix x_centered(
    center_observations ? n : 0,
    center_observations ? p : 0
  );
  Rcpp::NumericVector y_centered(n);
  if (!use_sufficient_statistics) {
    for (int i = 0; i < n; ++i) {
      y_centered[i] = center_observations ? y[i] - y_mean : y[i];
      if (center_observations) {
        for (int j = 0; j < p; ++j) {
          x_centered(i, j) = X(i, j) - x_mean[j];
        }
      }
    }
  }

  std::vector<double> x_squared(p, 0.0);
  for (int j = 0; j < p; ++j) {
    if (use_sufficient_statistics) {
      x_squared[j] = summary_XtX.diagonal(j);
    } else {
      const double* x_column = center_observations
        ? x_centered.begin() + static_cast<std::size_t>(n) * j
        : X.begin() + static_cast<std::size_t>(n) * j;
      for (int i = 0; i < n; ++i) {
        x_squared[j] += x_column[i] * x_column[i];
      }
    }
  }

  const int stored_rows = store_samples ? number_of_draws : 0;
  // Allocate posterior storage only for prior families used by this fit.
  Rcpp::NumericMatrix coefficient_samples(stored_rows, p);
  Rcpp::NumericVector intercept_samples(stored_rows);
  Rcpp::NumericVector residual_var_samples(stored_rows);
  Rcpp::NumericMatrix block_pve_samples(
    compute_pve ? stored_rows : 0,
    compute_pve ? number_of_blocks : 0
  );
  Rcpp::NumericVector total_pve_samples(
    compute_pve ? stored_rows : 0
  );
  Rcpp::NumericVector cross_block_pve_samples(
    compute_pve ? stored_rows : 0
  );
  Rcpp::NumericMatrix normal_var_samples(
    stored_rows, has_normal ? number_of_blocks : 0
  );
  Rcpp::IntegerMatrix inclusion_samples(
    stored_rows,
    model_size[bayeslinreg::prior_model_index(PriorModel::SpikeSlab)]
  );
  Rcpp::NumericMatrix pi_samples(
    stored_rows, has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericMatrix slab_var_samples(
    stored_rows, has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericMatrix local_var_samples(
    stored_rows,
    model_size[bayeslinreg::prior_model_index(PriorModel::GlobalLocal)]
  );
  Rcpp::NumericMatrix tau_sq_samples(
    stored_rows, has_global_local ? number_of_blocks : 0
  );
  Rcpp::IntegerMatrix multi_component_samples(
    stored_rows,
    model_size[bayeslinreg::prior_model_index(PriorModel::SpikeMultiSlab)]
  );
  Rcpp::List multi_pi_samples(number_of_blocks);
  Rcpp::NumericMatrix multi_var_samples(
    stored_rows, has_spike_multi_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector coefficient_mean(p);
  Rcpp::NumericVector coefficient_m2(p);
  std::vector<double> coefficient_delta(p, 0.0);
  // Only within-prior-block covariance matrices are returned. Accumulate the
  // matching block cross-products rather than a global p-by-p matrix.
  Rcpp::List coefficient_cov_m2(number_of_blocks);
  for (int block = 0; block < number_of_blocks; ++block) {
    if (!store_samples && store_coefficient_cov) {
      const int block_size = static_cast<int>(block_predictors[block].size());
      coefficient_cov_m2[block] = Rcpp::NumericMatrix(
        block_size, block_size
      );
    } else {
      coefficient_cov_m2[block] = R_NilValue;
    }
  }
  double intercept_mean = 0.0;
  double intercept_m2 = 0.0;
  double residual_var_mean = 0.0;
  double residual_var_m2 = 0.0;
  Rcpp::NumericVector block_pve_mean(
    compute_pve ? number_of_blocks : 0
  );
  Rcpp::NumericVector block_pve_m2(
    compute_pve ? number_of_blocks : 0
  );
  double total_pve_mean = 0.0;
  double total_pve_m2 = 0.0;
  double cross_block_pve_mean = 0.0;
  double cross_block_pve_m2 = 0.0;
  Rcpp::NumericVector normal_var_mean(
    has_normal ? number_of_blocks : 0
  );
  Rcpp::NumericVector normal_var_m2(
    has_normal ? number_of_blocks : 0
  );
  Rcpp::NumericVector inclusion_sum(
    model_size[bayeslinreg::prior_model_index(PriorModel::SpikeSlab)]
  );
  Rcpp::NumericVector pi_mean(
    has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector pi_m2(
    has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector slab_var_mean(
    has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector slab_var_m2(
    has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector local_var_mean(
    model_size[bayeslinreg::prior_model_index(PriorModel::GlobalLocal)]
  );
  Rcpp::NumericVector local_var_m2(
    model_size[bayeslinreg::prior_model_index(PriorModel::GlobalLocal)]
  );
  Rcpp::NumericVector tau_sq_mean(
    has_global_local ? number_of_blocks : 0
  );
  Rcpp::NumericVector tau_sq_m2(
    has_global_local ? number_of_blocks : 0
  );
  Rcpp::List multi_component_sum(number_of_blocks);
  Rcpp::List multi_pi_mean(number_of_blocks);
  Rcpp::List multi_pi_m2(number_of_blocks);
  Rcpp::NumericVector multi_var_mean(
    has_spike_multi_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector multi_var_m2(
    has_spike_multi_slab ? number_of_blocks : 0
  );
  std::vector<double> coefficient(p, 0.0);
  std::vector<int> inclusion(
    model_size[bayeslinreg::prior_model_index(PriorModel::SpikeSlab)], 1
  );
  std::vector<int> multi_component(
    model_size[bayeslinreg::prior_model_index(PriorModel::SpikeMultiSlab)], 0
  );
  std::vector<double> local_var(
    model_size[bayeslinreg::prior_model_index(PriorModel::GlobalLocal)], 1.0
  );
  std::vector<double> local_aux(
    model_size[bayeslinreg::prior_model_index(PriorModel::GlobalLocal)], 1.0
  );
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
  for (int block = 0; block < number_of_blocks; ++block) {
    if (prior_models[block] == PriorModel::Normal) {
      normal_var[block] = Rcpp::NumericVector::is_na(fixed_var[block])
        ? normal_scale[block] / (normal_shape[block] + 1.0)
        : fixed_var[block];
    }
    if (prior_models[block] == PriorModel::SpikeSlab) {
      pi[block] = pi_alpha[block] / (pi_alpha[block] + pi_beta[block]);
      slab_var[block] = Rcpp::NumericVector::is_na(fixed_var[block])
        ? spike_var_scale[block] / (spike_var_shape[block] + 1.0)
        : fixed_var[block];
    }
    if (prior_models[block] == PriorModel::GlobalLocal) {
      tau_sq[block] = Rcpp::NumericVector::is_na(fixed_global_var[block])
        ? global_scale[block] * global_scale[block]
        : fixed_global_var[block];
    }
    if (prior_models[block] == PriorModel::SpikeMultiSlab) {
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
      multi_var[block] = Rcpp::NumericVector::is_na(fixed_var[block])
        ? multi_var_scale[block] / (multi_var_shape[block] + 1.0)
        : fixed_var[block];
      if (store_samples) {
        multi_pi_samples[block] = Rcpp::NumericMatrix(
          stored_rows, alpha_values.size()
        );
      } else {
        const std::vector<int>& predictors = block_predictors[block];
        multi_component_sum[block] = Rcpp::NumericMatrix(
          predictors.size(), alpha_values.size()
        );
        multi_pi_mean[block] = Rcpp::NumericVector(alpha_values.size());
        multi_pi_m2[block] = Rcpp::NumericVector(alpha_values.size());
      }
    } else {
      multi_pi_samples[block] = R_NilValue;
      multi_component_sum[block] = R_NilValue;
      multi_pi_mean[block] = R_NilValue;
      multi_pi_m2[block] = R_NilValue;
    }
  }
  const double posterior_shape =
    residual_shape + static_cast<double>(likelihood_df) / 2.0;
  int retained_index = 0;
  int next_progress_percent = 10;
  int last_reported_iteration = 0;
  const int residual_refresh_interval = 100;
  const double residual_sse_relative_tolerance =
    std::sqrt(std::numeric_limits<double>::epsilon());
  auto reconstruct_sufficient_state = [&](const bool fitted_is_current,
                                            const int iteration) {
    if (!fitted_is_current) {
      summary_XtX.multiply(coefficient, fitted_crossproduct);
      center_dot = summary_XtX.center_dot(coefficient);
      for (int j = 0; j < p; ++j) {
        corrected_rhs[j] = summary_Xty[j] - fitted_crossproduct[j];
      }
    }
    double linear = 0.0;
    double quadratic = 0.0;
    for (int j = 0; j < p; ++j) {
      linear += coefficient[j] * summary_Xty[j];
      quadratic += coefficient[j] * summary_XtX.centered_fitted(
        fitted_crossproduct[j], j, center_dot
      );
    }
    const double reconstructed_sse =
      summary_yty - 2.0 * linear + quadratic;
    const double scale = std::max({
      1.0, std::abs(summary_yty), 2.0 * std::abs(linear),
      std::abs(quadratic)
    });
    const double tolerance = residual_sse_relative_tolerance * scale;
    if (!std::isfinite(reconstructed_sse)) {
      Rcpp::stop(
        "The reconstructed residual SSE is non-finite at Gibbs iteration %d.",
        iteration
      );
    }
    if (reconstructed_sse < -tolerance) {
      Rcpp::stop(
        "The reconstructed residual SSE is materially negative at Gibbs "
        "iteration %d (SSE %.6g; tolerance %.6g). The supplied sufficient "
        "statistics are incompatible.",
        iteration, reconstructed_sse, tolerance
      );
    }
    residual_sse = std::max(0.0, reconstructed_sse);
  };
  std::vector<BlockRng> block_rng;
  RSessionRng serial_rng;
  MixtureWorkspace serial_mixture_workspace;
  std::vector<MixtureWorkspace> parallel_mixture_workspace;
  if (nthreads > 1) {
    if constexpr (is_parallel_block_matrix<SummaryMatrix>::value) {
      const std::uint64_t seed_high = static_cast<std::uint64_t>(
        R::runif(0.0, 4294967296.0)
      );
      const std::uint64_t seed_low = static_cast<std::uint64_t>(
        R::runif(0.0, 4294967296.0)
      );
      const std::uint64_t base_seed = (seed_high << 32) | seed_low;
      block_rng.reserve(summary_XtX.block_count());
      parallel_mixture_workspace.resize(summary_XtX.block_count());
      for (int block = 0; block < summary_XtX.block_count(); ++block) {
        block_rng.emplace_back(
          splitmix64(base_seed + 2 * static_cast<std::uint64_t>(block)),
          splitmix64(base_seed + 2 * static_cast<std::uint64_t>(block) + 1)
        );
      }
    }
  }

  for (int iteration = 1; iteration <= iterations; ++iteration) {
    // Update each coefficient from its univariate conditional normal.
    bool parallel_sweep = false;
    if constexpr (is_parallel_block_matrix<SummaryMatrix>::value) {
      if (nthreads > 1) {
        parallel_coefficient_sweep(
          summary_XtX, block_id, prior_models, model_local_index, x_squared,
          residual_var, normal_var, pi, slab_var, tau_sq, local_var,
          multi_gamma, multi_pi, multi_var, learn_residual_var, coefficient,
          corrected_rhs, inclusion, multi_component, block_rng,
          parallel_mixture_workspace, residual_sse, nthreads
        );
        parallel_sweep = true;
      }
    }
    if (!parallel_sweep) {
      for (int j = 0; j < p; ++j) {
        const int block = block_id[j] - 1;
        const PriorModel model = prior_models[block];
        const int local_index = model_local_index[j];
        const double old_coefficient = coefficient[j];
        double conditional_numerator = 0.0;
        double partial_rhs = 0.0;
        if (use_sufficient_statistics) {
          partial_rhs = summary_XtX.corrected_value(
            corrected_rhs, j, center_dot
          ) + x_squared[j] * old_coefficient;
          conditional_numerator = partial_rhs;
        } else {
          const double* x_column_begin = center_observations
            ? x_centered.begin() + static_cast<std::size_t>(n) * j
            : X.begin() + static_cast<std::size_t>(n) * j;
          const Eigen::Map<const Eigen::VectorXd> x_column(
            x_column_begin, n
          );
          const Eigen::Map<const Eigen::VectorXd> residual(
            residuals.data(), n
          );
          conditional_numerator = x_column.dot(residual) +
            x_squared[j] * old_coefficient;
        }
        coefficient[j] = draw_coefficient_coordinate(
          model, block, local_index, x_squared[j], residual_var,
          conditional_numerator, normal_var, pi, slab_var, tau_sq, local_var,
          multi_gamma, multi_pi, multi_var, inclusion, multi_component,
          serial_rng, serial_mixture_workspace
        );
        const double coefficient_change = coefficient[j] - old_coefficient;
        if (use_sufficient_statistics) {
          if (coefficient_change != 0.0) {
            summary_XtX.update(corrected_rhs, j, coefficient_change);
            summary_XtX.update_center_dot(center_dot, j, coefficient_change);
            if (learn_residual_var) {
              residual_sse += -2.0 * coefficient_change * partial_rhs +
                (coefficient[j] * coefficient[j] -
                 old_coefficient * old_coefficient) * x_squared[j];
            }
          }
        } else if (coefficient_change != 0.0) {
          const double* x_column_begin = center_observations
            ? x_centered.begin() + static_cast<std::size_t>(n) * j
            : X.begin() + static_cast<std::size_t>(n) * j;
          const Eigen::Map<const Eigen::VectorXd> x_column(
            x_column_begin, n
          );
          Eigen::Map<Eigen::VectorXd> residual(residuals.data(), n);
          residual.noalias() -= coefficient_change * x_column;
        }
      }
    }

    // Streaming triangular blocks update only coordinates that remain in the
    // ascending scan. Reconstruct the full state once after every sweep so it
    // is ready for retained-draw summaries and the next iteration.
    bool streaming_state_reconstructed = false;
    if constexpr (has_streaming_triangular_blocks<SummaryMatrix>::value) {
      if (summary_XtX.has_streaming_blocks()) {
        summary_XtX.multiply(coefficient, fitted_crossproduct);
        center_dot = summary_XtX.center_dot(coefficient);
        for (int j = 0; j < p; ++j) {
          corrected_rhs[j] = summary_Xty[j] - fitted_crossproduct[j];
        }
        streaming_state_reconstructed = true;
      }
    }

    // Reconstruct accumulated state periodically to limit floating-point drift.
    if (iteration % residual_refresh_interval == 0) {
      if (use_sufficient_statistics) {
        if (learn_residual_var) {
          reconstruct_sufficient_state(
            streaming_state_reconstructed, iteration
          );
        } else if (!streaming_state_reconstructed) {
          summary_XtX.multiply(coefficient, fitted_crossproduct);
          center_dot = summary_XtX.center_dot(coefficient);
          for (int j = 0; j < p; ++j) {
            corrected_rhs[j] = summary_Xty[j] - fitted_crossproduct[j];
          }
        }
      } else {
        const double* design_begin = center_observations
          ? x_centered.begin()
          : X.begin();
        const Eigen::Map<const Eigen::MatrixXd> design(design_begin, n, p);
        const Eigen::Map<const Eigen::VectorXd> coefficient_vector(
          coefficient.data(), p
        );
        const Eigen::Map<const Eigen::VectorXd> response(y_centered.begin(), n);
        Eigen::Map<Eigen::VectorXd> residual(residuals.data(), n);
        residual.noalias() = response - design * coefficient_vector;
      }
    }

    if (has_normal) {
      for (int block = 0; block < number_of_blocks; ++block) {
        if (prior_models[block] != PriorModel::Normal) {
          continue;
        }
        if (!Rcpp::NumericVector::is_na(fixed_var[block])) {
          continue;
        }
        const std::vector<int>& predictors = block_predictors[block];
        double coefficient_sum_of_squares = 0.0;
        for (std::size_t index = 0; index < predictors.size(); ++index) {
          const int j = predictors[index];
          coefficient_sum_of_squares += coefficient[j] * coefficient[j];
        }
        const double normal_posterior_scale =
          normal_scale[block] + 0.5 * coefficient_sum_of_squares;
        normal_var[block] = 1.0 / R::rgamma(
          normal_shape[block] + 0.5 * predictors.size(),
          1.0 / normal_posterior_scale
        );
      }
    }

    if (has_spike_slab) {
      for (int block = 0; block < number_of_blocks; ++block) {
        if (prior_models[block] != PriorModel::SpikeSlab) {
          continue;
        }
        const std::vector<int>& predictors = block_predictors[block];
        const bool learn_slab_var =
          Rcpp::NumericVector::is_na(fixed_var[block]);
        int number_included = 0;
        double included_sum_of_squares = 0.0;
        for (std::size_t index = 0; index < predictors.size(); ++index) {
          const int j = predictors[index];
          const int local_index = model_local_index[j];
          number_included += inclusion[local_index];
          if (learn_slab_var && inclusion[local_index] == 1) {
            included_sum_of_squares += coefficient[j] * coefficient[j];
          }
        }
        pi[block] = R::rbeta(
          pi_alpha[block] + number_included,
          pi_beta[block] + predictors.size() - number_included
        );
        if (learn_slab_var) {
          const double slab_posterior_scale =
            spike_var_scale[block] + 0.5 * included_sum_of_squares;
          slab_var[block] = 1.0 / R::rgamma(
            spike_var_shape[block] + 0.5 * number_included,
            1.0 / slab_posterior_scale
          );
        }
      }
    }

    if (has_spike_multi_slab) {
      for (int block = 0; block < number_of_blocks; ++block) {
        if (prior_models[block] != PriorModel::SpikeMultiSlab) {
          continue;
        }
        const int component_count = multi_gamma[block].size();
        const std::vector<int>& predictors = block_predictors[block];
        const bool learn_multi_var =
          Rcpp::NumericVector::is_na(fixed_var[block]);
        std::vector<int> counts(component_count, 0);
        int number_nonzero = 0;
        double scaled_sum_of_squares = 0.0;
        for (std::size_t index = 0; index < predictors.size(); ++index) {
          const int j = predictors[index];
          const int component = multi_component[model_local_index[j]];
          ++counts[component];
          if (component > 0) {
            ++number_nonzero;
            if (learn_multi_var) {
              scaled_sum_of_squares += coefficient[j] * coefficient[j] /
                multi_gamma[block][component];
            }
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
        if (learn_multi_var) {
          const double posterior_scale = multi_var_scale[block] +
            0.5 * scaled_sum_of_squares;
          multi_var[block] = 1.0 / R::rgamma(
            multi_var_shape[block] + 0.5 * number_nonzero,
            1.0 / posterior_scale
          );
        }
      }
    }

    if (has_global_local) {
      bool parallel_local_variance = false;
      if constexpr (is_parallel_block_matrix<SummaryMatrix>::value) {
        if (nthreads > 1) {
          parallel_local_variance_sweep(
            summary_XtX, block_id, prior_models, model_local_index,
            coefficient, tau_sq, local_a, local_aux, local_var, block_rng,
            nthreads
          );
          parallel_local_variance = true;
        }
      }
      for (int block = 0; block < number_of_blocks; ++block) {
        if (prior_models[block] != PriorModel::GlobalLocal) {
          continue;
        }
        const std::vector<int>& predictors = block_predictors[block];
        for (std::size_t index = 0; index < predictors.size(); ++index) {
          const int j = predictors[index];
          const int local_index = model_local_index[j];
          const double raw_chi =
            coefficient[j] * coefficient[j] / tau_sq[block];
          const double chi = std::max(
            raw_chi,
            std::numeric_limits<double>::min()
          );
          if (!parallel_local_variance) {
            local_var[local_index] = draw_gig(
              local_a[block] - 0.5,
              chi,
              2.0 * local_aux[local_index]
            );
          }
          local_aux[local_index] = R::rgamma(
            local_a[block] + local_b[block],
            1.0 / (1.0 + local_var[local_index])
          );
          local_aux[local_index] = std::max(
            local_aux[local_index], std::numeric_limits<double>::min()
          );
        }

        if (Rcpp::NumericVector::is_na(fixed_global_var[block])) {
          double tau_rate = 1.0 / global_aux[block];
          for (std::size_t index = 0; index < predictors.size(); ++index) {
            const int j = predictors[index];
            tau_rate += coefficient[j] * coefficient[j] /
              (2.0 * local_var[model_local_index[j]]);
          }
          const double tau_precision = std::max(R::rgamma(
            (static_cast<double>(predictors.size()) + 1.0) / 2.0,
            1.0 / tau_rate
          ), std::numeric_limits<double>::min());
          tau_sq[block] = 1.0 / tau_precision;
          const double global_aux_rate =
            1.0 / (global_scale[block] * global_scale[block]) +
            1.0 / tau_sq[block];
          const double global_aux_precision = std::max(
            R::rgamma(1.0, 1.0 / global_aux_rate),
            std::numeric_limits<double>::min()
          );
          global_aux[block] = 1.0 / global_aux_precision;
        }
      }
    }

    if (learn_residual_var) {
      double sum_squared_residuals = residual_sse;
      if (use_sufficient_statistics) {
        const double preliminary_tolerance = residual_sse_relative_tolerance *
          std::max(1.0, std::abs(summary_yty));
        if (residual_sse < -preliminary_tolerance) {
          reconstruct_sufficient_state(
            streaming_state_reconstructed, iteration
          );
        }
        sum_squared_residuals = residual_sse;
      } else {
        const Eigen::Map<const Eigen::VectorXd> residual(
          residuals.data(), n
        );
        sum_squared_residuals =
          residual_sse_offset + residual.squaredNorm();
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
      std::vector<double> block_pve(number_of_blocks, 0.0);
      double total_pve = 0.0;
      double cross_block_pve = 0.0;
      if (compute_pve) {
        std::vector<double> standalone_sum_squares(number_of_blocks, 0.0);
        std::vector<double> allocated_sum_squares(number_of_blocks, 0.0);
        double total_sum_squares = 0.0;
        if (use_sufficient_statistics) {
          std::vector<double> centered_fitted(p, 0.0);
          for (int j = 0; j < p; ++j) {
            centered_fitted[j] = summary_Xty[j] -
              summary_XtX.corrected_value(corrected_rhs, j, center_dot);
            total_sum_squares += coefficient[j] * centered_fitted[j];
          }
          total_sum_squares = std::max(0.0, total_sum_squares);
          for (int block = 0; block < number_of_blocks; ++block) {
            standalone_sum_squares[block] = std::max(
              0.0,
              summary_XtX.block_quadratic(
                coefficient, block_predictors[block], block_id, block
              )
            );
            const std::vector<int>& predictors = block_predictors[block];
            for (std::size_t index = 0; index < predictors.size(); ++index) {
              const int j = predictors[index];
              allocated_sum_squares[block] +=
                coefficient[j] * centered_fitted[j];
            }
          }
        } else {
          std::vector<double> total_fitted(n, 0.0);
          for (int i = 0; i < n; ++i) {
            total_fitted[i] = y_centered[i] - residuals[i];
            total_sum_squares += total_fitted[i] * total_fitted[i];
          }
          std::vector<double> block_fitted(n, 0.0);
          for (int block = 0; block < number_of_blocks; ++block) {
            std::fill(block_fitted.begin(), block_fitted.end(), 0.0);
            const std::vector<int>& predictors = block_predictors[block];
            for (std::size_t index = 0; index < predictors.size(); ++index) {
              const int j = predictors[index];
              const double* x_column = center_observations
                ? x_centered.begin() + static_cast<std::size_t>(n) * j
                : X.begin() + static_cast<std::size_t>(n) * j;
              for (int i = 0; i < n; ++i) {
                block_fitted[i] += x_column[i] * coefficient[j];
              }
            }
            for (int i = 0; i < n; ++i) {
              standalone_sum_squares[block] +=
                block_fitted[i] * block_fitted[i];
              allocated_sum_squares[block] +=
                block_fitted[i] * total_fitted[i];
            }
          }
        }
        const double variance_df = static_cast<double>(likelihood_df);
        const double total_signal_variance = total_sum_squares / variance_df;
        const double pve_denominator = total_signal_variance + residual_var;
        double standalone_total = 0.0;
        for (int block = 0; block < number_of_blocks; ++block) {
          standalone_total += standalone_sum_squares[block];
          const double block_sum_squares = pve_type == PveType::Standalone
            ? standalone_sum_squares[block]
            : allocated_sum_squares[block];
          block_pve[block] =
            block_sum_squares / variance_df / pve_denominator;
        }
        total_pve = total_signal_variance / pve_denominator;
        cross_block_pve =
          (total_sum_squares - standalone_total) /
            variance_df / pve_denominator;
      }
      double intercept_location = intercept_y_mean;
      for (int j = 0; j < p; ++j) {
        const int block = block_id[j] - 1;
        const PriorModel model = prior_models[block];
        const int local_index = model_local_index[j];
        if (store_samples) {
          coefficient_samples(retained_index, j) = coefficient[j];
          if (model == PriorModel::SpikeSlab) {
            inclusion_samples(retained_index, local_index) =
              inclusion[local_index];
          }
          if (model == PriorModel::GlobalLocal) {
            local_var_samples(retained_index, local_index) =
              local_var[local_index];
          }
          if (model == PriorModel::SpikeMultiSlab) {
            multi_component_samples(retained_index, local_index) =
              multi_component[local_index] + 1;
          }
        }
        intercept_location -= intercept_x_mean[j] * coefficient[j];
      }
      const double intercept_draw = fit_intercept
        ? R::rnorm(
            intercept_location,
            std::sqrt(residual_var / effective_n)
          )
        : 0.0;
      if (store_samples) {
        intercept_samples[retained_index] = intercept_draw;
        residual_var_samples[retained_index] = residual_var;
        if (compute_pve) {
          for (int block = 0; block < number_of_blocks; ++block) {
            block_pve_samples(retained_index, block) = block_pve[block];
          }
          total_pve_samples[retained_index] = total_pve;
          cross_block_pve_samples[retained_index] = cross_block_pve;
        }
        if (has_normal) {
          for (int block = 0; block < number_of_blocks; ++block) {
            if (prior_models[block] == PriorModel::Normal) {
              normal_var_samples(retained_index, block) = normal_var[block];
            }
          }
        }
        if (has_spike_slab) {
          for (int block = 0; block < number_of_blocks; ++block) {
            if (prior_models[block] == PriorModel::SpikeSlab) {
              pi_samples(retained_index, block) = pi[block];
              slab_var_samples(retained_index, block) = slab_var[block];
            }
          }
        }
        if (has_global_local) {
          for (int block = 0; block < number_of_blocks; ++block) {
            if (prior_models[block] == PriorModel::GlobalLocal) {
              tau_sq_samples(retained_index, block) = tau_sq[block];
            }
          }
        }
        if (has_spike_multi_slab) {
          for (int block = 0; block < number_of_blocks; ++block) {
            if (prior_models[block] != PriorModel::SpikeMultiSlab) {
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
        const int online_count = retained_index + 1;
        for (int j = 0; j < p; ++j) {
          const double delta = coefficient[j] - coefficient_mean[j];
          coefficient_delta[j] = delta;
          coefficient_mean[j] += delta / static_cast<double>(online_count);
          coefficient_m2[j] +=
            delta * (coefficient[j] - coefficient_mean[j]);
          const int block = block_id[j] - 1;
          const PriorModel model = prior_models[block];
          const int local_index = model_local_index[j];
          if (model == PriorModel::SpikeSlab) {
            inclusion_sum[local_index] += inclusion[local_index];
          } else if (model == PriorModel::GlobalLocal) {
            double mean = local_var_mean[local_index];
            double m2 = local_var_m2[local_index];
            update_online_moment(
              local_var[local_index], online_count, mean, m2
            );
            local_var_mean[local_index] = mean;
            local_var_m2[local_index] = m2;
          } else if (model == PriorModel::SpikeMultiSlab) {
            Rcpp::NumericMatrix block_component_sum =
              multi_component_sum[block];
            block_component_sum(
              block_local_index[j], multi_component[local_index]
            ) += 1.0;
          }
        }
        if (store_coefficient_cov) {
          const double covariance_factor =
            static_cast<double>(online_count - 1) / online_count;
          for (int block = 0; block < number_of_blocks; ++block) {
            const std::vector<int>& predictors = block_predictors[block];
            Rcpp::NumericMatrix covariance_m2 = coefficient_cov_m2[block];
            const int block_size = static_cast<int>(predictors.size());
            for (int column = 0; column < block_size; ++column) {
              const double column_delta =
                coefficient_delta[predictors[column]];
              double* destination = covariance_m2.begin() +
                static_cast<std::size_t>(block_size) * column;
              for (int row = 0; row < block_size; ++row) {
                destination[row] += covariance_factor *
                  coefficient_delta[predictors[row]] * column_delta;
              }
            }
          }
        }
        update_online_moment(
          intercept_draw, online_count, intercept_mean, intercept_m2
        );
        update_online_moment(
          residual_var, online_count, residual_var_mean, residual_var_m2
        );
        if (compute_pve) {
          for (int block = 0; block < number_of_blocks; ++block) {
            double mean = block_pve_mean[block];
            double m2 = block_pve_m2[block];
            update_online_moment(block_pve[block], online_count, mean, m2);
            block_pve_mean[block] = mean;
            block_pve_m2[block] = m2;
          }
          update_online_moment(
            total_pve, online_count, total_pve_mean, total_pve_m2
          );
          update_online_moment(
            cross_block_pve, online_count,
            cross_block_pve_mean, cross_block_pve_m2
          );
        }
        for (int block = 0; block < number_of_blocks; ++block) {
          const PriorModel model = prior_models[block];
          if (model == PriorModel::Normal) {
            double mean = normal_var_mean[block];
            double m2 = normal_var_m2[block];
            update_online_moment(normal_var[block], online_count, mean, m2);
            normal_var_mean[block] = mean;
            normal_var_m2[block] = m2;
          } else if (model == PriorModel::SpikeSlab) {
            double pi_block_mean = pi_mean[block];
            double pi_block_m2 = pi_m2[block];
            update_online_moment(
              pi[block], online_count, pi_block_mean, pi_block_m2
            );
            pi_mean[block] = pi_block_mean;
            pi_m2[block] = pi_block_m2;
            double variance_mean = slab_var_mean[block];
            double variance_m2 = slab_var_m2[block];
            update_online_moment(
              slab_var[block], online_count, variance_mean, variance_m2
            );
            slab_var_mean[block] = variance_mean;
            slab_var_m2[block] = variance_m2;
          } else if (model == PriorModel::GlobalLocal) {
            double mean = tau_sq_mean[block];
            double m2 = tau_sq_m2[block];
            update_online_moment(tau_sq[block], online_count, mean, m2);
            tau_sq_mean[block] = mean;
            tau_sq_m2[block] = m2;
          } else if (model == PriorModel::SpikeMultiSlab) {
            Rcpp::NumericVector block_pi_mean = multi_pi_mean[block];
            Rcpp::NumericVector block_pi_m2 = multi_pi_m2[block];
            for (int component = 0;
                 component < static_cast<int>(multi_pi[block].size());
                 ++component) {
              double mean = block_pi_mean[component];
              double m2 = block_pi_m2[component];
              update_online_moment(
                multi_pi[block][component], online_count, mean, m2
              );
              block_pi_mean[component] = mean;
              block_pi_m2[component] = m2;
            }
            double mean = multi_var_mean[block];
            double m2 = multi_var_m2[block];
            update_online_moment(multi_var[block], online_count, mean, m2);
            multi_var_mean[block] = mean;
            multi_var_m2[block] = m2;
          }
        }
      }
      ++retained_index;
    }

    if (next_progress_percent <= 100) {
      int threshold = static_cast<int>(std::ceil(
        static_cast<double>(iterations) * next_progress_percent / 100.0
      ));
      if (iteration >= threshold) {
        do {
          next_progress_percent += 10;
          if (next_progress_percent > 100) {
            break;
          }
          threshold = static_cast<int>(std::ceil(
            static_cast<double>(iterations) * next_progress_percent / 100.0
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
      Rcpp::Named("multi_var_samples") = multi_var_samples,
      Rcpp::Named("block_pve_samples") = block_pve_samples,
      Rcpp::Named("total_pve_samples") = total_pve_samples,
      Rcpp::Named("cross_block_pve_samples") = cross_block_pve_samples
    );
  }
  Rcpp::List summaries = Rcpp::List::create(
    Rcpp::Named("number_of_draws") = number_of_draws,
    Rcpp::Named("coefficient_mean") = coefficient_mean,
    Rcpp::Named("coefficient_m2") = coefficient_m2,
    Rcpp::Named("intercept_mean") = intercept_mean,
    Rcpp::Named("intercept_m2") = intercept_m2,
    Rcpp::Named("residual_var_mean") = residual_var_mean,
    Rcpp::Named("residual_var_m2") = residual_var_m2,
    Rcpp::Named("normal_var_mean") = normal_var_mean,
    Rcpp::Named("normal_var_m2") = normal_var_m2,
    Rcpp::Named("inclusion_sum") = inclusion_sum,
    Rcpp::Named("pi_mean") = pi_mean,
    Rcpp::Named("pi_m2") = pi_m2,
    Rcpp::Named("slab_var_mean") = slab_var_mean,
    Rcpp::Named("slab_var_m2") = slab_var_m2,
    Rcpp::Named("local_var_mean") = local_var_mean,
    Rcpp::Named("local_var_m2") = local_var_m2,
    Rcpp::Named("tau_sq_mean") = tau_sq_mean,
    Rcpp::Named("tau_sq_m2") = tau_sq_m2,
    Rcpp::Named("multi_component_sum") = multi_component_sum,
    Rcpp::Named("multi_pi_mean") = multi_pi_mean,
    Rcpp::Named("multi_pi_m2") = multi_pi_m2,
    Rcpp::Named("multi_var_mean") = multi_var_mean,
    Rcpp::Named("multi_var_m2") = multi_var_m2,
    Rcpp::Named("block_pve_mean") = block_pve_mean,
    Rcpp::Named("block_pve_m2") = block_pve_m2,
    Rcpp::Named("total_pve_mean") = total_pve_mean,
    Rcpp::Named("total_pve_m2") = total_pve_m2,
    Rcpp::Named("cross_block_pve_mean") = cross_block_pve_mean,
    Rcpp::Named("cross_block_pve_m2") = cross_block_pve_m2
  );
  if (store_coefficient_cov) {
    summaries["coefficient_cov_m2"] = coefficient_cov_m2;
  }
  return summaries;
}

}  // namespace bayeslinreg

#endif
