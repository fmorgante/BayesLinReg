#ifndef BAYESLINREG_PARALLEL_BLOCKS_H
#define BAYESLINREG_PARALLEL_BLOCKS_H

#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>
#include "summary_matrices.h"

namespace bayeslinreg {

class BlockMultiplyWorker : public RcppParallel::Worker {
 public:
  BlockMultiplyWorker(
      const BlockSummaryMatrix& matrix,
      const std::vector<double>& coefficient,
      std::vector<double>& fitted)
    : matrix_(matrix), coefficient_(coefficient), fitted_(fitted) {}

  void operator()(const std::size_t begin, const std::size_t end) {
    for (std::size_t block = begin; block < end; ++block) {
      matrix_.multiply_block(
        static_cast<int>(block), coefficient_, fitted_
      );
    }
  }

 private:
  const BlockSummaryMatrix& matrix_;
  const std::vector<double>& coefficient_;
  std::vector<double>& fitted_;
};

class EigenBlockMultiplyWorker : public RcppParallel::Worker {
 public:
  EigenBlockMultiplyWorker(
      const EigenBlockSummaryMatrix& matrix,
      const std::vector<double>& coefficient,
      std::vector<double>& fitted)
    : matrix_(matrix), coefficient_(coefficient), fitted_(fitted) {}

  void operator()(const std::size_t begin, const std::size_t end) {
    for (std::size_t block = begin; block < end; ++block) {
      matrix_.multiply_block(static_cast<int>(block), coefficient_, fitted_);
    }
  }

 private:
  const EigenBlockSummaryMatrix& matrix_;
  const std::vector<double>& coefficient_;
  std::vector<double>& fitted_;
};

inline void BlockSummaryMatrix::multiply(
    const std::vector<double>& coefficient,
    std::vector<double>& fitted) const {
  std::fill(fitted.begin(), fitted.end(), 0.0);
  if (nthreads_ <= 1 || blocks_.size() <= 1) {
    for (int block = 0; block < block_count(); ++block) {
      multiply_block(block, coefficient, fitted);
    }
    return;
  }
  BlockMultiplyWorker worker(*this, coefficient, fitted);
  RcppParallel::parallelFor(
    0, blocks_.size(), worker, 1, nthreads_
  );
}

inline void EigenBlockSummaryMatrix::multiply(
    const std::vector<double>& coefficient,
    std::vector<double>& fitted) const {
  std::fill(fitted.begin(), fitted.end(), 0.0);
  if (nthreads_ <= 1 || blocks_.size() <= 1) {
    for (int block = 0; block < block_count(); ++block) {
      multiply_block(block, coefficient, fitted);
    }
    return;
  }
  EigenBlockMultiplyWorker worker(*this, coefficient, fitted);
  RcppParallel::parallelFor(0, blocks_.size(), worker, 1, nthreads_);
}

inline std::uint64_t splitmix64(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

class BlockRng {
 public:
  BlockRng() : state_(0), increment_(1) {}

  BlockRng(const std::uint64_t state, const std::uint64_t sequence) {
    seed(state, sequence);
  }

  void seed(const std::uint64_t state, const std::uint64_t sequence) {
    state_ = 0;
    increment_ = (sequence << 1u) | 1u;
    next_uint32();
    state_ += state;
    next_uint32();
  }

  double uniform() {
    const std::uint64_t high = static_cast<std::uint64_t>(next_uint32());
    const std::uint64_t low = static_cast<std::uint64_t>(next_uint32());
    const std::uint64_t bits = (high << 32) | low;
    return (static_cast<double>(bits >> 11) + 0.5) *
      (1.0 / 9007199254740992.0);
  }

  double normal() {
    const double radius = std::sqrt(-2.0 * std::log(uniform()));
    const double angle = 2.0 * std::acos(-1.0) * uniform();
    return radius * std::cos(angle);
  }

 private:
  std::uint32_t next_uint32() {
    const std::uint64_t old_state = state_;
    state_ = old_state * UINT64_C(6364136223846793005) + increment_;
    const std::uint32_t shifted = static_cast<std::uint32_t>(
      ((old_state >> 18u) ^ old_state) >> 27u
    );
    const std::uint32_t rotation = static_cast<std::uint32_t>(old_state >> 59u);
    return (shifted >> rotation) |
      (shifted << ((-rotation) & 31));
  }

  std::uint64_t state_;
  std::uint64_t increment_;
};

template <typename BlockMatrix>
class BlockCoefficientWorker : public RcppParallel::Worker {
 public:
  BlockCoefficientWorker(
      const BlockMatrix& matrix,
      const int* block_id,
      const std::vector<PriorModel>& block_model,
      const std::vector<int>& model_local_index,
      const std::vector<double>& x_squared,
      const double residual_var,
      const std::vector<double>& normal_var,
      const std::vector<double>& pi,
      const std::vector<double>& slab_var,
      const std::vector<double>& tau_sq,
      const std::vector<double>& local_var,
      const std::vector< std::vector<double> >& multi_gamma,
      const std::vector< std::vector<double> >& multi_pi,
      const std::vector<double>& multi_var,
      const bool learn_residual_var,
      std::vector<double>& coefficient,
      std::vector<double>& corrected_rhs,
      std::vector<int>& inclusion,
      std::vector<int>& multi_component,
      std::vector<BlockRng>& rng,
      std::vector<double>& residual_sse_change)
    : matrix_(matrix), block_id_(block_id), block_model_(block_model),
      model_local_index_(model_local_index), x_squared_(x_squared),
      residual_var_(residual_var), normal_var_(normal_var), pi_(pi),
      slab_var_(slab_var), tau_sq_(tau_sq), local_var_(local_var),
      multi_gamma_(multi_gamma), multi_pi_(multi_pi), multi_var_(multi_var),
      learn_residual_var_(learn_residual_var), coefficient_(coefficient),
      corrected_rhs_(corrected_rhs), inclusion_(inclusion),
      multi_component_(multi_component), rng_(rng),
      residual_sse_change_(residual_sse_change) {}

  void operator()(const std::size_t begin, const std::size_t end) {
    for (std::size_t gram_block = begin; gram_block < end; ++gram_block) {
      double sse_change = 0.0;
      BlockRng& generator = rng_[gram_block];
      const std::vector<int>& predictors =
        matrix_.block_predictors(static_cast<int>(gram_block));
      for (const int j : predictors) {
        const int prior_block = block_id_[j] - 1;
        const PriorModel model = block_model_[prior_block];
        const int local_index = model_local_index_[j];
        const double old_coefficient = coefficient_[j];
        const double partial_rhs = matrix_.corrected_value(
          corrected_rhs_, j, 0.0
        ) +
          x_squared_[j] * old_coefficient;

        if (model == PriorModel::SpikeMultiSlab) {
          const int component_count = multi_gamma_[prior_block].size();
          std::vector<double> weights(component_count, 0.0);
          std::vector<double> conditional_vars(component_count, 0.0);
          std::vector<double> conditional_means(component_count, 0.0);
          double maximum_log_weight = -std::numeric_limits<double>::infinity();
          for (int component = 0; component < component_count; ++component) {
            weights[component] = std::log(std::max(
              multi_pi_[prior_block][component],
              std::numeric_limits<double>::min()
            ));
            if (component > 0) {
              const double prior_var =
                multi_gamma_[prior_block][component] * multi_var_[prior_block];
              conditional_vars[component] = 1.0 / (
                x_squared_[j] / residual_var_ + 1.0 / prior_var
              );
              conditional_means[component] = conditional_vars[component] *
                partial_rhs / residual_var_;
              weights[component] +=
                0.5 * std::log(conditional_vars[component] / prior_var) +
                conditional_means[component] * conditional_means[component] /
                  (2.0 * conditional_vars[component]);
            }
            maximum_log_weight = std::max(
              maximum_log_weight, weights[component]
            );
          }
          double weight_total = 0.0;
          for (int component = 0; component < component_count; ++component) {
            weights[component] = std::exp(
              weights[component] - maximum_log_weight
            );
            weight_total += weights[component];
          }
          const double threshold = generator.uniform() * weight_total;
          double cumulative_weight = 0.0;
          int selected_component = component_count - 1;
          for (int component = 0; component < component_count; ++component) {
            cumulative_weight += weights[component];
            if (threshold <= cumulative_weight) {
              selected_component = component;
              break;
            }
          }
          multi_component_[local_index] = selected_component;
          coefficient_[j] = selected_component == 0
            ? 0.0
            : conditional_means[selected_component] +
                std::sqrt(conditional_vars[selected_component]) *
                  generator.normal();
        } else {
          const double prior_precision = model == PriorModel::Fixed
            ? 0.0
            : (model == PriorModel::GlobalLocal
                ? 1.0 / tau_sq_[prior_block] / local_var_[local_index]
                : (model == PriorModel::SpikeSlab
                    ? 1.0 / slab_var_[prior_block]
                    : 1.0 / normal_var_[prior_block]));
          const double conditional_var = 1.0 / (
            x_squared_[j] / residual_var_ + prior_precision
          );
          const double conditional_mean =
            conditional_var * partial_rhs / residual_var_;
          if (model == PriorModel::SpikeSlab) {
            const double epsilon = std::numeric_limits<double>::epsilon();
            const double bounded_pi = std::min(
              std::max(pi_[prior_block], epsilon), 1.0 - epsilon
            );
            const double log_inclusion_odds =
              std::log(bounded_pi) - std::log1p(-bounded_pi) +
              0.5 * std::log(conditional_var / slab_var_[prior_block]) +
              conditional_mean * conditional_mean / (2.0 * conditional_var);
            const double inclusion_probability = log_inclusion_odds >= 0.0
              ? 1.0 / (1.0 + std::exp(-log_inclusion_odds))
              : std::exp(log_inclusion_odds) /
                  (1.0 + std::exp(log_inclusion_odds));
            inclusion_[local_index] =
              generator.uniform() < inclusion_probability ? 1 : 0;
          }
          coefficient_[j] =
            model != PriorModel::SpikeSlab || inclusion_[local_index] == 1
            ? conditional_mean + std::sqrt(conditional_var) * generator.normal()
            : 0.0;
        }

        const double coefficient_change = coefficient_[j] - old_coefficient;
        if (coefficient_change != 0.0) {
          matrix_.update(corrected_rhs_, j, coefficient_change);
          if (learn_residual_var_) {
            sse_change += -2.0 * coefficient_change * partial_rhs +
              (coefficient_[j] * coefficient_[j] -
               old_coefficient * old_coefficient) * x_squared_[j];
          }
        }
      }
      residual_sse_change_[gram_block] = sse_change;
    }
  }

 private:
  const BlockMatrix& matrix_;
  const int* block_id_;
  const std::vector<PriorModel>& block_model_;
  const std::vector<int>& model_local_index_;
  const std::vector<double>& x_squared_;
  double residual_var_;
  const std::vector<double>& normal_var_;
  const std::vector<double>& pi_;
  const std::vector<double>& slab_var_;
  const std::vector<double>& tau_sq_;
  const std::vector<double>& local_var_;
  const std::vector< std::vector<double> >& multi_gamma_;
  const std::vector< std::vector<double> >& multi_pi_;
  const std::vector<double>& multi_var_;
  bool learn_residual_var_;
  std::vector<double>& coefficient_;
  std::vector<double>& corrected_rhs_;
  std::vector<int>& inclusion_;
  std::vector<int>& multi_component_;
  std::vector<BlockRng>& rng_;
  std::vector<double>& residual_sse_change_;
};

template <typename BlockMatrix>
void parallel_coefficient_sweep(
    const BlockMatrix& matrix,
    const Rcpp::IntegerVector& block_id,
    const std::vector<PriorModel>& block_model,
    const std::vector<int>& model_local_index,
    const std::vector<double>& x_squared,
    const double residual_var,
    const std::vector<double>& normal_var,
    const std::vector<double>& pi,
    const std::vector<double>& slab_var,
    const std::vector<double>& tau_sq,
    const std::vector<double>& local_var,
    const std::vector< std::vector<double> >& multi_gamma,
    const std::vector< std::vector<double> >& multi_pi,
    const std::vector<double>& multi_var,
    const bool learn_residual_var,
    std::vector<double>& coefficient,
    std::vector<double>& corrected_rhs,
    std::vector<int>& inclusion,
    std::vector<int>& multi_component,
    std::vector<BlockRng>& rng,
    double& residual_sse,
    const int nthreads) {
  std::vector<double> residual_sse_change(matrix.block_count(), 0.0);
  BlockCoefficientWorker<BlockMatrix> worker(
    matrix, block_id.begin(), block_model, model_local_index,
    x_squared, residual_var, normal_var, pi, slab_var, tau_sq, local_var,
    multi_gamma, multi_pi, multi_var, learn_residual_var, coefficient,
    corrected_rhs, inclusion, multi_component, rng, residual_sse_change
  );
  RcppParallel::parallelFor(
    0, matrix.block_count(), worker, 1, nthreads
  );
  if (learn_residual_var) {
    for (const double change : residual_sse_change) residual_sse += change;
  }
}

template <typename SummaryMatrix>
struct is_parallel_block_matrix : std::false_type {};

template <>
struct is_parallel_block_matrix<BlockSummaryMatrix> : std::true_type {};

template <>
struct is_parallel_block_matrix<EigenBlockSummaryMatrix> : std::true_type {};

template <typename SummaryMatrix>
struct has_streaming_triangular_blocks : std::false_type {};

template <>
struct has_streaming_triangular_blocks<BlockSummaryMatrix> : std::true_type {};

}  // namespace bayeslinreg

#endif
