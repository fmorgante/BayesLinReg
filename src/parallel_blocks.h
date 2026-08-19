#ifndef BAYESLINREG_PARALLEL_BLOCKS_H
#define BAYESLINREG_PARALLEL_BLOCKS_H

#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>
#include "coefficient_updates.h"
#include "gig_sampler.h"
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

class LDBlockMultiplyWorker : public RcppParallel::Worker {
 public:
  LDBlockMultiplyWorker(
      const LDSummaryMatrix& matrix,
      const std::vector<double>& coefficient,
      std::vector<double>& fitted)
    : matrix_(matrix), coefficient_(coefficient), fitted_(fitted) {}

  void operator()(const std::size_t begin, const std::size_t end) {
    for (std::size_t block = begin; block < end; ++block) {
      matrix_.multiply_block(static_cast<int>(block), coefficient_, fitted_);
    }
  }

 private:
  const LDSummaryMatrix& matrix_;
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

inline void LDSummaryMatrix::multiply(
    const std::vector<double>& coefficient,
    std::vector<double>& fitted) const {
  std::fill(fitted.begin(), fitted.end(), 0.0);
  if (nthreads_ <= 1 || blocks_.size() <= 1) {
    for (int block = 0; block < block_count(); ++block) {
      multiply_block(block, coefficient, fitted);
    }
    return;
  }
  LDBlockMultiplyWorker worker(*this, coefficient, fitted);
  RcppParallel::parallelFor(0, blocks_.size(), worker, 1, nthreads_);
}

// Reusable retained-draw PVE storage. Each LD block owns one output row, so
// workers never contend for writes. Reduction is deliberately performed later
// in LD-block order to keep results independent of thread scheduling.
struct LDPveWorkspace {
  int ld_blocks = 0;
  int prior_blocks = 0;
  std::vector<double> total;
  std::vector<double> contribution_scale;
  std::vector<double> standalone;
  std::vector<double> allocated;

  void ensure_size(const int new_ld_blocks, const int new_prior_blocks) {
    if (ld_blocks == new_ld_blocks && prior_blocks == new_prior_blocks) return;
    ld_blocks = new_ld_blocks;
    prior_blocks = new_prior_blocks;
    total.resize(ld_blocks);
    contribution_scale.resize(ld_blocks);
    const std::size_t block_entries =
      static_cast<std::size_t>(ld_blocks) * prior_blocks;
    standalone.resize(block_entries);
    allocated.resize(block_entries);
  }
};

inline void LDSummaryMatrix::pve_block_quadratics(
    const int ld_block,
    const std::vector<double>& coefficient,
    const int* prior_block,
    const int number_of_prior_blocks,
    double& total,
    double& contribution_scale,
    double* standalone,
    double* allocated) const {
  total = 0.0;
  contribution_scale = 0.0;
  std::fill(standalone, standalone + number_of_prior_blocks, 0.0);
  std::fill(allocated, allocated + number_of_prior_blocks, 0.0);

  const Block& block = blocks_[ld_block];
  const bool single_prior_block = number_of_prior_blocks == 1;
  for (int column = 0; column < block.size; ++column) {
    const int global_column = block.global[column];
    const double column_coefficient = coefficient[global_column];
    const double diagonal_contribution =
      column_coefficient * diagonal(global_column) * column_coefficient;
    total += diagonal_contribution;
    contribution_scale += std::abs(diagonal_contribution);

    if (!single_prior_block) {
      const int column_prior_block = prior_block[global_column] - 1;
      standalone[column_prior_block] += diagonal_contribution;
      allocated[column_prior_block] += diagonal_contribution;
    }

    for (int position = block.indptr[column];
         position < block.indptr[column + 1]; ++position) {
      const int row = block.type == 0
        ? column + 1 + position - block.indptr[column]
        : block.row_index[position];
      const int global_row = block.global[row];
      const double edge_contribution = column_coefficient *
        scale_[global_column] * off_diagonal_scale_ * block.data[position] *
        scale_[global_row] * coefficient[global_row];
      total += 2.0 * edge_contribution;
      contribution_scale += 2.0 * std::abs(edge_contribution);

      if (!single_prior_block) {
        const int column_prior_block = prior_block[global_column] - 1;
        const int row_prior_block = prior_block[global_row] - 1;
        if (column_prior_block == row_prior_block) {
          standalone[column_prior_block] += 2.0 * edge_contribution;
          allocated[column_prior_block] += 2.0 * edge_contribution;
        } else {
          allocated[column_prior_block] += edge_contribution;
          allocated[row_prior_block] += edge_contribution;
        }
      }
    }
  }

  if (single_prior_block) {
    standalone[0] = total;
    allocated[0] = total;
  }
}

class LDPveWorker : public RcppParallel::Worker {
 public:
  LDPveWorker(
      const LDSummaryMatrix& matrix,
      const std::vector<double>& coefficient,
      const int* prior_block,
      const int number_of_prior_blocks,
      LDPveWorkspace& workspace)
    : matrix_(matrix), coefficient_(coefficient), prior_block_(prior_block),
      number_of_prior_blocks_(number_of_prior_blocks), workspace_(workspace) {}

  void operator()(const std::size_t begin, const std::size_t end) {
    for (std::size_t ld_block = begin; ld_block < end; ++ld_block) {
      const std::size_t offset = ld_block * number_of_prior_blocks_;
      matrix_.pve_block_quadratics(
        static_cast<int>(ld_block), coefficient_, prior_block_,
        number_of_prior_blocks_, workspace_.total[ld_block],
        workspace_.contribution_scale[ld_block],
        workspace_.standalone.data() + offset,
        workspace_.allocated.data() + offset
      );
    }
  }

 private:
  const LDSummaryMatrix& matrix_;
  const std::vector<double>& coefficient_;
  const int* prior_block_;
  int number_of_prior_blocks_;
  LDPveWorkspace& workspace_;
};

inline void parallel_ld_pve_quadratics(
    const LDSummaryMatrix& matrix,
    const std::vector<double>& coefficient,
    const Rcpp::IntegerVector& prior_block,
    const int number_of_prior_blocks,
    LDPveWorkspace& workspace,
    std::vector<double>& standalone,
    std::vector<double>& allocated,
    double& total,
    double& contribution_scale,
    const int nthreads) {
  const int ld_blocks = matrix.block_count();
  workspace.ensure_size(ld_blocks, number_of_prior_blocks);
  LDPveWorker worker(
    matrix, coefficient, prior_block.begin(), number_of_prior_blocks, workspace
  );
  if (nthreads > 1 && ld_blocks > 1) {
    RcppParallel::parallelFor(0, ld_blocks, worker, 1, nthreads);
  } else {
    worker(0, ld_blocks);
  }

  total = 0.0;
  contribution_scale = 0.0;
  std::fill(standalone.begin(), standalone.end(), 0.0);
  std::fill(allocated.begin(), allocated.end(), 0.0);
  for (int ld_block = 0; ld_block < ld_blocks; ++ld_block) {
    total += workspace.total[ld_block];
    contribution_scale += workspace.contribution_scale[ld_block];
    const std::size_t offset =
      static_cast<std::size_t>(ld_block) * number_of_prior_blocks;
    for (int block = 0; block < number_of_prior_blocks; ++block) {
      standalone[block] += workspace.standalone[offset + block];
      allocated[block] += workspace.allocated[offset + block];
    }
  }
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

  int bernoulli(const double probability) {
    return uniform() < probability ? 1 : 0;
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
      std::vector<MixtureWorkspace>& mixture_workspace,
      std::vector<double>& residual_sse_change)
    : matrix_(matrix), block_id_(block_id), block_model_(block_model),
      model_local_index_(model_local_index), x_squared_(x_squared),
      residual_var_(residual_var), normal_var_(normal_var), pi_(pi),
      slab_var_(slab_var), tau_sq_(tau_sq), local_var_(local_var),
      multi_gamma_(multi_gamma), multi_pi_(multi_pi), multi_var_(multi_var),
      learn_residual_var_(learn_residual_var), coefficient_(coefficient),
      corrected_rhs_(corrected_rhs), inclusion_(inclusion),
      multi_component_(multi_component), rng_(rng),
      mixture_workspace_(mixture_workspace),
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

        MixtureWorkspace& workspace = mixture_workspace_[gram_block];
        coefficient_[j] = draw_coefficient_coordinate(
          model, prior_block, local_index, x_squared_[j], residual_var_,
          partial_rhs, normal_var_, pi_, slab_var_, tau_sq_, local_var_,
          multi_gamma_, multi_pi_, multi_var_, inclusion_, multi_component_,
          generator, workspace
        );

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
  std::vector<MixtureWorkspace>& mixture_workspace_;
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
    std::vector<MixtureWorkspace>& mixture_workspace,
    double& residual_sse,
    const int nthreads) {
  std::vector<double> residual_sse_change(matrix.block_count(), 0.0);
  BlockCoefficientWorker<BlockMatrix> worker(
    matrix, block_id.begin(), block_model, model_local_index,
    x_squared, residual_var, normal_var, pi, slab_var, tau_sq, local_var,
    multi_gamma, multi_pi, multi_var, learn_residual_var, coefficient,
    corrected_rhs, inclusion, multi_component, rng, mixture_workspace,
    residual_sse_change
  );
  RcppParallel::parallelFor(
    0, matrix.block_count(), worker, 1, nthreads
  );
  if (learn_residual_var) {
    for (const double change : residual_sse_change) residual_sse += change;
  }
}

template <typename BlockMatrix>
class BlockLocalVarianceWorker : public RcppParallel::Worker {
 public:
  BlockLocalVarianceWorker(
      const BlockMatrix& matrix,
      const int* block_id,
      const std::vector<PriorModel>& block_model,
      const std::vector<int>& model_local_index,
      const std::vector<double>& coefficient,
      const std::vector<double>& tau_sq,
      const double* local_a,
      const std::vector<double>& local_aux,
      std::vector<double>& local_var,
      std::vector<BlockRng>& rng)
    : matrix_(matrix), block_id_(block_id), block_model_(block_model),
      model_local_index_(model_local_index), coefficient_(coefficient),
      tau_sq_(tau_sq), local_a_(local_a), local_aux_(local_aux),
      local_var_(local_var), rng_(rng) {}

  void operator()(const std::size_t begin, const std::size_t end) {
    for (std::size_t gram_block = begin; gram_block < end; ++gram_block) {
      BlockRng& generator = rng_[gram_block];
      const std::vector<int>& predictors =
        matrix_.block_predictors(static_cast<int>(gram_block));
      for (const int j : predictors) {
        const int prior_block = block_id_[j] - 1;
        if (block_model_[prior_block] != PriorModel::GlobalLocal) continue;
        const int local_index = model_local_index_[j];
        const double raw_chi =
          coefficient_[j] * coefficient_[j] / tau_sq_[prior_block];
        local_var_[local_index] = sample_gig(
          local_a_[prior_block] - 0.5,
          std::max(raw_chi, std::numeric_limits<double>::min()),
          2.0 * std::max(
            local_aux_[local_index], std::numeric_limits<double>::min()
          ),
          generator
        );
      }
    }
  }

 private:
  const BlockMatrix& matrix_;
  const int* block_id_;
  const std::vector<PriorModel>& block_model_;
  const std::vector<int>& model_local_index_;
  const std::vector<double>& coefficient_;
  const std::vector<double>& tau_sq_;
  const double* local_a_;
  const std::vector<double>& local_aux_;
  std::vector<double>& local_var_;
  std::vector<BlockRng>& rng_;
};

template <typename BlockMatrix>
void parallel_local_variance_sweep(
    const BlockMatrix& matrix,
    const Rcpp::IntegerVector& block_id,
    const std::vector<PriorModel>& block_model,
    const std::vector<int>& model_local_index,
    const std::vector<double>& coefficient,
    const std::vector<double>& tau_sq,
    const Rcpp::NumericVector& local_a,
    const std::vector<double>& local_aux,
    std::vector<double>& local_var,
    std::vector<BlockRng>& rng,
    const int nthreads) {
  BlockLocalVarianceWorker<BlockMatrix> worker(
    matrix, block_id.begin(), block_model, model_local_index, coefficient,
    tau_sq, local_a.begin(), local_aux, local_var, rng
  );
  RcppParallel::parallelFor(
    0, matrix.block_count(), worker, 1, nthreads
  );
}

template <typename SummaryMatrix>
struct is_parallel_block_matrix : std::false_type {};

template <>
struct is_parallel_block_matrix<BlockSummaryMatrix> : std::true_type {};

template <>
struct is_parallel_block_matrix<EigenBlockSummaryMatrix> : std::true_type {};

template <>
struct is_parallel_block_matrix<LDSummaryMatrix> : std::true_type {};

template <typename SummaryMatrix>
struct has_streaming_triangular_blocks : std::false_type {};

template <>
struct has_streaming_triangular_blocks<BlockSummaryMatrix> : std::true_type {};

template <>
struct has_streaming_triangular_blocks<LDSummaryMatrix> : std::true_type {};

template <typename SummaryMatrix>
struct has_streaming_serial_sweep_order : std::false_type {};

template <>
struct has_streaming_serial_sweep_order<BlockSummaryMatrix> : std::true_type {};

template <typename SummaryMatrix>
struct has_parallel_ld_pve : std::false_type {};

template <>
struct has_parallel_ld_pve<LDSummaryMatrix> : std::true_type {};

}  // namespace bayeslinreg

#endif
