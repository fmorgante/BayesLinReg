#ifndef BAYESLINREG_SUMMARY_MATRICES_H
#define BAYESLINREG_SUMMARY_MATRICES_H

#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>
#include "sampler_types.h"

namespace bayeslinreg {

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
    const Eigen::Map<const Eigen::MatrixXd> matrix(
      matrix_.begin(), p, p
    );
    const Eigen::Map<const Eigen::VectorXd> coefficient_vector(
      coefficient.data(), p
    );
    Eigen::Map<Eigen::VectorXd> fitted_vector(fitted.data(), p);
    fitted_vector.noalias() =
      matrix.selfadjointView<Eigen::Upper>() * coefficient_vector;
  }

  double center_dot(const std::vector<double>&) const { return 0.0; }

  double centered_fitted(
      const double fitted,
      const int,
      const double) const {
    return fitted;
  }

  void update_center_dot(double&, const int, const double) const {}

  double block_quadratic(
      const std::vector<double>& coefficient,
      const std::vector<int>& predictors,
      const Rcpp::IntegerVector&,
      const int) const {
    double result = 0.0;
    for (std::size_t column = 0; column < predictors.size(); ++column) {
      const int j = predictors[column];
      for (std::size_t row = 0; row < predictors.size(); ++row) {
        const int k = predictors[row];
        result += coefficient[k] * matrix_(k, j) * coefficient[j];
      }
    }
    return result;
  }

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

  double block_quadratic(
      const std::vector<double>& coefficient,
      const std::vector<int>& predictors,
      const Rcpp::IntegerVector& block_id,
      const int block) const {
    double sparse_result = 0.0;
    double block_center_dot = 0.0;
    for (std::size_t index = 0; index < predictors.size(); ++index) {
      const int j = predictors[index];
      block_center_dot += center_[j] * coefficient[j];
      for (Eigen::MappedSparseMatrix<double>::InnerIterator entry(matrix_, j);
           entry; ++entry) {
        if (block_id[entry.row()] - 1 == block) {
          sparse_result += coefficient[entry.row()] * entry.value() *
            coefficient[j];
        }
      }
    }
    return sparse_result - block_center_dot * block_center_dot;
  }

 private:
  const Eigen::MappedSparseMatrix<double>& matrix_;
  const Rcpp::NumericVector& center_;
};

// A block-diagonal sufficient-statistics matrix. Each block may be dense,
// general sparse, or lower-triangular symmetric sparse. The triangular
// representation streams updates only to coordinates that have not yet been
// visited in the current ascending Gibbs sweep, avoiding a reverse index.
class BlockSummaryMatrix {
 public:
  BlockSummaryMatrix(
      const Rcpp::List& matrices,
      const Rcpp::List& indices,
      const Rcpp::IntegerVector& types,
      const Rcpp::NumericVector& center,
      const int nthreads = 1)
    : center_(center), p_(center.size()),
      global_block_(p_, -1), global_local_(p_, -1), diagonal_(p_, 0.0),
      nthreads_(nthreads), has_streaming_blocks_(false) {
    if (matrices.size() < 1 || matrices.size() != indices.size() ||
        matrices.size() != types.size()) {
      Rcpp::stop("Invalid block-diagonal `XtX` representation.");
    }
    blocks_.reserve(matrices.size());
    for (int block_index = 0; block_index < matrices.size(); ++block_index) {
      Block block;
      block.type = bayeslinreg::gram_storage_from_code(types[block_index]);
      const Rcpp::IntegerVector mapping = indices[block_index];
      block.size = mapping.size();
      block.global.resize(block.size);
      for (int local = 0; local < block.size; ++local) {
        const int global = mapping[local] - 1;
        if (global < 0 || global >= p_ || global_block_[global] >= 0) {
          Rcpp::stop("Gram-block predictor indices must partition predictors.");
        }
        block.global[local] = global;
        global_block_[global] = block_index;
        global_local_[global] = local;
      }

      if (block.type == GramStorage::Dense) {
        const Rcpp::NumericMatrix matrix = matrices[block_index];
        if (matrix.nrow() != block.size || matrix.ncol() != block.size) {
          Rcpp::stop("Dense Gram-block dimensions do not match their indices.");
        }
        block.dense = matrix.begin();
        for (int local = 0; local < block.size; ++local) {
          diagonal_[block.global[local]] =
            block.dense[local + static_cast<std::size_t>(block.size) * local];
        }
      } else if (block.type == GramStorage::SparseGeneral ||
                 block.type == GramStorage::SparseLowerTriangular) {
        const Rcpp::S4 matrix = matrices[block_index];
        const Rcpp::IntegerVector dimensions = matrix.slot("Dim");
        if (dimensions[0] != block.size || dimensions[1] != block.size) {
          Rcpp::stop("Sparse Gram-block dimensions do not match their indices.");
        }
        const Rcpp::IntegerVector column_pointer = matrix.slot("p");
        const Rcpp::IntegerVector row_index = matrix.slot("i");
        const Rcpp::NumericVector values = matrix.slot("x");
        block.column_pointer = column_pointer.begin();
        block.row_index = row_index.begin();
        block.values = values.begin();
        block.nonzeros = values.size();
        for (int column = 0; column < block.size; ++column) {
          for (int position = block.column_pointer[column];
               position < block.column_pointer[column + 1]; ++position) {
            if (block.row_index[position] == column) {
              diagonal_[block.global[column]] = block.values[position];
            }
          }
        }
        if (block.type == GramStorage::SparseLowerTriangular) {
          has_streaming_blocks_ = true;
        }
      }
      blocks_.push_back(std::move(block));
    }
    for (int global = 0; global < p_; ++global) {
      if (global_block_[global] < 0) {
        Rcpp::stop("Gram-block predictor indices must cover every predictor.");
      }
    }
  }

  int cols() const { return p_; }

  bool has_streaming_blocks() const { return has_streaming_blocks_; }

  // Lower-triangular sparse blocks can update only block-local coordinates
  // that have not yet been visited.  When ETA order differs from the source
  // Gram order, retain the global scan slots but place each streaming block's
  // predictors in its local triangular order within those slots.
  std::vector<int> serial_sweep_order() const {
    std::vector<int> order(p_);
    for (int j = 0; j < p_; ++j) order[j] = j;
    for (const Block& block : blocks_) {
      if (block.type != GramStorage::SparseLowerTriangular) continue;
      std::vector<int> slots = block.global;
      std::sort(slots.begin(), slots.end());
      for (int local = 0; local < block.size; ++local) {
        order[slots[local]] = block.global[local];
      }
    }
    return order;
  }

  double diagonal(const int j) const {
    return diagonal_[j] - center_[j] * center_[j];
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
    const Block& block = blocks_[global_block_[j]];
    const int local = global_local_[j];
    if (block.type == GramStorage::Dense) {
      const double* column = block.dense +
        static_cast<std::size_t>(block.size) * local;
      for (int row = 0; row < block.size; ++row) {
        rhs[block.global[row]] -= column[row] * change;
      }
      return;
    }
    for (int position = block.column_pointer[local];
         position < block.column_pointer[local + 1]; ++position) {
      rhs[block.global[block.row_index[position]]] -=
        block.values[position] * change;
    }
  }

  void multiply(
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const;

  int block_count() const {
    return static_cast<int>(blocks_.size());
  }

  const std::vector<int>& block_predictors(const int block) const {
    return blocks_[block].global;
  }

  void multiply_block(
      const int block_index,
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const {
    const Block& block = blocks_[block_index];
    if (block.type == GramStorage::Dense) {
      for (int column = 0; column < block.size; ++column) {
        const double beta = coefficient[block.global[column]];
        const double* values = block.dense +
          static_cast<std::size_t>(block.size) * column;
        for (int row = 0; row < block.size; ++row) {
          fitted[block.global[row]] += values[row] * beta;
        }
      }
    } else if (block.type == GramStorage::SparseGeneral) {
      for (int column = 0; column < block.size; ++column) {
        const double beta = coefficient[block.global[column]];
        for (int position = block.column_pointer[column];
             position < block.column_pointer[column + 1]; ++position) {
          fitted[block.global[block.row_index[position]]] +=
            block.values[position] * beta;
        }
      }
    } else {
      for (int column = 0; column < block.size; ++column) {
        const int global_column = block.global[column];
        for (int position = block.column_pointer[column];
             position < block.column_pointer[column + 1]; ++position) {
          const int row = block.row_index[position];
          const int global_row = block.global[row];
          const double value = block.values[position];
          fitted[global_row] += value * coefficient[global_column];
          if (row != column) {
            fitted[global_column] += value * coefficient[global_row];
          }
        }
      }
    }
  }

  double center_dot(const std::vector<double>& coefficient) const {
    double result = 0.0;
    for (int j = 0; j < p_; ++j) result += center_[j] * coefficient[j];
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

  double block_quadratic(
      const std::vector<double>& coefficient,
      const std::vector<int>& predictors,
      const Rcpp::IntegerVector& block_id,
      const int prior_block) const {
    double result = 0.0;
    double block_center_dot = 0.0;
    for (const int j : predictors) {
      block_center_dot += center_[j] * coefficient[j];
      const Block& block = blocks_[global_block_[j]];
      const int local = global_local_[j];
      if (block.type == GramStorage::Dense) {
        const double* column = block.dense +
          static_cast<std::size_t>(block.size) * local;
        for (int row = 0; row < block.size; ++row) {
          const int global_row = block.global[row];
          if (block_id[global_row] - 1 == prior_block) {
            result += coefficient[global_row] * column[row] * coefficient[j];
          }
        }
      } else if (block.type == GramStorage::SparseGeneral) {
        for (int position = block.column_pointer[local];
             position < block.column_pointer[local + 1]; ++position) {
          const int global_row = block.global[block.row_index[position]];
          if (block_id[global_row] - 1 == prior_block) {
            result += coefficient[global_row] * block.values[position] *
              coefficient[j];
          }
        }
      } else {
        for (int position = block.column_pointer[local];
             position < block.column_pointer[local + 1]; ++position) {
          const int row = block.row_index[position];
          const int global_row = block.global[row];
          if (block_id[global_row] - 1 == prior_block) {
            const double contribution = coefficient[global_row] *
              block.values[position] * coefficient[j];
            result += row == local ? contribution : 2.0 * contribution;
          }
        }
      }
    }
    return result - block_center_dot * block_center_dot;
  }

 private:
  struct Block {
    GramStorage type = GramStorage::Dense;
    int size = 0;
    int nonzeros = 0;
    const double* dense = NULL;
    const int* column_pointer = NULL;
    const int* row_index = NULL;
    const double* values = NULL;
    std::vector<int> global;
  };

  const Rcpp::NumericVector& center_;
  int p_;
  std::vector<Block> blocks_;
  std::vector<int> global_block_;
  std::vector<int> global_local_;
  std::vector<double> diagonal_;
  int nthreads_;
  bool has_streaming_blocks_;
};

// A block-diagonal LD correlation matrix representing the working Gram matrix
// as diag(scale) * R * diag(scale). Each block stores only strict-lower
// correlations; the unit diagonal of R is implicit. Interval blocks omit row
// indices, while indexed blocks retain them for irregular sparse patterns.
class LDSummaryMatrix {
 public:
  LDSummaryMatrix(
      const Rcpp::List& input_blocks,
      const Rcpp::List& indices,
      const Rcpp::NumericVector& scale,
      const double ld_shrink,
      const int nthreads = 1)
    : scale_(scale), p_(scale.size()), global_block_(p_, -1),
      global_local_(p_, -1), off_diagonal_scale_(1.0 - ld_shrink),
      nthreads_(nthreads) {
    if (!std::isfinite(ld_shrink) || ld_shrink < 0.0 || ld_shrink >= 1.0) {
      Rcpp::stop("`ld_shrink` must be finite and in [0, 1).");
    }
    if (input_blocks.size() < 1 || input_blocks.size() != indices.size()) {
      Rcpp::stop("Invalid LD block representation.");
    }
    blocks_.reserve(input_blocks.size());
    for (int block_index = 0; block_index < input_blocks.size(); ++block_index) {
      const Rcpp::List input = input_blocks[block_index];
      const Rcpp::IntegerVector mapping = indices[block_index];
      Block block;
      block.type = Rcpp::as<int>(input["type"]);
      block.size = Rcpp::as<int>(input["size"]);
      const Rcpp::NumericVector data = input["data"];
      const Rcpp::IntegerVector indptr = input["indptr"];
      const Rcpp::IntegerVector row_index = input["row_index"];
      if ((block.type != 0 && block.type != 1) || block.size < 1 ||
          mapping.size() != block.size || indptr.size() != block.size + 1 ||
          indptr[0] != 0 || indptr[block.size] != data.size() ||
          (block.type == 1 && row_index.size() != data.size())) {
        Rcpp::stop("Invalid compressed LD block.");
      }
      block.data = data.begin();
      block.indptr = indptr.begin();
      block.row_index = row_index.begin();
      block.global.resize(block.size);
      for (int local = 0; local < block.size; ++local) {
        const int global = mapping[local] - 1;
        if (global < 0 || global >= p_ || global_block_[global] >= 0) {
          Rcpp::stop("LD block indices must partition predictors.");
        }
        block.global[local] = global;
        global_block_[global] = block_index;
        global_local_[global] = local;
        if (!std::isfinite(scale_[global]) || scale_[global] <= 0.0) {
          Rcpp::stop("LD predictor scales must be positive and finite.");
        }
      }
      for (int column = 0; column < block.size; ++column) {
        if (block.indptr[column] > block.indptr[column + 1]) {
          Rcpp::stop("Invalid LD `indptr` array.");
        }
        const int count = block.indptr[column + 1] - block.indptr[column];
        if (block.type == 0 && column + count >= block.size) {
          Rcpp::stop("LD interval extends beyond its block.");
        }
        if (block.type == 1) {
          for (int position = block.indptr[column];
               position < block.indptr[column + 1]; ++position) {
            if (block.row_index[position] <= column ||
                block.row_index[position] >= block.size) {
              Rcpp::stop("Indexed LD entries must be strict-lower triangular.");
            }
          }
        }
      }
      blocks_.push_back(std::move(block));
    }
    for (int global = 0; global < p_; ++global) {
      if (global_block_[global] < 0) {
        Rcpp::stop("LD block indices must cover every predictor.");
      }
    }
  }

  int cols() const { return p_; }
  int block_count() const { return static_cast<int>(blocks_.size()); }
  bool has_streaming_blocks() const { return true; }

  const std::vector<int>& block_predictors(const int block) const {
    return blocks_[block].global;
  }

  double diagonal(const int j) const {
    return scale_[j] * scale_[j];
  }

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
    rhs[j] -= diagonal(j) * change;
    const Block& block = blocks_[global_block_[j]];
    const int column = global_local_[j];
    const double column_scale = scale_[j];
    for (int position = block.indptr[column];
         position < block.indptr[column + 1]; ++position) {
      const int row = block.type == 0
        ? column + 1 + position - block.indptr[column]
        : block.row_index[position];
      const int global_row = block.global[row];
      rhs[global_row] -=
        column_scale * off_diagonal_scale_ * block.data[position] *
        scale_[global_row] * change;
    }
  }

  void multiply(
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const;

  void multiply_block(
      const int block_index,
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const {
    const Block& block = blocks_[block_index];
    for (int column = 0; column < block.size; ++column) {
      const int global_column = block.global[column];
      fitted[global_column] +=
        diagonal(global_column) * coefficient[global_column];
      for (int position = block.indptr[column];
           position < block.indptr[column + 1]; ++position) {
        const int row = block.type == 0
          ? column + 1 + position - block.indptr[column]
          : block.row_index[position];
        const int global_row = block.global[row];
        const double value = scale_[global_column] * off_diagonal_scale_ *
          block.data[position] * scale_[global_row];
        fitted[global_row] += value * coefficient[global_column];
        fitted[global_column] += value * coefficient[global_row];
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

  double block_quadratic(
      const std::vector<double>& coefficient,
      const std::vector<int>& predictors,
      const Rcpp::IntegerVector& block_id,
      const int prior_block) const {
    double result = 0.0;
    for (const int j : predictors) {
      result += coefficient[j] * diagonal(j) * coefficient[j];
      const Block& block = blocks_[global_block_[j]];
      const int column = global_local_[j];
      for (int position = block.indptr[column];
           position < block.indptr[column + 1]; ++position) {
        const int row = block.type == 0
          ? column + 1 + position - block.indptr[column]
          : block.row_index[position];
        const int global_row = block.global[row];
        if (block_id[global_row] - 1 == prior_block) {
          result += 2.0 * coefficient[j] * scale_[j] *
            off_diagonal_scale_ * block.data[position] * scale_[global_row] *
            coefficient[global_row];
        }
      }
    }
    return result;
  }

 private:
  struct Block {
    int type = 0;
    int size = 0;
    const double* data = NULL;
    const int* indptr = NULL;
    const int* row_index = NULL;
    std::vector<int> global;
  };

  const Rcpp::NumericVector& scale_;
  int p_;
  std::vector<Block> blocks_;
  std::vector<int> global_block_;
  std::vector<int> global_local_;
  double off_diagonal_scale_;
  int nthreads_;
};

// A block-diagonal low-rank matrix represented as G_b = Q_b' Q_b. Besides
// avoiding p_b-by-p_b Gram matrices, each block keeps its transformed residual
// w_b - Q_b beta_b, so a coefficient update costs O(q_b).
class EigenBlockSummaryMatrix {
 public:
  EigenBlockSummaryMatrix(
      const Rcpp::List& designs,
      const Rcpp::List& responses,
      const Rcpp::List& indices,
      const int p,
      const int nthreads = 1)
    : p_(p), global_block_(p, -1), global_local_(p, -1),
      diagonal_(p, 0.0), nthreads_(nthreads) {
    if (designs.size() < 1 || designs.size() != responses.size() ||
        designs.size() != indices.size()) {
      Rcpp::stop("Invalid block eigen representation.");
    }
    blocks_.reserve(designs.size());
    for (int block_index = 0; block_index < designs.size(); ++block_index) {
      Block block;
      const Rcpp::NumericMatrix design = designs[block_index];
      const Rcpp::NumericVector response = responses[block_index];
      const Rcpp::IntegerVector mapping = indices[block_index];
      block.rows = design.nrow();
      block.cols = design.ncol();
      if (block.rows < 1 || block.cols < 1 || response.size() != block.rows ||
          mapping.size() != block.cols) {
        Rcpp::stop("Eigen-block dimensions are inconsistent.");
      }
      block.design = design.begin();
      block.response.assign(response.begin(), response.end());
      block.residual = block.response;
      block.transformed_fitted.resize(block.rows, 0.0);
      block.global.resize(block.cols);
      for (int local = 0; local < block.cols; ++local) {
        const int global = mapping[local] - 1;
        if (global < 0 || global >= p_ || global_block_[global] >= 0) {
          Rcpp::stop("Eigen-block predictor indices must partition predictors.");
        }
        block.global[local] = global;
        global_block_[global] = block_index;
        global_local_[global] = local;
        const double* column = block.design +
          static_cast<std::size_t>(block.rows) * local;
        for (int row = 0; row < block.rows; ++row) {
          diagonal_[global] += column[row] * column[row];
        }
      }
      blocks_.push_back(std::move(block));
    }
    for (int global = 0; global < p_; ++global) {
      if (global_block_[global] < 0) {
        Rcpp::stop("Eigen-block predictor indices must cover every predictor.");
      }
    }
    pve_fitted_.reserve(blocks_.size());
    for (const Block& block : blocks_) {
      pve_fitted_.emplace_back(block.rows, 0.0);
    }
    pve_touched_.assign(blocks_.size(), 0);
    pve_touched_blocks_.reserve(blocks_.size());
  }

  int cols() const { return p_; }
  int block_count() const { return static_cast<int>(blocks_.size()); }
  double diagonal(const int j) const { return diagonal_[j]; }

  const std::vector<int>& block_predictors(const int block) const {
    return blocks_[block].global;
  }

  double corrected_value(
      const std::vector<double>&,
      const int j,
      const double) const {
    const Block& block = blocks_[global_block_[j]];
    const int local = global_local_[j];
    const double* column = block.design +
      static_cast<std::size_t>(block.rows) * local;
    const Eigen::Map<const Eigen::VectorXd> design_column(
      column, block.rows
    );
    const Eigen::Map<const Eigen::VectorXd> residual(
      block.residual.data(), block.rows
    );
    return design_column.dot(residual);
  }

  void update(
      std::vector<double>&,
      const int j,
      const double change) const {
    Block& block = blocks_[global_block_[j]];
    const int local = global_local_[j];
    const double* column = block.design +
      static_cast<std::size_t>(block.rows) * local;
    const Eigen::Map<const Eigen::VectorXd> design_column(
      column, block.rows
    );
    Eigen::Map<Eigen::VectorXd> residual(
      block.residual.data(), block.rows
    );
    residual.noalias() -= change * design_column;
  }

  void multiply_block(
      const int block_index,
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const {
    Block& block = blocks_[block_index];
    std::vector<double>& transformed_fitted = block.transformed_fitted;
    std::fill(transformed_fitted.begin(), transformed_fitted.end(), 0.0);
    for (int local = 0; local < block.cols; ++local) {
      const double beta = coefficient[block.global[local]];
      const double* column = block.design +
        static_cast<std::size_t>(block.rows) * local;
      for (int row = 0; row < block.rows; ++row) {
        transformed_fitted[row] += column[row] * beta;
      }
    }
    for (int row = 0; row < block.rows; ++row) {
      block.residual[row] = block.response[row] - transformed_fitted[row];
    }
    for (int local = 0; local < block.cols; ++local) {
      const double* column = block.design +
        static_cast<std::size_t>(block.rows) * local;
      double value = 0.0;
      for (int row = 0; row < block.rows; ++row) {
        value += column[row] * transformed_fitted[row];
      }
      fitted[block.global[local]] = value;
    }
  }

  void multiply(
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const;

  double center_dot(const std::vector<double>&) const { return 0.0; }
  double centered_fitted(const double fitted, const int, const double) const {
    return fitted;
  }
  void update_center_dot(double&, const int, const double) const {}

  double block_quadratic(
      const std::vector<double>& coefficient,
      const std::vector<int>& predictors,
      const Rcpp::IntegerVector&,
      const int) const {
    pve_touched_blocks_.clear();
    for (const int j : predictors) {
      const int block_index = global_block_[j];
      const Block& block = blocks_[block_index];
      const int local = global_local_[j];
      if (pve_touched_[block_index] == 0) {
        std::fill(
          pve_fitted_[block_index].begin(),
          pve_fitted_[block_index].end(),
          0.0
        );
        pve_touched_[block_index] = 1;
        pve_touched_blocks_.push_back(block_index);
      }
      const double* column = block.design +
        static_cast<std::size_t>(block.rows) * local;
      Eigen::Map<Eigen::VectorXd> fitted(
        pve_fitted_[block_index].data(), block.rows
      );
      const Eigen::Map<const Eigen::VectorXd> design_column(
        column, block.rows
      );
      fitted.noalias() += coefficient[j] * design_column;
    }
    double result = 0.0;
    for (const int block_index : pve_touched_blocks_) {
      const std::vector<double>& values = pve_fitted_[block_index];
      const Eigen::Map<const Eigen::VectorXd> fitted(
        values.data(), values.size()
      );
      result += fitted.squaredNorm();
      pve_touched_[block_index] = 0;
    }
    return result;
  }

 private:
  struct Block {
    int rows = 0;
    int cols = 0;
    const double* design = NULL;
    std::vector<int> global;
    std::vector<double> response;
    mutable std::vector<double> residual;
    mutable std::vector<double> transformed_fitted;
  };

  int p_;
  mutable std::vector<Block> blocks_;
  std::vector<int> global_block_;
  std::vector<int> global_local_;
  std::vector<double> diagonal_;
  int nthreads_;
  mutable std::vector< std::vector<double> > pve_fitted_;
  mutable std::vector<unsigned char> pve_touched_;
  mutable std::vector<int> pve_touched_blocks_;
};

}  // namespace bayeslinreg

#endif
