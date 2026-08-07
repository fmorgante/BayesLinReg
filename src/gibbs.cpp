#include <RcppEigen.h>
#include <RcppParallel.h>
#include <R_ext/Rdynload.h>
#include <GIGrvg.h>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <type_traits>
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
// general sparse, or symmetric sparse. Symmetric sparse blocks retain one
// triangle and build a reverse index that refers back to the original values.
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
      nthreads_(nthreads) {
    if (matrices.size() < 1 || matrices.size() != indices.size() ||
        matrices.size() != types.size()) {
      Rcpp::stop("Invalid block-diagonal `XtX` representation.");
    }
    blocks_.reserve(matrices.size());
    for (int block_index = 0; block_index < matrices.size(); ++block_index) {
      Block block;
      block.type = types[block_index];
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

      if (block.type == 0) {
        const Rcpp::NumericMatrix matrix = matrices[block_index];
        if (matrix.nrow() != block.size || matrix.ncol() != block.size) {
          Rcpp::stop("Dense Gram-block dimensions do not match their indices.");
        }
        block.dense = matrix.begin();
        for (int local = 0; local < block.size; ++local) {
          diagonal_[block.global[local]] =
            block.dense[local + static_cast<std::size_t>(block.size) * local];
        }
      } else if (block.type == 1 || block.type == 2) {
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
        if (block.type == 2) {
          build_reverse_index(block);
        }
      } else {
        Rcpp::stop("Unknown Gram-block storage type.");
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
    if (block.type == 0) {
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
    if (block.type == 2) {
      for (int reverse = block.reverse_pointer[local];
           reverse < block.reverse_pointer[local + 1]; ++reverse) {
        rhs[block.global[block.reverse_neighbor[reverse]]] -=
          block.values[block.reverse_position[reverse]] * change;
      }
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
    if (block.type == 0) {
      for (int column = 0; column < block.size; ++column) {
        const double beta = coefficient[block.global[column]];
        const double* values = block.dense +
          static_cast<std::size_t>(block.size) * column;
        for (int row = 0; row < block.size; ++row) {
          fitted[block.global[row]] += values[row] * beta;
        }
      }
    } else if (block.type == 1) {
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
      if (block.type == 0) {
        const double* column = block.dense +
          static_cast<std::size_t>(block.size) * local;
        for (int row = 0; row < block.size; ++row) {
          const int global_row = block.global[row];
          if (block_id[global_row] - 1 == prior_block) {
            result += coefficient[global_row] * column[row] * coefficient[j];
          }
        }
      } else if (block.type == 1) {
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
    int type = -1;
    int size = 0;
    int nonzeros = 0;
    const double* dense = NULL;
    const int* column_pointer = NULL;
    const int* row_index = NULL;
    const double* values = NULL;
    std::vector<int> global;
    std::vector<int> reverse_pointer;
    std::vector<int> reverse_neighbor;
    std::vector<int> reverse_position;
  };

  static void build_reverse_index(Block& block) {
    block.reverse_pointer.assign(block.size + 1, 0);
    for (int column = 0; column < block.size; ++column) {
      for (int position = block.column_pointer[column];
           position < block.column_pointer[column + 1]; ++position) {
        const int row = block.row_index[position];
        if (row != column) ++block.reverse_pointer[row + 1];
      }
    }
    for (int row = 0; row < block.size; ++row) {
      block.reverse_pointer[row + 1] += block.reverse_pointer[row];
    }
    const int off_diagonal = block.reverse_pointer[block.size];
    block.reverse_neighbor.resize(off_diagonal);
    block.reverse_position.resize(off_diagonal);
    std::vector<int> next = block.reverse_pointer;
    for (int column = 0; column < block.size; ++column) {
      for (int position = block.column_pointer[column];
           position < block.column_pointer[column + 1]; ++position) {
        const int row = block.row_index[position];
        if (row == column) continue;
        const int destination = next[row]++;
        block.reverse_neighbor[destination] = column;
        block.reverse_position[destination] = position;
      }
    }
  }

  const Rcpp::NumericVector& center_;
  int p_;
  std::vector<Block> blocks_;
  std::vector<int> global_block_;
  std::vector<int> global_local_;
  std::vector<double> diagonal_;
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
    double result = 0.0;
    for (int row = 0; row < block.rows; ++row) {
      result += column[row] * block.residual[row];
    }
    return result;
  }

  void update(
      std::vector<double>&,
      const int j,
      const double change) const {
    Block& block = blocks_[global_block_[j]];
    const int local = global_local_[j];
    const double* column = block.design +
      static_cast<std::size_t>(block.rows) * local;
    for (int row = 0; row < block.rows; ++row) {
      block.residual[row] -= column[row] * change;
    }
  }

  void multiply_block(
      const int block_index,
      const std::vector<double>& coefficient,
      std::vector<double>& fitted) const {
    Block& block = blocks_[block_index];
    std::vector<double> transformed_fitted(block.rows, 0.0);
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

void BlockSummaryMatrix::multiply(
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

void EigenBlockSummaryMatrix::multiply(
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

std::uint64_t splitmix64(std::uint64_t value) {
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
      const int* block_model,
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
        const int model = block_model_[prior_block];
        const int local_index = model_local_index_[j];
        const double old_coefficient = coefficient_[j];
        const double partial_rhs = matrix_.corrected_value(
          corrected_rhs_, j, 0.0
        ) +
          x_squared_[j] * old_coefficient;

        if (model == 3) {
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
          const double prior_precision = model == 4
            ? 0.0
            : (model == 2
                ? 1.0 / tau_sq_[prior_block] / local_var_[local_index]
                : (model == 1
                    ? 1.0 / slab_var_[prior_block]
                    : 1.0 / normal_var_[prior_block]));
          const double conditional_var = 1.0 / (
            x_squared_[j] / residual_var_ + prior_precision
          );
          const double conditional_mean =
            conditional_var * partial_rhs / residual_var_;
          if (model == 1) {
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
          coefficient_[j] = model != 1 || inclusion_[local_index] == 1
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
  const int* block_model_;
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
    const Rcpp::IntegerVector& block_model,
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
    matrix, block_id.begin(), block_model.begin(), model_local_index,
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

}  // namespace

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]

// [[Rcpp::export]]
Rcpp::List blm_build_info_cpp() {
#ifdef EIGEN_USE_BLAS
  const bool eigen_blas = true;
#else
  const bool eigen_blas = false;
#endif
  return Rcpp::List::create(
    Rcpp::Named("eigen_blas") = eigen_blas,
    Rcpp::Named("dense_kernel") = eigen_blas ? "external BLAS" : "Eigen"
  );
}

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
  std::vector< std::vector<int> > block_predictors(number_of_blocks);
  std::vector<int> model_local_index(p, -1);
  std::vector<int> block_local_index(p, -1);
  int model_size[5] = {0, 0, 0, 0, 0};
  bool has_normal = false;
  bool has_spike_slab = false;
  bool has_global_local = false;
  bool has_spike_multi_slab = false;
  for (int block = 0; block < number_of_blocks; ++block) {
    has_normal = has_normal || block_model[block] == 0;
    has_spike_slab = has_spike_slab || block_model[block] == 1;
    has_global_local = has_global_local || block_model[block] == 2;
    has_spike_multi_slab =
      has_spike_multi_slab || block_model[block] == 3;
  }
  for (int j = 0; j < p; ++j) {
    const int block = block_id[j] - 1;
    const int model = block_model[block];
    block_local_index[j] =
      static_cast<int>(block_predictors[block].size());
    block_predictors[block].push_back(j);
    model_local_index[j] = model_size[model];
    ++model_size[model];
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
    stored_rows, model_size[1]
  );
  Rcpp::NumericMatrix pi_samples(
    stored_rows, has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericMatrix slab_var_samples(
    stored_rows, has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericMatrix local_var_samples(
    stored_rows, model_size[2]
  );
  Rcpp::NumericMatrix tau_sq_samples(
    stored_rows, has_global_local ? number_of_blocks : 0
  );
  Rcpp::IntegerMatrix multi_component_samples(
    stored_rows, model_size[3]
  );
  Rcpp::List multi_pi_samples(number_of_blocks);
  Rcpp::NumericMatrix multi_var_samples(
    stored_rows, has_spike_multi_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector coefficient_sum(p);
  Rcpp::NumericVector coefficient_sum_sq(p);
  // Stored draws are used to compute their covariance after sampling.
  const int covariance_dimension =
    !store_samples && store_coefficient_cov ? p : 0;
  Rcpp::NumericMatrix coefficient_crossprod(
    covariance_dimension, covariance_dimension
  );
  double intercept_sum = 0.0;
  double intercept_sum_sq = 0.0;
  double residual_var_sum = 0.0;
  double residual_var_sum_sq = 0.0;
  Rcpp::NumericVector block_pve_sum(
    compute_pve ? number_of_blocks : 0
  );
  Rcpp::NumericVector block_pve_sum_sq(
    compute_pve ? number_of_blocks : 0
  );
  double total_pve_sum = 0.0;
  double total_pve_sum_sq = 0.0;
  double cross_block_pve_sum = 0.0;
  double cross_block_pve_sum_sq = 0.0;
  Rcpp::NumericVector normal_var_sum(
    has_normal ? number_of_blocks : 0
  );
  Rcpp::NumericVector normal_var_sum_sq(
    has_normal ? number_of_blocks : 0
  );
  Rcpp::NumericVector inclusion_sum(model_size[1]);
  Rcpp::NumericVector pi_sum(
    has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector pi_sum_sq(
    has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector slab_var_sum(
    has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector slab_var_sum_sq(
    has_spike_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector local_var_sum(model_size[2]);
  Rcpp::NumericVector local_var_sum_sq(model_size[2]);
  Rcpp::NumericVector tau_sq_sum(
    has_global_local ? number_of_blocks : 0
  );
  Rcpp::NumericVector tau_sq_sum_sq(
    has_global_local ? number_of_blocks : 0
  );
  Rcpp::List multi_component_sum(number_of_blocks);
  Rcpp::List multi_pi_sum(number_of_blocks);
  Rcpp::List multi_pi_sum_sq(number_of_blocks);
  Rcpp::NumericVector multi_var_sum(
    has_spike_multi_slab ? number_of_blocks : 0
  );
  Rcpp::NumericVector multi_var_sum_sq(
    has_spike_multi_slab ? number_of_blocks : 0
  );
  std::vector<double> coefficient(p, 0.0);
  std::vector<int> inclusion(model_size[1], 1);
  std::vector<int> multi_component(model_size[3], 0);
  std::vector<double> local_var(model_size[2], 1.0);
  std::vector<double> local_aux(model_size[2], 1.0);
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
    if (block_model[block] == 0) {
      normal_var[block] =
        normal_scale[block] / (normal_shape[block] + 1.0);
    }
    if (block_model[block] == 1) {
      pi[block] = pi_alpha[block] / (pi_alpha[block] + pi_beta[block]);
      slab_var[block] =
        spike_var_scale[block] / (spike_var_shape[block] + 1.0);
    }
    if (block_model[block] == 2) {
      tau_sq[block] = global_scale[block] * global_scale[block];
    }
    if (block_model[block] == 3) {
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
        const std::vector<int>& predictors = block_predictors[block];
        multi_component_sum[block] = Rcpp::NumericMatrix(
          predictors.size(), alpha_values.size()
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
  const int residual_refresh_interval = 100;
  std::vector<BlockRng> block_rng;
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
          summary_XtX, block_id, block_model, model_local_index, x_squared,
          residual_var, normal_var, pi, slab_var, tau_sq, local_var,
          multi_gamma, multi_pi, multi_var, learn_residual_var, coefficient,
          corrected_rhs, inclusion, multi_component, block_rng, residual_sse,
          nthreads
        );
        parallel_sweep = true;
      }
    }
    if (!parallel_sweep) {
      for (int j = 0; j < p; ++j) {
      const int block = block_id[j] - 1;
      const int model = block_model[block];
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
        multi_component[local_index] = selected_component;
        coefficient[j] = selected_component == 0
          ? 0.0
          : R::rnorm(
              conditional_means[selected_component],
              std::sqrt(conditional_vars[selected_component])
            );
      } else {
        const double prior_precision = model == 4
          ? 0.0
          : (model == 2
              ? 1.0 / tau_sq[block] / local_var[local_index]
              : (model == 1
                  ? 1.0 / slab_var[block]
                  : 1.0 / normal_var[block]));
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
          inclusion[local_index] = static_cast<int>(
            R::rbinom(1.0, inclusion_probability)
          );
        }
        if (model != 1 || inclusion[local_index] == 1) {
          coefficient[j] = R::rnorm(
            conditional_mean,
            std::sqrt(conditional_var)
          );
        } else {
          coefficient[j] = 0.0;
        }
      }
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

    // Reconstruct accumulated state periodically to limit floating-point drift.
    if (iteration % residual_refresh_interval == 0) {
      if (use_sufficient_statistics) {
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
        if (block_model[block] != 0) {
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
        if (block_model[block] != 1) {
          continue;
        }
        const std::vector<int>& predictors = block_predictors[block];
        int number_included = 0;
        double included_sum_of_squares = 0.0;
        for (std::size_t index = 0; index < predictors.size(); ++index) {
          const int j = predictors[index];
          const int local_index = model_local_index[j];
          number_included += inclusion[local_index];
          if (inclusion[local_index] == 1) {
            included_sum_of_squares += coefficient[j] * coefficient[j];
          }
        }
        pi[block] = R::rbeta(
          pi_alpha[block] + number_included,
          pi_beta[block] + predictors.size() - number_included
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
        const std::vector<int>& predictors = block_predictors[block];
        std::vector<int> counts(component_count, 0);
        int number_nonzero = 0;
        double scaled_sum_of_squares = 0.0;
        for (std::size_t index = 0; index < predictors.size(); ++index) {
          const int j = predictors[index];
          const int component = multi_component[model_local_index[j]];
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
          local_var[local_index] = draw_gig(
            local_a[block] - 0.5,
            chi,
            2.0 * local_aux[local_index]
          );
          local_aux[local_index] = R::rgamma(
            local_a[block] + local_b[block],
            1.0 / (1.0 + local_var[local_index])
          );
        }

        double tau_rate = 1.0 / global_aux[block];
        for (std::size_t index = 0; index < predictors.size(); ++index) {
          const int j = predictors[index];
          tau_rate += coefficient[j] * coefficient[j] /
            (2.0 * local_var[model_local_index[j]]);
        }
        tau_sq[block] = 1.0 / R::rgamma(
          (static_cast<double>(predictors.size()) + 1.0) / 2.0,
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
        const double variance_df = static_cast<double>(
          effective_n - (fit_intercept ? 1 : 0)
        );
        const double total_signal_variance = total_sum_squares / variance_df;
        const double pve_denominator = total_signal_variance + residual_var;
        double standalone_total = 0.0;
        for (int block = 0; block < number_of_blocks; ++block) {
          standalone_total += standalone_sum_squares[block];
          const double block_sum_squares = pve_type_code == 0
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
      double intercept_mean = intercept_y_mean;
      for (int j = 0; j < p; ++j) {
        const int block = block_id[j] - 1;
        const int model = block_model[block];
        const int local_index = model_local_index[j];
        if (store_samples) {
          coefficient_samples(retained_index, j) = coefficient[j];
          if (model == 1) {
            inclusion_samples(retained_index, local_index) =
              inclusion[local_index];
          }
          if (model == 2) {
            local_var_samples(retained_index, local_index) =
              local_var[local_index];
          }
          if (model == 3) {
            multi_component_samples(retained_index, local_index) =
              multi_component[local_index] + 1;
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
        if (compute_pve) {
          for (int block = 0; block < number_of_blocks; ++block) {
            block_pve_samples(retained_index, block) = block_pve[block];
          }
          total_pve_samples[retained_index] = total_pve;
          cross_block_pve_samples[retained_index] = cross_block_pve;
        }
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
          const int block = block_id[j] - 1;
          const int model = block_model[block];
          const int local_index = model_local_index[j];
          if (model == 1) {
            inclusion_sum[local_index] += inclusion[local_index];
          } else if (model == 2) {
            local_var_sum[local_index] += local_var[local_index];
            local_var_sum_sq[local_index] +=
              local_var[local_index] * local_var[local_index];
          } else if (model == 3) {
            Rcpp::NumericMatrix block_component_sum =
              multi_component_sum[block];
            block_component_sum(
              block_local_index[j], multi_component[local_index]
            ) += 1.0;
          }
        }
        if (store_coefficient_cov) {
          Eigen::Map<Eigen::MatrixXd> crossproduct(
            coefficient_crossprod.begin(), p, p
          );
          const Eigen::Map<const Eigen::VectorXd> coefficient_vector(
            coefficient.data(), p
          );
          crossproduct.noalias() +=
            coefficient_vector * coefficient_vector.transpose();
        }
        intercept_sum += intercept_draw;
        intercept_sum_sq += intercept_draw * intercept_draw;
        residual_var_sum += residual_var;
        residual_var_sum_sq += residual_var * residual_var;
        if (compute_pve) {
          for (int block = 0; block < number_of_blocks; ++block) {
            block_pve_sum[block] += block_pve[block];
            block_pve_sum_sq[block] += block_pve[block] * block_pve[block];
          }
          total_pve_sum += total_pve;
          total_pve_sum_sq += total_pve * total_pve;
          cross_block_pve_sum += cross_block_pve;
          cross_block_pve_sum_sq +=
            cross_block_pve * cross_block_pve;
        }
        for (int block = 0; block < number_of_blocks; ++block) {
          const int model = block_model[block];
          if (model == 0) {
            normal_var_sum[block] += normal_var[block];
            normal_var_sum_sq[block] +=
              normal_var[block] * normal_var[block];
          } else if (model == 1) {
            pi_sum[block] += pi[block];
            pi_sum_sq[block] += pi[block] * pi[block];
            slab_var_sum[block] += slab_var[block];
            slab_var_sum_sq[block] += slab_var[block] * slab_var[block];
          } else if (model == 2) {
            tau_sq_sum[block] += tau_sq[block];
            tau_sq_sum_sq[block] += tau_sq[block] * tau_sq[block];
          } else if (model == 3) {
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
      Rcpp::Named("multi_var_samples") = multi_var_samples,
      Rcpp::Named("block_pve_samples") = block_pve_samples,
      Rcpp::Named("total_pve_samples") = total_pve_samples,
      Rcpp::Named("cross_block_pve_samples") = cross_block_pve_samples
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
    Rcpp::Named("multi_var_sum_sq") = multi_var_sum_sq,
    Rcpp::Named("block_pve_sum") = block_pve_sum,
    Rcpp::Named("block_pve_sum_sq") = block_pve_sum_sq,
    Rcpp::Named("total_pve_sum") = total_pve_sum,
    Rcpp::Named("total_pve_sum_sq") = total_pve_sum_sq,
    Rcpp::Named("cross_block_pve_sum") = cross_block_pve_sum,
    Rcpp::Named("cross_block_pve_sum_sq") = cross_block_pve_sum_sq
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
    const double summary_yty,
    const bool center_observations,
    const double residual_sse_offset,
    const bool compute_pve,
    const int pve_type_code) {
  const DenseSummaryMatrix summary_matrix(summary_XtX);
  return blm_gibbs_core(
    y, X, residual_shape, residual_scale, iterations, burnin, thin,
    progress_callback, block_id, block_model, normal_shape, normal_scale,
    pi_alpha, pi_beta, spike_var_shape, spike_var_scale, global_scale,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, fit_intercept, intercept_x_mean,
    intercept_y_mean, use_sufficient_statistics, summary_matrix, summary_Xty,
    summary_yty, center_observations, residual_sse_offset, compute_pve,
    pve_type_code, 1
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
    const double summary_yty,
    const bool compute_pve,
    const int pve_type_code) {
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
    intercept_y_mean, true, summary_matrix, summary_Xty, summary_yty,
    true, 0.0, compute_pve, pve_type_code, 1
  );
}

// [[Rcpp::export]]
Rcpp::List blm_gibbs_block_rcpp_cpp(
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
    const Rcpp::List& summary_XtX,
    const Rcpp::List& summary_indices,
    const Rcpp::IntegerVector& summary_types,
    const Rcpp::NumericVector& summary_center,
    const Rcpp::NumericVector& summary_Xty,
    const double summary_yty,
    const bool compute_pve,
    const int pve_type_code,
    const int nthreads) {
  if (nthreads < 1) {
    Rcpp::stop("`nthreads` must be a positive integer.");
  }
  if (nthreads > 1) {
    for (int j = 0; j < summary_center.size(); ++j) {
      if (summary_center[j] != 0.0) {
        Rcpp::stop(
          "Threaded Gram-block sampling requires zero predictor means."
        );
      }
    }
  }
  const Rcpp::NumericVector empty_y;
  const Rcpp::NumericMatrix empty_X(0, 0);
  const BlockSummaryMatrix summary_matrix(
    summary_XtX, summary_indices, summary_types, summary_center, nthreads
  );
  return blm_gibbs_core(
    empty_y, empty_X, residual_shape, residual_scale, iterations, burnin, thin,
    progress_callback, block_id, block_model, normal_shape, normal_scale,
    pi_alpha, pi_beta, spike_var_shape, spike_var_scale, global_scale,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, fit_intercept, intercept_x_mean,
    intercept_y_mean, true, summary_matrix, summary_Xty, summary_yty,
    true, 0.0, compute_pve, pve_type_code, nthreads
  );
}

// [[Rcpp::export]]
Rcpp::List blm_gibbs_eigen_block_rcpp_cpp(
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
    const Rcpp::List& transformed_X,
    const Rcpp::List& transformed_y,
    const Rcpp::List& eigen_indices,
    const Rcpp::NumericVector& summary_Xty,
    const double summary_yty,
    const bool compute_pve,
    const int pve_type_code,
    const int nthreads) {
  if (nthreads < 1) {
    Rcpp::stop("`nthreads` must be a positive integer.");
  }
  const Rcpp::NumericVector empty_y;
  const Rcpp::NumericMatrix empty_X(0, 0);
  const EigenBlockSummaryMatrix summary_matrix(
    transformed_X, transformed_y, eigen_indices, summary_Xty.size(), nthreads
  );
  return blm_gibbs_core(
    empty_y, empty_X, residual_shape, residual_scale, iterations, burnin, thin,
    progress_callback, block_id, block_model, normal_shape, normal_scale,
    pi_alpha, pi_beta, spike_var_shape, spike_var_scale, global_scale,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, fit_intercept, intercept_x_mean,
    intercept_y_mean, true, summary_matrix, summary_Xty, summary_yty,
    true, 0.0, compute_pve, pve_type_code, nthreads
  );
}
