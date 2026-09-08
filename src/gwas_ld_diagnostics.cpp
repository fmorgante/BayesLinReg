#include <RcppEigen.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using RcppParallel::RVector;
using RcppParallel::Worker;

struct DiagnosticBlockView {
  DiagnosticBlockView(
      const Rcpp::List& block,
      const int offset)
      : data(Rcpp::as<Rcpp::NumericVector>(block["data"])),
        indptr(Rcpp::as<Rcpp::IntegerVector>(block["indptr"])),
        row_index(Rcpp::as<Rcpp::IntegerVector>(block["row_index"])),
        storage_type(Rcpp::as<int>(block["type"])),
        size(Rcpp::as<int>(block["size"])),
        offset(offset) {}

  RVector<double> data;
  RVector<int> indptr;
  RVector<int> row_index;
  int storage_type;
  int size;
  int offset;
};

class GwasLdDiagnosticWorker : public Worker {
 public:
  GwasLdDiagnosticWorker(
      const Rcpp::List& blocks,
      const Rcpp::IntegerVector& block_offset,
      const Rcpp::NumericVector& z,
      const Rcpp::IntegerVector& window_block,
      const Rcpp::IntegerVector& core_start,
      const Rcpp::IntegerVector& core_end,
      const Rcpp::IntegerVector& expanded_start,
      const Rcpp::IntegerVector& expanded_end,
      const Rcpp::IntegerVector& group_offset,
      const Rcpp::IntegerVector& group,
      const double ld_shrink,
      const double variance_floor,
      Rcpp::NumericVector& predicted,
      Rcpp::NumericVector& conditional_variance,
      Rcpp::NumericVector& tagging,
      Rcpp::NumericVector& conditional_z,
      Rcpp::NumericVector& statistic,
      Rcpp::NumericVector& flip_log_likelihood_ratio,
      Rcpp::IntegerVector& predictors_used,
      Rcpp::IntegerVector& failed)
      : z_(z),
        window_block_(window_block),
        core_start_(core_start),
        core_end_(core_end),
        expanded_start_(expanded_start),
        expanded_end_(expanded_end),
        group_offset_(group_offset),
        group_(group),
        ld_scale_(1.0 - ld_shrink),
        variance_floor_(variance_floor),
        predicted_(predicted),
        conditional_variance_(conditional_variance),
        tagging_(tagging),
        conditional_z_(conditional_z),
        statistic_(statistic),
        flip_log_likelihood_ratio_(flip_log_likelihood_ratio),
        predictors_used_(predictors_used),
        failed_(failed) {
    blocks_.reserve(blocks.size());
    for (int index = 0; index < blocks.size(); ++index) {
      blocks_.emplace_back(
        Rcpp::as<Rcpp::List>(blocks[index]), block_offset[index]
      );
    }
  }

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t window = begin; window < end; ++window) {
      diagnose_window(window);
    }
  }

 private:
  std::vector<DiagnosticBlockView> blocks_;
  RVector<double> z_;
  RVector<int> window_block_;
  RVector<int> core_start_;
  RVector<int> core_end_;
  RVector<int> expanded_start_;
  RVector<int> expanded_end_;
  RVector<int> group_offset_;
  RVector<int> group_;
  double ld_scale_;
  double variance_floor_;
  RVector<double> predicted_;
  RVector<double> conditional_variance_;
  RVector<double> tagging_;
  RVector<double> conditional_z_;
  RVector<double> statistic_;
  RVector<double> flip_log_likelihood_ratio_;
  RVector<int> predictors_used_;
  RVector<int> failed_;

  void diagnose_window(const std::size_t window) {
    const DiagnosticBlockView& block = blocks_[window_block_[window]];
    const int core_first = core_start_[window];
    const int core_last = core_end_[window];
    const int expanded_first = expanded_start_[window];
    const int expanded_last = expanded_end_[window];
    const int expanded_size = expanded_last - expanded_first + 1;
    const int group_first = group_offset_[window];

    std::vector<int> predictor_one;
    std::vector<int> predictor_two;
    std::vector<int> target_one;
    std::vector<int> target_two;
    predictor_one.reserve((expanded_size + 1) / 2);
    predictor_two.reserve((expanded_size + 1) / 2);
    target_one.reserve((core_last - core_first + 2) / 2);
    target_two.reserve((core_last - core_first + 2) / 2);
    std::vector<int> predictor_one_position(expanded_size, -1);
    std::vector<int> predictor_two_position(expanded_size, -1);
    std::vector<int> target_one_position(expanded_size, -1);
    std::vector<int> target_two_position(expanded_size, -1);

    for (int local = 0; local < expanded_size; ++local) {
      const int global = expanded_first + local;
      const int assignment = group_[group_first + local];
      if (assignment == 1) {
        predictor_one_position[local] = predictor_one.size();
        predictor_one.push_back(global);
        if (global >= core_first && global <= core_last) {
          target_one_position[local] = target_one.size();
          target_one.push_back(global);
        }
      } else {
        predictor_two_position[local] = predictor_two.size();
        predictor_two.push_back(global);
        if (global >= core_first && global <= core_last) {
          target_two_position[local] = target_two.size();
          target_two.push_back(global);
        }
      }
    }

    MatrixXd covariance_one = MatrixXd::Identity(
      predictor_one.size(), predictor_one.size()
    );
    MatrixXd covariance_two = MatrixXd::Identity(
      predictor_two.size(), predictor_two.size()
    );
    MatrixXd cross_one_two = MatrixXd::Zero(
      target_one.size(), predictor_two.size()
    );
    MatrixXd cross_two_one = MatrixXd::Zero(
      target_two.size(), predictor_one.size()
    );

    for (int column = expanded_first; column <= expanded_last; ++column) {
      const int column_local = column - expanded_first;
      const std::size_t stored_first =
        static_cast<std::size_t>(block.indptr[column]);
      const std::size_t stored_end =
        static_cast<std::size_t>(block.indptr[column + 1]);
      for (std::size_t position = stored_first;
           position < stored_end; ++position) {
        const int row = block.storage_type == 0
          ? column + 1 + static_cast<int>(position - stored_first)
          : block.row_index[position];
        if (row > expanded_last) break;
        if (row < expanded_first) continue;
        const int row_local = row - expanded_first;
        const double value = ld_scale_ * block.data[position];
        const int column_group = group_[group_first + column_local];
        const int row_group = group_[group_first + row_local];
        if (column_group == 1 && row_group == 1) {
          const int i = predictor_one_position[column_local];
          const int j = predictor_one_position[row_local];
          covariance_one(i, j) = value;
          covariance_one(j, i) = value;
        } else if (column_group == 2 && row_group == 2) {
          const int i = predictor_two_position[column_local];
          const int j = predictor_two_position[row_local];
          covariance_two(i, j) = value;
          covariance_two(j, i) = value;
        } else if (column_group == 1) {
          const int target_column = target_one_position[column_local];
          const int target_row = target_two_position[row_local];
          if (target_column >= 0) {
            cross_one_two(
              target_column, predictor_two_position[row_local]
            ) = value;
          }
          if (target_row >= 0) {
            cross_two_one(
              target_row, predictor_one_position[column_local]
            ) = value;
          }
        } else {
          const int target_column = target_two_position[column_local];
          const int target_row = target_one_position[row_local];
          if (target_column >= 0) {
            cross_two_one(
              target_column, predictor_one_position[row_local]
            ) = value;
          }
          if (target_row >= 0) {
            cross_one_two(
              target_row, predictor_two_position[column_local]
            ) = value;
          }
        }
      }
    }

    bool successful = true;
    successful = diagnose_direction(
      covariance_two, cross_one_two, predictor_two, target_one, block.offset
    ) && successful;
    successful = diagnose_direction(
      covariance_one, cross_two_one, predictor_one, target_two, block.offset
    ) && successful;
    failed_[window] = successful ? 0 : 1;
  }

  bool diagnose_direction(
      const MatrixXd& covariance,
      const MatrixXd& cross,
      const std::vector<int>& predictor,
      const std::vector<int>& target,
      const int block_offset) {
    if (target.empty() || predictor.empty()) return true;
    Eigen::LLT<MatrixXd> factor(covariance);
    if (factor.info() != Eigen::Success) return false;

    VectorXd predictor_z(predictor.size());
    for (std::size_t i = 0; i < predictor.size(); ++i) {
      predictor_z[static_cast<Eigen::Index>(i)] =
        z_[block_offset + predictor[i]];
    }
    const VectorXd weights = factor.solve(predictor_z);
    if (factor.info() != Eigen::Success || !weights.allFinite()) return false;
    const VectorXd fitted = cross * weights;

    // If covariance = L L', then h_j is the squared column norm of
    // L^{-1} cross'. The multiple-RHS triangular solve is BLAS-eligible when
    // the package is compiled with EIGEN_USE_BLAS.
    const MatrixXd whitened = factor.matrixL().solve(cross.transpose());
    const VectorXd target_tagging = whitened.colwise().squaredNorm();
    const double tolerance = std::sqrt(
      std::numeric_limits<double>::epsilon()
    );
    bool successful = true;
    for (std::size_t i = 0; i < target.size(); ++i) {
      double h = target_tagging[static_cast<Eigen::Index>(i)];
      if (!std::isfinite(h) || h < -tolerance || h > 1.0 + tolerance ||
          !std::isfinite(fitted[static_cast<Eigen::Index>(i)])) {
        successful = false;
        continue;
      }
      h = std::max(0.0, std::min(1.0, h));
      const double variance = std::max(variance_floor_, 1.0 - h);
      const int global = block_offset + target[i];
      const double residual = z_[global] - fitted[static_cast<Eigen::Index>(i)];
      const double standardized = residual / std::sqrt(variance);
      predicted_[global] = fitted[static_cast<Eigen::Index>(i)];
      conditional_variance_[global] = variance;
      tagging_[global] = h;
      conditional_z_[global] = standardized;
      statistic_[global] = standardized * standardized;
      flip_log_likelihood_ratio_[global] =
        -2.0 * z_[global] * fitted[static_cast<Eigen::Index>(i)] / variance;
      predictors_used_[global] = predictor.size();
    }
    return successful;
  }
};

}  // namespace

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]

// [[Rcpp::export]]
Rcpp::List diagnose_gwas_ld_cpp(
    const Rcpp::List& blocks,
    const Rcpp::IntegerVector& block_offset,
    const Rcpp::NumericVector& z,
    const Rcpp::IntegerVector& window_block,
    const Rcpp::IntegerVector& core_start,
    const Rcpp::IntegerVector& core_end,
    const Rcpp::IntegerVector& expanded_start,
    const Rcpp::IntegerVector& expanded_end,
    const Rcpp::IntegerVector& group_offset,
    const Rcpp::IntegerVector& group,
    const double ld_shrink,
    const double conditional_variance_floor,
    const int nthreads) {
  const int block_count = blocks.size();
  const int windows = core_start.size();
  if (block_count < 1 || block_offset.size() != block_count ||
      window_block.size() != windows || core_end.size() != windows ||
      expanded_start.size() != windows || expanded_end.size() != windows ||
      group_offset.size() != windows + 1 || nthreads < 1 ||
      group_offset[0] != 0 || group_offset[windows] != group.size()) {
    Rcpp::stop("Invalid native GWAS-LD diagnostic inputs.");
  }
  int expected_offset = 0;
  for (int index = 0; index < block_count; ++index) {
    const Rcpp::List block = blocks[index];
    const int size = Rcpp::as<int>(block["size"]);
    if (block_offset[index] != expected_offset || size < 1) {
      Rcpp::stop("Invalid native GWAS-LD diagnostic block offsets.");
    }
    expected_offset += size;
  }
  if (expected_offset != z.size()) {
    Rcpp::stop("Native GWAS-LD diagnostic blocks do not cover `z`.");
  }
  for (int window = 0; window < windows; ++window) {
    const int block_index = window_block[window];
    if (block_index < 0 || block_index >= block_count) {
      Rcpp::stop("Invalid native GWAS-LD diagnostic window block.");
    }
    const int size = Rcpp::as<int>(
      Rcpp::as<Rcpp::List>(blocks[block_index])["size"]
    );
    if (core_start[window] < 0 || core_end[window] < core_start[window] ||
        core_end[window] >= size ||
        expanded_start[window] < 0 ||
        expanded_start[window] > core_start[window] ||
        expanded_end[window] < core_end[window] ||
        expanded_end[window] >= size ||
        group_offset[window] > group_offset[window + 1] ||
        group_offset[window + 1] - group_offset[window] !=
          expanded_end[window] - expanded_start[window] + 1) {
      Rcpp::stop("Invalid native GWAS-LD diagnostic window bounds.");
    }
  }
  Rcpp::NumericVector predicted(z.size(), NA_REAL);
  Rcpp::NumericVector conditional_variance(z.size(), NA_REAL);
  Rcpp::NumericVector tagging(z.size(), NA_REAL);
  Rcpp::NumericVector conditional_z(z.size(), NA_REAL);
  Rcpp::NumericVector statistic(z.size(), NA_REAL);
  Rcpp::NumericVector flip_log_likelihood_ratio(z.size(), NA_REAL);
  Rcpp::IntegerVector predictors_used(z.size());
  Rcpp::IntegerVector failed(windows);

  GwasLdDiagnosticWorker worker(
    blocks, block_offset, z, window_block, core_start, core_end,
    expanded_start, expanded_end, group_offset, group, ld_shrink,
    conditional_variance_floor, predicted, conditional_variance, tagging,
    conditional_z, statistic, flip_log_likelihood_ratio, predictors_used,
    failed
  );
  if (nthreads > 1 && windows > 1) {
    RcppParallel::parallelFor(0, windows, worker, 1, nthreads);
  } else {
    worker(0, windows);
  }
  return Rcpp::List::create(
    Rcpp::Named("predicted_z") = predicted,
    Rcpp::Named("conditional_variance") = conditional_variance,
    Rcpp::Named("tagging") = tagging,
    Rcpp::Named("conditional_z") = conditional_z,
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("flip_log_likelihood_ratio") =
      flip_log_likelihood_ratio,
    Rcpp::Named("predictors_used") = predictors_used,
    Rcpp::Named("windows") = windows,
    Rcpp::Named("window_failed") = failed
  );
}
