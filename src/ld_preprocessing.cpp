#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

inline double symmetric_value(
    const double* matrix,
    const int p,
    const int row,
    const int column) {
  return 0.5 * matrix[row + static_cast<std::size_t>(p) * column] +
    0.5 * matrix[column + static_cast<std::size_t>(p) * row];
}

std::string storage_limit_message(
    const std::string& label,
    const int first,
    const int last,
    const std::uint64_t required,
    const std::uint64_t limit) {
  std::ostringstream message;
  message << "Dense LD block `" << label << "` contains a connected range "
          << (first + 1) << "-" << (last + 1) << " requiring at least "
          << required << " stored lower-triangular values, exceeding the "
          << "per-block limit of " << limit
          << ". Split or sparsify this block.";
  return message.str();
}

}  // namespace

// Validate, find exact contiguous zero-LD boundaries, and compress a dense LD
// matrix without allocating full-size logical, transpose, triplet, or child
// matrices. The storage limit is exposed only to make the overflow guard
// testable without constructing an enormous matrix.
// [[Rcpp::export]]
Rcpp::List compress_dense_ld_blocks_cpp(
    const Rcpp::NumericMatrix& matrix,
    const std::string& label,
    const double max_stored_values) {
  const int p = matrix.nrow();
  if (p < 1 || matrix.ncol() != p) {
    Rcpp::stop("Dense LD input must be a nonempty square matrix.");
  }
  if (!std::isfinite(max_stored_values) || max_stored_values < 1.0 ||
      max_stored_values > static_cast<double>(
        std::numeric_limits<int>::max()
      ) || std::floor(max_stored_values) != max_stored_values) {
    Rcpp::stop("Invalid dense LD compressed-storage limit.");
  }
  const std::uint64_t storage_limit =
    static_cast<std::uint64_t>(max_stored_values);
  const double* values = matrix.begin();
  std::vector<int> nonzero_count(p, 0);
  std::vector<int> last_nonzero_row(p);
  std::vector<int> boundaries;
  boundaries.reserve(p);

  double maximum_raw_value = 0.0;
  double maximum_asymmetry = 0.0;
  double maximum_correlation = 1.0;
  double maximum_diagonal_deviation = 0.0;
  for (int column = 0; column < p; ++column) {
    const double diagonal = values[
      column + static_cast<std::size_t>(p) * column
    ];
    if (!std::isfinite(diagonal)) {
      Rcpp::stop("`R[[\"%s\"]]` must be a finite numeric square matrix.",
                 label.c_str());
    }
    maximum_raw_value = std::max(maximum_raw_value, std::abs(diagonal));
    maximum_diagonal_deviation = std::max(
      maximum_diagonal_deviation, std::abs(diagonal - 1.0)
    );
    last_nonzero_row[column] = column;
  }

  int furthest_connected_row = -1;
  for (int column = 0; column < p - 1; ++column) {
    if ((column & 255) == 0) Rcpp::checkUserInterrupt();
    const double* lower_column =
      values + static_cast<std::size_t>(p) * column;
    for (int row = column + 1; row < p; ++row) {
      const double lower = lower_column[row];
      const double upper = values[
        column + static_cast<std::size_t>(p) * row
      ];
      if (!std::isfinite(lower) || !std::isfinite(upper)) {
        Rcpp::stop("`R[[\"%s\"]]` must be a finite numeric square matrix.",
                   label.c_str());
      }
      maximum_raw_value = std::max(
        maximum_raw_value, std::max(std::abs(lower), std::abs(upper))
      );
      maximum_asymmetry = std::max(
        maximum_asymmetry, std::abs(lower - upper)
      );
      const double correlation = 0.5 * lower + 0.5 * upper;
      maximum_correlation = std::max(
        maximum_correlation, std::abs(correlation)
      );
      if (correlation != 0.0) {
        ++nonzero_count[column];
        last_nonzero_row[column] = row;
      }
    }
    furthest_connected_row = std::max(
      furthest_connected_row, last_nonzero_row[column]
    );
    if (furthest_connected_row <= column) boundaries.push_back(column);
  }
  if (boundaries.empty() || boundaries.back() != p - 1) {
    boundaries.push_back(p - 1);
  }

  const double symmetry_tolerance =
    std::sqrt(std::numeric_limits<double>::epsilon()) *
      std::max(1.0, maximum_raw_value);
  if (maximum_asymmetry > symmetry_tolerance) {
    Rcpp::stop("`R[[\"%s\"]]` must be symmetric.", label.c_str());
  }
  if (maximum_diagonal_deviation > 1e-6) {
    Rcpp::stop("LD block `%s` must have a unit diagonal.", label.c_str());
  }
  if (maximum_correlation >
      1.0 + std::sqrt(std::numeric_limits<double>::epsilon())) {
    Rcpp::stop("LD block `%s` contains a correlation outside [-1, 1].",
               label.c_str());
  }

  Rcpp::List output(boundaries.size());
  int first = 0;
  for (std::size_t block_index = 0;
       block_index < boundaries.size(); ++block_index) {
    const int last = boundaries[block_index];
    const int size = last - first + 1;
    std::uint64_t indexed_values = 0;
    std::uint64_t interval_values = 0;
    for (int column = first; column <= last; ++column) {
      indexed_values += static_cast<std::uint64_t>(nonzero_count[column]);
      interval_values += static_cast<std::uint64_t>(
        last_nonzero_row[column] - column
      );
    }
    const bool indexed_valid = indexed_values <= storage_limit;
    const bool interval_valid = interval_values <= storage_limit;
    if (!indexed_valid && !interval_valid) {
      Rcpp::stop(storage_limit_message(
        label, first, last,
        std::min(indexed_values, interval_values), storage_limit
      ));
    }
    const long double pointer_bytes =
      4.0L * static_cast<long double>(size + 1);
    const long double indexed_bytes =
      12.0L * static_cast<long double>(indexed_values) + pointer_bytes;
    const long double interval_bytes =
      8.0L * static_cast<long double>(interval_values) + pointer_bytes;
    const bool use_interval = interval_valid &&
      (!indexed_valid || interval_bytes <= indexed_bytes);
    const int stored_values = static_cast<int>(
      use_interval ? interval_values : indexed_values
    );

    Rcpp::NumericVector data(stored_values);
    Rcpp::IntegerVector indptr(size + 1);
    Rcpp::IntegerVector row_index = use_interval
      ? Rcpp::IntegerVector(0)
      : Rcpp::IntegerVector(stored_values);
    int offset = 0;
    for (int column = first; column <= last; ++column) {
      const int local_column = column - first;
      indptr[local_column] = offset;
      if (use_interval) {
        const int final_row = last_nonzero_row[column];
        for (int row = column + 1; row <= final_row; ++row) {
          data[offset++] = symmetric_value(values, p, row, column);
        }
      } else {
        for (int row = column + 1; row <= last; ++row) {
          const double correlation = symmetric_value(values, p, row, column);
          if (correlation != 0.0) {
            data[offset] = correlation;
            row_index[offset] = row - first;
            ++offset;
          }
        }
      }
    }
    indptr[size] = offset;
    output[block_index] = Rcpp::List::create(
      Rcpp::Named("size") = size,
      Rcpp::Named("type") = use_interval ? 0 : 1,
      Rcpp::Named("storage") = use_interval
        ? "interval_triangular"
        : "indexed_triangular",
      Rcpp::Named("data") = data,
      Rcpp::Named("indptr") = indptr,
      Rcpp::Named("row_index") = row_index
    );
    first = last + 1;
  }
  return output;
}
