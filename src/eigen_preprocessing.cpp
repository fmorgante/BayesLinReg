#include <RcppEigen.h>
#include <cmath>
#include <cstddef>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

// Compute quantities derived from one eigen block without materializing
// squared eigenvectors or other p-by-q temporary matrices.
// [[Rcpp::export]]
Rcpp::List prepare_eigen_statistics_cpp(
    const Rcpp::NumericMatrix& eigenvectors,
    const Rcpp::NumericVector& eigenvalues,
    const Rcpp::NumericVector& crossproduct) {
  const int p = eigenvectors.nrow();
  const int q = eigenvectors.ncol();
  if (p < 1 || q < 1 || eigenvalues.size() != q ||
      crossproduct.size() != p) {
    Rcpp::stop("Invalid eigen preprocessing dimensions.");
  }

  const Eigen::Map<const Eigen::MatrixXd> vectors(
    eigenvectors.begin(), p, q
  );
  const Eigen::Map<const Eigen::VectorXd> crossproduct_vector(
    crossproduct.begin(), p
  );
  Eigen::VectorXd projected(q);
  projected.noalias() = vectors.transpose() * crossproduct_vector;

  Rcpp::NumericVector transformed_response(q);
  Rcpp::NumericVector projected_crossproduct(p);
  Rcpp::NumericVector diagonal(p);
  Eigen::Map<Eigen::VectorXd> projected_crossproduct_vector(
    projected_crossproduct.begin(), p
  );
  projected_crossproduct_vector.noalias() = vectors * projected;

  for (int component = 0; component < q; ++component) {
    const double eigenvalue = eigenvalues[component];
    if (!std::isfinite(eigenvalue) || eigenvalue <= 0.0) {
      Rcpp::stop("Eigen preprocessing requires positive finite eigenvalues.");
    }
    transformed_response[component] =
      projected[component] / std::sqrt(eigenvalue);
    const double* vector_column =
      eigenvectors.begin() + static_cast<std::size_t>(p) * component;
    for (int predictor = 0; predictor < p; ++predictor) {
      const double value = vector_column[predictor];
      diagonal[predictor] += eigenvalue * value * value;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("transformed_response") = transformed_response,
    Rcpp::Named("projected_crossproduct") = projected_crossproduct,
    Rcpp::Named("diagonal") = diagonal
  );
}

// Write Q = Lambda^(1/2) U' directly in final sampler order and scale,
// avoiding transpose, multiplication, sweep, and subsetting temporaries in R.
// An empty source_order selects the original predictor order.
// [[Rcpp::export]]
Rcpp::NumericMatrix build_scaled_eigen_factor_cpp(
    const Rcpp::NumericMatrix& eigenvectors,
    const Rcpp::NumericVector& eigenvalues,
    const Rcpp::NumericVector& predictor_scale,
    const Rcpp::IntegerVector& source_order) {
  const int p = eigenvectors.nrow();
  const int q = eigenvectors.ncol();
  const int output_predictors = source_order.size() == 0
    ? p
    : source_order.size();
  if (p < 1 || q < 1 || eigenvalues.size() != q ||
      predictor_scale.size() != output_predictors) {
    Rcpp::stop("Invalid scaled eigen-factor dimensions.");
  }

  std::vector<double> square_root_eigenvalue(q);
  for (int component = 0; component < q; ++component) {
    const double eigenvalue = eigenvalues[component];
    if (!std::isfinite(eigenvalue) || eigenvalue <= 0.0) {
      Rcpp::stop("Eigen-factor construction requires positive eigenvalues.");
    }
    square_root_eigenvalue[component] = std::sqrt(eigenvalue);
  }

  Rcpp::NumericMatrix factor(q, output_predictors);
  for (int output = 0; output < output_predictors; ++output) {
    const int source = source_order.size() == 0
      ? output
      : source_order[output] - 1;
    const double scale = predictor_scale[output];
    if (source < 0 || source >= p || !std::isfinite(scale) || scale <= 0.0) {
      Rcpp::stop("Invalid predictor order or scale in eigen preprocessing.");
    }
    double* factor_column =
      factor.begin() + static_cast<std::size_t>(q) * output;
    for (int component = 0; component < q; ++component) {
      factor_column[component] = square_root_eigenvalue[component] *
        eigenvectors.begin()[
          source + static_cast<std::size_t>(p) * component
        ] / scale;
    }
  }
  return factor;
}
