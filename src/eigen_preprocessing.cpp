#include <RcppEigen.h>
#include <algorithm>
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

// Prepare a low-rank likelihood from a general factor G = Q'Q. This is used
// when predictor-specific scaling means that eigenvectors of the LD
// correlation matrix are no longer eigenvectors of the working Gram matrix.
// The transformed response is the least-squares solution of Q'w = s, and the
// returned cross-product is its projection Q'w.
// [[Rcpp::export]]
Rcpp::List prepare_eigen_factor_statistics_cpp(
    const Rcpp::NumericMatrix& factor,
    const Rcpp::NumericVector& crossproduct) {
  const int q = factor.nrow();
  const int p = factor.ncol();
  if (q < 1 || p < 1 || q > p || crossproduct.size() != p) {
    Rcpp::stop("Invalid eigen-factor preprocessing dimensions.");
  }

  const Eigen::Map<const Eigen::MatrixXd> Q(factor.begin(), q, p);
  const Eigen::Map<const Eigen::VectorXd> s(crossproduct.begin(), p);
  Eigen::MatrixXd row_gram(q, q);
  row_gram.noalias() = Q * Q.transpose();
  Eigen::LDLT<Eigen::MatrixXd> decomposition(row_gram);
  if (decomposition.info() != Eigen::Success || !decomposition.isPositive()) {
    Rcpp::stop("The scaled eigen factor is numerically rank deficient.");
  }

  Eigen::VectorXd right_hand_side(q);
  right_hand_side.noalias() = Q * s;
  Eigen::VectorXd response = decomposition.solve(right_hand_side);
  if (decomposition.info() != Eigen::Success || !response.allFinite()) {
    Rcpp::stop("The scaled eigen-factor projection could not be solved.");
  }

  Rcpp::NumericVector transformed_response(q);
  Rcpp::NumericVector projected_crossproduct(p);
  Rcpp::NumericVector diagonal(p);
  std::copy(
    response.data(), response.data() + q, transformed_response.begin()
  );
  Eigen::Map<Eigen::VectorXd> projected(
    projected_crossproduct.begin(), p
  );
  projected.noalias() = Q.transpose() * response;
  for (int predictor = 0; predictor < p; ++predictor) {
    diagonal[predictor] = Q.col(predictor).squaredNorm();
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
