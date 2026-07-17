#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <vector>

// [[Rcpp::export]]
Rcpp::List blm_gibbs_rcpp_cpp(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& X,
    const Rcpp::NumericVector& prior_var,
    const double residual_shape,
    const double residual_scale,
    const int iterations,
    const int burnin,
    const int thin,
    const bool verbose) {
  Rcpp::RNGScope scope;

  const int n = y.size();
  const int p = X.ncol();
  const int number_of_draws = (iterations - burnin - 1) / thin + 1;

  std::vector<double> x_mean(p, 0.0);
  double y_mean = 0.0;
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

  Rcpp::NumericMatrix x_centered(n, p);
  Rcpp::NumericVector y_centered(n);
  for (int i = 0; i < n; ++i) {
    y_centered[i] = y[i] - y_mean;
    for (int j = 0; j < p; ++j) {
      x_centered(i, j) = X(i, j) - x_mean[j];
    }
  }

  std::vector<double> x_squared(p, 0.0);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) {
      x_squared[j] += x_centered(i, j) * x_centered(i, j);
    }
  }

  Rcpp::NumericMatrix coefficient_samples(number_of_draws, p);
  Rcpp::NumericVector intercept_samples(number_of_draws);
  Rcpp::NumericVector residual_var_samples(number_of_draws);
  std::vector<double> coefficient(p, 0.0);
  std::vector<double> residuals(n);
  for (int i = 0; i < n; ++i) {
    residuals[i] = y_centered[i];
  }

  double residual_var = residual_scale / (residual_shape + 1.0);
  const double posterior_shape =
    residual_shape + static_cast<double>(n - 1) / 2.0;
  int retained_index = 0;

  for (int iteration = 1; iteration <= iterations; ++iteration) {
    if (verbose) {
      Rcpp::Rcout << "\rGibbs iteration " << iteration << "/"
                  << iterations << std::flush;
    }

    // Update each coefficient from its univariate conditional normal.
    for (int j = 0; j < p; ++j) {
      double conditional_numerator = 0.0;
      for (int i = 0; i < n; ++i) {
        residuals[i] += x_centered(i, j) * coefficient[j];
        conditional_numerator += x_centered(i, j) * residuals[i];
      }
      const double conditional_var = 1.0 / (
        x_squared[j] / residual_var + 1.0 / prior_var[j]
      );
      const double conditional_mean =
        conditional_var * conditional_numerator / residual_var;
      coefficient[j] = R::rnorm(
        conditional_mean,
        std::sqrt(conditional_var)
      );
      for (int i = 0; i < n; ++i) {
        residuals[i] -= x_centered(i, j) * coefficient[j];
      }
    }

    double sum_squared_residuals = 0.0;
    for (int i = 0; i < n; ++i) {
      sum_squared_residuals += residuals[i] * residuals[i];
    }
    const double posterior_scale =
      residual_scale + 0.5 * sum_squared_residuals;
    residual_var = 1.0 / R::rgamma(posterior_shape, 1.0 / posterior_scale);

    if (iteration > burnin &&
        (iteration - burnin - 1) % thin == 0) {
      double intercept_mean = y_mean;
      for (int j = 0; j < p; ++j) {
        coefficient_samples(retained_index, j) = coefficient[j];
        intercept_mean -= x_mean[j] * coefficient[j];
      }
      intercept_samples[retained_index] = R::rnorm(
        intercept_mean,
        std::sqrt(residual_var / n)
      );
      residual_var_samples[retained_index] = residual_var;
      ++retained_index;
    }

    if (iteration % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  if (verbose) {
    Rcpp::Rcout << std::endl;
  }

  return Rcpp::List::create(
    Rcpp::Named("coefficient_samples") = coefficient_samples,
    Rcpp::Named("intercept_samples") = intercept_samples,
    Rcpp::Named("residual_var_samples") = residual_var_samples
  );
}
