#include <RcppEigen.h>
#include <RcppParallel.h>
#include <R_ext/Rdynload.h>
#include <GIGrvg.h>
#include <cmath>
#include "gibbs_core.h"
#include "gig_sampler.h"

namespace bayeslinreg {

typedef SEXP (*gig_sampler_type)(int, double, double, double);

double draw_gig_gigrvg(
    const double lambda,
    const double chi,
    const double psi) {
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

class RUniformRng {
 public:
  double uniform() {
    return ::unif_rand();
  }
};

double draw_gig(const double lambda, const double chi, const double psi) {
  RUniformRng rng;
  try {
    return sample_gig(lambda, chi, psi, rng);
  } catch (const std::exception& error) {
    Rcpp::stop("%s", error.what());
  }
  return NA_REAL;
}

}  // namespace bayeslinreg

using bayeslinreg::BlockSummaryMatrix;
using bayeslinreg::DenseSummaryMatrix;
using bayeslinreg::EigenBlockSummaryMatrix;
using bayeslinreg::LDSummaryMatrix;
using bayeslinreg::SparseSummaryMatrix;
using bayeslinreg::blm_gibbs_core;
using bayeslinreg::draw_gig;

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
    draws[index] = bayeslinreg::draw_gig_gigrvg(
      lambda, chi, psi
    );
  }
  return draws;
}

// [[Rcpp::export]]
Rcpp::NumericVector draw_gig_native_rcpp_cpp(
    const int n,
    const double lambda,
    const double chi,
    const double psi) {
  Rcpp::RNGScope scope;
  if (n < 1) {
    Rcpp::stop("GIG sample size must be a positive integer.");
  }
  Rcpp::NumericVector draws(n);
  bayeslinreg::RUniformRng rng;
  try {
    for (int index = 0; index < n; ++index) {
      draws[index] = bayeslinreg::sample_gig(lambda, chi, psi, rng);
    }
  } catch (const std::exception& error) {
    Rcpp::stop("%s", error.what());
  }
  return draws;
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
    fixed_global_var,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, fixed_var, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, likelihood_df, fit_intercept,
    intercept_x_mean,
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
    fixed_global_var,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, fixed_var, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, likelihood_df, fit_intercept,
    intercept_x_mean,
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
    fixed_global_var,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, fixed_var, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, likelihood_df, fit_intercept,
    intercept_x_mean,
    intercept_y_mean, true, summary_matrix, summary_Xty, summary_yty,
    true, 0.0, compute_pve, pve_type_code, nthreads
  );
}

// [[Rcpp::export]]
Rcpp::List blm_gibbs_ld_rcpp_cpp(
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
    const Rcpp::List& ld_blocks,
    const Rcpp::List& ld_indices,
    const Rcpp::NumericVector& ld_scale,
    const double ld_shrink,
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
  const LDSummaryMatrix summary_matrix(
    ld_blocks, ld_indices, ld_scale, ld_shrink, nthreads
  );
  return blm_gibbs_core(
    empty_y, empty_X, residual_shape, residual_scale, iterations, burnin, thin,
    progress_callback, block_id, block_model, normal_shape, normal_scale,
    pi_alpha, pi_beta, spike_var_shape, spike_var_scale, global_scale,
    fixed_global_var,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, fixed_var, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, likelihood_df, fit_intercept,
    intercept_x_mean,
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
    fixed_global_var,
    local_a, local_b, multi_gamma_list, multi_pi_alpha_list, multi_var_shape,
    multi_var_scale, fixed_var, learn_residual_var, fixed_residual_var, store_samples,
    store_coefficient_cov, effective_n, likelihood_df, fit_intercept,
    intercept_x_mean,
    intercept_y_mean, true, summary_matrix, summary_Xty, summary_yty,
    true, 0.0, compute_pve, pve_type_code, nthreads
  );
}
