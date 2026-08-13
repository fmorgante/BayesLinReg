# BayesLinReg

This R package implements Bayesian linear regression with independent normal,
spike-and-slab, BayesR-style spike-and-multiple-slab, and beta-prime
global-local priors on the regression
coefficients. The global-local family defaults to the Strawderman-Berger prior
and includes the horseshoe as a special case. Multiple regression fits can
combine multiple predictor blocks through a BGLR-style `ETA` interface. Each
block controls its own standardization, prior family, and prior parameters,
while coefficients are always returned on their original scale. Gibbs sampling
is available in R and Rcpp with optional parallel chains.

> [!CAUTION]
> **DISCLAIMER:** This package was created by OpenAI Codex, supervised by Fabio Morgante. It has not been reviewed and tested carefully.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("fmorgante/BayesLinReg")
```

Then load the package with:

```r
library(BayesLinReg)
```

Source installations use Eigen's native dense kernels by default. Users with
an optimized BLAS can opt into Eigen's external BLAS delegation at install
time:

```r
install.packages("remotes")
remotes::install_github("fmorgante/BayesLinReg", configure.args="--enable-eigen-blas")
```

or

```sh
R CMD INSTALL --configure-args="--enable-eigen-blas" BayesLinReg
```

Run `blm_build_info()` to report the compiled backend. See `INSTALL` for the
equivalent `install.packages()` and environment-variable forms.

## Examples

### Raw data

`blm()` always receives predictors and coefficient priors through `ETA`. The
following example fits a ten-predictor model with normal coefficient priors and
known residual variance:

```r
set.seed(123)
n <- 50
p <- 10
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("x", seq_len(p))

beta <- c(2, -1.5, rep(0, p - 2))
y <- drop(1 + X %*% beta + rnorm(n))

fit <- blm(
  y,
  ETA = list(X = X, model = "Normal", var_shape = 2, var_scale = 10),
  residual_var = 1,
  store_samples = TRUE
)

fit$ETA$ETA1$coefficient_mean
fit$intercept_mean
```

Posterior summaries and uncertainty-aware predictions use familiar S3
methods:

```r
summary(fit)
summary(fit, coefficients = "top", max_coefficients = 20)
summary(
  fit,
  coefficients = "top",
  max_coefficients = 20,
  rank_by = "inclusion_probability"
)
predict(fit, X[1:5, , drop = FALSE], interval = "credible")
predict(fit, X[1:5, , drop = FALSE], interval = "prediction", level = 0.9)
```

The default summary is block-first. It reports all coefficients automatically
only for fits with at most 20 predictors; larger coefficient tables are
omitted unless `coefficients = "top"` or `coefficients = "all"` is requested.
Quantiles are computed only for reported coefficients.

Credible intervals describe uncertainty in the conditional mean. Prediction
intervals additionally include residual variation. Both require
`store_samples = TRUE`; point predictions and mean/SD summaries remain
available for memory-efficient online fits. Interval predictions are processed
in chunks of at most 1,000 observations by default; use `chunk_size` in
`predict()` to change that bound.

The available models are `"Normal"`, `"SpikeSlab"`, `"SpikeMultiSlab"`,
`"GlobalLocal"`, and `"Fixed"`. A `"Fixed"` block gives its coefficients a
flat prior while retaining the package's separate intercept. Fixed predictors
are checked jointly for full column rank after centering, including predictors
split across multiple fixed blocks.
For every model, `residual_var` may be supplied as a fixed value or learned
from `residual_shape` and `residual_scale`.
For mixed priors, use a named `ETA` list whose blocks specify their own
predictors, model, standardization, and prior parameters.

The `Normal`, `SpikeSlab`, and `SpikeMultiSlab` models can hold their shared
coefficient or slab variance fixed with `var`. For example:

```r
fixed_slab_fit <- blm(
  y,
  ETA = list(X = X, model = "SpikeSlab", var = 0.5),
  residual_var = 1
)
```

For `SpikeMultiSlab`, component variances are `gamma * var`. A fixed `var` is
mutually exclusive with `var_shape`, `var_scale`, and `expected_pve`.

For `GlobalLocal`, use `global_var` to fix the global variance while retaining
sampled local variances:

```r
fixed_global_fit <- blm(
  y,
  ETA = list(X = X, model = "GlobalLocal", global_var = 0.01),
  residual_var = 1
)
```

`global_var` is mutually exclusive with `global_scale`, `expected_nonzero`,
`reference_residual_var`, and `expected_pve`.

For a high-dimensional global-local block, the half-Cauchy global-scale prior
can be calibrated from an expected number of nonzero coefficients and a
reference residual variance:

```r
fit_global_local <- blm(
  y,
  ETA = list(
    X = X,
    model = "GlobalLocal",
    local_shape = c(a = 0.5, b = 0.5),
    expected_nonzero = 2,
    reference_residual_var = 1
  ),
  residual_shape = 2,
  residual_scale = 1
)
```

For a block with `p` predictors, this sets `global_scale` to
`expected_nonzero / (p - expected_nonzero) *
sqrt(reference_residual_var / n)`. Supply either the two calibration fields or
an explicit `global_scale`, not both.

All penalized coefficient-prior models can instead calibrate their scale from an
expected proportion of response variance explained:

```r
fit_pve <- blm(
  y,
  ETA = list(
    X = X,
    model = "SpikeSlab",
    pi = c(a = 1, b = 9),
    var_shape = 3,
    expected_pve = 0.3
  ),
  residual_shape = 6
)
```

By default, `blm()` uses `var(y)` as the reference response variance; supply
`reference_response_var` to use an external value. For Normal, SpikeSlab, and
SpikeMultiSlab, `expected_pve` determines `var_scale` while `var_shape`
continues to control prior concentration. For GlobalLocal, combine
`expected_pve` with `expected_nonzero`; the expected PVE values across all
blocks determine the reference residual variance used by the sparsity-based
global-scale calibration. This does not directly moment-match the
GlobalLocal signal variance because its beta-prime local prior need not have a
finite variance moment.

When every block supplies `expected_pve`, omitting `residual_scale` calibrates
the residual inverse-gamma prior as
`(residual_shape - 1) * (1 - sum(expected_pve)) *
reference_response_var`. Thus its prior mean is the response variance not
assigned to the predictor blocks. This requires `residual_shape > 1`. Supply
`residual_scale` explicitly to use an independently specified residual prior.
Fits containing a `"Fixed"` block always require an explicit `residual_scale`
when residual variance is learned, because fixed effects do not have an
`expected_pve` prior target.

By default, `blm()` computes posterior summaries online without retaining
individual draws or full coefficient covariance matrices. Request either form
of posterior storage explicitly when needed:

```r
fit_summary <- blm(
  y,
  ETA = list(X = X, model = "Normal"),
  residual_var = 1,
  store_samples = TRUE,
  store_coefficient_cov = TRUE
)
```

Individual draws and convergence diagnostics are unavailable for the default
summary-only fit. Every `ETA` block always returns a named `coefficient_var`
vector. Set `store_coefficient_cov = TRUE` to request its full
`coefficient_cov` matrix. With online summaries, covariance storage and work
are quadratic within each `ETA` block.

### Posterior variance explained

Set `compute_pve = TRUE` to calculate PVE only at retained posterior draws:

```r
fit_with_pve <- blm(
  y,
  ETA = list(
    first = list(X = X[, 1, drop = FALSE], model = "Normal"),
    second = list(X = X[, -1, drop = FALSE], model = "SpikeSlab")
  ),
  residual_var = 1,
  compute_pve = TRUE,
  pve_type = "standalone"
)
```

Each block returns `pve_mean`, `pve_var`, and, when posterior draws are stored,
`pve_samples`. The fit also returns total and cross-block PVE. Standalone block
PVE measures the variance of that block's linear predictor; because correlated
blocks have cross-covariance, standalone values need not sum to total PVE.
`pve_type = "allocated"` distributes those covariance terms across blocks so
their PVE values sum to total PVE, although an allocated block value can be
negative. `cross_block_pve` is reported under either definition.

PVE calculation is disabled by default and therefore has no default runtime
cost. When enabled, it adds matrix-vector work only after burn-in and at draws
retained by `thin`; there is no separate PVE thinning parameter.

### Sufficient statistics

Use `blm_ss()` when the original response and predictor matrix are unavailable:

```r
fit_ss <- blm_ss(
  n = nrow(X),
  XtX = crossprod(X),
  Xty = crossprod(X, y),
  yty = sum(y^2),
  X_means = colMeans(X),
  y_mean = mean(y),
  ETA = list(model = "Normal"),
  residual_var = 1
)
```

`n`, `XtX`, and `Xty` are required. Learning the residual variance additionally
requires `yty`. Supply `X_means` and `y_mean` together to fit an intercept;
otherwise `blm_ss()` fits a no-intercept model and warns. For multiple prior
blocks, each `ETA` block uses `indices` to select a disjoint set of columns from
`XtX`.

`XtX` may also be a compressed sparse `dgCMatrix` or `dsCMatrix`. Sparse input
requires `version = "Rcpp"`. Full eigenvalue-based validation is optional through
`check_psd = TRUE` and is disabled by default to avoid its cubic initialization
cost; requesting it for sparse input temporarily constructs a dense matrix.

For exactly block-diagonal cross-products, `XtX` may instead be a list of
dense or sparse Gram matrices. Predictor order is their concatenated order,
and omitted cross-block entries are assumed to be zero. Gram blocks do not
need to align with `ETA` prior blocks: a prior block may span several Gram
blocks, and shared prior parameters are updated using all of its predictors.

```r
fit_block_ss <- blm_ss(
  n = n,
  XtX = list(region_1 = XtX_region_1, region_2 = XtX_region_2),
  Xty = c(Xty_region_1, Xty_region_2),
  yty = yty,
  ETA = list(model = "SpikeMultiSlab"),
  residual_var = 1
)
```

For symmetric sparse matrices supplied directly or in a list,
`XtX_storage = "speed"` expands both triangles, while `"memory"` chooses the
smaller exact representation.
The default `"auto"` uses `XtX_memory_limit` to retain triangular storage when
expansion would exceed the requested internal-memory budget. The selected
representation and estimated bytes are returned in `fit$XtX_storage`. The
memory representation stores a lower triangle and
streams updates to not-yet-visited coefficients, then reconstructs the complete
right-hand-side state once per Gibbs sweep; it does not allocate a reverse
adjacency index.

When the working predictor means are zero, separate Gram blocks can be updated
within one chain using `nthreads`. This requires the Rcpp sampler and one chain:

```r
fit_threaded_ss <- blm_ss(
  n = n,
  XtX = list(region_1 = XtX_region_1, region_2 = XtX_region_2),
  Xty = c(Xty_region_1, Xty_region_2),
  yty = yty,
  X_means = numeric(length(c(Xty_region_1, Xty_region_2))),
  y_mean = 0,
  ETA = list(model = "SpikeMultiSlab"),
  nthreads = 4,
  nchains = 1
)
```

Threading is across Gram blocks; coefficient updates within a block remain
sequential. Prior blocks may still cross Gram-block boundaries. Threaded runs
are reproducible for fixed inputs, seed, and thread count, but do not reproduce
the serial sampler draw-for-draw.

### GWAS summary statistics

`blm_gwas()` fits directly from an R data frame of additive quantitative-trait
GWAS results and a reusable LD object, without materializing `XtX`. The GWAS
object must contain `CHR`, `ID`, `POS`, `A1`, `A0`, `N`, `BETA`, and `SE`.

```r
ld <- as_blm_ld(
  R = list(chr1 = R_chr1, chr2 = R_chr2),
  variants = list(chr1 = variants_chr1, chr2 = variants_chr2)
)

fit_gwas <- blm_gwas(
  gwas = gwas_results,
  ld = ld,
  ETA = list(model = "SpikeMultiSlab"),
  scale = "standardized",
  residual_var = 1,
  ld_shrink = 0.01
)
```

`R` contains signed correlations, not squared correlations. List elements are
treated as exactly independent. `as_blm_ld()` also detects exact contiguous
sub-blocks within each element and stores one triangle with an implicit unit
diagonal. LD blocks control computation and may be chromosomes or smaller
regions; `ETA` blocks independently control coefficient priors and may cross
LD-block boundaries.

GWAS variants are matched by `ID`, position, and alleles. Reversed alleles are
handled by changing effect orientation, while unresolved strand-ambiguous or
incompatible variants are excluded. Returned coefficients are effects per
input GWAS `A1` allele. `residual_df_gwas` controls the finite-sample marginal-
regression conversion and defaults to `N - 2`; it is distinct from the fitted
model's `residual_var`, `residual_shape`, and `residual_scale` arguments.
Character variant IDs should be used for multi-block `ETA` definitions when
harmonization may exclude variants. Numeric indices are rejected after an LD
variant is excluded because their original positional meaning is ambiguous.

The Gibbs sampler retains the LD object's predictor order even when `ETA`
blocks cross LD blocks or list predictors in a different order. Consequently,
arbitrary prior layouts do not require reconstructing or recompressing LD.
Saved `blm_ld` objects carry an internal format version; recreate an object
with the current `as_blm_ld()` if a later package version reports that its
format is unsupported.

When LD comes from an external reference panel, fixing `residual_var` is
recommended because the reconstructed `XtX`, `Xty`, and `yty` need not be
mutually compatible. For a standardized quantitative trait,
`residual_var = 1` is a conservative choice. `ld_shrink` applies
`(1 - ld_shrink) * R + ld_shrink * I` inside the native LD operator without
copying or modifying the stored LD object. Use `diagnose_blm_ld()` to inspect
block definiteness and `regularize_blm_ld()` to create a reusable regularized
LD object when structural repair is needed. Shrinkage of a block that is too
large for eigendiagnostics regularizes it but does not certify that it is PSD.

With `scale = "original"`, `reference_response_var` is required. Without it,
the default `scale = "auto"` constructs a standardized working likelihood.
GWAS statistics do not identify the phenotype mean, so `blm_gwas()` fits no
intercept and predictions are genetic scores.

For users who need explicit sufficient statistics,

`compute_ss_from_gwas()` converts marginal ordinary least-squares effect
estimates, their standard errors, a common sample size, and signed LD
correlations into centered cross-products for `blm_ss()` using the
finite-sample SuSiE-RSS transformation. `LD` can be one matrix or a list of
matrices when cross-block LD is assumed to be exactly zero:

```r
ss <- compute_ss_from_gwas(
  beta = marginal_beta,
  se = marginal_se,
  LD = LD_correlations,
  n = gwas_sample_size,
  response_var = phenotype_variance
)

fit_gwas <- blm_ss(
  n = ss$n,
  XtX = ss$XtX,
  Xty = ss$Xty,
  yty = ss$yty,
  X_means = ss$X_means,
  y_mean = ss$y_mean,
  ETA = list(model = "SpikeMultiSlab"),
  residual_shape = 2,
  residual_scale = ss$reference_response_var
)
```

For block LD, use `LD = list(region_1 = LD_1, region_2 = LD_2)`. The returned
`ss$XtX` is a matching list and can be passed directly to `blm_ss()`; its prior
blocks may cross LD-block boundaries. Symmetric sparse `dsCMatrix` inputs keep
their one-triangle storage through this conversion. With `output = "eigen"`,
each scaled LD block is instead decomposed independently and the raw
eigenvector, eigenvalue, and tolerance outputs retain the same block names.
After selecting positive eigenpairs and calculating retained fractions, these
lists can be passed directly to `blm_ss_eigen()`.

Omitting `response_var` constructs the statistics on the standardized
predictor and response scale. With an in-sample LD matrix and compatible OLS
statistics, the reconstructed cross-products recover the individual-data
likelihood. With reference LD or other GWAS models, they are approximate
working cross-products. The function performs basic input validation but does
not diagnose general GWAS/LD mismatch.

Set `output = "eigen"` to compute and return the complete eigendecomposition as
`XtX_eigenvectors_raw` and `XtX_eigenvalues_raw`. No eigenpairs are filtered or
modified. Negative eigenvalues produce a warning; the user must select strictly
positive eigenpairs and calculate `XtX_prop_var` before calling
`blm_ss_eigen()`.

### Eigen sufficient statistics

`blm_ss_eigen()` accepts precomputed positive eigenpairs of the centered
predictor cross-product. This avoids storing or traversing a dense `p` by `p`
matrix when the retained rank is much smaller than `p`:

```r
X_centered <- sweep(X, 2, colMeans(X), FUN = "-")
decomposition <- eigen(crossprod(X_centered), symmetric = TRUE)
keep <- decomposition$values >
  sqrt(.Machine$double.eps) * max(decomposition$values)

fit_eigen <- blm_ss_eigen(
  n = nrow(X),
  XtX_eigenvectors = decomposition$vectors[, keep, drop = FALSE],
  XtX_eigenvalues = decomposition$values[keep],
  XtX_prop_var = 1,
  Xty = crossprod(X, y),
  yty = sum(y^2),
  X_means = colMeans(X),
  y_mean = mean(y),
  ETA = list(model = "Normal"),
  residual_var = 1
)
```

`XtX_prop_var = 1` declares that every positive eigenpair is supplied and gives
an exact re-expression of the sufficient-statistic likelihood. Supplying a
truncated eigenspace with `XtX_prop_var < 1` gives an explicitly approximate
posterior. The eigendecomposition is intentionally not computed inside
`blm_ss_eigen()`, so it can be precomputed once and reused across fits.

For an exactly block-diagonal centered cross-product, the eigenvectors and
eigenvalues may instead be matching lists. The sampler retains separate
low-rank transformed designs and residuals, avoiding a global mostly-zero
pseudo-design. `XtX_prop_var` may be a scalar or one value per block, and ETA
prior blocks may cross eigen-block boundaries. With list input, `nthreads > 1`
updates independent eigen blocks concurrently within one chain; coefficient
updates remain sequential within each block. Nonzero `X_means` are supported
because the supplied eigenpairs describe the already-centered cross-product.
