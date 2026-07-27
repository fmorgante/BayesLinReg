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

It was created by OpenAI Codex, supervised by Fabio Morgante. It has not been reviewed and tested carefully.

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
  residual_var = 1
)

fit$ETA$ETA1$coefficient_mean
fit$intercept_mean
```

The available models are `"Normal"`, `"SpikeSlab"`, `"SpikeMultiSlab"`, and
`"GlobalLocal"`.
For every model, `residual_var` may be supplied as a fixed value or learned
from `residual_shape` and `residual_scale`.
For mixed priors, use a named `ETA` list whose blocks specify their own
predictors, model, standardization, and prior parameters.

For a high-dimensional global-local block, the half-Cauchy global-scale prior
can be calibrated from an expected number of nonzero coefficients and a
reference residual standard deviation:

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

All coefficient-prior models can instead calibrate their scale from an
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
  residual_shape = 6,
  residual_scale = 5
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

By default, `blm()` returns the retained posterior draws. For large models,
use `store_samples = FALSE` to compute posterior summaries online and keep the
fitted object smaller:

```r
fit_summary <- blm(
  y,
  ETA = list(X = X, model = "Normal"),
  residual_var = 1,
  store_samples = FALSE,
  store_coefficient_cov = FALSE
)
```

Individual draws and convergence diagnostics are unavailable for a
summary-only fit. Every `ETA` block always returns a named `coefficient_var`
vector. Set `store_coefficient_cov = FALSE` to omit its full
`coefficient_cov` matrix; this also avoids the quadratic-size covariance
accumulator when `store_samples = FALSE`.

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
requires`version = "Rcpp"`. Full eigenvalue-based validation is optional through
`check_psd = TRUE` and is disabled by default to avoid its cubic initialization
cost; requesting it for sparse input temporarily constructs a dense matrix.
