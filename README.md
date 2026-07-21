# BayesLinReg

This R package implements Bayesian linear regression with independent normal,
spike-and-slab, and beta-prime global-local priors on the regression
coefficients. The global-local family defaults to the Strawderman-Berger prior
and includes the horseshoe as a special case. Multiple regression fits can
combine multiple predictor blocks through a BGLR-style `ETA` interface. Each
block controls its own standardization, prior family, and prior parameters,
while coefficients are always returned on their original scale. Gibbs sampling
is available in R and Rcpp with optional parallel chains.

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

## Example

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

The available models are `"Normal"`, `"SpikeSlab"`, and `"GlobalLocal"`.
For mixed priors, use a named `ETA` list whose blocks specify their own
predictors, model, standardization, and prior parameters.

It was created by OpenAI Codex, supervised by Fabio Morgante. It has not been reviewed and tested carefully.
