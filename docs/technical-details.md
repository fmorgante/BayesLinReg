# BayesLinReg statistical models and implementation

This document describes the statistical parameterizations and computational
algorithms implemented in BayesLinReg. It complements the R function
documentation. The equations below describe the current sampler rather than
generic versions of similarly named models.

## 1. Statistical models

### 1.1 Likelihood, blocks, and working scale

Let the response be $y\in\mathbb{R}^n$, and partition the predictors into
blocks $X_1,\ldots,X_B$. The model is

$$
y = \mu\mathbf{1}+\sum_{b=1}^B X_b\beta_b+e,
\qquad e\sim N(0,\sigma_e^2I_n).
$$

The intercept is integrated out of the coefficient updates by centering the
response and every predictor. Within each block, `standardize = TRUE` (the
default) also divides predictor $j$ by its sample standard deviation
$s_j$. If $\theta_j$ denotes the coefficient sampled on this working
scale, the reported coefficient and intercept are

$$
\beta_j=\frac{\theta_j}{s_j},\qquad
\mu=\bar y-\sum_j\bar x_j\beta_j.
$$

More precisely, with the flat intercept prior implicit in centering,

$$
\mu\mid\beta,\sigma_e^2,y,X\sim
N(\bar y-\sum_j\bar x_j\beta_j,\frac{\sigma_e^2}{n}).
$$

At each retained coefficient and residual-variance draw, the package also draws
the intercept from this conditional normal distribution. Reported intercept
uncertainty and prediction intervals therefore include the conditional
$\sigma_e^2/n$ contribution. A sufficient-statistic fit omits the intercept
entirely when `X_means` and `y_mean` are not supplied.

When `standardize = FALSE`, $s_j=1$. Priors described below apply to the
working-scale coefficients $\theta_j$. Different blocks may use different
prior families, and every block has its own hyperparameters.

Define the centered working design by $Z$, its columns by $z_j$, and

$$
d_j=z_j'z_j,\qquad
r_{-j}=y_c-\sum_{k\ne j}z_k\theta_k,\qquad
t_j=z_j'r_{-j}.
$$

For any normal prior $\theta_j\sim N(0,v_j)$, the coefficient update is

$$
V_j=(\frac{d_j}{\sigma_e^2}+\frac{1}{v_j})^{-1},
\qquad
m_j=V_j\frac{t_j}{\sigma_e^2},
\qquad
\theta_j\mid-\sim N(m_j,V_j).
$$

For a flat prior, set $1/v_j=0$.

Throughout this document, $x\sim\mathrm{IG}(a,b)$ means the
inverse-gamma distribution with density proportional to
$x^{-(a+1)}\exp(-b/x)$. Thus, the code samples it as $x=1/G$, where
$G\sim\mathrm{Gamma}(a,\text{rate}=b)$. Coefficient-prior variances
are not multiplied by $\sigma_e^2$.

### 1.2 Fixed

For a `Fixed` block,

$$
p(\theta_b)\propto 1.
$$

Each coefficient uses the common normal update with prior precision zero:

$$
V_j=\frac{\sigma_e^2}{d_j},\qquad m_j=\frac{t_j}{d_j}.
$$

All fixed predictors across all blocks must be jointly full column rank on the
centered working scale. The intercept remains separate, so a constant predictor
must not be supplied as a fixed effect.

### 1.3 Normal

For a `Normal` block with $p_b$ predictors,

$$
\theta_j\mid\sigma_{\beta,b}^2
  \stackrel{\mathrm{ind}}{\sim}N(0,\sigma_{\beta,b}^2),
\qquad
\sigma_{\beta,b}^2\sim\mathrm{IG}(a_b,b_b).
$$

The user-facing parameters are `var_shape` $=a_b$ and `var_scale`
$=b_b$. Their defaults are 2 and 1. Coefficients use the common normal
update with $v_j=\sigma_{\beta,b}^2$. The shared variance update is

$$
\sigma_{\beta,b}^2\mid-\sim\mathrm{IG}(
  a_b+\frac{p_b}{2},
  b_b+\frac{1}{2}\sum_{j\in b}\theta_j^2
).
$$

Alternatively, positive `var` fixes $\sigma_{\beta,b}^2$ and skips this
inverse-gamma update.

### 1.4 SpikeSlab

For a `SpikeSlab` block,

$$
\delta_j\mid\pi_b\sim\mathrm{Bernoulli}(\pi_b),
\qquad
\theta_j\mid\delta_j,\sigma_{\beta,b}^2\sim
\begin{cases}
\delta_0,&\delta_j=0,\\
N(0,\sigma_{\beta,b}^2),&\delta_j=1,
\end{cases}
$$

$$
\pi_b\sim\mathrm{Beta}(a_{\pi,b},b_{\pi,b}),
\qquad
\sigma_{\beta,b}^2\sim\mathrm{IG}(a_b,b_b).
$$

The defaults are `pi = c(a = 1, b = 1)`, `var_shape = 2`, and
`var_scale = 1`. Alternatively, positive `var` fixes
$\sigma_{\beta,b}^2$ and skips its inverse-gamma update. Let $V_j,m_j$ be the slab conditional parameters from the
common normal update with $v_j=\sigma_{\beta,b}^2$. The posterior log odds
used by the implementation are

$$
\ell_j=log\frac{\pi_b}{1-\pi_b}
+\frac12\log\frac{V_j}{\sigma_{\beta,b}^2}
+\frac{m_j^2}{2V_j}.
$$

Then

$$
\delta_j\mid-\sim\mathrm{Bernoulli}\{\mathrm{logit}^{-1}(\ell_j)\},
$$

and $\theta_j=0$ if $\delta_j=0$, otherwise
$\theta_j\sim N(m_j,V_j)$. If $m_b=\sum_{j\in b}\delta_j$,

$$
\pi_b\mid-\sim
\mathrm{Beta}(a_{\pi,b}+m_b,b_{\pi,b}+p_b-m_b),
$$

$$
\sigma_{\beta,b}^2\mid-\sim\mathrm{IG}(
  a_b+\frac{m_b}{2},
  b_b+\frac12\sum_{j:\delta_j=1}\theta_j^2
).
$$

### 1.5 SpikeMultiSlab

`SpikeMultiSlab` is a BayesR-style point-normal mixture with a learned shared
variance. Let $c_j\in\{0,\ldots,K-1\}$ be the component assignment,
$\gamma_0=0$, and $0<\gamma_1<\cdots<\gamma_{K-1}$. The hierarchy is

$$
c_j\mid\boldsymbol\pi_b\sim\mathrm{Categorical}(\boldsymbol\pi_b),
$$

$$
\theta_j\mid c_j,\sigma_{\beta,b}^2\sim
\begin{cases}
\delta_0,&c_j=0,\\
N(0,\gamma_{c_j}\sigma_{\beta,b}^2),&c_j>0,
\end{cases}
$$

$$
\boldsymbol\pi_b\sim\mathrm{Dirichlet}(\boldsymbol\alpha_b),
\qquad
\sigma_{\beta,b}^2\sim\mathrm{IG}(a_b,b_b).
$$

The defaults are

```r
gamma = c(0, 0.01, 0.1, 1)
alpha = rep(1, length(gamma))
var_shape = 2
var_scale = 1
```

Alternatively, positive `var` fixes $\sigma_{\beta,b}^2$; the component
variances are then `gamma * var`, and the inverse-gamma update is skipped.
For all three variance-mixture models, `var` is mutually exclusive with
`var_shape`, `var_scale`, and `expected_pve`.

For component $c>0$, define

$$
V_{jc}=\{\frac{d_j}{\sigma_e^2}+
  \frac{1}{\gamma_c\sigma_{\beta,b}^2}\}^{-1},
\qquad m_{jc}=V_{jc}\frac{t_j}{\sigma_e^2}.
$$

The unnormalized log component weights are

$$
L_{j0}=\log\pi_{b0},
$$

$$
L_{jc}=\log\pi_{bc}
+\frac12\log\frac{V_{jc}}{\gamma_c\sigma_{\beta,b}^2}
+\frac{m_{jc}^2}{2V_{jc}},\qquad c>0.
$$

The code normalizes these weights with a log-sum-exp calculation, samples
$c_j$, and sets $\theta_j=0$ for the spike or draws
$N(m_{jc},V_{jc})$ for a slab. If $n_{bc}$ is the number of indices $j$ in
block $b$ for which $c_j=c$,

$$
\boldsymbol\pi_b\mid-\sim
\mathrm{Dirichlet}(\alpha_{b0}+n_{b0},\ldots,
                         \alpha_{b,K-1}+n_{b,K-1}),
$$

$$
\sigma_{\beta,b}^2\mid-\sim\mathrm{IG}(
  a_b+\frac12\sum_{c>0}n_{bc},
  b_b+\frac12\sum_{j:c_j>0}\frac{\theta_j^2}{\gamma_{c_j}}
).
$$

### 1.6 GlobalLocal

For a `GlobalLocal` block,

$$
\theta_j\mid\tau_b^2,\psi_j\sim N(0,\tau_b^2\psi_j),
\qquad \psi_j\sim\mathrm{BetaPrime}(a_b,b_b),
\qquad \tau_b\sim C^+(0,s_b).
$$

Here `local_shape = c(a, b)` and `global_scale` $=s_b$. The default
`local_shape = c(1, 0.5)` gives the Strawderman--Berger local prior;
`c(0.5, 0.5)` gives the horseshoe local prior.

Alternatively, positive `global_var` fixes $\tau_b^2$ while the local
variances $\psi_j$ remain sampled. In that case the $\tau_b^2$ and $\zeta_b$
updates below are skipped. `global_var` is mutually exclusive with
`global_scale`, `expected_nonzero`, `reference_residual_var`, and
`expected_pve`.

The beta-prime distribution is represented through

$$
\psi_j\mid\xi_j\sim\mathrm{Gamma}(a_b,\text{rate}=\xi_j),
\qquad \xi_j\sim\mathrm{Gamma}(b_b,\text{rate}=1).
$$

The half-Cauchy global scale is represented through

$$
\tau_b^2\mid\zeta_b\sim\mathrm{IG}(1/2,1/\zeta_b),
\qquad \zeta_b\sim\mathrm{IG}(1/2,1/s_b^2).
$$

Coefficient updates use $v_j=\tau_b^2\psi_j$. The remaining full
conditionals are

$$
\psi_j\mid-\sim\mathrm{GIG}(
  a_b-\frac12,\frac{\theta_j^2}{\tau_b^2},2\xi_j
),
$$

where the GIG density is proportional to
$x^{\lambda-1}\exp\{-(\chi/x+\omega x)/2\}$, and

$$
\xi_j\mid-\sim
\mathrm{Gamma}(a_b+b_b,\text{rate}=1+\psi_j),
$$

$$
\tau_b^2\mid-\sim\mathrm{IG}(
  \frac{p_b+1}{2},
  \frac{1}{\zeta_b}+\frac12\sum_{j\in b}\frac{\theta_j^2}{\psi_j}
),
$$

$$
\zeta_b\mid-\sim\mathrm{IG}(
  1,\frac{1}{s_b^2}+\frac{1}{\tau_b^2}
).
$$

The Rcpp sampler uses an independent implementation of the Hörmann--Leydold
(2014) GIG rejection algorithms. Its small-parameter envelope is evaluated on
the log scale, and positive gamma auxiliaries are bounded at the smallest
positive normal double if a valid draw underflows numerically. For
block-parallel sufficient-statistic fits,
conditionally independent local-variance updates run concurrently across
computational blocks using per-block C++ random-number streams. Local auxiliary
and shared global-variance updates remain serial after the workers join. The
`GIGrvg` implementation remains an internal statistical and performance
reference.

### 1.7 Residual variance

The residual variance can be fixed by supplying `residual_var`. Otherwise,

$$
\sigma_e^2\sim\mathrm{IG}(a_e,b_e),
$$

where `residual_shape` $=a_e$ and `residual_scale` $=b_e$. Let $d_L$ denote
the likelihood dimension recorded as `likelihood_df`. For ordinary raw-data
fits and sufficient statistics with an identified intercept, $d_L=n-1$; for
uncentered sufficient statistics without an intercept, $d_L=n$. Its update is

$$
\sigma_e^2\mid-\sim\mathrm{IG}(
  a_e+\frac{d_L}{2},
  b_e+\frac12\|y_c-Z\theta\|^2
).
$$

For a truncated eigen representation, the residual sum of squares also contains
the stored orthogonal offset described in Section 2.4.

### 1.8 Prior calibration with `expected_pve`

Let $V_y$ be `reference_response_var`, let $h_b$ be a block's
`expected_pve`, let $G_b=Z_b'Z_b$ be its working-scale Gram block, and
define

$$
V_{g,b}=h_bV_y,\qquad
D_b=\frac{\mathrm{tr}(G_b)}{d_L}
   =\sum_{j\in b}\frac{Z_j'Z_j}{d_L}.
$$

Thus the calibration uses the same `likelihood_df` denominator as posterior
PVE. Predictor standardization itself continues to use sample standard
deviations: each standardized column has $Z_j'Z_j=n-1$, so
$D_b=p_b(n-1)/d_L$, reducing to $p_b$ when $d_L=n-1$. The calibration matches
the prior mean block signal variance $E(f_b'f_b/d_L)$ under independent,
zero-mean coefficient priors; it does not claim that the prior mean of the
random PVE ratio equals $h_b$.

For `Normal`,

$$
b_b=(a_b-1)\frac{V_{g,b}}{D_b}.
$$

For `SpikeSlab`, with
$q_b=E(\pi_b)=a_{\pi,b}/(a_{\pi,b}+b_{\pi,b})$,

$$
b_b=(a_b-1)\frac{V_{g,b}}{q_bD_b}.
$$

For `SpikeMultiSlab`, define

$$
\bar\gamma_b=\sum_c
\frac{\alpha_{bc}}{\sum_k\alpha_{bk}}\gamma_c.
$$

Then

$$
b_b=(a_b-1)\frac{V_{g,b}}{D_b\bar\gamma_b}.
$$

These three calibrations require `var_shape > 1`. For `GlobalLocal`, the
beta-prime local prior need not possess a finite signal-variance moment, so
`expected_pve` does not directly match $V_{g,b}$. Instead, if all penalized
blocks supply `expected_pve`, define $H=\sum_bh_b$ and

$$
\sigma_0^2=(1-H)V_y,
\qquad
s_b=\frac{p_{0b}}{p_b-p_{0b}}\sqrt{\frac{\sigma_0^2}{n}},
$$

where $p_{0b}$ is `expected_nonzero`. This is the expected-sparsity scale
heuristic used by the package.

If all blocks supply `expected_pve`, `residual_scale` may be omitted. The
package then sets

$$
b_e=(a_e-1)(1-H)V_y,
$$

so the inverse-gamma prior mean is $(1-H)V_y$. This requires
`residual_shape > 1`.

### 1.9 Posterior PVE

At a retained draw, let

$$
f_b=Z_b\theta_b,\qquad f=\sum_bf_b,
\qquad d=d_L=\texttt{likelihood\_df}.
$$

By default, $d=n-1$ with an intercept and $d=n$ without one; summary-data
interfaces can instead supply an explicit likelihood dimension.

The total signal variance and total PVE are

$$
V_f=\frac{f'f}{d},
\qquad
\mathrm{PVE}_{\mathrm{total}}=
\frac{V_f}{V_f+\sigma_e^2}.
$$

Standalone block PVE is

$$
\mathrm{PVE}^{\mathrm{standalone}}_b=
\frac{f_b'f_b/d}{V_f+\sigma_e^2}.
$$

Because correlated blocks can have nonzero cross-products, these values need
not sum to total PVE. The reported cross-block PVE is

$$
\frac{f'f-\sum_bf_b'f_b}{d(V_f+\sigma_e^2)}.
$$

Allocated block PVE instead uses

$$
\mathrm{PVE}^{\mathrm{allocated}}_b=
\frac{f_b'f/d}{V_f+\sigma_e^2}.
$$

Allocated values sum to total PVE but can be negative. PVE is evaluated only
for retained draws.

### 1.10 Prediction

For new predictors $x_{\ast}$, each retained draw gives the conditional mean

$$
\eta_{\ast}^{(s)}=\mu^{(s)}+x_{\ast}'\beta^{(s)}.
$$

Point prediction uses the posterior mean of this quantity. A credible interval
uses empirical quantiles of $\{\eta_{\ast}^{(s)}\}$. A prediction interval uses
quantiles of the equal-weight normal mixture

$$
\frac1S\sum_{s=1}^S N(\eta_{\ast}^{(s)},\sigma_e^{2(s)}).
$$

The mixture quantiles are found deterministically by root-finding; no additional
posterior-predictive simulation is performed. Each $\mu^{(s)}$ is a retained
conditional intercept draw, so the credible distribution includes intercept
uncertainty. The predictive mixture additionally includes the new observation's
residual variation through $\sigma_e^{2(s)}$.

## 2. Implementation details

### 2.1 Gibbs schedule and random numbers

One iteration performs the following operations:

1. Update coefficients one at a time within the current working design.
2. Update block-specific prior hyperparameters.
3. Update the residual variance when it is learned.
4. Store the draw or update online posterior accumulators when the iteration is
   retained.

Burn-in draws and iterations skipped by `thin` are not stored and do not incur
PVE calculations. Random draws use R's random-number state, so `set.seed()`
controls reproducibility. Independent chains use parallel-safe streams managed
through `future`. Within-chain threaded block sampling uses deterministic
block-specific streams for a fixed seed and thread count.
Serial and within-chain block-parallel sweeps call the same coefficient-
conditional implementation through RNG-specific adapters. This keeps the
Normal, SpikeSlab, GlobalLocal, SpikeMultiSlab, and Fixed update formulas in
one place while preserving R's serial random-number sequence.

### 2.2 Individual-level fitting with `blm()`

`blm()` centers and optionally standardizes each predictor block, concatenates
the working matrices, and maintains the residual vector

$$
r=y_c-Z\theta.
$$

After a coefficient change $\Delta_j$, it performs

$$
r\leftarrow r-z_j\Delta_j.
$$

Thus, one complete dense coefficient sweep costs $O(np)$. The residual is
reconstructed periodically from $y_c-Z\theta$ to limit accumulated
floating-point drift. Dense matrix-vector reconstruction uses Eigen and can use
external BLAS when the package is configured accordingly.

### 2.3 Direct sufficient-statistic fitting with `blm_ss()`

Let the supplied uncentered statistics be

$$
G=X'X,\qquad g=X'y,\qquad y_2=y'y.
$$

When means are supplied, the centered statistics are

$$
G_c=G-n\bar x\bar x',\qquad
g_c=g-n\bar x\bar y,\qquad
y_{2c}=y_2-n\bar y^2.
$$

For sparse input, the dense rank-one centering correction is tracked separately
and is not inserted into the sparse matrix. The sampler maintains

$$
q=g-G\theta.
$$

The coefficient numerator is recovered as

$$
t_j=(g_c-G_c\theta)_j+(G_c)_{jj}\theta_j
$$

with the centering correction applied separately. After a coefficient change,

$$
q\leftarrow q-G_{\cdot j}\Delta_j.
$$

The residual sum of squares is maintained incrementally as

$$
\mathrm{SSE}=y_{2c}-2\theta'g_c+\theta'G_c\theta
$$

and periodically reconstructed to control numerical drift. Direct sufficient-
statistic sampling does not construct individual-level pseudo-observations.

Supported Gram representations are:

- A dense base-R matrix. Dense reconstruction uses an Eigen self-adjoint
  matrix-vector product.
- A general sparse `dgCMatrix`, traversed through `RcppEigen`.
- A list of exactly independent dense or sparse Gram blocks. Omitted
  cross-block entries are assumed to be zero.
- A symmetric sparse `dsCMatrix`, supplied directly or within list input. The
  speed representation expands both triangles. The memory representation uses
  the block kernel with one or more blocks, keeps a lower triangle, updates only
  not-yet-visited coordinates during the ascending Gibbs scan, and reconstructs
  the complete $q$ vector once per sweep. It therefore avoids a reverse
  adjacency index while preserving an exact serial Gibbs transition.

With zero working predictor means, exactly independent Gram blocks can be
sampled concurrently through `RcppParallel`. Updates within a Gram block remain
sequential. Prior blocks need not coincide with Gram blocks; shared prior
hyperparameters are updated after the coefficient sweep.

The optional `check_psd = TRUE` validation materializes the centered Gram
matrix and performs a full eigendecomposition, costing $O(p^3)$ time and
$O(p^2)$ memory. Its default is `FALSE`.

### 2.4 Low-rank sufficient-statistic fitting with `blm_ss_eigen()`

Suppose the centered Gram matrix is represented by retained eigenpairs

$$
G_c\approx U\Lambda U'.
$$

The implementation constructs the reduced algebraic design and response

$$
Q=\Lambda^{1/2}U',\qquad
w=\Lambda^{-1/2}U'g_c.
$$

Initialization computes $U'g_c$, $UU'g_c$, and the diagonal of
$U\Lambda U'$ in C++. It then writes the scaled factor $Q$ directly in final
sampler order. This avoids materializing squared eigenvectors, an explicit
transpose, and additional scaling or reordering copies of the $p$-by-$q$
representation.

Then

$$
Q'Q=U\Lambda U',\qquad Q'w=UU'g_c,
$$

and coefficient sampling uses the reduced residual $w-Q\theta$. A coordinate
update costs $O(q)$, where $q$ is retained rank, and a sweep costs
$O(pq)$. The coordinate kernel maps each contiguous column of $Q$ and the
block residual into Eigen vectors, using vectorized dot-product and scaled
residual-update operations. List input applies this representation
independently to exact eigen blocks and can process independent blocks
concurrently. Periodic state reconstruction reuses one transformed fitted-
value workspace per eigen block, so it does not allocate a new length-$q_b$
vector on every reconstruction.

When `XtX_prop_var = 1`, the supplied eigenvectors must span `Xty`, and the
representation is treated as exact. Values below one explicitly define an
approximate posterior based on the retained eigenspace. When `yty` is supplied,
the implementation records

$$
c=y_{2c}-w'w\ge0
$$

and uses

$$
\mathrm{SSE}=c+\|w-Q\theta\|^2.
$$

Block-level eigen PVE reuses allocated work vectors and clears only blocks
touched by the requested prior block.

### 2.5 LD-native GWAS fitting with `blm_gwas()`

`as_blm_ld(R, variants)` converts one signed correlation matrix, or a named
list of exactly independent matrices, into a reusable LD object. Every variant
table contains `CHR`, `ID`, `POS`, `A1`, and `A0`, where `A1` is the dosage
allele used to calculate the correlations. The constructor stores one strict
lower triangle and treats the unit diagonal as implicit. Contiguous rows use
an index-free `data` and `indptr` representation; irregular sparse patterns
retain row indices.

Within every supplied matrix, exact contiguous block-diagonal structure is
detected without thresholding. The supplied list element remains a parent
reporting group, while maximal exact sub-blocks become independent
computational blocks. Correlations omitted between list elements or detected
sub-blocks are treated as exactly zero.

`blm_gwas()` accepts an in-memory table with columns `CHR`, `ID`, `POS`, `A1`,
`A0`, `N`, `BETA`, and `SE`. It matches variants to the LD object, validates
position and alleles, changes the sign of marginal effects when allele dosage
orientation is reversed, and restores coefficients to the input GWAS `A1`
orientation on output. Unresolved ambiguous or incompatible variants are not
included. Retained variants currently require a common `N`.
Character-indexed `ETA` blocks are filtered by variant ID. Numeric indices are
rejected when an LD variant is excluded because positions before and after
harmonization do not have an unambiguous common meaning.

Writing $\nu_{\mathrm{GWAS}}$ for `residual_df_gwas`, define

$$
z_j=\frac{\widehat\beta_j}{s_j},\qquad
a_j=\frac{n-1}{z_j^2+\nu_{\mathrm{GWAS}}}.
$$

The default $\nu_{\mathrm{GWAS}}=n-2$ corresponds to an intercept and one
tested predictor. This input is separate from `residual_var`,
`residual_shape`, and `residual_scale`, which describe residual variation in
the fitted joint Bayesian regression. It is also separate from the joint
model's `likelihood_df`. Because the reconstructed GWAS statistics are
centered, `blm_gwas()` uses `likelihood_df = n - 1` for the residual-variance
update and PVE normalization even though the phenotype mean is unavailable and
no intercept is returned.

On the standardized working scale,

$$
G=(n-1)R,\qquad
g_j=\sqrt{n-1}\sqrt{a_j}z_j,\qquad
y_2=n-1.
$$

For `ld_shrink = lambda`, the LD-native operator substitutes

$$
R_\lambda=(1-\lambda)R+\lambda I.
$$

Only a scalar multiplier is passed to the native operator: stored
off-diagonal values are not copied or modified, and the implicit diagonal
remains one. This regularization can reduce sensitivity to external-reference
LD but does not make reconstructed summary statistics exact.

On the original scale, `reference_response_var` supplies $V_y$ and

$$
d_j=\frac{V_ya_j}{s_j^2},\qquad
G=D^{1/2}RD^{1/2},\qquad
g_j=d_j\widehat\beta_j,\qquad
y_2=(n-1)V_y.
$$

The LD kernel applies $D^{1/2}RD^{1/2}\theta$ directly and never constructs
the full Gram matrix. During an ascending Gibbs sweep, strict-lower entries
propagate coefficient changes only to coordinates that have not yet been
visited. The complete right-hand-side state is reconstructed once after each
sweep. Independent exact LD sub-blocks may be processed concurrently, while
updates remain sequential inside a connected block. For $m_b$ stored values
in block $b$, a sweep costs $O(p+\sum_b m_b)$.

For posterior PVE, one fused pass over each compressed LD block accumulates the
total quadratic form, every `ETA` block's standalone quadratic form, and its
allocated contribution. Independent LD blocks are processed concurrently and
then reduced in fixed LD-block order, so results do not depend on thread
scheduling. Reusable workspace costs $O(BK)$ doubles for $B$ LD blocks and $K$
`ETA` blocks; with the usual one or two `ETA` blocks this is small. The PVE
pass remains $O(p+\sum_b m_b)$ per retained draw and is never run during
burn-in or skipped iterations.

Sampler coordinates always retain LD order. Prior-block membership, predictor
scaling, and posterior extraction are stored as separate mappings, so `ETA`
blocks can cross LD blocks or request a different output order without
rebuilding the compressed correlations. Partial harmonization filters the
compressed triangle directly and recomputes exact contiguous sub-blocks from
the retained triplets, avoiding an intermediate `sparseMatrix`. Reusable LD
objects carry a validated internal format version. Indexed triangles are
checked for strict-lower, sorted, unique row indices within every column, and
regularization reports are checked against current block names, parents, and
sizes, so incompatible serialized objects fail before entering compiled code.
Reports created before the provenance and final-floor fields were introduced
remain valid when their shared structural and numerical fields are consistent.

The returned `ld_harmonization` vector separates GWAS-only and LD-only entries
from location-mismatched, allele-mismatched, and unresolved ambiguous matched
variants. These categories are mutually exclusive. Its aggregate `excluded`
count is measured in input-table entries: an unmatched entry contributes one,
whereas a matched variant pair rejected from both inputs contributes two.

GWAS summary statistics do not identify a phenotype mean, so this interface
fits no intercept. Reference-panel LD, meta-analysis statistics, mixed-model
statistics, or covariate-adjusted marginal estimates generally define an
approximate working likelihood rather than exact sufficient statistics.
Because learning the residual variance requires a compatible `G`, `g`, and
`y2`, a fixed `residual_var` is recommended with external LD. On the
standardized scale, `residual_var = 1` is a conservative robust choice.

`diagnose_blm_ld()` performs blockwise eigenvalue diagnostics subject to
dimension and estimated-memory limits. `regularize_blm_ld()` returns a new LD
object and an audit report. It can eigen-repair manageable blocks, shrink
off-diagonal correlations without densification, or choose between these
policies automatically. Large blocks that receive shrinkage without an
eigendecomposition are marked as uncertified in the report. Structural PSD
repair and per-fit `ld_shrink` have
different roles: a positive-definite matrix may still be mismatched to the
GWAS population.

After eigen repair restores a unit diagonal, the implementation applies the
smallest additional identity shrinkage needed to enforce the requested final
eigenvalue floor. A fit preserves the input object's original regularization
report unchanged. Its separate `ld_regularization_block_map` links each fitted
post-harmonization block to the source report row, gives fitted and source
predictor counts, and indicates whether subsetting occurred. The fit's
`ld_block_table` describes the resulting computational blocks.

A harmonized block obtained only by deleting matching rows and columns is a
principal submatrix of its repaired source block. Cauchy interlacing therefore
guarantees that its minimum eigenvalue is at least the source block's reported
post-repair minimum eigenvalue. The source value remains a valid lower bound,
but it is not presented as the fitted block's exact eigenvalue.

### 2.6 GWAS reconstruction with `compute_ss_from_gwas()`

For marginal estimates $\widehat\beta_j$, standard errors $s_j$, common
sample size $n$, and signed LD correlation matrix $R$, define

$$
z_j=\frac{\widehat\beta_j}{s_j},\qquad
a_j=\frac{n-1}{z_j^2+\mathrm{residual\_df}}.
$$

On the standardized response and predictor scale, the function returns

$$
G=(n-1)R,\qquad
g_j=\sqrt{n-1}\sqrt{a_j}z_j,\qquad
y_2=n-1.
$$

With response variance $V_y$, the original-scale reconstruction uses

$$
G_{jj}=\frac{V_ya_j}{s_j^2},\qquad
G=D^{1/2}RD^{1/2},\qquad
g_j=G_{jj}\widehat\beta_j,\qquad
y_2=(n-1)V_y,
$$

where $D=\mathrm{diag}(G_{11},\ldots,G_{pp})$. These are working
sufficient statistics when LD and the marginal statistics are not derived from
the same individual-level sample. The function performs structural validation
but intentionally does not perform general summary-statistic/LD mismatch
diagnostics. Diagonal deviations within the unit-diagonal validation tolerance
are normalized to exactly one before reconstruction.

### 2.7 Numerical safeguards

The package:

- rejects constant predictors in individual, direct sufficient-statistic,
  eigen, and GWAS interfaces;
- jointly rank-checks all `Fixed` predictors;
- validates symmetry, finite values, dimensions, names, and block partitions;
- periodically reconstructs residual or cross-product state;
- uses log-scale component probabilities and stable logistic calculations;
- bounds tiny GIG parameters and gamma auxiliaries away from zero;
- rejects execution-control values that exceed the C++ integer range;
- clamps residual SSE and PVE quadratic forms only when their negative values
  are within scale-aware floating-point tolerances, and stops when a
  reconstructed residual SSE is materially negative; and
- optionally checks positive semidefiniteness and joint compatibility of
  sufficient statistics.

These safeguards do not make an externally estimated LD matrix compatible with
GWAS statistics; that remains an input-modeling responsibility.

### 2.8 Posterior storage and summaries

The public fitting functions default to `store_samples = FALSE` and
`store_coefficient_cov = FALSE`, so their default output contains online
posterior means and marginal variances without individual draws or full
coefficient covariance matrices.

With `store_samples = TRUE`, retained coefficient, variance, inclusion,
component, PVE, intercept, and residual-variance draws are stored as applicable.
With `store_samples = FALSE`, means and centered second moments are accumulated
online with Welford updates. Multiple chains are combined with Chan's parallel
moment formulas, including the between-chain contribution to coefficient
covariances. This avoids the cancellation inherent in subtracting squared raw
sums when a posterior mean is large.
Convergence conversion excludes fixed residual and coefficient-prior
variances, which are stored as constant values for posterior reporting but are
not sampled parameters. If a fit has no remaining sampled scalar or PVE
quantity, convergence assessment reports that no diagnostic target is
available.

`store_coefficient_cov = TRUE` additionally accumulates a centered coefficient
second-moment matrix for every retained draw within each `ETA` block, requiring
$O(\sum_b p_b^2)$ storage and work. Turning it off retains only marginal
coefficient variances.
The scalable `summary()` method reports block-level summaries by default and
requires explicit selection before producing large coefficient tables.

All returned coefficients and their covariance summaries are transformed back
to the original predictor scales. Prediction therefore accepts predictors on
the same scales and in the same block structure as those supplied at fitting.
Interval predictions process new observations in reusable chunks, bounding the
temporary fitted-draw matrix by `chunk_size` times the number of retained draws.

## 3. Source map

The main implementations are located in:

- `R/blm.R`: individual-level interface and user documentation.
- `R/blm_ss.R`: direct sufficient-statistic validation and storage planning.
- `R/blm_ss_eigen.R`: low-rank transformation and validation.
- `R/blm_gwas.R`: GWAS validation, harmonization, scaling, and fitting.
- `R/ld_matrix.R`: LD construction, exact sub-block detection, and compressed
  storage.
- `R/compute_ss_from_gwas.R`: GWAS/LD reconstruction.
- `R/fit_preparation.R`: shared block layouts, prior arguments, and sampler
  execution used by all fitting interfaces.
- `R/prior_specification.R`: prior validation, normalization, calibration, and
  predictor preparation.
- `R/mcmc_validation.R`: shared bounded-integer, MCMC, chain, thread, and PVE
  control validation.
- `R/gibbs_r.R`: reference R Gibbs sampler.
- `R/sampler_interface.R`: Rcpp dispatch, progress handling, and chain
  combination.
- `R/posterior_conversion.R`: retained-draw conversion and convergence helper
  calculations.
- `src/sampler_types.h`: validated internal prior, Gram-storage, and PVE enums.
- `src/coefficient_updates.h`: shared coefficient conditional draws, mixture
  workspace, and serial R-session RNG adapter.
- `src/summary_matrices.h`: dense, sparse, block, eigen, and LD-native matrix
  backends.
- `src/parallel_blocks.h`: independent-block workers and per-block random-number
  streams that use the shared coefficient-update kernel.
- `src/gibbs_core.h`: templated coefficient, hyperparameter, PVE, and posterior
  accumulation kernel.
- `src/gibbs.cpp`: exported Rcpp entry points and backend dispatch.
- `src/eigen_preprocessing.cpp`: memory-efficient eigen sufficient-statistic
  preprocessing and scaled-factor construction.
- `R/fit_result.R`: posterior assembly and scale restoration.
- `R/fit_methods.R`: coefficient extraction, prediction, and summaries.
- `inst/benchmarks/benchmark_ld_preprocessing.R`: configurable LD construction,
  compressed-filtering, and `ETA`-ordering benchmarks.
- `inst/benchmarks/benchmark_gig_sampler.R`: compares an independently
  implemented Hörmann--Leydold GIG rejection sampler with the current GIGrvg
  reference backend.
