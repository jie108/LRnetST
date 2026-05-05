# LRnetSTv2 Implementation Notes

Updated: 2026-05-04

`LRnetSTv2` is a local-only package derived from `LRnetST`. The root `.gitignore`
excludes the original package folder from version control:

```text
LRnetST/
```

## Purpose

The package implements hill-climbing DAG learning for spatial transcriptomics data
(zero-inflated, mixed continuous/binary nodes). The design follows the same
improvement strategy used in `dagbagMv2` relative to `dagbagM`:

- replace the single overloaded `delta` cache with separate add/delete/reverse caches
  plus version stamps;
- fix the `acyclicUpdate` double-`if` bug in the reverse block;
- replace global `srand/rand` with a local `mt19937` generator;
- precompute row-selection indicators (`indRow`) once before the HC loop;
- replace full-OLS/full-logistic list returns with loglik-only internal functions;
- pre-allocate a shared Eigen workspace to avoid per-call heap allocation in the
  hot loop.

The R-facing API is preserved except that `hcSC_boot_parallel` is merged into
`hcSC_boot(n.thread=)` using the `future` backend.

---

## Key Structural Difference from dagbagM: Row Subsetting

LRnetST introduces a spatial-transcriptomics-specific modification to the BIC
score computation for continuous nodes. Ordinary continuous regression uses all
`n` samples; here, each continuous node `i` (a log-normalised count) is fitted on
only the **nonzero** rows:

```r
indRow_i <- Y[, i] != 0
```

The motivation is biological: the zero-inflation is assumed to reflect whether
the gene is expressed at all (captured separately by a binary indicator node),
so only the expressed cells should inform inference about expression level.

The row-selection set is fixed for each node throughout the HC run because the
data matrix `Y` does not change during search. In the original `hcSC1` it was
recomputed every step; `LRnetSTv2` precomputes it once.

---

## Main C++ Changes

File:

```text
LRnetSTv2/src/dag_package_SC.cpp
```

### 1. Separate Caches with Version Stamps

The original `hcSC1` used a single `delta(p, p)` matrix. Its semantics were
overloaded: `delta(j, i)` recorded the score change for deletion when `j -> i`
was present, or for addition when it was absent. Reverse deltas were the sum of
two entries.

`LRnetSTv2` uses three separate cache arrays plus a per-node version counter,
matching the `dagbagMv2` design:

```cpp
std::vector<OneCache> addCache(p * p);
std::vector<OneCache> deleteCache(p * p);
std::vector<OneCache> reverseCache(p * p);
std::vector<int>      version(p, 0);
```

Each `OneCache` entry stores the cached score and delta along with the version
stamp(s) at which they were computed:

```cpp
struct OneCache {
  double value, scoreA, scoreB;
  int versionA, versionB;
  OneCache() : value(kInf), scoreA(kInf), scoreB(kInf), versionA(-1), versionB(-1) {}
};
```

A cache entry is valid iff its stored version equals the current `version[node]`.
Default stamp `-1` guarantees a miss on first access (`version` initialises to 0).
When an accepted operation changes a node's parent set, that node's version is
incremented. Add and delete increment `version[to]`. Reverse increments both
`version[to]` and `version[from]`.

The reverse cache reuses `deleteCache[idx].scoreA` for the "node `to` after
deletion" component, since both the delete and the reverse-from-`to`-perspective
require the identical computation.

### 2. Fixed acyclicUpdate Double-if Bug

In the original `acyclicUpdate` function, the `lastOper == 3` (reverse) block
contained two independent `if` statements for the same candidate type:

```cpp
// Original — two separate ifs can both fire
if(oper[2]==1 && acyStatus(i,j) && fromDe[i] && toAn[j]){
    acyStatus(i,j) = false;
}
if(oper[2]==1 && (!acyStatus(i,j)) && toDe[i] && fromAn[j]){
    acyStatus(i,j) = acyclicCheck(...);  // re-check can overwrite the false just set
}
```

If both condition sets hold simultaneously (the candidate pair lies in both the
"new path created" region and the "old path freed" region of the reversed edge),
the first block sets the entry to `false`, and the second block immediately
re-checks and may set it back to `true`. This can produce an incorrectly `true`
acyclicity status.

`LRnetSTv2` uses `if/else if`, matching the `dagbagMv2` fix:

```cpp
if (acyStatus(i,j) == true && last.fromDe[i] && last.toAn[j]) {
    // New edge may create cycle; if also in freed-path region, BFS is authoritative.
    if (last.toDe[i] && last.fromAn[j])
        acyStatus(i,j) = acyclic_add(i, j, parents) ? true : false;
    else
        acyStatus(i,j) = false;
} else if (acyStatus(i,j) == false && last.toDe[i] && last.fromAn[j]) {
    // Old path was freed; re-check.
    acyStatus(i,j) = acyclic_add(i, j, parents) ? true : false;
}
```

The same pattern is applied to the reverse-candidate (`operType == 3`) sub-block.

### 3. mt19937 Replaces srand/rand

The original tied tie-breaking to the global C `rand()` via `srand(seed)`. This
is non-reentrant and low quality. `LRnetSTv2` seeds a local generator per HC run:

```cpp
std::mt19937 rng(static_cast<unsigned int>(seed));
std::uniform_int_distribution<int> coin(0, 1);
```

Restarts space seeds by 101 to ensure diversity:

```cpp
const int iseed = static_cast<int>(
    static_cast<unsigned int>(seed) +
    static_cast<unsigned int>(i) * 101U);
```

### 4. Precomputed indRow (Row-Selection Indicators)

The original `hcSC1` called `indRow_i = (Y_i != 0)` on every step of the outer
node loop, repeating `n` comparisons per node per step even though `Y` is
constant during HC.

`LRnetSTv2` precomputes before entering the HC loop:

```cpp
std::vector<std::vector<char>> indRow(p, std::vector<char>(n, 1));
std::vector<int> n_sel(p, n);
for (int i = 0; i < p; ++i) {
    if (types[i] == "c") {
        int cnt = 0;
        for (int r = 0; r < n; ++r) {
            indRow[i][r] = (Y(r, i) != 0.0) ? 1 : 0;
            cnt += indRow[i][r];
        }
        n_sel[i] = cnt;
    }
}
```

`indRow[i]` uses `char` (not `LogicalVector`) to avoid Rcpp overhead in the inner
loop. `n_sel[i]` is the effective sample size used in the BIC penalty for
continuous node `i`.

### 5. Loglik-Only Score Functions

The original `fastLM2_` returned a full R list (coefficients, fitted values,
residuals, loglikelihood) from which only `loglikelihood` was used. Similarly,
`fastLR2_` returned a full logistic regression list.

`LRnetSTv2` replaces these with direct scalar-returning functions:

```cpp
double continuous_loglik(const ConstMatRef& X, const ConstVecRef& Y);
double logistic_loglik(const ConstMatRef& X, const ConstVecRef& Y);
```

`continuous_loglik` uses `rankUpdate` to fill only the lower triangle of `X'X`
(~2× fewer FLOPs for symmetric gram matrices), as in `dagbagMv2`:

```cpp
MatrixXd XtX(q, q);
XtX.setZero();
XtX.selfadjointView<Lower>().rankUpdate(X.adjoint());
LLT<MatrixXd> llt(XtX);
```

### 6. Pre-Allocated Workspace

The original code allocated new `Rcpp::NumericMatrix` objects for every design
matrix in the HC inner loop via `matColSelAddOne`, `matRowSelInd`, `addCol`, and
`delCol`. These are hot-path allocations.

`LRnetSTv2` pre-allocates a single workspace and fills it in-place:

```cpp
MatrixXd ws(n, p + 1);   // shared across all score computations
VectorXd y_buf(n);        // row-selected response buffer
```

For continuous nodes, `fill_design_sel` writes selected rows directly into the
top of `ws`. `continuous_loglik` then receives an explicit copy of
`ws.topLeftCorner(n_sel, q)` because the block's outer stride (= `n`) differs
from `n_sel`, which prevents `Ref<const MatrixXd>` from binding without a copy:

```cpp
const MatrixXd X(ws.topLeftCorner(n_sel, q));  // small copy, negligible vs Cholesky
const double loglik = continuous_loglik(X, y_buf.head(n_sel));
```

For binary nodes, `fill_design` writes all `n` rows and `ws.leftCols(q)` binds
directly.

---

## Bug Analysis: Constant Columns After Row Subsetting

### The Setup

`SC_prepare` creates a paired node structure: for each gene `i`, the data matrix
contains a continuous column `logCount_i` (columns `1..p`) and a binary indicator
column `S_i = (logCount_i > 0)` (columns `p+1..2p`). It whitelists the edge
`S_i -> logCount_i` to ensure the indicator is always a parent of the logcount.

After row selection on `logCount_i != 0`, the indicator column `S_i` is
**constant** in the selected rows (always 1, because `S_i = 1 ↔ logCount_i > 0`
by construction). With the intercept also all-ones, the design matrix has two
identical columns:

```
matRowSelInd([S_i, ones], indRow_i) = [[1, 1],
                                        [1, 1],
                                        ...  ]   (n_sel × 2, rank 1)
```

### Original LRnetST Behaviour (Near-Singular but Functionally Correct)

The original `fastLM2_` calls Eigen's `LLT` on `X'X` without checking
`llt.info()`. For the all-ones design matrix, `X'X = [[n_sel, n_sel],[n_sel, n_sel]]`.
The Cholesky gives `L[1,1] = sqrt(n_sel - n_sel)`. In floating-point arithmetic
this may be exactly `0` or a tiny residual, but not `NaN`.

When `L[1,1] ≈ 0`, the triangular solve distributes the fit across `beta[0]` and
`beta[1]` with large but opposite-sign values. Their sum satisfies:

```
beta[0] + beta[1] ≈ mean(Y_sel)
```

so `fitted ≈ mean(Y_sel) * ones`. The computed loglikelihood is therefore the
correct **intercept-only** log-likelihood. The BIC then counts `q = 2` (both
columns), giving:

```
BIC_orig = -2 * logL_intercept + log(n_sel) * (2 - 1)
         = BIC_intercept_only + log(n_sel)
```

This `+log(n_sel)` offset is present in **both** the current and candidate
scores whenever a logcount node is being scored. It cancels exactly in every
delta comparison:

```
delta = BIC_new - BIC_old
      = (correct_BIC_new + log(n_sel)) - (correct_BIC_old + log(n_sel))
      = correct_delta
```

**Verified empirically**: running the same data through original `LRnetST` and
`LRnetSTv2` produces identical adjacency matrices. The per-node BIC absolute
scores differ by exactly `log(n_sel_i)` for continuous nodes and 0 for binary
nodes, confirming the cancellation.

### Why It Is Still Fragile

1. **Relies on undefined / implementation-specific Eigen behaviour.** The
   triangular solve with a zero (or near-zero) diagonal element is not
   guaranteed to produce the intercept-only fit. Different Eigen versions,
   compiler flags (`-O2` vs `-O3`, FMA), or target architectures may produce
   `NaN`, `Inf`, or a subtly wrong residual sum.

2. **Breaks silently if the indicator is not exactly collinear.** If floating-
   point scaling causes `S_i` to differ from `ones` by even a single ULP in any
   selected row, the Cholesky succeeds fully and `beta[0]`, `beta[1]` are both
   estimated, distorting the BIC penalty.

3. **The `+log(n_sel)` offset inflates absolute BIC values.** Tools that use
   node scores directly (score comparison across models, information criteria for
   model selection) would receive biased values.

### LRnetSTv2 Fix

`fill_design_sel` silently drops any column that is constant across the selected
rows before passing the design matrix to `continuous_loglik`:

```cpp
auto fill_col_if_nonconstant = [&](int src) {
    double first_val = 0.0;
    bool first_set = false, is_constant = true;
    for (int r = 0; r < n; ++r) {
        if (!indRow[r]) continue;
        const double v = Y(r, src);
        ws(r_out++, col) = v;
        if (!first_set) { first_val = v; first_set = true; }
        else if (v != first_val) is_constant = false;
    }
    if (!is_constant) ++col;
};
```

This makes `XtX` always full-rank for continuous nodes: `XtX` is (at most)
`q × q` where `q` = non-constant parents + 1 (intercept). The `+log(n_sel)`
overcount is removed, so absolute BIC values are correct.

---

## R Layer Changes

Files:

```text
LRnetSTv2/R/hcSC.R
LRnetSTv2/R/auxiliary.R
LRnetSTv2/R/score_shd_new.R
```

### hcSC.R

**Input validation.** Both `hcSC` and `hcSC_boot` validate `Y`, `nodeType`,
`whiteList`, and `blackList` at the R level before the C++ call. Error messages
identify the parameter and the size mismatch.

**L2 scaling loop.** `1:p` replaced by `seq_len(p)` throughout.

**Seed scheme.** The original used `set.seed(i*1001+seed)` for bootstrap
resampling and `i*11+seed` for the HC. These multiplicative schemes can overflow
for large `i`. `LRnetSTv2` generates all seeds upfront:

```r
if (seed > .Machine$integer.max - 2L * n.boot)
    stop("seed is too large relative to n.boot")
boot_seeds <- seed + seq_len(n.boot)         # R RNG per bootstrap resample
hc_seeds   <- seed + n.boot + seq_len(n.boot)  # C++ HC per bootstrap resample
```

The ranges `[seed+1, seed+n.boot]` and `[seed+n.boot+1, seed+2*n.boot]` are
disjoint, ensuring the bootstrap sampling and HC tie-breaking use independent
seeds.

**Parallel bootstrap via `future`.** The `hcSC_boot_parallel` function has been
merged into `hcSC_boot(n.thread=)`. When `n.thread > 1`, a `future::multisession`
plan is activated for the duration of the bootstrap loop and restored on exit:

```r
oplan <- future::plan(future::multisession, workers = n.thread)
on.exit(future::plan(oplan), add = TRUE)
futs <- lapply(seq_len(n.boot), function(b) future::future(run_one(b), substitute=FALSE))
results <- lapply(futs, future::value)
```

`substitute = FALSE` ensures `run_one(b)` is evaluated eagerly in the calling
frame so that `b` is captured by value, not by reference.

**`boot_dense`.** Unchanged. Rejection-sampling ensures every bootstrap column
has at least `bootDensityThre` fraction of nonzero entries, preserving the
zero-inflation structure needed by the continuous node regressions.

### auxiliary.R

`1:(n.parents[i] - 1)` in `moral_graph` and `1:(n.par - 1)` in `vstructures`
replaced by `seq_len(...)` to avoid the length-zero error when a node has exactly
one parent.

### score_shd_new.R

`library(igraph)` was called inside `score_shd` to check whitelist acyclicity.
This is an impure side-effect in a scoring function and depends on `igraph` being
installed. Replaced by a comment that the caller is responsible for providing an
acyclic whitelist; acyclicity of the output is guaranteed by the greedy
`edgeOnLoop` check performed when inserting each candidate edge.

---

## DESCRIPTION and NAMESPACE

| Field | LRnetST (original) | LRnetSTv2 |
|---|---|---|
| Version | 2.0 | 2.1 |
| Imports | `Rcpp` | `Rcpp`, `future`, `stats` |
| LinkingTo | `Rcpp, RcppEigen, RcppNumerical, doParallel, foreach` | `Rcpp, RcppEigen, RcppNumerical` |
| Exported | `hcSC, hcSC_boot, hcSC_boot_parallel, SC_prepare, score_shd, cal_order, moral_graph, vstructures, skeleton, compare.vstructures` | same minus `hcSC_boot_parallel` |

`doParallel` and `foreach` are removed from `LinkingTo` (they are R packages, not
C++ header libraries; listing them there was incorrect). Parallelism is now handled
entirely in R via `future`.

---

## Summary of Changes

| # | Category | Original LRnetST | LRnetSTv2 |
|---|---|---|---|
| 1 | Cache | Single overloaded `delta(p,p)` | Separate `addCache`, `deleteCache`, `reverseCache` with version stamps |
| 2 | Acyclicity cache | Double-if in reverse block (can overwrite false→true) | `if/else if` prevents spurious re-check |
| 3 | RNG | Global `srand(seed)` / `rand()` | Per-run `std::mt19937`, seeded at construction |
| 4 | `indRow` | Recomputed each step inner loop | Precomputed once before HC loop; `std::vector<char>` |
| 5 | Score functions | `fastLM2_` / `fastLR2_` return full OLS/logistic lists | `continuous_loglik` / `logistic_loglik` return scalar |
| 6 | `AtA` proxy | `selfadjointView` proxy return (potential UB on some compilers) | Removed; `rankUpdate` inlined with explicit `MatrixXd XtX` |
| 7 | Workspace | Per-call `Rcpp::NumericMatrix` heap allocation in hot loop | Pre-allocated `MatrixXd ws(n, p+1)` reused each call |
| 8 | Constant columns | `S_i` column becomes all-ones after row selection; near-singular XtX tolerated silently | `fill_design_sel` drops constant columns; XtX always full-rank |
| 9 | Seed scheme | `i*1001+seed` (R) / `i*11+seed` (C++), overflow-prone | `seed + seq_len(n.boot)`, overflow guard, disjoint ranges |
| 10 | Parallelism | `foreach + doParallel` in separate `hcSC_boot_parallel` | `future::multisession` via `hcSC_boot(n.thread=)` |
| 11 | `seq_len` | `1:(n-1)` in `moral_graph`, `vstructures` | `seq_len(n-1)` |
| 12 | `library()` in function | `library(igraph)` inside `score_shd` | Removed |
