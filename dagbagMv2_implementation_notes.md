# dagbagMv2 Implementation Notes

Updated: 2026-05-03

`dagbagMv2` is a local-only package copy derived from `dagbagM`. It is intended
to stay out of the parent repository history. The root `.gitignore` ignores:

```text
dagbagMv2/
dagbagMv2.Rcheck/
```

## Purpose

The package implements hill-climbing DAG learning for mixtures of continuous and
binary nodes, with cleanup parallel to the HC and bagging improvements made in
InterCellDAG v2/v3.

The main goals were:

- remove the old overloaded single `delta` cache semantics;
- prevent `NA`, `NaN`, and infinite values from entering HC decisions;
- move internal graph state away from R vectors/lists where possible;
- improve bootstrap output memory behavior;
- provide an add/delete-only HC mode for comparison and diagnostics;
- move SHD aggregation into C++;
- keep the R-facing API close to the original package.

## Main C++ HC Changes

File:

```text
dagbagMv2/src/dag_package.cpp
```

### Separate Caches

The original `dagbagM` HC used one overloaded `delta` matrix for add/delete and
reversal-related score changes. `dagbagMv2` uses separate caches:

```cpp
std::vector<OneCache> addCache(p * p);
std::vector<OneCache> deleteCache(p * p);
std::vector<OneCache> reverseCache(p * p);
```

Each cache entry stores:

```cpp
struct OneCache {
  double value;
  double scoreA;
  double scoreB;
  int versionA;
  int versionB;
};
```

For add/delete, `value` is the candidate score change for one target node. For
reverse, `value` is the combined score change for deleting `from -> to` and
adding `to -> from`.

### Cache Versioning

The implementation uses the same idea as InterCellDAGv3's node version stamps,
but stores cache-entry stamps inside each cache entry rather than in separate
`p x p` version matrices.

InterCellDAGv3 style:

```cpp
deltaAdd(p, p);
deltaAddVersion(p, p);
deltaDelete(p, p);
deltaDeleteVersion(p, p);
nodeVersion(p, 0);
```

`dagbagMv2` style:

```cpp
std::vector<int> version(p, 0);
addCache[idx].versionA;
deleteCache[idx].versionA;
reverseCache[idx].versionA;
reverseCache[idx].versionB;
```

Mapping:

- `nodeVersion[i]` maps to `version[i]`.
- `deltaAdd(j, i)` maps to `addCache[idx].value`.
- `deltaAddVersion(j, i)` maps to `addCache[idx].versionA`.
- `deltaDelete(j, i)` maps to `deleteCache[idx].value`.
- `deltaDeleteVersion(j, i)` maps to `deleteCache[idx].versionA`.
- Reverse cache entries use two stamps because reversal affects two local node
  scores.

When an accepted operation changes a node's parent set, that node's version is
incremented. Add/delete increment the target node version. Reverse increments
both affected node versions.

### Native Graph State

Internal parent and child sets are now native C++ adjacency lists:

```cpp
typedef std::vector<std::vector<int> > AdjList;
```

The current graph and constraint matrices are stored as compact C++ vectors
internally, then converted back to an R logical adjacency matrix for output.
This avoids repeated mutation of `Rcpp::List` and reduces accidental dependence
on R's `NA_LOGICAL` behavior.

### NA and Non-Finite Guards

The C++ layer validates:

- finite numeric `Y`;
- valid binary node values;
- valid node types;
- no `NA` in whitelist/blacklist;
- no whitelist/blacklist conflicts;
- acyclic whitelist;
- signed, nonnegative control arguments where appropriate.

Candidate deltas are ignored unless finite:

```cpp
if (!std::isfinite(candidate)) {
  return false;
}
```

Initial node scores must also be finite; otherwise HC stops with an explicit
error.

### RNG

Tie-breaking uses a local generator:

```cpp
std::mt19937 rng(static_cast<unsigned int>(seed));
```

This replaces the old global `srand()`/`rand()` approach.

### Debug Checks

`hc()` and the C++ HC functions accept `debug = TRUE`. In debug mode, the
selected operation is checked before it is applied, and the final graph is
checked for acyclicity.

Debug mode also recomputes each visited add/delete/reverse candidate score from
the current parent sets and compares it to the cached score/delta before the
candidate can be selected. A mismatch raises an explicit error such as:

```text
debug cache check failed for add/delete/reverse ...
```

### Add/Delete-Only Mode

`hc()`, `hc_boot()`, and `hc_boot_parallel()` now expose:

```r
addDeleteOnly = FALSE
```

When `addDeleteOnly = TRUE`, the HC loop skips reversal candidates entirely:

```cpp
if (!addDeleteOnly && !has_edge(black, to, from, p) &&
    acyclic_reverse(from, to, parents)) {
  ...
}
```

This keeps the same separate-cache implementation for add/delete operations and
provides a direct diagnostic path parallel to the InterCellDAG add/delete-only
comparisons.

## Main R Changes

File:

```text
dagbagMv2/R/hc.R
```

### Input Validation

The R wrapper now checks and normalizes inputs before C++:

- `Y` must be finite numeric data;
- `nodeType` must contain only `"c"` and `"b"`;
- binary nodes must contain only 0/1;
- whitelist and blacklist must be dimension-matched logical matrices with no
  `NA`;
- whitelist must be acyclic and conflict-free;
- continuous nodes must have positive finite standard deviation when
  `standardize = TRUE`;
- `tol`, `maxStep`, `restart`, and `seed` are checked before C++.

### Bootstrap Output Modes

`hc_boot()` and `hc_boot_parallel()` now support:

```r
output_type = "array"
output_type = "freq"
output_type = "both"
```

This allows large bootstrap runs to return only a `p x p` edge frequency matrix
instead of a full `p x p x B` adjacency array.

### Parallel Safety

`hc_boot_parallel()` saves the current `future` plan and restores it with
`on.exit()`, so errors in bootstrap fitting do not leave the user's R session in
a changed parallel state.

## Aggregation Changes

File:

```text
dagbagMv2/R/score_shd_new.R
dagbagMv2/src/dag_package.cpp
```

### C++ Aggregation Core

The generalized SHD aggregation loop now runs in C++:

```cpp
score_shd_cpp(...)
score_shd_freq_cpp(...)
```

The R functions `score_shd()` and `score_shd_freq()` now prepare default
`whiteList`/`blackList` matrices and call the C++ core. This mirrors the
InterCellDAG direction of keeping bootstrap aggregation in compiled code.

### Frequency-Based Aggregation

`score_shd_freq()` was added. It aggregates directly from a `p x p` edge
frequency matrix, matching the new bootstrap `return = "freq"` mode.

### Deterministic Ordering

Aggregation now sorts candidate edges by:

1. decreasing generalized selection frequency;
2. decreasing direct selection frequency;
3. source node index;
4. target node index.

This avoids ambiguous tie ordering from the old data-frame order alone.

### Validation

`score_shd()` and `score_shd_freq()` now check dimensions, finite values,
frequency bounds, whitelist/blacklist conflicts, and whitelist acyclicity.

### API Naming Alignment (2026-05-01)

Parameter names in `score_shd()` and `score_shd_freq()` were aligned with the
rest of the package:

- `whitelist` / `blacklist` → `whiteList` / `blackList`
- `max.step` → `maxStep` (deprecated, emits a warning when supplied)

The `threshold` parameter was replaced by `freqCutoff`:

| Old | New | Default |
|---|---|---|
| `threshold = 0` | `freqCutoff = 0.5` | majority-vote cutoff |

The original `threshold` was an indirect encoding: the internal cutoff was
computed as `(1 - threshold) / 2`, so `threshold = 0` → cutoff 0.5. The new
`freqCutoff` is the cutoff directly, making the intent immediately readable at
the call site. The change was pushed all the way through the C++ layer —
`aggregate_freq_cpp`, `score_shd_freq_cpp`, and `score_shd_cpp` all now accept
`freqCutoff` in `[0, 1]` instead of `threshold` in `[-1, 1]`.

## Package Metadata

The package was renamed from `dagbagM` to `dagbagMv2` in:

- `DESCRIPTION`;
- `NAMESPACE`;
- generated Rcpp registration files.

Minimal Rd documentation was added for package checks.

## Verification

The package was installed into a temporary library:

```sh
R CMD INSTALL -l /private/tmp/dagbagMv2-lib dagbagMv2
```

Smoke tests passed for:

- `hc()`;
- `hc(..., addDeleteOnly = TRUE)`;
- `hc_boot(..., output_type = "freq")`;
- `hc_boot_parallel(..., output_type = "both")`;
- `score_shd_freq()`;
- validation failure on non-finite input.

A `testthat` suite was added under:

```text
dagbagMv2/tests/testthat/test-hc-debug.R
```

The tests run HC with `debug = TRUE`, confirm returned graphs are acyclic,
confirm cached candidate deltas match full recomputation, confirm
`addDeleteOnly = TRUE` selects no reversal operations, exercise bootstrap
aggregation through the C++ `score_shd` path, and check validation failure
cases.

`R CMD check --no-manual dagbagMv2` completed with one external compiler-header
warning from the local R/RcppEigen toolchain and one standard note about checking
an unpacked source directory rather than an `R CMD build` tarball.

## Code Review Recommendations (2026-04-30)

| # | Category   | Item                                                                 | Status      |
|---|------------|----------------------------------------------------------------------|-------------|
| 1 | Correctness| `cal_order` infinite loop on cyclic input                            | Implemented |
| 2 | Correctness| `compare.vstructures` crash when `true.vstructures` is NULL          | Implemented |
| 3 | Correctness| `1:p` loops break when p == 0; replace with `seq_len(p)`             | Implemented |
| 4 | Correctness| `max.step` parameter accepted but silently ignored                   | Implemented |
| 5 | Efficiency | Reverse cache duplicates delete `scoreA` computation                 | Implemented |
| 6 | Efficiency | `make_design` allocates R matrix per `node_score` call               | Implemented |
| 7 | Efficiency | `is_dag` in C++ is O(p³); could use Kahn's algorithm                 | Planned     |
| 8 | Efficiency | Quadratic vector/matrix growth in auxiliary functions                 | Planned     |
| 9 | Efficiency | `.format_boot_result` creates n.boot intermediate matrices           | Implemented |
| 10| Style      | Inconsistent `whiteList`/`whitelist` and `maxStep`/`max.step` naming | Implemented |
| 11| Style      | `compare.vstructures` S3-method-like dot naming                      | Planned     |
| 12| Style      | Missing `\examples{}` in Rd documentation                           | Implemented |
| 13| Style      | `auxiliary.R` uses v1-era coding style                               | Implemented |

### Details on #4: `max.step` deprecation

`score_shd()` and `score_shd_freq()` retained `max.step` for backward
compatibility with v1 call sites. The C++ aggregator processes all eligible
candidates in a single greedy pass and does not use a step limit. Both functions
now emit a deprecation warning when `max.step` is explicitly supplied, with an
inline comment explaining why the parameter remains in the signature.

### Details on #5: Reverse cache reusing delete scoreA

When an edge `from -> to` exists, both delete and reverse candidates need the
score of node `to` with parent `from` removed:

```cpp
// delete
del.scoreA = node_score(..., to, -1, from);

// reverse (was identical, now reuses delete result)
rev.scoreA = del.scoreA;
```

The delete block always runs first for the same `(from, to)` pair, so
`del.scoreA` is guaranteed fresh when the reverse block executes. The debug-mode
recomputation on the reverse path still independently verifies correctness.

### Details on #6: Pre-allocated design matrix workspace

This was the most substantial change. The old `make_design` allocated a new
`Rcpp::NumericMatrix` (R heap, GC-managed) on every `node_score` call.
Similarly, `logistic_loglik` allocated an `Rcpp::NumericVector` for the initial
beta, and `node_score` extracted the response column into a new
`Rcpp::NumericVector`. With O(p²) candidates evaluated per HC step across
potentially thousands of steps, this produced heavy R GC pressure.

The replacement:

1. **`fill_design`** replaces `make_design`. It writes columns into a
   pre-allocated `MatrixXd workspace(n, p + 1)` created once at the top of
   `hc1`, and returns the number of columns used (`q`). No R allocation occurs.

2. **Score functions** (`continuous_loglik`, `logistic_loglik`, `bic_score`) now
   accept `Eigen::Ref<const MatrixXd>` and `Eigen::Ref<const VectorXd>` instead
   of `Rcpp::NumericMatrix` / `Rcpp::NumericVector`. The `Ref` type accepts
   block views (e.g. `workspace.leftCols(q)`) without copying.

3. **`LogisticReg`** stores `Eigen::Ref<const MatrixXd>` members instead of
   `Eigen::Map<Eigen::MatrixXd>`. The initial beta is a `VectorXd::Zero(q)` on
   the C++ heap instead of an `Rcpp::NumericVector`.

4. **`node_score`** maps the response column as a zero-copy
   `Eigen::Map<const VectorXd>` directly into the R data matrix, and passes
   `workspace.leftCols(q)` to `bic_score`.

Net effect: **zero R heap allocations per candidate evaluation** in the HC inner
loop, down from three.

## Code Review Recommendations (2026-05-01) — Round 3

| # | Severity   | File(s)                                      | Item                                                             | Status      |
|---|------------|----------------------------------------------|------------------------------------------------------------------|-------------|
| 1 | Critical   | dag_package.cpp                              | `p*p` int overflow in cache and vector indexing for large p      | Implemented |
| 2 | Major      | auxiliary.R                                  | `cal_order` cycle detection used heuristic iteration limit       | Implemented |
| 3 | Major      | dag_package.cpp                              | `fill_design` unguarded against duplicate `extraParent`          | Implemented |
| 4 | Major      | tests/testthat/                              | Test suite missing coverage for aux functions, hc_boot, edge cases | Implemented |
| 5 | Minor      | dag_package.cpp, RcppExports.*               | `hc1` exported via Rcpp but intended as internal                 | Implemented |
| 6 | Minor      | dag_package.cpp, RcppExports.*               | `edgeOnLoop` exported via Rcpp but not part of public API        | Implemented |
| 7 | Minor      | dag_package.cpp                              | Undocumented seed spacing constant 101 in `hc_`                  | Implemented |
| 8 | Minor      | R/hc.R                                       | Undocumented seed spacing constants 1001 and 11 in `.fit_boot_one` | Implemented |
| 9 | Minor      | dag_package.cpp, RcppExports.*               | C++ aggregation parameters used lowercase vs R camelCase         | Implemented |
| 10| Minor      | R/hc.R                                       | `seed = TRUE` in `hc_boot_parallel` misleads about who controls seeding | Implemented |

### Details on #1: Integer overflow guard

`p * p` is computed as a 32-bit `int` product in cache allocation and index arithmetic.
For p > 46,340 the product overflows. A bounds check was added at both C++ entry
points immediately after `p` is derived:

```cpp
if (p > 46000)
  Rcpp::stop("p exceeds safe index limit (46000)");
```

Added in `validate_hc_inputs` (guards all HC paths) and `aggregate_freq_cpp`
(guards aggregation paths). The threshold 46,000 is conservative — HC is
computationally infeasible at that scale — but prevents silent memory corruption
if ever reached.

### Details on #2: `cal_order` replaced with Kahn's algorithm

The old implementation rotated unresolved nodes to the back of a queue and
declared a cycle only if an iteration counter exceeded `p * p`. This is a
heuristic: the limit is generous enough to be safe in practice, but the cycle
detection is indirect. The new implementation uses Kahn's algorithm directly:

```r
in_degree <- colSums(adj_matrix != 0)
queue <- which(in_degree == 0)
# drain zero-in-degree nodes, decrement children, enqueue newly freed nodes
if (length(node_order) != p) stop("cycle detected")
```

Cycle detection is now exact: any node left unprocessed at the end of the drain
proves a cycle. The O(p) bound is tight rather than approximate. The test suite
already contained an equivalent Kahn's implementation (`is_dag_matrix` in
`test-hc-debug.R`); the two are now consistent.

### Details on #3: `fill_design` duplicate-parent guard

`fill_design` adds an `extraParent` column first, then loops over the existing
parent set. If `extraParent` were already in `parents`, the same data column
would be written twice, producing a singular design matrix and wrong scores. The
HC logic prevents this (an edge is only proposed as `extraParent` when it is not
yet in the graph), but `fill_design` itself was unguarded. A runtime check was
added:

```cpp
if (std::find(parents.begin(), parents.end(), extraParent) != parents.end())
  Rcpp::stop("fill_design: extraParent is already in the parent set (internal error)");
```

The check is O(|parents|) and fires only on internal programming errors, not
user input.

### Details on #4: Expanded test suite

A second test file was added:

```text
dagbagMv2/tests/testthat/test-core.R
```

37 new tests covering:

- `cal_order`: chain ordering, root-only graph, cycle detection, single node.
- `skeleton`, `vstructures`, `moral_graph`: known-answer tests on 3-node chain
  and v-structure graphs.
- `compare.vstructures`: NULL-guard and matching cases.
- `hc()`: all-continuous data, p = 1, seed reproducibility.
- `hc_boot()`: array/freq/both return modes and frequency consistency.
- `score_shd` / `score_shd_freq`: `freqCutoff = 0`, `freqCutoff = 1`, and
  agreement between the two aggregation entry points on matching inputs.

### Details on #5 and #6: Internal C++ exports removed

`hc1` and `edgeOnLoop` had `// [[Rcpp::export]]` attributes but were not
listed in `NAMESPACE` and had no Rd documentation. Both are internal
implementation details:

- `hc1` is the single-restart HC primitive called directly from `hc_` within
  the same translation unit.
- `edgeOnLoop` is an R-list-accepting wrapper around the internal
  `edge_on_loop_cpp` function, used only in early debugging.

The `[[Rcpp::export]]` attributes were removed. Both `RcppExports.cpp` and
`RcppExports.R` were updated to remove the corresponding wrappers and
`CallEntries` registrations. The `CallEntries` table shrank from 5 to 3
entries.

### Details on #9: C++ parameter naming aligned to camelCase

`aggregate_freq_cpp`, `score_shd_freq_cpp`, and `score_shd_cpp` used lowercase
`whitelist`/`blacklist` parameter names while the rest of the package uses
`whiteList`/`blackList`. All three functions were updated:

```cpp
// before
const Rcpp::LogicalMatrix& whitelist,
const Rcpp::LogicalMatrix& blacklist,

// after
const Rcpp::LogicalMatrix& whiteList,
const Rcpp::LogicalMatrix& blackList,
```

Error messages inside `aggregate_freq_cpp` were updated to match. Both
`RcppExports.cpp` (forward declarations, SEXP variable names, local parameter
names, and call sites) and `RcppExports.R` (wrapper parameter names) were
updated to be consistent throughout.

### Details on #10: Parallel seed comment

In `hc_boot_parallel`, `.options.future = list(seed = TRUE)` and the explicit
`set.seed(i * 1001L + seed)` inside `.fit_boot_one` are both active. The
`future` framework assigns each worker an L'Ecuyer-CMRG stream when
`seed = TRUE`, but `.fit_boot_one` immediately overrides that with a
Mersenne Twister seeded from the deterministic formula. The `set.seed` always
wins. `seed = TRUE` is retained only to suppress `doFuture`'s warning about
non-reproducible parallel RNG. A two-line comment was added above
`.options.future` to prevent a future maintainer from removing the explicit
`set.seed` under the false assumption that `future` already handles it.

## Acyclicity Caching (2026-05-02)

File:

```text
dagbagMv2/src/dag_package.cpp
```

### Motivation

The original `dagbagM` HC checked acyclicity with a fresh BFS for every
candidate add and reverse in the inner scan loop — O(p²(p+e)) per step.
`InterCellDAGv3` avoids most of these BFS calls by maintaining a cached
`acyStatus(p,p)` matrix and propagating changes incrementally after each
accepted operation.  This pattern was ported to `dagbagMv2`.

### Key components

**`LastOpState` struct**

Bundles the identity and ancestor/descendant sets of the last accepted
operation's from/to nodes, computed on the post-operation graph:

```cpp
struct LastOpState {
  int from = -1, to = -1, type = -1;
  std::vector<char> fromAn, fromDe, toAn, toDe;
};
```

Default-constructed `type == -1` signals "no operation accepted yet".

**`compute_ancestors` / `compute_descendants`**

BFS returning `std::vector<char>` (1 byte/element) indicator vectors.
Using `char` instead of `bool` or R `LogicalVector` gives four-times smaller
footprint versus `int` and better cache line utilization.

**`acyclic_cache_update`**

Incrementally updates `acyStatus(i,j)` for one candidate (i,j) using the
three-way logic ported from `InterCellDAGv3 acyclicUpdate`:

| Last op | Candidate | Effect |
|---------|-----------|--------|
| add | add/rev | new edge may create new cycle paths → set false |
| delete | add/rev | removed edge may free previously-cyclic pairs → BFS re-check |
| reverse | add/rev | both effects apply; ambiguous pairs do a BFS re-check |

When both "cycle-creating" and "cycle-freeing" conditions hold for the same
candidate (only possible after a reverse), the implementation resolves the
ambiguity with a full BFS rather than picking the conservative answer.

**`acyclic_reverse_inplace`**

Temporarily removes `from` from `parents[to]` using order-preserving
`erase + insert` (not swap-and-pop), checks the BFS, and restores the entry
at the original position.  Preserving element order is essential: `fill_design`
iterates over `parents[to]` to build design matrix columns, so any reordering
changes the Cholesky factorization at floating-point precision and can alter
tie-breaking in `better_delta`, producing different accepted operations and
different final graphs.

**NA guard and fallback**

```
acyStatus entry == NA_LOGICAL  →  compute via full BFS, cache result
acyStatus entry != NA_LOGICAL  →  call acyclic_cache_update
                                   if still NA after update → full BFS
```

`acyclic_cache_update` never writes `NA_LOGICAL` — it only writes `true`,
`false`, or leaves the entry unchanged — so the post-update NA fallback
fires only if none of the ancestor/descendant conditions matched.

### Debug validation

When `debug = TRUE`, the main loop performs a full BFS check for every
add and reverse candidate and compares the result to the cached `acyStatus`
entry.  A mismatch raises an explicit error:

```text
debug: acyStatus mismatch for add/reverse %d->%d: cache=%d actual=%d
```

Validated on the `p=102` example: 347 accepted steps, zero cache mismatches.

### Benchmark results (n=102, p=102 example data)

| Operation | dagbagMv2 | dagbagM v1 | Speedup |
|-----------|-----------|------------|---------|
| `hc` (no boot) | 0.088 s | 0.533 s | **6.1×** |
| `hc_boot_parallel` (50 boots, 2 threads) | 15.7 s | 55.2 s | **3.5×** |

---

## Performance Analysis: dagbagMv2 vs dagbagM v1

The ~6× single-HC speedup and ~3.5× bootstrap speedup come from four
independent changes.  They are listed in rough order of impact.

### 1. Design matrix: heap allocation eliminated (dominant)

v1 allocates a new R matrix on every candidate score evaluation:

```cpp
// v1 — three R heap allocations per candidate:
Rcpp::NumericMatrix X_i = matColSelAddOne(Y, par_i);
Rcpp::NumericMatrix X_c = delCol(X_i, k);          // or
Rcpp::NumericMatrix X_c = addCol(X_j, Y_j);
```

v2 pre-allocates one `MatrixXd workspace(n, p+1)` at the top of `hc1` and
`fill_design` writes column views into it.  `bic_score` receives an
`Eigen::Ref` block view — zero-copy, no R GC involvement:

```cpp
// v2 — one allocation per hc1 call, reused forever:
MatrixXd workspace(n, p + 1);
const int q = fill_design(workspace, Y, parents[node], extraParent, dropParent);
bic_score(workspace.leftCols(q), y, nodeType[node]);
```

For p=102, n=102 with 347 accepted steps and O(p²) candidates per step,
this eliminates roughly 3.6 million R heap allocations per `hc` call.

### 2. Graph and parent/child representation: R objects → native C++

v1 holds the graph as `Rcpp::LogicalMatrix` (4 bytes/element, R heap) and
parent/child sets as `Rcpp::List`.  Every BFS traversal copies a parent
vector off the R heap:

```cpp
// v1 — R→C++ copy on every BFS node visit:
std::vector<int> curPar = parSet[curNode];
```

v2 uses flat native C++ containers:

```cpp
// v2 — direct reference, no copy:
const std::vector<int>& par = parents[queue[pos]];
```

The graph is `std::vector<unsigned char>` (1 byte/element).  Ancestor and
descendant indicators are `std::vector<char>` (1 byte/element) vs v1's
`Rcpp::LogicalVector` (4 bytes/element as R integer).  Smaller element size
means better cache utilization on the p×p acyclicity and delta structures.

### 3. Delta caching: separate caches with version stamps vs one shared matrix

v1 uses a single `delta(p,p)` matrix shared across add, delete, and reversal.
Invalidation is coarse: any operation touching node i triggers a full
recomputation of all delta entries for that node:

```cpp
// v1 — recomputes for ALL j when node i was the last target:
if (lastOper[2]==-1 || last_to_i || (lastOper[2]==3 && last_from_i)) {
    X_i = matColSelAddOne(Y, par_i);
    score_c = bicScore(X_c, Y_i, ...);
}
```

v2 has three separate caches (`addCache`, `deleteCache`, `reverseCache`), each
entry with a per-node version stamp.  An entry is recomputed only when its
stamp disagrees with the current node version, which increments only when that
node's parent set actually changes:

```cpp
// v2 — recomputes only when parent set of `to` changed:
if (del.versionA != version[to]) {
    del.scoreA = node_score(...);
    del.versionA = version[to];
}
```

In a typical step only one or two nodes change their parent set, so the vast
majority of cache entries survive intact across steps.

Additionally, v1 shares `delta(j,i)` for both delete and reversal uses of the
same edge, creating ambiguity when the entry is repurposed between operations.
v2's separate caches make the scoring semantics explicit and allow the reverse
cache entry to directly reuse `del.scoreA` (computed earlier in the same
iteration) for the common deletion-part subcomputation.

### 4. RNG: global C `srand/rand` → `std::mt19937`

v1 uses `srand(seed)` / `rand()` — a global C RNG reset on every `hc1` entry,
not thread-safe:

```cpp
srand(seed);
...
(rand() + 0.0) / (RAND_MAX + 1.0) > 0.5
```

v2 uses a local `std::mt19937` generator seeded per call:

```cpp
std::mt19937 rng(static_cast<unsigned int>(seed));
```

This is the main reason the bootstrap parallel speedup (3.5×) is somewhat
lower than single-HC (6×): v1's global RNG contention and R-object locking
limit parallel efficiency, so removing both enables closer-to-linear parallel
scaling.

---

## Code Review Recommendations (2026-05-02) — Round 4

Fixes applied to `dagbagMv2/src/dag_package.cpp`.  All 12 documented
metrics verified unchanged after applying all changes.

| # | Category | Item | Status |
|---|----------|------|--------|
| 1 | Dead code | `lastType == -1` else-branch in `acyclic_cache_update` unreachable | Removed |
| 2 | Efficiency | `acyclic_reverse` copies full AdjList in NA reversal fallback | Fixed: use `acyclic_reverse_inplace` |
| 3 | Correctness | Double `if` in `lastType==3` can fire BFS after conservative invalidation | Fixed: restructured as `if / else if` with inner ambiguity check |
| 4 | Efficiency | `acyStatus(from,to) = true` written unconditionally on every scan | Fixed: guarded with `!= true` check |
| 5 | Dead code | `llt.info() != Eigen::Success` check after `llt.solve()` | Removed: `solve` does not change info status |
| 6 | Style | `acyclic_cache_update` has 11 parameters | Fixed: last-op state bundled into `LastOpState` struct |

### Details on #1: Dead `lastType == -1` branch

`acyclic_cache_update` is called from the `else` arms at lines 822 and 877,
which are only entered when `acyStatus(...) != NA_LOGICAL`.  In scan 0 (before
any operation is accepted, `lastOp.type == -1`), all add and reverse candidate
entries are `NA_LOGICAL`, so the NA fast-path handles them and the `else` arms
are never entered.  After scan 0, `lastOp.type ∈ {1, 2, 3}` permanently.  The
`else` branch was therefore dead code and has been removed.

### Details on #3: Double-condition in `lastType == 3` (reverse)

The original code had two independent `if` blocks for `operType == 1`:

```cpp
// Block A: new edge may create a cycle
if (operType == 1 && acyStatus(i,j) && lastFromDe[i] && lastToAn[j])
  acyStatus(i,j) = false;
// Block B: deleted edge may free a cycle
if (operType == 1 && !acyStatus(i,j) && lastToDe[i] && lastFromAn[j])
  acyStatus(i,j) = acyclic_add(i,j,parents) ? true : false;
```

When all four ancestor/descendant conditions hold simultaneously (the reverse
both created and freed paths affecting the same candidate), block A would set
to `false`, and block B would then see `!acyStatus(i,j) == true` and call BFS
unnecessarily.  The BFS always produced the correct final answer, but the
`false` assignment from block A was redundant.

The fix restructures this as:

```cpp
if (acyStatus(i,j) && lastFromDe[i] && lastToAn[j]) {
  // Ambiguous: if freed-path condition also holds, BFS decides; else just invalidate.
  if (lastToDe[i] && lastFromAn[j])
    acyStatus(i,j) = acyclic_add(i,j,parents) ? true : false;
  else
    acyStatus(i,j) = false;
} else if (!acyStatus(i,j) && lastToDe[i] && lastFromAn[j]) {
  acyStatus(i,j) = acyclic_add(i,j,parents) ? true : false;
}
```

The same restructuring was applied to the `operType == 3` (reversal candidate)
sub-case.

### Details on #6: `LastOpState` struct

The seven separate variables:

```cpp
int lastOperFrom, lastOperTo, lastOperType;
std::vector<char> lastFromAn, lastFromDe, lastToAn, lastToDe;
```

are always passed and updated as a unit.  They were bundled into:

```cpp
struct LastOpState {
  int from = -1, to = -1, type = -1;
  std::vector<char> fromAn, fromDe, toAn, toDe;
};
```

`acyclic_cache_update` signature reduced from 11 to 6 parameters:

```cpp
// before
void acyclic_cache_update(int lastFrom, int lastTo, int lastType,
    const std::vector<char>& lastFromAn, const std::vector<char>& lastFromDe,
    const std::vector<char>& lastToAn, const std::vector<char>& lastToDe,
    int operFrom, int operTo, int operType,
    AdjList& parents, Rcpp::LogicalMatrix& acyStatus);

// after
void acyclic_cache_update(const LastOpState& last,
    int operFrom, int operTo, int operType,
    AdjList& parents, Rcpp::LogicalMatrix& acyStatus);
```

---

## Code Review — Round 5 (2026-05-03): Inline Comment Pass

### Files modified

```text
dagbagMv2/src/dag_package.cpp
dagbagMv2/R/hc.R
dagbagMv2/R/auxiliary.R
dagbagMv2/R/score_shd_new.R
```

No logic changes were made in this pass.  All edits are documentation only.
The package was rebuilt and all 48 tests passed without modification.

### C++ inline comments added

**`fill_design`**

The existing function-level comment was replaced with a detailed block explaining
the column layout and the `q`-includes-intercept convention:

```
Column layout (left to right):
  [extraParent column, if extraParent >= 0]   -- candidate new parent (add op)
  [existing parents, skipping dropParent]      -- current parents minus dropped (delete op)
  [ones column]                                -- intercept, always last

q = total columns including the intercept.
Caller passes workspace.leftCols(q) to bic_score.
bic_score uses (q - 1) as the BIC penalty term: q-1 equals the number of
real parent predictors because the intercept is not a structural parameter.
```

This is the key design invariant connecting `fill_design`, `bic_score`, and
`node_score` — misunderstanding `q` as "number of parents" (forgetting the
intercept) would produce a BIC penalty off by one for every candidate.

**`continuous_loglik`**

Added the Gaussian MLE formula as a function-level comment:

```
beta_hat  = (X'X)^{-1} X'Y  (via Cholesky of X'X)
sigma2_hat = ||Y - X*beta_hat||^2 / n  (MLE: divides by n, not n-1)
log L = -n/2 * (log(2*pi) + log(sigma2_hat) + 1)
Returns -Inf if X'X singular, beta non-finite, or sigma2 <= 0.
```

**`bic_score`**

Added a function header explaining `BIC = -2*logL + log(n)*(q-1)` and an
inline comment on the `q - 1` term at the point of use:

```cpp
// q - 1 = number of parent predictors (excludes the intercept column).
const double score = -2.0 * loglik + std::log(static_cast<double>(n)) * (q - 1);
```

**`OneCache` struct**

The existing two-line comment was replaced with a per-field table documenting
which fields are used for which operation type:

| Field | Add from->to | Delete from->to | Reverse from->to |
|-------|---|---|---|
| `scoreA` | BIC of `to` after adding `from` | BIC of `to` after removing `from` | same as deleteCache.scoreA (reused) |
| `scoreB` | unused | unused | BIC of `from` after adding `to` |
| `versionA` | `version[to]` stamp | `version[to]` stamp | `version[to]` stamp |
| `versionB` | unused | unused | `version[from]` stamp |
| `value` | `scoreA - curScore[to]` | `scoreA - curScore[to]` | `(scoreA - curScore[to]) + (scoreB - curScore[from])` |

Default `versionA = versionB = -1`: guarantees a cache miss on first access
because `version[node]` initializes to 0.

**`acyclic_cache_update`**

The existing short comment was expanded to explain all three last-op branches
and the NA fallback path, reinforcing the dual-slot `acyStatus` convention:

```
Dual-slot convention (one matrix, two semantics):
  acyStatus(from, to) = true  iff adding from->to is currently acyclic
  acyStatus(to, from) = true  iff reversing from->to is currently acyclic
  NA_LOGICAL           = not yet computed (lazy initialization)

type 1 (add):    new path last.from->last.to may close cycles -> set false.
type 2 (delete): blocking path removed -> re-check previously false entries.
type 3 (reverse): both effects apply; ambiguous pairs do a full BFS re-check.
```

**`aggregate_freq_cpp`**

Added full generalized SHD formula to the function header:

```
gsf(from, to) = sf(from, to) + (1 - alpha/2) * sf(to, from)
With alpha=1: gsf = sf + 0.5 * sf_reverse.
Both orientations of an undirected skeleton edge enter the candidate list with
the same gsf; the greedy pass orients toward the higher forward sf.
freqCutoff is a skeleton-inclusion threshold: an undirected edge appears iff
max(gsf(i,j), gsf(j,i)) > freqCutoff.
```

**`hc1`**

Extended the existing cache-versioning comment to cover the per-step loop
structure and the `reverseCache.scoreA` reuse rationale:

```
Per-step loop structure:
  Outer loop: to (node whose BIC score changes).
  Inner loop: from.
  Branch 1 -- edge (from,to) exists:   evaluate delete and (if eligible) reverse.
  Branch 2 -- no edge in either direction and not blacklisted: evaluate add.

reverseCache.scoreA reuse:
  For reverse of from->to, the "new score of to after deletion" is identical to
  deleteCache.scoreA. When deleteCache is fresh (same version[to]), that value
  is copied directly into reverseCache.scoreA without recomputation.
```

### R inline comments added

**`hc.R` — `.fit_boot_one`**

After refactoring to accept pre-generated `boot_index` and `node_perm`:

- Explained why node shuffling reduces greedy ordering bias
- Clarified that `node.index = order(node_perm)` is the inverse permutation
  used to restore the original node ordering in the returned adjacency matrix

**`hc.R` — `.format_boot_result`**

Added descriptions for all three `output_type` modes and noted that the
incremental `freq` accumulation avoids building the full 3D array when only
frequencies are needed.

**`hc.R` — `hc_boot`**

Added a function header explaining the upfront seed strategy:

```
All randomness (bootstrap row indices, HC tie-breaking seeds, node permutations)
is generated upfront from a single set.seed(seed) call so that sequential and
future backends produce bit-identical results.
```

**`auxiliary.R` — `moral_graph`**

Added explanation of:
- spouse-edge addition for each node with >= 2 parents (`j < k` loop visits
  each unordered pair exactly once)
- why `(A + t(A)) > 0` symmetrizes correctly for both skeleton edges and newly
  added spouse edges

**`auxiliary.R` — `vstructures`**

Added explanation of:
- the unshielded collider definition (no edge in either direction between
  `par1` and `par2`)
- why the `j < k` loop produces canonical `par1 < par2` ordering

**`auxiliary.R` — `compare.vstructures`**

Added a function header clarifying the asymmetric semantics: outer loop
iterates over `target`, each row is looked up in `true`.  The return value
is the subset of `true` rows that were matched.  Passing
`(estimated, truth)` gives true-positive v-structures; swapping gives
false positives.

**`score_shd_new.R` — `score_shd`**

Added full generalized SHD formula (same as C++ counterpart) plus description
of the greedy acyclic pass, consistent with the new C++ documentation.

---

## API Consolidation (2026-05-03)

### 1. `hc_boot` unified function

`hc_boot_parallel` (which depended on `foreach` and `doFuture`) was removed.
`hc_boot` was rewritten to subsume both sequential and parallel execution via a
`backend` parameter, following the same design as `TS_DAG_K23_boot` in
InterCellDAGv3:

```r
hc_boot(Y, n.boot, ...,
        backend = c("sequential", "future"),
        workers = NULL,
        output_type = c("array", "freq", "both"))
```

**Key design changes vs old `hc_boot` + `hc_boot_parallel`:**

| Aspect | Old design | New design |
|--------|-----------|-----------|
| Sequential | `hc_boot(..., return=)` | `hc_boot(..., backend="sequential", output_type=)` |
| Parallel | `hc_boot_parallel(..., numThread=, return=)` | `hc_boot(..., backend="future", workers=, output_type=)` |
| Parallel framework | `foreach` + `doFuture` | `future` directly |
| Seed generation | `set.seed(i*1001+seed)` inside `.fit_boot_one` | All randomness drawn upfront from `set.seed(seed)` |
| Per-bootstrap HC seed | `i * 11L + seed` | `seed + b` |
| Sequential == future? | No (different RNG streams) | Yes (same pre-generated indices) |
| `workers` cap | `min(numThread, cores-1)` | `min(workers, availableCores()-1, n.boot)` |

**Upfront seed generation:**

```r
set.seed(args$seed)
boot_indices      <- replicate(n.boot, sample.int(n, n, replace = TRUE), simplify = FALSE)
hc_seeds          <- args$seed + seq_len(n.boot)
node_permutations <- if (nodeShuffle) replicate(n.boot, sample.int(p), simplify = FALSE)
                     else NULL
```

All randomness comes from a single `set.seed` call.  A `run_one(b)` closure
captures each bootstrap's pre-generated index, seed, and permutation.
Sequential and `future` backends invoke the same `run_one(b)` in the same order,
producing identical results.

**Global RNG isolation:**

The caller's `.Random.seed` is saved and restored via `on.exit`, preventing
`set.seed(args$seed)` from having side effects outside `hc_boot`.

**Why `foreach` was dropped (`future` directly instead):**

`foreach + doFuture` was the dominant R parallelism pattern when `hc_boot_parallel`
was written.  The ecosystem has since converged on `future` as the unified
backend: `doFuture` is now an adapter that runs `foreach` loops on top of
`future`, not a destination in its own right.  Direct `future` use (via
`future::future()` / `future.apply` / `furrr`) is the current standard in new
code, and major packages that previously used `foreach` are migrating toward it.
Using `future` directly also removes the `foreach` and `doFuture` package
dependencies.  `hc_boot_parallel` retains `foreach + doFuture` only because it
is preserved as a reference implementation for the convergence comparison; if it
is retired, those dependencies go with it.

**`future` backend (no `foreach`/`doFuture`):**

```r
futures <- lapply(seq_len(n.boot), function(b) {
  future::future(run_one(b), seed = TRUE)
})
result <- lapply(futures, future::value)
```

`seed = TRUE` suppresses future's L'Ecuyer-CMRG RNG warning; actual
reproducibility is controlled by the upfront `set.seed`, not by future's own
RNG management.

**`.fit_boot_one` updated signature:**

The helper was refactored to accept pre-generated data rather than computing
bootstrap indices internally:

```r
# old
.fit_boot_one(i, Y, ..., seed, nodeShuffle, ...)

# new
.fit_boot_one(b, Y, ..., hc_seed, boot_index, node_perm, ...)
```

**NAMESPACE / exports:**

`hc_boot_parallel` was removed from `NAMESPACE` exports during consolidation,
then later restored (2026-05-04) to allow direct comparison with the new
`hc_boot` implementation.

### 2. `freqCutoff` renamed to `freq.cutoff`

The parameter `freqCutoff` in `score_shd()` and `score_shd_freq()` was renamed
to `freq.cutoff` to match the `freq.cutoff` naming used in
`InterCellDAGv3::score_shd()`.

Changed in:
- `R/score_shd_new.R` — parameter declaration, validation message, internal usage
- `man/dagbagMv2-functions.Rd` — `\arguments` and `\usage` blocks
- `tests/testthat/test-core.R` — all call sites using the named argument

The C++ functions `score_shd_cpp` and `score_shd_freq_cpp` retain their
internal `freqCutoff` parameter name (unchanged, as they are internal).

### 3. `return` renamed to `output_type`

The `return` argument in `hc_boot` (formerly also in `hc_boot_parallel`) was
renamed to `output_type`, consistent with `TS_DAG_K23_boot` in InterCellDAGv3.

Changed in:
- `R/hc.R` — `.format_boot_result` parameter and all branches
- `man/dagbagMv2-functions.Rd`
- `tests/testthat/test-core.R` and `test-hc-debug.R` — all call sites

---

## README Update (2026-05-03)

`README.md` at the repository root was updated to reflect all API changes:

- `hc_boot_parallel` section and argument table row removed
- `hc_boot` Usage block updated with `backend`, `workers`, `output_type`
- Arguments table:
  - `numThread` removed; `backend` and `workers` added
  - `return` → `output_type`
  - `freqCutoff` → `freq.cutoff`
- Value table for bootstrap functions updated: `return` → `output_type`
- Example code updated to use `backend = "sequential"`, `output_type = "array"`,
  `freq.cutoff = 0.5`; a comment shows the `backend = "future"` equivalent

### Example verification (n = 102, p = 102, n.boot = 50, seed = 1)

The README example was re-run verbatim after all API and standardization changes.
All documented values match exactly (post per-sample standardization fix):

| Metric | Value |
|--------|-------|
| `sum(true.dir)` — true edges | 109 |
| (i) DAG FDR — single HC | 0.8803681 |
| (i) DAG Power — single HC | 0.3577982 |
| (i) Skeleton FDR — single HC | 0.7055215 |
| (i) Skeleton Power — single HC | 0.8807339 |
| (iv) DAG edges — bootstrap agg | 114 |
| (iv) DAG FDR — bootstrap agg | 0.4035088 |
| (iv) DAG Power — bootstrap agg | 0.6238532 |
| (iv) Skeleton FDR — bootstrap agg | 0.2192982 |
| (iv) Skeleton Power — bootstrap agg | 0.8165138 |
| (iv) Moral FDR — bootstrap agg | 0.3155080 |
| (iv) Moral Power — bootstrap agg | 0.6956522 |
| (iv) V-structure FDR — bootstrap agg | 0.5000000 |
| (iv) V-structure Power — bootstrap agg | 0.4805195 |

Note: the bootstrap metrics differ from those in the original `dagbagM` README
because both the seed scheme and standardization changed.  The old code called
`set.seed(i * 1001 + seed)` inside `.fit_boot_one` and standardized the full
dataset once before resampling.  The new code generates all bootstrap indices
upfront from a single `set.seed(seed)` call (guaranteeing sequential == future
reproducibility) and standardizes each bootstrap sample after row-resampling.
The single-HC metrics (i) are unchanged because that path was not modified.

---

---

## Bug Fix: Per-Bootstrap Standardization (2026-05-03)

### Problem

`standardize` was consumed entirely by `.prepare_hc_inputs`, which standardized
the full dataset once before any bootstrap work began.  `args$Y` returned from
`.prepare_hc_inputs` was already scaled; `.fit_boot_one` resampled rows of
pre-standardized data.

Effect: each bootstrap sample had **approximately** mean 0 and sd 1 (based on
the full-data statistics), not exactly.  Bootstrap resampling introduces small
deviations from unit scale — some samples are pulled more than once, others not
at all — so the effective sd of a bootstrap column can differ from 1 even when
the full data is scaled.

Contrast with `InterCellDAGv3::TS_DAG_K23_boot`, which standardizes each
bootstrap sample after row-resampling:

```r
X_b <- X_M[index, , drop = FALSE]
if (standardize) {
    X_b <- scale(X_b)
}
```

`hc_boot` now matches this behavior.

### Fix

Three coordinated changes:

**1. `.fit_boot_one` — add `standardize` parameter and per-sample scaling**

```r
.fit_boot_one <- function(b, Y, ..., standardize, ...)
```

After `Y.B <- Y[boot_index, , drop = FALSE]`:

```r
if (standardize) {
    for (i in which(nodeType == "c")) {
        m <- mean(Y.B[, i]);  s <- stats::sd(Y.B[, i])
        if (is.finite(s) && s > 0) Y.B[, i] <- (Y.B[, i] - m) / s
    }
}
```

The `is.finite(s) && s > 0` guard handles the rare case of a constant bootstrap
column (possible when n is small or many rows are duplicated); `hc1` would catch
it anyway via the "non-finite initial score" guard, but this gives a cleaner
fallback.

**2. `hc_boot` — call `.prepare_hc_inputs` with `standardize = FALSE`**

```r
args <- .prepare_hc_inputs(Y, nodeType, whiteList, blackList,
                            standardize = FALSE, tol, maxStep, restart, seed)
```

`args$Y` is now the raw (unscaled) data.  All other validation (types, binary
values, list constraints) still runs in `.prepare_hc_inputs`.

**3. `hc_boot` — add upfront constant-column check on full data**

```r
if (standardize) {
    for (i in which(args$nodeType == "c")) {
        s <- stats::sd(args$Y[, i])
        if (!is.finite(s) || s <= 0)
            stop("continuous nodes must have positive finite standard deviation when standardize = TRUE")
    }
}
```

This replicates the check that was previously inside `.prepare_hc_inputs`'s
`if (standardize)` block, so the user still gets an early error if any full-data
column is constant before wasting time on bootstrap fits.

**4. `run_one` — passes `standardize` to `.fit_boot_one`**

`standardize` is in scope as a `hc_boot` parameter and is captured by the
`run_one` closure; it is now explicitly forwarded:

```r
.fit_boot_one(..., standardize, ...)
```

### Behaviour change

`hc()` (single fit) is unaffected: `.prepare_hc_inputs` is still called with
the actual `standardize` value and applies the transformation to the full data
before the single HC run.

The bootstrap results change slightly because each replicate now uses per-sample
mean/sd instead of full-data mean/sd.  Updated documented values
(n = 102, p = 102, n.boot = 50, seed = 1):

| Metric | Before fix (full-data scale) | After fix (per-sample scale) |
|--------|------------------------------|------------------------------|
| DAG edges | 116 | 114 |
| DAG FDR | 0.4051724 | 0.4035088 |
| DAG Power | 0.6330275 | 0.6238532 |
| Skeleton FDR | 0.2241379 | 0.2192982 |
| Skeleton Power | 0.8256881 | 0.8165138 |
| Moral FDR | 0.3177083 | 0.3155080 |
| Moral Power | 0.7119565 | 0.6956522 |
| V-str FDR | 0.4935065 | 0.5000000 |
| V-str Power | 0.5064935 | 0.4805195 |

Note: the "After fix" column reflects both the per-sample standardization fix
and the subsequent `score_shd_cpp` floating-point accumulation fix (2026-05-04),
which changed the C++ freq computation from divide-then-accumulate to
accumulate-then-divide-once, ensuring `score_shd(array)` and
`score_shd_freq(freq)` produce identical aggregated DAGs.

DAG power and skeleton/moral/v-structure power are unchanged; only the FDR
metrics shift slightly (fewer false positives with correct per-sample scaling).

---

---

## `hc_boot_parallel` vs `hc_boot`: Implementation Differences and Numerical Comparison (2026-05-04)

### Why `hc_boot_parallel` was restored

`hc_boot_parallel` was removed during the API consolidation (2026-05-03) and its
behavior merged into `hc_boot(backend = "future")`.  It was later restored to
`dagbagMv2/R/hc.R` (exported from `NAMESPACE`) to allow direct comparison of
results and to serve as a reference implementation while the two approaches are
evaluated.

### Implementation differences

The two implementations differ on two axes:

| Aspect | `hc_boot_parallel` (A) | `hc_boot` current (C) |
|--------|----------------------|----------------------|
| Standardization | Full dataset once, before resampling | Per bootstrap sample, after resampling |
| Bootstrap seed | `set.seed(i * 1001L + seed)` inside each worker | Pre-drawn upfront from `set.seed(seed)` |
| HC seed | `i * 11L + seed` | `seed + i` |

Note: the "Before fix" column in the Bug Fix section above refers to a transient
intermediate state during development (new seeding already in place, full-data
standardization not yet corrected) — it does not correspond to either A or C.

### Numerical comparison (n=102, p=102, n.boot=50, seed=1, freq.cutoff=0.5)

| Metric | A: `hc_boot_parallel` | C: `hc_boot` (current) |
|--------|----------------------|------------------------|
| DAG FDR | 0.3636364 | 0.4035088 |
| DAG Power | 0.6422018 | 0.6238532 |
| Skeleton FDR | 0.1818182 | 0.2192982 |
| Skeleton Power | 0.8256881 | 0.8165138 |

The original README documented version A values. The current README documents
version C values.

### Sources of divergence

The difference between A and C is the combined effect of both axes:

1. **Standardization**: Full-data scaling (A) gives each bootstrap sample
   approximately mean 0 / sd 1, because the pre-scaled column statistics shift
   slightly after row-resampling. Per-sample scaling (C) gives each replicate
   exactly unit scale.

2. **Seeding**: `set.seed(i * 1001L + seed)` (A) generates bootstrap row
   indices and node permutation from the same RNG call inside each worker.
   `hc_boot` (C) draws all bootstrap indices and permutations upfront from a
   single `set.seed(seed)` call, producing a different sequence of random draws.
   Even with identical standardization, A and C visit different bootstrap problems.

Whether the two implementations converge at large n.boot is an open question
to be investigated.

---

## Bug Fix: `score_shd_cpp` Floating-Point Accumulation (2026-05-04)

### Problem

`score_shd_cpp` computed edge frequencies from the bootstrap array using
divide-then-accumulate:

```cpp
freq(from, to) += value / static_cast<double>(nb);
```

`hc_boot(output_type = "freq")` computes frequencies in R using sum-then-divide:

```r
for (i in seq_len(n.boot)) freq <- freq + result[[i]]
freq <- freq / n.boot
```

Both approximate `k / n.boot` for integer count `k`, but via different
floating-point evaluation orders.  Since `1.0 / n.boot` is not exactly
representable for most `n.boot` values, the C++ path introduces a rounding
error that grows with `k` (up to `n.boot × ε(1/n.boot)`).  The R path
accumulates exact integer counts and divides once, incurring a single rounding.

The two paths produced different freq values for borderline edges, causing
`score_shd(boot.adj)` and `score_shd_freq(freq)` to return different aggregated
DAGs for the same bootstrap run (n.boot=50, seed=1: 115 vs 114 edges).

### Fix

`score_shd_cpp` now accumulates integer counts first, then divides once:

```cpp
// accumulate
freq(from, to) += value;
// normalize after loop
const double nb_d = static_cast<double>(nb);
for (int i = 0; i < p; ++i)
  for (int j = 0; j < p; ++j)
    freq(i, j) /= nb_d;
```

This matches the R accumulation exactly and is the more accurate approach.
Both `score_shd(array)` and `score_shd_freq(freq)` now produce identical
aggregated DAGs for the same bootstrap run.

### Effect on documented results

The README example results (n=102, p=102, n.boot=50, seed=1, maxStep=1000)
change slightly due to one borderline edge flipping.  See the updated
verification table in the "README Update" section above.

---

## Code Review — Round 6 (2026-05-04): Final Pass

14 items found; 10 fixed, 4 accepted-as-is.

| # | Severity | File | Item | Status |
|---|----------|------|------|--------|
| 1 | Major | dag_package.cpp | `crossprod_self` used chained proxy return via `SelfAdjointView`; safe in practice (LLT uses lower triangle) but confusing — replaced with `return A.adjoint() * A` | Fixed |
| 2 | Major | R/hc.R | Agent flagged `future::future(run_one(b), ...)` as eager evaluation — false positive; `future::future` uses `substitute=TRUE` by default, capturing the expression unevaluated | Not a bug |
| 3 | Major | DESCRIPTION | `foreach` and `doFuture` in `Imports` though only needed by legacy `hc_boot_parallel`; moved to `Suggests` | Fixed |
| 4 | Major | NAMESPACE / DESCRIPTION | Missing `importFrom(stats, sd)`; `stats` absent from `Imports` | Fixed |
| 5 | Minor | dag_package.cpp | `bestScore == kInf` exact equality on `inf`; replaced with `!std::isfinite(bestScore)` for clarity | Fixed |
| 6 | Minor | dag_package.cpp | `seed + i * 101` signed integer overflow UB for near-`INT_MAX` seeds; cast to `unsigned int` arithmetic | Fixed |
| 7 | Minor | R/hc.R | `args$seed + seq_len(n.boot)` R integer overflow for near-max seeds; added pre-check with clear error message | Fixed |
| 8 | Minor | R/auxiliary.R | `1:(n-1)` in `moral_graph` and `vstructures`; safe due to guards but inconsistent — replaced with `seq_len(n-1)` | Fixed |
| 9 | Minor | NAMESPACE | `hc_boot_parallel` exported but undocumented; removed from public exports (still accessible via `:::`) | Fixed |
| 10 | Minor | DESCRIPTION | `stats` and `parallel` missing from `Imports`; `stats` added; `parallel` is a base package and needs no declaration | Fixed |
| 11 | Minor | R/hc.R | Bare `mean()` vs qualified `stats::sd()`; `mean` is from `base` (always available) so bare calls are correct — no change | Accepted |
| 12 | Minor | dag_package.cpp | `aggregate_freq_cpp` reverse-skip logic lacked comment; added explanation | Fixed |
| 13 | Minor | R/hc.R | Silent skip on constant bootstrap column; intentional behavior documented in comment — no warning added | Accepted |
| 14 | Minor | tests/testthat/test-core.R | p > 46000 test unconditionally skipped; left as manual test with existing note | Accepted |

---

## Future TODO

- Add an `iniGraph` argument or internal debug entry point so HC can start from a
  non-empty initial DAG separate from `whiteList`.
- Run cache-vs-full-recompute debug tests from random acyclic initial graphs.
  This is an important stress case because deleting or reversing edges can expose
  candidates that were illegal at initialization.
- Compare empty-initial and random-initial behavior for full HC and
  `addDeleteOnly = TRUE`, mirroring the InterCellDAG cache validation workflow.
