## LRnetSTv2: Learning Directed Acyclic Graphs for Ligands and Receptors based on Spatial Transcriptomics Data

<img src="Fig1A_new.png" width="700" align="center">

<img src="Fig1B_new.png" width="700" align="center">

- [Reference](#Reference)
- [Overview](#Overview)
- [Installation](#Installation)
- [Usage](#Usage)
- [Arguments](#Arguments)
- [Value](#Value)
- [Examples](#Examples)
- [Contributions](#Contributions)

## Reference

Shrabanti Chowdhury, Sammy Ferri-Borgogno, Peng Yang, Wenyi Wang, Jie Peng, Samuel C Mok, Pei Wang.
Learning directed acyclic graphs for ligands and receptors based on spatially resolved transcriptomic data of ovarian cancer.
Briefings in Bioinformatics, Volume 26, Issue 2, March 2025, https://doi.org/10.1093/bib/bbaf085

## Overview

LRnetSTv2 is an improved implementation of LRnetST for learning directed acyclic graphs (DAGs) from spatial transcriptomics (ST) data. The data contain mixed continuous (log-transformed count) and binary (on/off indicator) nodes, with zero-inflation modelled explicitly via paired node structures. For a detailed account of the C++ and R-layer changes relative to LRnetST, see the implementation notes (`LRnetSTv2_implementation_notes.md`).

## Installation

```r
library(devtools)
install_github("jie108/LRnetST", subdir="LRnetSTv2")
```

or alternatively

```r
install.packages("remotes")
remotes::install_github("jie108/LRnetST", subdir="LRnetSTv2")
```

## Usage

```
SC_prepare: Prepare a log-count ST matrix for hcSC / hcSC_boot by creating paired
  (logCount, indicator) node structure with whitelisted indicator → logcount edges.

LRnetSTv2::SC_prepare(logCount)


hcSC: Learn a DAG from ST data (no bootstrap) by hill climbing for mixtures of
  continuous and binary variables.

LRnetSTv2::hcSC(Y, nodeType, whiteList, blackList, scale, tol, maxStep, restart, seed, verbose)


hcSC_boot: Learn a DAG from each bootstrap resample of the ST data; supports
  parallel execution via the future backend.

LRnetSTv2::hcSC_boot(Y, n.boot, nodeType, whiteList, blackList, scale, tol, maxStep,
                     restart, seed, nodeShuffle, bootDensityThre,
                     backend, workers, verbose, output_type)


score_shd: Aggregate a 3D bootstrap array of DAGs by minimising generalised structural
  Hamming distance (gSHD).

LRnetSTv2::score_shd(boot.adj, alpha, freq.cutoff, whiteList, blackList, max.step, verbose)


score_shd_freq: Aggregate a bootstrap frequency matrix by minimising gSHD.
  Pairs with hcSC_boot(output_type = "freq") to avoid storing the full 3D array.

LRnetSTv2::score_shd_freq(freq, alpha, freq.cutoff, whiteList, blackList, max.step, verbose)
```

## Arguments

### Arguments for `SC_prepare`

| Parameter  | Description |
| :--------- | :---------- |
| logCount   | An n by p matrix of log-transformed count values (e.g. log2(count + 1)). Rows are spots/cells, columns are genes. |

### Arguments for `hcSC` and `hcSC_boot`

| Parameter | Default | Description |
| :-------- | :-----: | :---------- |
| Y | | An n by p data matrix: n – sample size, p – number of variables. When using the SC workflow, pass `prep$Y` from `SC_prepare`. |
| n.boot *(hcSC_boot only)* | 100 | Number of bootstrap resamples. |
| nodeType | NULL | A character vector of length p specifying node type: `"c"` for continuous, `"b"` for binary. When using `SC_prepare`, pass `prep$nodeType`. Defaults to all `"c"` when NULL. |
| whiteList | NULL | A p by p logical matrix; entry `[i,j] = TRUE` forces edge i → j into every learned DAG. When using `SC_prepare`, pass `prep$whiteList` (indicator_i → logcount_i edges). |
| blackList | NULL | A p by p logical matrix; entry `[i,j] = TRUE` forbids edge i → j. When using `SC_prepare`, pass `prep$blackList` (logcount_i → indicator_i edges). Diagonal is always blacklisted. |
| scale | TRUE | Logical: L2-normalise each continuous column so that `‖Y[,i]‖²/n = 1` (zero pattern is preserved). |
| tol | 1e-6 | Minimum BIC improvement required to accept a hill-climbing step. |
| maxStep | 2000 | Maximum number of hill-climbing steps per restart. |
| restart | 1 | Number of random restarts. The best-scoring DAG across all restarts is returned. |
| seed | 1 | Integer seed for the mt19937 random number generator (used for restart tie-breaking and, in `hcSC_boot`, for bootstrap resampling). |
| nodeShuffle *(hcSC_boot only)* | FALSE | Logical: randomly permute the node ordering before each bootstrap DAG search. |
| bootDensityThre *(hcSC_boot only)* | 0.1 | Minimum column-wise nonzero fraction allowed in any bootstrap resample (rejection sampling). Must be strictly between 0 and the lowest observed column nonzero rate. |
| backend *(hcSC_boot only)* | `"sequential"` | Execution backend: `"sequential"` (default, runs in the current session) or `"future"` (launches a `future::multisession` plan). |
| workers *(hcSC_boot only)* | NULL | Number of parallel workers for `backend = "future"`. `NULL` uses all available cores minus one, capped at `n.boot`. |
| output_type *(hcSC_boot only)* | `"array"` | Output format: `"array"` (p × p × n.boot integer array), `"freq"` (p × p frequency matrix), or `"both"` (a list with both). Use `"freq"` or `"both"` with `score_shd_freq` to avoid holding the full array in memory. |
| verbose | FALSE | Logical: print step-by-step information. |

### Arguments for `score_shd` and `score_shd_freq`

| Parameter | Default | Description |
| :-------- | :-----: | :---------- |
| boot.adj *(score_shd only)* | | A p by p by B numeric array of bootstrap adjacency matrices. Typically `hcSC_boot(output_type = "array")` output. |
| freq *(score_shd_freq only)* | | A p by p numeric matrix of bootstrap edge frequencies (values in [0, 1]). Typically `hcSC_boot(output_type = "freq")` output. |
| alpha | 1 | Generalised SHD weight: `gSF(i,j) = SF(i,j) + (1 − α/2)·SF(j,i)`. Larger α penalises undirected edges more and produces sparser aggregated DAGs. |
| freq.cutoff | 0.5 | Minimum generalised score required for a candidate edge to be considered. Edges with `gSF(i,j) ≤ freq.cutoff` are excluded. |
| whiteList | NULL | A p by p logical or 0/1 matrix: entry `[i,j] = TRUE/1` forces edge i → j into the aggregated DAG. **Pass the same `whiteList` used in `hcSC_boot`** so that forced edges are initialised correctly before the greedy aggregation. |
| blackList | NULL | A p by p logical or 0/1 matrix: entry `[i,j] = TRUE/1` forbids edge i → j. **Pass the same `blackList` used in `hcSC_boot`** to keep the aggregated DAG consistent with the per-replicate constraints. |
| max.step | NULL | Deprecated; emits a warning and has no effect. |
| verbose | FALSE | Logical: print edge-addition information. |

## Value

### Value for `SC_prepare`

A list with four components:

| Object    | Description |
| :-------- | :---------- |
| Y         | n by 2p matrix: `cbind(logCount, logCount > 0)`. Columns 1..p are continuous; columns p+1..2p are binary indicators. |
| nodeType  | Character vector of length 2p: `"c"` for columns 1..p, `"b"` for columns p+1..2p. |
| whiteList | 2p by 2p logical matrix with `whiteList[p+i, i] = TRUE` (indicator_i → logcount_i). |
| blackList | 2p by 2p logical matrix with `blackList[i, p+i] = TRUE` (logcount_i → indicator_i) and `TRUE` on the diagonal. |

### Value for `hcSC`

A list with four components:

| Object     | Description |
| :--------- | :---------- |
| adjacency  | p by p 0/1 integer matrix: the adjacency matrix of the learned DAG (`adj[i,j] = 1` means edge i → j). |
| score      | Numeric vector of BIC scores at each accepted hill-climbing step. |
| operations | Integer matrix recording the operation (1 = add, 2 = delete, 3 = reverse) and the two nodes involved at each step. |
| deltaMin   | Numeric vector of the minimum candidate score change evaluated at each step. |

### Value for `hcSC_boot`

Depends on `output_type`:

| `output_type` | Return value |
| :------------ | :----------- |
| `"array"` | p × p × n.boot integer array; slice `[,,b]` is the adjacency matrix from resample b. |
| `"freq"` | p × p numeric matrix of edge frequencies (proportion of resamples containing each edge). |
| `"both"` | A list with `$adjacency` (array) and `$freq` (frequency matrix). |

### Value for `score_shd` and `score_shd_freq`

A p by p 0/1 integer matrix: the adjacency matrix of the aggregated DAG.

## Examples

### Example 1: All-continuous nodes (built-in example dataset)

```r
rm(list=ls())
library(LRnetSTv2)
data(example)

Y.n      <- example$Y        # 102 × 102 log-count matrix
true.dir <- example$true.dir # 102 × 102 true DAG (109 edges)
p <- ncol(Y.n)   # 102
n <- nrow(Y.n)   # 102

# (i) Single-run DAG learning
res <- LRnetSTv2::hcSC(
  Y       = Y.n,
  nodeType = rep("c", p),
  scale   = TRUE, maxStep = 1000, tol = 1e-6,
  restart = 1, seed = 1, verbose = FALSE
)
adj.single <- res$adjacency   # 322 edges

#  Evaluation
## DAG
sum(adj.single == 1 & true.dir == 0) / sum(adj.single == 1)   # FDR:  0.8803681
sum(adj.single == 1 & true.dir == 1) / sum(true.dir == 1)  # Power: 0.3577982

## Skeleton
true.ske    <- LRnetSTv2::skeleton(true.dir)
adj.single.ske <- LRnetSTv2::skeleton(adj.single)
sum(adj.single.ske == 1 & true.ske == 0) / sum(adj.single.ske == 1)  # FDR: 0.7055215  
sum(adj.single.ske == 1 & true.ske == 1) / sum(true.ske == 1)     # Power: 0.8807339


# (ii) Bootstrap DAG learning (sequential; use backend = "future" for parallel)
boot.adj <- LRnetSTv2::hcSC_boot(
  Y       = Y.n, n.boot = 100,
  nodeType = rep("c", p),
  scale   = TRUE, tol = 1e-6, maxStep = 1000,
  restart = 1, seed = 1, nodeShuffle = TRUE,
  bootDensityThre = 0.1, 
  backend = "future", workers = 5,
  output_type = "array",
  verbose = FALSE
)

# (iii) Bootstrap aggregation
adj.bag <- LRnetSTv2::score_shd(boot.adj, alpha = 1, freq.cutoff = 0.5)

# (iv) Evaluation
## DAG
sum(adj.bag == 1 & true.dir == 0) / sum(adj.bag == 1)   # FDR:    0.3557692
sum(adj.bag == 1 & true.dir == 1) / sum(true.dir == 1)  # Power: 0.6146789

## Skeleton
true.ske    <- LRnetSTv2::skeleton(true.dir)
adj.bag.ske <- LRnetSTv2::skeleton(adj.bag)
sum(adj.bag.ske == 1 & true.ske == 0) / sum(adj.bag.ske == 1)  # FDR:   0.1346154
sum(adj.bag.ske == 1 & true.ske == 1) / sum(true.ske == 1)     # Power: 0.8256881

## Moral graph
true.moral    <- LRnetSTv2::moral_graph(true.dir)
adj.bag.moral <- LRnetSTv2::moral_graph(adj.bag)
sum(adj.bag.moral == 1 & true.moral == 0) / sum(adj.bag.moral == 1)  # FDR:   0.1666667
sum(adj.bag.moral == 1 & true.moral == 1) / sum(true.moral == 1)     # Power: 0.7065217

## V-structures
true.vstr    <- LRnetSTv2::vstructures(true.dir)
adj.bag.vstr <- LRnetSTv2::vstructures(adj.bag)
vstr.corr    <- LRnetSTv2::compare.vstructures(adj.bag.vstr, true.vstr)
1 - nrow(vstr.corr) / nrow(adj.bag.vstr)  # FDR:   0.3076923
nrow(vstr.corr) / nrow(true.vstr)          # Power: 0.4675325
```

### Example 2: SC workflow with zero-inflated data

This example demonstrates the recommended workflow for ST count data. Starting from
the built-in `example$Y` (all-continuous), we simulate zero-inflation by independently
masking each entry with a Bernoulli(0.75) draw. The observed data are `(Z, U)` where
`Z = Y * U` (zero-inflated log-counts) and `U` (binary on/off indicators). The true
DAG has a 2p × 2p block structure: the original Z→Z edges plus a fixed U_j → Z_j
edge for every gene j.

```r
rm(list=ls())
library(LRnetSTv2)
data(example)

Y.n      <- example$Y        # 102 × 102 log-count matrix (no zeros)
true.dir <- example$true.dir # 102 × 102 true DAG (109 edges)
p <- ncol(Y.n)
n <- nrow(Y.n)

# (i) Simulate zero-inflated SC data
set.seed(42)
U <- matrix(rbinom(n * p, size = 1, prob = 0.75), nrow = n, ncol = p)
Z <- Y.n * U          # Z_ij = Y_ij * U_ij; ~25% zeros

# Observed data: n × 2p matrix (Z columns first, then U columns)
Y.sc      <- cbind(Z, U)
node.type <- c(rep("c", p), rep("b", p))

# (ii) Construct white/blacklists: U_j -> Z_j (white), Z_j -> U_j (black)
whiteList <- matrix(FALSE, 2*p, 2*p)
blackList <- matrix(FALSE, 2*p, 2*p)
diag(blackList) <- TRUE
for (j in seq_len(p)) {
  whiteList[p + j, j] <- TRUE   # U_j -> Z_j always included
  blackList[j, p + j] <- TRUE   # Z_j -> U_j always excluded
}

# (iii) True 2p × 2p DAG
true.dir.2p <- matrix(0L, 2*p, 2*p)
true.dir.2p[1:p, 1:p] <- true.dir          # Z -> Z block (109 edges)
for (j in seq_len(p)) true.dir.2p[p+j, j] <- 1L  # U_j -> Z_j (102 edges)
# total: 211 true edges

# (iv) Bootstrap DAG learning
boot.sc <- LRnetSTv2::hcSC_boot(
  Y               = Y.sc, n.boot = 100,
  nodeType        = node.type,
  whiteList       = whiteList, blackList = blackList,
  scale           = TRUE, tol = 1e-6, maxStep = 1000,
  restart         = 1, seed = 1, nodeShuffle = TRUE,
  bootDensityThre = 0.05, 
  output_type = "array",
  backend = "future", workers = 5, verbose = FALSE
)

# (v) Bootstrap aggregation
adj.bag.2p <- LRnetSTv2::score_shd(boot.sc, alpha = 1, freq.cutoff = 0.5,
                                    whiteList = whiteList, blackList=blackList)

# (vi) Evaluation on the full 2p × 2p model
## DAG
sum(adj.bag.2p == 1 & true.dir.2p == 0) / sum(adj.bag.2p == 1)   # FDR:   0.2470588
sum(adj.bag.2p == 1 & true.dir.2p == 1) / sum(true.dir.2p == 1)  # Power: 0.6066351

## Skeleton
true.ske.2p    <- LRnetSTv2::skeleton(true.dir.2p)
adj.bag.ske.2p <- LRnetSTv2::skeleton(adj.bag.2p)
sum(adj.bag.ske.2p == 1 & true.ske.2p == 0) / sum(adj.bag.ske.2p == 1)  # FDR:   0.08235294
sum(adj.bag.ske.2p == 1 & true.ske.2p == 1) / sum(true.ske.2p == 1)     # Power: 0.7393365

## Moral graph
true.moral.2p    <- LRnetSTv2::moral_graph(true.dir.2p)
adj.bag.moral.2p <- LRnetSTv2::moral_graph(adj.bag.2p)
sum(adj.bag.moral.2p == 1 & true.moral.2p == 0) / sum(adj.bag.moral.2p == 1)  # FDR:   0.2336066
sum(adj.bag.moral.2p == 1 & true.moral.2p == 1) / sum(true.moral.2p == 1)     # Power: 0.4734177

## V-structures
true.vstr.2p    <- LRnetSTv2::vstructures(true.dir.2p)  # 186 true v-structures
adj.bag.vstr.2p <- LRnetSTv2::vstructures(adj.bag.2p)
vstr.corr.2p    <- LRnetSTv2::compare.vstructures(adj.bag.vstr.2p, true.vstr.2p)
1 - nrow(vstr.corr.2p) / nrow(adj.bag.vstr.2p)  # FDR:   0.5945946
nrow(vstr.corr.2p) / nrow(true.vstr.2p)          # Power: 0.1612903
```

## Contributions

If you find bugs or have suggestions, please email the maintainer at <jiepeng108@gmail.com>. Contributions (via pull requests or otherwise) are welcome.
