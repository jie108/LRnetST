\name{LRnetST-functions}
\alias{SC_prepare}
\alias{hcSC}
\alias{hcSC_boot}
\alias{score_shd}
\alias{score_shd_freq}
\alias{cal_order}
\alias{moral_graph}
\alias{vstructures}
\alias{skeleton}
\alias{compare.vstructures}
\title{LRnetST user functions}
\description{
User-facing helpers for preparing spatial transcriptomics data, learning DAGs
with hill climbing, bootstrap DAG learning, bootstrap aggregation, and graph
summaries.
}
\usage{
SC_prepare(logCount)

hcSC(Y, nodeType = NULL, whiteList = NULL, blackList = NULL,
     scale = TRUE, tol = 1e-6, maxStep = 2000L,
     restart = 1L, seed = 1L, verbose = FALSE)

hcSC_boot(Y, n.boot = 100L, nodeType = NULL, whiteList = NULL,
          blackList = NULL, scale = TRUE, tol = 1e-6,
          maxStep = 2000L, restart = 1L, seed = 1L,
          nodeShuffle = FALSE, bootDensityThre = 0.1,
          backend = c("sequential", "future"),
          workers = NULL, verbose = FALSE,
          output_type = c("array", "freq", "both"))

score_shd(boot.adj, alpha = 1, freq.cutoff = 0.5,
          whiteList = NULL, blackList = NULL,
          max.step = NULL, verbose = FALSE)

score_shd_freq(freq, alpha = 1, freq.cutoff = 0.5,
               whiteList = NULL, blackList = NULL,
               max.step = NULL, verbose = FALSE)

cal_order(adj_matrix)

moral_graph(adj.matrix)

vstructures(adj.matrix)

skeleton(adj.matrix)

compare.vstructures(target.vstructures, true.vstructures)
}
\arguments{
\item{logCount}{Numeric log-count matrix with rows as spots/cells and columns as genes.}
\item{Y}{Numeric data matrix. For the spatial transcriptomics workflow, use
\code{SC_prepare(logCount)$Y}.}
\item{nodeType}{Character vector with entries \code{"c"} for continuous nodes and
\code{"b"} for binary nodes. Defaults to all continuous when \code{NULL}.}
\item{whiteList}{Logical or 0/1 matrix of required directed edges.}
\item{blackList}{Logical or 0/1 matrix of forbidden directed edges.}
\item{scale}{Whether to L2-normalise continuous columns before fitting.}
\item{tol}{Minimum BIC improvement threshold.}
\item{maxStep}{Maximum hill-climbing steps per restart.}
\item{restart}{Number of hill-climbing restarts.}
\item{seed}{Nonnegative integer seed for bootstrap sampling, node shuffling, and
hill-climbing tie-breaking.}
\item{verbose}{Whether to print progress.}
\item{n.boot}{Number of bootstrap resamples.}
\item{nodeShuffle}{Whether to randomly permute node order within each bootstrap fit.}
\item{bootDensityThre}{Minimum allowed column-wise nonzero fraction in each
bootstrap resample.}
\item{backend}{Execution backend: \code{"sequential"} runs replicates in the
current session; \code{"future"} dispatches them with \code{future::multisession}.}
\item{workers}{Number of parallel workers for \code{backend = "future"}.
\code{NULL} uses all available cores minus one, capped at \code{n.boot}.}
\item{output_type}{Bootstrap output mode: \code{"array"} returns a
\code{p x p x B} adjacency array; \code{"freq"} returns an edge frequency matrix;
\code{"both"} returns a list with both.}
\item{boot.adj}{A \code{p by p by B} array of bootstrap DAG adjacency matrices.}
\item{freq}{A \code{p by p} edge frequency matrix with values in \code{[0, 1]}.}
\item{alpha}{Positive generalized SHD reversal weight parameter.}
\item{freq.cutoff}{Edge inclusion cutoff in \code{[0, 1]}.}
\item{max.step}{Deprecated legacy argument with no effect.}
\item{adj_matrix, adj.matrix}{Square adjacency matrix for a directed graph.}
\item{target.vstructures}{Candidate v-structure matrix with columns
\code{par1}, \code{child}, and \code{par2}.}
\item{true.vstructures}{Reference v-structure matrix with columns
\code{par1}, \code{child}, and \code{par2}.}
}
\value{
\code{SC_prepare} returns a list with \code{Y}, \code{nodeType},
\code{whiteList}, and \code{blackList}. \code{hcSC} returns a list with a learned
adjacency matrix and fit diagnostics. \code{hcSC_boot} returns an adjacency
array, a frequency matrix, or a list with both, depending on \code{output_type}.
\code{score_shd} and \code{score_shd_freq} return an integer adjacency matrix.
Graph summary helpers return a topological order, skeleton/moral adjacency
matrix, v-structure matrix, or matching v-structures.
}
\examples{
set.seed(1)
logCount <- matrix(rexp(60, rate = 1), nrow = 20, ncol = 3)
logCount[logCount < 0.8] <- 0
prep <- SC_prepare(logCount)

fit <- hcSC(prep$Y, prep$nodeType, prep$whiteList, prep$blackList,
            maxStep = 20, seed = 1)
fit$adjacency

boot_freq <- hcSC_boot(prep$Y, n.boot = 3, nodeType = prep$nodeType,
                       whiteList = prep$whiteList,
                       blackList = prep$blackList,
                       maxStep = 20, output_type = "freq")
score_shd_freq(boot_freq, whiteList = prep$whiteList,
               blackList = prep$blackList)
}
\keyword{models}
