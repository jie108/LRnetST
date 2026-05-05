\name{LRnetST-package}
\alias{LRnetST-package}
\alias{LRnetST}
\docType{package}
\title{LRnetST: DAG Learning for Spatial Transcriptomics Data}
\description{
  Improved implementation of LRnetST for learning directed acyclic graphs (DAGs)
  from spatial transcriptomics data containing mixed continuous (log-count) and
  binary (on/off indicator) nodes. Includes score caching, acyclicity cache
  correctness fixes, reproducible RNG, and future-based parallel bootstrap.
}
\details{
  The main workflow for ST data:
  \enumerate{
    \item Prepare data with \code{SC_prepare} to create paired (logCount, indicator)
          node structure.
    \item Learn DAGs from bootstrap resamples with \code{hcSC_boot}.
    \item Aggregate the bootstrap DAGs with \code{score_shd}.
  }
  For a single DAG without bootstrap, use \code{hcSC} directly.
}
\author{
Jie Peng

Maintainer: Jie Peng <jiepeng@ucdavis.edu>
}
\references{
  Chowdhury S, Ferri-Borgogno S, Yang P, Wang W, Peng J, Mok SC, Wang P.
  Learning directed acyclic graphs for ligands and receptors based on spatially
  resolved transcriptomic data of ovarian cancer.
  Briefings in Bioinformatics, 26(2), 2025. https://doi.org/10.1093/bib/bbaf085
}
\keyword{package}
