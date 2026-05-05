hcSC <- function(Y, nodeType = NULL, whiteList = NULL, blackList = NULL,
                 scale = TRUE, tol = 1e-6, maxStep = 2000L, restart = 10L,
                 seed = 1L, verbose = FALSE) {
  ## Y: n by p data matrix
  ## nodeType: length-p character vector, elements "c" (continuous) or "b" (binary)
  ## whiteList / blackList: p by p logical matrices
  ## scale: L2-normalise continuous columns so ||Y[,i]||^2/n = 1 (preserves zero pattern)
  ## tol: minimum score improvement for HC; maxStep: maximum HC steps
  ## restart: number of HC restarts (distinct tie-breaking seeds)
  ## seed: integer seed for tie-breaking RNG
  ## Returns: list from hcSC_ (adjacency matrix, scores, operations, deltaMin)

  p <- ncol(Y)
  n <- nrow(Y)

  if (is.null(whiteList)) {
    whiteList <- matrix(FALSE, p, p)
  } else {
    if (!is.matrix(whiteList) || !is.logical(whiteList) ||
        nrow(whiteList) != p || ncol(whiteList) != p)
      stop(paste("whiteList must be a", p, "by", p, "logical matrix"))
  }

  if (is.null(blackList)) {
    blackList <- matrix(FALSE, p, p)
    diag(blackList) <- TRUE
  } else {
    if (!is.matrix(blackList) || !is.logical(blackList) ||
        nrow(blackList) != p || ncol(blackList) != p)
      stop(paste("blackList must be a", p, "by", p, "logical matrix"))
  }

  if (is.null(nodeType)) {
    nodeType <- rep("c", p)
  } else {
    if (!is.character(nodeType) || length(nodeType) != p ||
        !all(nodeType %in% c("c", "b")))
      stop(paste("nodeType must be a length", p,
                 "character vector with elements 'c' or 'b'"))
  }

  if (scale) {
    for (i in seq_len(p)) {
      if (nodeType[i] == "c") {
        l2n <- sqrt(sum(Y[, i]^2) / n)
        if (l2n > 0) Y[, i] <- Y[, i] / l2n
      }
    }
  }

  hcSC_(Y, nodeType, whiteList, blackList, tol, maxStep, restart, seed, verbose)
}


hcSC_boot <- function(Y, n.boot = 1L, nodeType = NULL, whiteList = NULL,
                      blackList = NULL, scale = TRUE, tol = 1e-6,
                      maxStep = 2000L, restart = 10L, seed = 1L,
                      nodeShuffle = FALSE, bootDensityThre = 0.1,
                      n.thread = 1L, verbose = FALSE) {
  ## Y: n by p data matrix
  ## n.boot: number of bootstrap resamples
  ## scale: L2-normalise continuous columns (applied once, before bootstrap loop)
  ## bootDensityThre: minimum column-wise nonzero rate in any bootstrap resample
  ## n.thread: number of parallel workers (1 = sequential); uses future backend
  ## Returns: p by p by n.boot logical array of bootstrap adjacency matrices

  p <- ncol(Y)
  n <- nrow(Y)
  density.min <- min(apply(Y != 0, 2, mean))

  if (is.null(whiteList)) {
    whiteList <- matrix(FALSE, p, p)
  } else {
    if (!is.matrix(whiteList) || !is.logical(whiteList) ||
        nrow(whiteList) != p || ncol(whiteList) != p)
      stop(paste("whiteList must be a", p, "by", p, "logical matrix"))
  }

  if (is.null(blackList)) {
    blackList <- matrix(FALSE, p, p)
    diag(blackList) <- TRUE
  } else {
    if (!is.matrix(blackList) || !is.logical(blackList) ||
        nrow(blackList) != p || ncol(blackList) != p)
      stop(paste("blackList must be a", p, "by", p, "logical matrix"))
  }

  if (is.null(nodeType)) {
    nodeType <- rep("c", p)
  } else {
    if (!is.character(nodeType) || length(nodeType) != p ||
        !all(nodeType %in% c("c", "b")))
      stop(paste("nodeType must be a length", p,
                 "character vector with elements 'c' or 'b'"))
  }

  if (bootDensityThre <= 0 || bootDensityThre >= density.min)
    stop(paste("bootDensityThre must be strictly between 0 and",
               density.min, "(lowest column nonzero rate in data)"))

  ## Scale once on the original data; bootstrap resamples inherit the scaled columns.
  if (scale) {
    for (i in seq_len(p)) {
      if (nodeType[i] == "c") {
        l2n <- sqrt(sum(Y[, i]^2) / n)
        if (l2n > 0) Y[, i] <- Y[, i] / l2n
      }
    }
  }

  ## Generate all seeds upfront: boot b uses seed + b for R sampling,
  ## seed + n.boot + b for the C++ HC restarts.
  if (seed > .Machine$integer.max - 2L * n.boot)
    stop("seed is too large relative to n.boot; reduce seed or n.boot")
  boot_seeds  <- seed + seq_len(n.boot)       # R RNG per bootstrap
  hc_seeds    <- seed + n.boot + seq_len(n.boot)  # C++ HC per bootstrap

  run_one <- function(b) {
    set.seed(boot_seeds[[b]])
    Y.pick <- boot_dense(Y, bootDensityThre)

    if (nodeShuffle) {
      node.rand  <- sample(seq_len(p), p, replace = FALSE)
      node.index <- sort(node.rand, decreasing = FALSE, index.return = TRUE)$ix
    } else {
      node.rand  <- seq_len(p)
      node.index <- seq_len(p)
    }

    Y.B          <- Y.pick[, node.rand, drop = FALSE]
    node.B       <- nodeType[node.rand]
    whiteList.B  <- whiteList[node.rand, node.rand, drop = FALSE]
    blackList.B  <- blackList[node.rand, node.rand, drop = FALSE]

    res    <- hcSC_(Y.B, node.B, whiteList.B, blackList.B,
                    tol, maxStep, restart, hc_seeds[[b]], verbose)
    adjRes <- res$adjacency
    adjRes[node.index, node.index, drop = FALSE]
  }

  if (n.thread > 1L) {
    if (!requireNamespace("future", quietly = TRUE))
      stop("package 'future' is required for parallel bootstrap (n.thread > 1)")
    oplan <- future::plan(future::multisession, workers = n.thread)
    on.exit(future::plan(oplan), add = TRUE)

    futs <- vector("list", n.boot)
    for (b in seq_len(n.boot))
      futs[[b]] <- future::future(run_one(b), substitute = FALSE)
    results <- lapply(futs, future::value)
  } else {
    results <- vector("list", n.boot)
    for (b in seq_len(n.boot)) {
      if (verbose) message("bootstrap sample ", b, " ...")
      results[[b]] <- run_one(b)
    }
  }

  array(as.numeric(unlist(results)), dim = c(p, p, n.boot))
}


boot_dense <- function(Y, bootDensityThre) {
  ## Bootstrap resample of Y ensuring every column has at least bootDensityThre
  ## fraction of nonzero entries (rejection sampling).
  n <- nrow(Y)
  repeat {
    s.pick <- sample(seq_len(n), n, replace = TRUE)
    Y.pick <- Y[s.pick, , drop = FALSE]
    if (!any(apply(Y.pick != 0, 2, mean) < bootDensityThre)) break
  }
  Y.pick
}


SC_prepare <- function(logCount) {
  ## Input: log(count+1) transformed single-cell count matrix (n by p)
  ## Output: list(Y, nodeType, whiteList, blackList) ready for hcSC / hcSC_boot
  ## Data structure: Y = cbind(logCount, logCount > 0)
  ##   - columns 1..p  are continuous (logCount values)
  ##   - columns p+1..2p are binary  (on/off indicators)
  ##   - whiteList: (indicator_i -> logcount_i) always included
  ##   - blackList: (logcount_i -> indicator_i) always excluded

  n <- nrow(logCount)
  p <- ncol(logCount)
  S <- matrix(as.numeric(logCount != 0), n, p)
  Y <- cbind(logCount, S)

  node.type <- c(rep("c", p), rep("b", p))

  whiteList <- matrix(FALSE, 2 * p, 2 * p)
  blackList <- matrix(FALSE, 2 * p, 2 * p)
  diag(blackList) <- TRUE

  for (i in seq_len(p)) {
    whiteList[p + i, i] <- TRUE   # indicator_i -> logcount_i (whitelisted)
    blackList[i, p + i] <- TRUE   # logcount_i -> indicator_i (blacklisted)
  }

  list(Y = Y, nodeType = node.type,
       whiteList = whiteList, blackList = blackList)
}
