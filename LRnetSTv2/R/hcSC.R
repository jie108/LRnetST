.is_dag_matrix <- function(adj) {
  ## Kahn topological sort: returns TRUE iff adj is acyclic.
  p     <- nrow(adj)
  indeg <- colSums(adj != 0)
  queue <- which(indeg == 0)
  seen  <- 0L
  while (length(queue) > 0L) {
    node  <- queue[1L]; queue <- queue[-1L]; seen <- seen + 1L
    for (child in which(adj[node, ] != 0)) {
      indeg[child] <- indeg[child] - 1L
      if (indeg[child] == 0L) queue <- c(queue, child)
    }
  }
  seen == p
}

.validate_controls <- function(tol, maxStep, restart, seed) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("tol must be a finite nonnegative scalar")
  if (!is.numeric(maxStep) || length(maxStep) != 1L ||
      is.na(maxStep) || maxStep < 0)
    stop("maxStep must be a nonnegative integer")
  if (!is.numeric(restart) || length(restart) != 1L ||
      is.na(restart) || restart < 1)
    stop("restart must be a positive integer")
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0)
    stop("seed must be a nonnegative integer")
  list(maxStep = as.integer(maxStep),
       restart = as.integer(restart),
       seed    = as.integer(seed))
}

.validate_hc_inputs <- function(Y, nodeType, whiteList, blackList,
                                tol, maxStep, restart, seed) {
  ## Validate and normalise all inputs shared by hcSC() and hcSC_boot().
  ## Returns raw (unscaled) Y; L2 scaling is the caller's responsibility.
  if (!is.matrix(Y) && !is.data.frame(Y))
    stop("Y must be a numeric matrix or data frame")
  Y <- as.matrix(Y); storage.mode(Y) <- "double"
  if (!is.numeric(Y) || anyNA(Y) || any(!is.finite(Y)))
    stop("Y must contain only finite numeric values")
  p <- ncol(Y)
  if (p < 1L || nrow(Y) < 1L)
    stop("Y must have positive numbers of rows and columns")

  if (is.null(nodeType)) {
    nodeType <- rep("c", p)
  }
  if (!is.character(nodeType) || length(nodeType) != p ||
      anyNA(nodeType) || !all(nodeType %in% c("c", "b")))
    stop("nodeType must be a character vector of length ncol(Y) with entries 'c' or 'b'")
  for (i in which(nodeType == "b")) {
    if (!all(Y[, i] %in% c(0, 1)))
      stop("binary nodes must contain only 0/1 values")
  }

  if (is.null(whiteList)) {
    whiteList <- matrix(FALSE, p, p)
  }
  if (!is.matrix(whiteList) || !identical(dim(whiteList), c(p, p)) || anyNA(whiteList))
    stop("whiteList must be a p by p matrix with no NA values")
  whiteList <- whiteList != 0; dimnames(whiteList) <- NULL
  if (any(diag(whiteList)))
    stop("whiteList diagonal must be FALSE")
  if (any(whiteList & t(whiteList)))
    stop("whiteList cannot contain both directions of an edge")
  if (!.is_dag_matrix(whiteList))
    stop("whiteList must be acyclic")

  if (is.null(blackList)) {
    blackList <- matrix(FALSE, p, p)
    diag(blackList) <- TRUE
  }
  if (!is.matrix(blackList) || !identical(dim(blackList), c(p, p)) || anyNA(blackList))
    stop("blackList must be a p by p matrix with no NA values")
  blackList <- blackList != 0; dimnames(blackList) <- NULL
  diag(blackList) <- TRUE
  if (any(whiteList & blackList))
    stop("whiteList and blackList conflict")

  controls <- .validate_controls(tol, maxStep, restart, seed)
  list(Y         = Y,
       nodeType  = nodeType,
       whiteList = whiteList,
       blackList = blackList,
       maxStep   = controls$maxStep,
       restart   = controls$restart,
       seed      = controls$seed)
}

.validate_boot_density <- function(Y, bootDensityThre) {
  if (!is.numeric(bootDensityThre) || length(bootDensityThre) != 1L ||
      !is.finite(bootDensityThre) || bootDensityThre <= 0)
    stop("bootDensityThre must be a positive finite scalar")
  nonzero_rates <- apply(Y != 0, 2, mean)
  density_min   <- min(nonzero_rates)
  if (bootDensityThre >= density_min)
    stop(sprintf(
      paste0("bootDensityThre (%.4f) must be strictly less than the minimum",
             " column nonzero rate (%.4f)"),
      bootDensityThre, density_min))
}

.resolve_future_workers <- function(workers, n.boot, verbose) {
  available <- max(1L, as.integer(future::availableCores()) - 1L)
  requested <- if (is.null(workers)) available else as.integer(workers)
  capped    <- max(1L, min(requested, available, n.boot))
  if (verbose && capped < requested)
    message("Capping workers from ", requested, " to ", capped,
            " based on available cores and n.boot.")
  capped
}

.format_boot_result <- function(result, p, n.boot, output_type) {
  if (output_type == "array")
    return(array(as.numeric(unlist(result, use.names = FALSE)),
                 dim = c(p, p, n.boot)))
  freq <- matrix(0, p, p)
  for (i in seq_len(n.boot)) freq <- freq + result[[i]]
  freq <- freq / n.boot
  diag(freq) <- 0
  if (output_type == "freq") return(freq)
  list(
    adjacency = array(as.numeric(unlist(result, use.names = FALSE)),
                      dim = c(p, p, n.boot)),
    freq = freq
  )
}

hcSC <- function(Y, nodeType = NULL, whiteList = NULL, blackList = NULL,
                 scale = TRUE, tol = 1e-6, maxStep = 2000L,
                 restart = 1L, seed = 1L, verbose = FALSE) {
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
    stop("scale must be TRUE or FALSE")
  args <- .validate_hc_inputs(Y, nodeType, whiteList, blackList,
                               tol, maxStep, restart, seed)
  if (scale) {
    for (i in which(args$nodeType == "c")) {
      l2n <- sqrt(sum(args$Y[, i]^2) / nrow(args$Y))
      if (!is.finite(l2n) || l2n <= 0)
        stop("continuous nodes must have positive finite L2 norm when scale = TRUE")
      args$Y[, i] <- args$Y[, i] / l2n
    }
  }
  hcSC_(args$Y, args$nodeType, args$whiteList, args$blackList,
        tol, args$maxStep, args$restart, args$seed, verbose)
}


boot_dense <- function(Y, bootDensityThre) {
  ## Bootstrap resample of Y (with replacement) rejecting any sample where
  ## any column falls below bootDensityThre fraction of nonzero entries.
  n <- nrow(Y)
  repeat {
    s.pick <- sample(seq_len(n), n, replace = TRUE)
    Y.pick <- Y[s.pick, , drop = FALSE]
    if (!any(apply(Y.pick != 0, 2, mean) < bootDensityThre)) break
  }
  Y.pick
}

.fit_boot_one_SC <- function(b, Y, nodeType, whiteList, blackList,
                              tol, maxStep, restart, hc_seed,
                              boot_seed, node_perm, scale,
                              bootDensityThre, verbose) {
  ## Fit one bootstrap HC replicate.
  ## boot_seed controls R's RNG for the rejection-sampled bootstrap draw;
  ## hc_seed controls the C++ mt19937 for HC tie-breaking.
  ## L2 scaling is applied per bootstrap sample after row-resampling so that
  ## every replicate sees exactly unit-scale data (not approximately, which
  ## would result from scaling the full dataset once before resampling).
  set.seed(boot_seed)
  Y.B <- boot_dense(Y, bootDensityThre)

  if (scale) {
    n.B <- nrow(Y.B)
    for (i in which(nodeType == "c")) {
      l2n <- sqrt(sum(Y.B[, i]^2) / n.B)
      if (is.finite(l2n) && l2n > 0) Y.B[, i] <- Y.B[, i] / l2n
    }
  }

  p <- ncol(Y)
  if (!is.null(node_perm)) {
    node.index  <- order(node_perm)
    Y.B         <- Y.B[, node_perm, drop = FALSE]
    node.B      <- nodeType[node_perm]
    whiteList.B <- whiteList[node_perm, node_perm, drop = FALSE]
    blackList.B <- blackList[node_perm, node_perm, drop = FALSE]
  } else {
    node.index  <- seq_len(p)
    node.B      <- nodeType
    whiteList.B <- whiteList
    blackList.B <- blackList
  }

  res <- hcSC_(Y.B, node.B, whiteList.B, blackList.B,
               tol, maxStep, restart, hc_seed, verbose)
  res$adjacency[node.index, node.index, drop = FALSE]
}

hcSC_boot <- function(Y, n.boot = 100L, nodeType = NULL,
                      whiteList = NULL, blackList = NULL,
                      scale = TRUE, tol = 1e-6, maxStep = 2000L,
                      restart = 1L, seed = 1L,
                      nodeShuffle = FALSE, bootDensityThre = 0.1,
                      backend = c("sequential", "future"),
                      workers = NULL,
                      verbose = FALSE,
                      output_type = c("array", "freq", "both")) {
  ## bootstrap = "sequential": replicates run one by one in the current session.
  ## backend   = "future":     replicates dispatched as future::multisession.
  ##
  ## All randomness (bootstrap RNG seeds, HC tie-breaking seeds, node permutations)
  ## is generated upfront from a single set.seed(seed) call so that sequential
  ## and future backends produce bit-identical results.
  backend     <- match.arg(backend)
  output_type <- match.arg(output_type)

  if (!is.numeric(n.boot) || length(n.boot) != 1L ||
      is.na(n.boot) || n.boot < 1)
    stop("n.boot must be a positive integer")
  n.boot <- as.integer(n.boot)
  if (!is.logical(nodeShuffle) || length(nodeShuffle) != 1L ||
      is.na(nodeShuffle))
    stop("nodeShuffle must be TRUE or FALSE")
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
    stop("scale must be TRUE or FALSE")
  if (!is.null(workers) &&
      (!is.numeric(workers) || length(workers) != 1L ||
       is.na(workers) || workers < 1))
    stop("workers must be NULL or a positive integer")

  args <- .validate_hc_inputs(Y, nodeType, whiteList, blackList,
                               tol, maxStep, restart, seed)
  .validate_boot_density(args$Y, bootDensityThre)

  p <- ncol(args$Y)

  ## Save and restore the caller's global RNG state so set.seed() below does
  ## not have side effects outside this function.
  old_rng <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (is.null(old_rng)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_rng, envir = .GlobalEnv)
    }
  }, add = TRUE)

  ## Generate all per-replicate randomness upfront from one seed so that
  ## sequential and future backends visit the same bootstrap problems.
  set.seed(args$seed)
  boot_seeds <- sample.int(.Machine$integer.max, n.boot)
  hc_seeds   <- sample.int(.Machine$integer.max, n.boot)
  node_permutations <- if (nodeShuffle) {
    replicate(n.boot, sample.int(p), simplify = FALSE)
  } else {
    NULL
  }

  run_one <- function(b) {
    if (verbose) message("Processing bootstrap sample ", b)
    node_perm <- if (!is.null(node_permutations)) node_permutations[[b]] else NULL
    .fit_boot_one_SC(b, args$Y, args$nodeType, args$whiteList, args$blackList,
                     tol, args$maxStep, args$restart, hc_seeds[b],
                     boot_seeds[b], node_perm, scale, bootDensityThre, verbose)
  }

  if (backend == "sequential") {
    result <- vector("list", n.boot)
    for (b in seq_len(n.boot)) result[[b]] <- run_one(b)
  } else {
    if (!requireNamespace("future", quietly = TRUE))
      stop("The future package is required for backend = 'future'.")
    nw       <- .resolve_future_workers(workers, n.boot, verbose)
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = nw)
    futures <- lapply(seq_len(n.boot), function(b) {
      future::future(run_one(b), seed = TRUE)
    })
    result <- lapply(futures, future::value)
  }

  .format_boot_result(result, p, n.boot, output_type)
}


SC_prepare <- function(logCount) {
  ## Input:  n by p log-count matrix (log2(count + 1) or similar).
  ## Output: list(Y, nodeType, whiteList, blackList) for hcSC / hcSC_boot.
  ##   Y = cbind(logCount, logCount != 0): columns 1..p continuous, p+1..2p binary.
  ##   whiteList: indicator_i -> logcount_i for every i (always a parent).
  ##   blackList: logcount_i -> indicator_i forbidden; diagonal forbidden.
  if (!is.matrix(logCount) && !is.data.frame(logCount))
    stop("logCount must be a numeric matrix or data frame")
  logCount <- as.matrix(logCount); storage.mode(logCount) <- "double"
  if (anyNA(logCount) || any(!is.finite(logCount)))
    stop("logCount must contain only finite numeric values")
  if (nrow(logCount) < 1L || ncol(logCount) < 1L)
    stop("logCount must have positive numbers of rows and columns")
  if (any(logCount < 0))
    stop("logCount must be nonnegative (log-counts cannot be negative)")
  n <- nrow(logCount)
  p <- ncol(logCount)
  S <- matrix(as.numeric(logCount != 0), n, p)
  Y <- cbind(logCount, S)

  node.type <- c(rep("c", p), rep("b", p))

  whiteList <- matrix(FALSE, 2 * p, 2 * p)
  blackList <- matrix(FALSE, 2 * p, 2 * p)
  diag(blackList) <- TRUE
  for (i in seq_len(p)) {
    whiteList[p + i, i] <- TRUE
    blackList[i, p + i] <- TRUE
  }

  list(Y = Y, nodeType = node.type,
       whiteList = whiteList, blackList = blackList)
}
