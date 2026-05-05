.validate_graph_constraint <- function(x, p, name, default = FALSE) {
  if (is.null(x)) {
    x <- matrix(default, p, p)
  }
  if (!is.matrix(x) || !identical(dim(x), c(p, p)) || anyNA(x)) {
    stop(name, " must be a p by p matrix with no NA values")
  }
  x <- x != 0
  dimnames(x) <- NULL
  x
}

.prepare_score_constraints <- function(p, whiteList, blackList) {
  whiteList <- .validate_graph_constraint(whiteList, p, "whiteList")
  blackList <- .validate_graph_constraint(blackList, p, "blackList")
  diag(blackList) <- TRUE
  list(whiteList = whiteList, blackList = blackList)
}

score_shd <- function(boot.adj, alpha = 1, freq.cutoff = 0.5,
                      whiteList = NULL, blackList = NULL,
                      max.step = NULL, verbose = FALSE) {
  ## Aggregate a bootstrap DAG array into a consensus DAG using generalised SHD.
  ##
  ## For each ordered pair (i, j), the generalised score is:
  ##   gsf(i, j) = sf(i, j) + (1 - alpha/2) * sf(j, i)
  ## where sf(i, j) = fraction of bootstrap DAGs containing edge i->j.
  ## A pair enters the candidate set iff gsf(i, j) > freq.cutoff.
  ## Candidates are added greedily in decreasing gsf order, skipping any
  ## edge that would create a cycle or violate a list constraint.
  if (!is.null(max.step)) {
    warning("max.step is deprecated and has no effect", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a positive finite scalar")
  }
  if (!is.numeric(freq.cutoff) || length(freq.cutoff) != 1L ||
      !is.finite(freq.cutoff) || freq.cutoff < 0 || freq.cutoff > 1) {
    stop("freq.cutoff must be a finite scalar in [0, 1]")
  }
  if (is.null(dim(boot.adj)) || length(dim(boot.adj)) != 3L ||
      dim(boot.adj)[1] != dim(boot.adj)[2]) {
    stop("boot.adj must be a p by p by B array")
  }
  p <- dim(boot.adj)[1]
  constraints <- .prepare_score_constraints(p, whiteList, blackList)
  score_shd_cpp(boot.adj, alpha, freq.cutoff,
                constraints$whiteList, constraints$blackList, verbose)
}

score_shd_freq <- function(freq, alpha = 1, freq.cutoff = 0.5,
                           whiteList = NULL, blackList = NULL,
                           max.step = NULL, verbose = FALSE) {
  ## Frequency-matrix variant of score_shd; pairs with hcSC_boot(output_type="freq")
  ## so large bootstrap runs do not need to retain all individual adjacency arrays.
  if (!is.null(max.step)) {
    warning("max.step is deprecated and has no effect", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a positive finite scalar")
  }
  if (!is.numeric(freq.cutoff) || length(freq.cutoff) != 1L ||
      !is.finite(freq.cutoff) || freq.cutoff < 0 || freq.cutoff > 1) {
    stop("freq.cutoff must be a finite scalar in [0, 1]")
  }
  if (!is.matrix(freq) || nrow(freq) != ncol(freq)) {
    stop("freq must be a square matrix")
  }
  freq <- as.matrix(freq)
  storage.mode(freq) <- "double"
  if (any(!is.finite(freq)) || any(freq < 0) || any(freq > 1)) {
    stop("freq must contain finite values in [0, 1]")
  }
  p <- nrow(freq)
  constraints <- .prepare_score_constraints(p, whiteList, blackList)
  score_shd_freq_cpp(freq, alpha, freq.cutoff,
                     constraints$whiteList, constraints$blackList, verbose)
}
