score_shd <- function(boot.adj, alpha = 1, threshold = 0,
                      whitelist = NULL, blacklist = NULL,
                      max.step = NULL, verbose = FALSE) {
  ## Aggregate bootstrap DAGs via generalised SHD frequency.
  ##
  ## boot.adj: p by p by nb logical/integer array of bootstrap adjacency matrices
  ## alpha:    reversal weight in GSF: gsf(i,j) = sf(i,j) + (1 - alpha/2)*sf(j,i)
  ## threshold: freq.cut = (1 - threshold)/2; edges with GSF <= freq.cut are excluded
  ## whitelist/blacklist: p by p 0/1 matrices (whitelist edges always included,
  ##   blacklist edges never included); whitelist must be acyclic
  ## Returns: p by p 0/1 adjacency matrix

  p  <- dim(boot.adj)[1]
  nb <- dim(boot.adj)[3]

  if (is.null(whitelist))  whitelist  <- matrix(0, p, p)
  if (is.null(blacklist))  blacklist  <- matrix(0, p, p)

  res     <- whitelist
  par.set <- lapply(seq_len(p), function(i) which(res[, i] == 1) - 1L)

  sele.freq <- apply(boot.adj, c(1, 2), mean)
  gen.freq  <- sele.freq + (1 - alpha / 2) * t(sele.freq)

  freq.cut <- (1 - threshold) / 2

  indices   <- which(!is.na(sele.freq), arr.ind = TRUE)
  df        <- data.frame(i   = indices[, 1],
                          j   = indices[, 2],
                          SF  = sele.freq[indices],
                          GSF = gen.freq[indices])
  df_sorted   <- df[order(-df$GSF), ]
  df_filtered <- df_sorted[df_sorted$GSF > freq.cut, ]

  if (nrow(df_filtered) > 0) {
    for (k in seq_len(nrow(df_filtered))) {
      from.node <- df_filtered[k, "i"]
      to.node   <- df_filtered[k, "j"]

      if (verbose)
        message(k, "th operation: ", from.node, " -> ", to.node)

      if (!blacklist[from.node, to.node] && !whitelist[from.node, to.node]) {
        if (!edgeOnLoop(from.node - 1L, to.node - 1L, par.set)) {
          res[from.node, to.node]  <- 1
          par.set[[to.node]] <- c(par.set[[to.node]], from.node - 1L)
        } else if (verbose) {
          message("  not acyclic, skipped")
        }
      }
    }
  }

  res
}
